#include "messy_main_ppd_bi.inc"

MODULE messy_gmxe_si

!
! DESCRIPTION
! -----------
! GMXe INTERFACE LAYER FOR ECHAM5/MESSY
!
! AUTHOR
! ------
! Holger Tost, IPA, Johannes Gutenberg University, Mainz, Germany
! Swen Metzger, Max Planck Institute for Chemistry, Mainz, Germany
! questions/suggestions: holger.tost@mpic.de
!
! Copyright 2006-2009+. All rights reserved.
! LAST MODIFICATIONS - Dec 2008, Swen Metzger, Kirsty Pringle, Holger Tost
!                    - Aug 2010, Holger Tost
!                    - 2011,     Holger Tost
!                    - 2012,     Holger Tost
!                    - 2013,     Holger Tost

!*****************************************************************************

  USE messy_gmxe_mem
  USE messy_gmxe
  
!  USE messy_main_tracer,        ONLY: t_ident, t_med_aerosol, STRLEN_MEDIUM
  USE messy_main_tracer,        ONLY: t_ident
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  ! POINTER TO STREAM ELEMENTS
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_blather_bi,    ONLY: info_bi, warning_bi, error_bi
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle, mtend_register, &
                                      mtend_get_start_l, mtend_id_t, &
                                      mtend_id_q,  mtend_id_xl,      &
                                      mtend_id_xi, mtend_id_tracer,  &
                                      mtend_add_l
#endif

  IMPLICIT NONE
  PRIVATE

  ! SUBROUTINES
  PUBLIC :: gmxe_initialize              ! initialization
  PUBLIC :: gmxe_global_start            ! read input fields
  PUBLIC :: gmxe_new_tracer              ! define tracers
  PUBLIC :: gmxe_init_tracer             ! initialize tracers
  PUBLIC :: gmxe_init_memory             ! allocate memory
  PUBLIC :: gmxe_vdiff                   ! distribute online emissions
  PUBLIC :: gmxe_driver                  ! calculate gmxe-physics
  PUBLIC :: gmxe_radiation               ! calls gmxe_driver (set in gmxe.nml)
  PUBLIC :: gmxe_physc                   ! calls gmxe_driver (set in gmxe.nml)
  PUBLIC :: gmxe_free_memory             ! deallocate radius field
  PUBLIC :: gmxe_init_coupling           ! coupling to echam5

  INTRINSIC ABS, ASSOCIATED, ALLOCATED, MAX, MIN, TRIM

  REAL(dp), DIMENSION(:), POINTER         :: sigma_str => NULL()
  REAL(dp), DIMENSION(:), POINTER         :: crdiv_str => NULL()

  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: rwet_str  => NULL(), &
                                             rdry_str  => NULL(), &
                                             ddry_str  => NULL(), &
                                             anum_str  => NULL(), &
                                             ccn_str   => NULL(), &
                                             cph_str   => NULL()

  TYPE (PTR_3D_ARRAY), DIMENSION(:,:), ALLOCATABLE :: diagaer

  ! define radius and density 5D pointer
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: cph  => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: ccn  => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: anum => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: rwet => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: rdry => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: ddry => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: t_3D         => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: p_3D         => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: arho_3D      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: rh_3D        => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: sh_3D        => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zsat_3D      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zgh2o_3D     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zclh2o_3D    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zcih2o_3D    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zaopt_3D     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zaclc_3D     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zaclcac_3D   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: zaclcov_2D   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zcdnc_3D     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zicnc_3D     => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zcrain_3D    => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: zprec_2D    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zsnow_3D     => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5aclc_3D    => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5aclcac_3D  => NULL()
!!$  REAL(dp), DIMENSION(:,:),   POINTER :: e5aclcov_2D  => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5clh2o_3D   => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5cih2o_3D   => NULL()
!!$  REAL(dp), DIMENSION(:,:),   POINTER :: e5prec_2D    => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5snow_3D    => NULL()
!!$  REAL(dp), DIMENSION(:,:,:), POINTER :: e5cdnc_3D    => NULL()
! REAL(dp), DIMENSION(:,:,:), POINTER :: e5cinc_3D    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: val_strat_3D => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: zvap => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zvap2 => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zprod => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: zminc => NULL()

  ! switch for tracer initialisation
  LOGICAL :: tracer_init_required = .false.

  CHARACTER(LEN=100)                      :: ctracername =''  ! tracer longname
  CHARACTER(LEN=2)                        :: csubname   = ''  ! tracer subname
   
  !--- CPL namelist fields:

  LOGICAL :: l_calc_emis           = .FALSE.
  LOGICAL :: l_tendency            = .FALSE.

  CHARACTER (LEN=STRLEN_MEDIUM), PUBLIC ::  &
                                 Tropopchannel= '', &
                                 TropopIndex = '', &
                                 Pscchannel   = '', &
                                 Pscreg      = '', &
                                 phase       = ''
  CHARACTER (LEN=STRLEN_MEDIUM), PUBLIC ::  &               !! mz_kp_20080112
                                 driver_call
 
  !-----
  ! 3D-field (from channel "hetchem") of Stratosphere region indicators
  !-----
  REAL(dp), DIMENSION(:,:,:), POINTER :: flt_stratreg => NULL(), &
                                         flt_pscreg0  => NULL(), &
                                         flt_pscreg   => NULL(), &
                                         flt_pscreg1  => NULL()
  !-----
  ! 2d-field for tropopause index (details via namelist)
  !-----
  REAL(dp), DIMENSION(:,:), POINTER :: tp_i0 => NULL()

  ! info about PSC channel / region (if available, ierr_psc = 0)
  INTEGER, SAVE :: ierr_psc,ierr_tropop

  INTEGER, SAVE :: idt_H2O

  LOGICAL,  DIMENSION(:),     ALLOCATABLE :: lrhoa            ! true if air density is needed for unit conversion
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zconvert,      & ! factor to convert unit including air density
                                             xconvert         ! factor to convert unit excluding air density
  
!-------------------------------------------------------------------------------------------
!   for emissions

  TYPE emspec
    CHARACTER(LEN=STRLEN_MEDIUM) :: name      ! tracer name
    INTEGER                      :: trac_idx  ! tracer index
    REAL(dp)                     :: molarmass ! molar mass
    REAL(dp)                     :: frac      ! fraction of the total emission flux used
    INTEGER                      :: mode      ! mode of the target species
    LOGICAL                      :: L_numb    ! species is a number concentration and not a mass
  END TYPE emspec

  TYPE emflux
    ! name for each flux (fix, helps for identification)
    CHARACTER(LEN=STRLEN_MEDIUM)        :: name
    ! 3D-array for the mass flux
    REAL(dp), DIMENSION(:,:,:), POINTER :: flux          => NULL()
    ! 2D-array for the mass flux
    REAL(dp), DIMENSION(:,:),   POINTER :: flux_2D       => NULL()
    ! 3D-array for the corresponding number flux (if it exists)
    REAL(dp), DIMENSION(:,:,:), POINTER :: nflux         => NULL()
    ! 2D-array for the corresponding number flux (if it exists)
    REAL(dp), DIMENSION(:,:),   POINTER :: nflux_2D      => NULL()
    ! 3D-array for the V(ertical)IND(ex) for NxD emissions
    REAL(dp), DIMENSION(:,:,:), POINTER :: vind          => NULL()
    ! density = native density of the emission flux
    REAL(dp)                            :: density
    ! total_frac = total scaling factor for the emission flux
    ! controlled via the parameters.inc file
    REAL(dp)                            :: total_frac
    ! num_spec_emis = number of species which get a tendency from this flux
    INTEGER                             :: num_spec_emis
    ! dim = used dimension of the emission flux array (2D,3D)
    INTEGER                             :: dim
    ! dim_orig = native dimension of the emission flux array (2D,3D)
    INTEGER                             :: dim_orig
    ! NxD = logical whether a 3D flux originates from a NxD flux
    LOGICAL                             :: NxD
    ! mode = native mode of the flux
    INTEGER                             :: mode
    ! scal_fac = scaling factor if for a flux another flux is used and scaled 
    REAL(dp)                            :: scal_fac                   
    ! diameter = value for the aerosol diameter associated with this flux
    !            used in case of determining the number from the mass flux
    REAL(dp)                            :: diameter
    ! fac_num_emis = conversion factor (depending on the emission flux) to convert
    !                mass mean to count median (also including density if required)
    REAL(dp)                            :: fac_num_emis
    ! unit = unit of the emission flux -> determines conversion of emission
    CHARACTER(LEN=STRLEN_MEDIUM)        :: unit
    ! flux_name = name of the corresponding channel element
    CHARACTER(LEN=STRLEN_OBJECT)        :: flux_name
    ! nflux_name = name of the corresponding number flux channel element
    CHARACTER(LEN=STRLEN_OBJECT)        :: nflux_name
    ! channel_name = name of the corresponding channel
    CHARACTER(LEN=STRLEN_CHANNEL)       :: channel_name
    TYPE(emspec), DIMENSION(:), POINTER :: specs         => NULL()
  END TYPE emflux

  TYPE(emflux), DIMENSION(:),   POINTER, SAVE :: emis_flux_array => NULL()
  TYPE(emflux), SAVE                          :: emis_flux_list(100)
  INTEGER, SAVE                               :: num_fluxes

  ! EMIS_CASK = Character array which contains for each of maximum 500 emission 
  !              fields the necessary information for emission assignment:
  !          1 = name of the emission object (for identification)
  !          2 = channel / channel name of emission flux
  !          3 = channel object name of mass emission flux
  !          4 = corresponding name emission flux (if exists)
  !          5 = list of tracers which should get a value from this emission
  !              ";" separated list of tracers (fullname, CASE SENSITIVE)
  !          6 = list of scaling factors for each tracer
  !              ";" separated list of REAL values
  ! mz_dk_20120203+
  CHARACTER(LEN=STRLEN_ULONG), DIMENSION(100,7), SAVE :: EMIS_CASK 
  ! mz_dk_20120203-
  
  !mz_sb_20161020+ coupling with oracle
  REAL(DP), POINTER :: nomode   => NULL() ! number of modes in ORACLE ! mz_se 20170313 pointed to NULL()
  REAL(DP), POINTER :: NfPOA    => NULL()  
  REAL(DP), POINTER :: NbbPOA   => NULL()  
  REAL(DP), POINTER :: NfSOAsv  => NULL() 
  REAL(DP), POINTER :: NbbSOAsv => NULL()    
  REAL(DP), POINTER :: NfSOAiv  => NULL()  
  REAL(DP), POINTER :: NbbSOAiv => NULL()   
  REAL(DP), POINTER :: NSOAv    => NULL()     
  REAL(dp), DIMENSION(:), POINTER :: tomode => NULL() ! type of modes in ORACLE 
  REAL(dp), DIMENSION(:), POINTER :: kMOM   => NULL() ! number of MOM species in each volatility bin
  !mz_sb_20161020- coupling with oracle

  NAMELIST /CPL/ l_calc_emis,   l_tendency,    &
                 Tropopchannel,  TropopIndex,                      & 
                 Pscchannel,   Pscreg,  phase,                     &
                 driver_call,                                     &
                 EMIS_CASK

  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: photol(1) = 'jval_gp'

  ! coupling to main data
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grmass     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: grvol      => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: rhum_3d    => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: press_3d   => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: pressi_3d  => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: deltaz     => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: aclc       => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc      => NULL()
#if defined(COSMO) || defined(MESSYDWARF)
  REAL(dp), POINTER, DIMENSION(:,:,:) :: vervel_3D  => NULL()
#endif
  REAL(dp), POINTER, DIMENSION(:,:)   :: vervel     => NULL()

#ifdef MESSYTENDENCY
  INTEGER                             :: my_handle_td ! thermodynamics
  INTEGER                             :: my_handle_em ! emissions
#endif

CONTAINS
!==============================================================================

  SUBROUTINE gmxe_initialize
  !
    USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, p_pe
    USE messy_gmxe_aerchem,      ONLY: init_aerchem
    USE messy_gmxe_soa,          ONLY: set_soa_params, l_gasprec, l_oxi, &
                                       ngas_max_SOA => ngas_max, noxi_MAX, &
                                       excl_str_soa, exclude_max_soa,    &
                                       noxi, ngas_soa => ngas, nsoa
    USE messy_gmxe_isorropia2,   ONLY: init_species_isorropia
    USE messy_gmxe_eqsam4clim,   ONLY: init_species_eqsam4clim


    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'gmxe_initialize' ! name of subroutine
    INTEGER                         :: iou    ! I/O unit
    INTEGER                         :: status ! error status
    INTEGER                         :: jm, jc   

    IF (p_parallel_io) &
    CALL start_message_bi(modstr,'INITIALIZATION',substr)

    ! Initialise namelist variables
    L_GASPREC(:,:) = .FALSE.
    L_OXI(:)       = .FALSE.
    
    !--- Read namelist and control variables:

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL gmxe_read_nml_ctrl(status, iou)
       IF (status /= 0)  CALL error_bi('error in gmxe_read_nml_ctrl', substr)
    END IF
    
    !--- Broadcast over processors:

    CALL p_bcast (lgmxe,      p_io)
    CALL p_bcast (loutput,    p_io)
    CALL p_bcast (lmass_diag, p_io)
    CALL p_bcast (lstrat,     p_io)
    CALL p_bcast (lpsc,       p_io)
    CALL p_bcast (lnucl,      p_io)
    CALL p_bcast (lcoat,      p_io)
    CALL p_bcast (lsize,      p_io)
    CALL p_bcast (lcoag,      p_io)
    CALL p_bcast (lcond,      p_io)
    CALL p_bcast (ladyn,      p_io)
    CALL p_bcast (lah2o,      p_io)
    CALL p_bcast (lgh2o,      p_io)
    CALL p_bcast (lshum,      p_io)
    CALL p_bcast (lclc,       p_io)
    CALL p_bcast (lacc,       p_io)
    CALL p_bcast (lclwc,      p_io)
    CALL p_bcast (lciwc,      p_io)
    CALL p_bcast (lcdnc,      p_io)
    CALL p_bcast (licnc,      p_io)
    CALL p_bcast (lgas,       p_io)
    CALL p_bcast (laerosol,   p_io)
    CALL p_bcast (lnumber,    p_io)
    CALL p_bcast (lwetrad,    p_io)
    CALL p_bcast (ldryrad,    p_io)
    CALL p_bcast (ldrydens,   p_io)
    CALL p_bcast (nlowermode, p_io)
    CALL p_bcast (nuppermode, p_io)
    CALL p_bcast (neqm,       p_io)
    CALL p_bcast (nnucl,      p_io)

    CALL p_bcast (ldry,       p_io)
    CALL p_bcast (lhyster,    p_io)
    
    CALL p_bcast (l_aerchem,   p_io)
    CALL p_bcast (l_oc_aging,  p_io)
    CALL p_bcast (l_soa,       p_io)
    CALL p_bcast (l_oracle,    p_io) !mz_sb_20161020
    ! for discretisation
    CALL p_bcast (nmod,        p_io)
    CALL p_bcast (sigma_nml,   p_io)
    CALL p_bcast (crdiv_nml,   p_io)
    CALL p_bcast (cmodes_nml,  p_io)
    ! for species
    DO jm=1,ngas_max
      CALL p_bcast (cask_gases(jm,1),  p_io)
      CALL p_bcast (cask_gases(jm,2),  p_io)
    END DO
    DO jm=1,nanions_max
      CALL p_bcast (cask_anions(jm,1),  p_io)
      CALL p_bcast (cask_anions(jm,2),  p_io)
    END DO
    DO jm=1,ncations_max
      CALL p_bcast (cask_cations(jm,1),  p_io)
      CALL p_bcast (cask_cations(jm,2),  p_io)
    END DO
    DO jm=1,nsolutes_max
      CALL p_bcast (cask_solutes(jm,1),  p_io)
      CALL p_bcast (cask_solutes(jm,2),  p_io)
    END DO
    DO jm=1,nbulk_max
      CALL p_bcast (cask_bulk(jm,1),  p_io)
      CALL p_bcast (cask_bulk(jm,2),  p_io)
    END DO
    ! for aerchem
    CALL p_bcast (umode,       p_io)
    CALL p_bcast (lmode,       p_io)
    ! for passive tracers
    CALL p_bcast (l_passive_aer, p_io) ! mz_dk_20120119
    CALL p_bcast (num_pa,     p_io)    ! mz_dk_20120119
    CALL p_bcast (pamode1,    p_io)    ! mz_dk_20120119
    CALL p_bcast (pamode2,    p_io)    ! mz_dk_20120119
    ! for SOA
    
    CALL p_bcast (lmode_SOA, p_io) 
    CALL p_bcast (umode_SOA, p_io) 
    DO jm=1,ngas_max_SOA
       DO jc=1,2
          CALL p_bcast (L_gasprec(jm,jc),p_io)
       END DO
    END DO
    DO jm=1,noxi_max
       CALL p_bcast (L_OXI(jm), p_io)  
    ENDDO
    DO jm=1,exclude_max_soa
       call p_bcast(excl_str_soa(jm), p_io)
    END DO

   !--- Read CPL namelist
    IF (p_parallel_io) THEN
      EMIS_CASK(:,:) =''
      iou = find_next_free_unit(100,200)
      CALL gmxe_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi("error in coupling namelist", substr)
    END IF

    !--- Broadcast over processors:

    CALL p_bcast(l_calc_emis,   p_io)
    CALL p_bcast(l_tendency,    p_io)
    
    CALL p_bcast(Tropopchannel,  p_io)
    CALL p_bcast(TropopIndex,   p_io)
    CALL p_bcast(Pscchannel,     p_io)
    CALL p_bcast(Pscreg,        p_io)
    CALL p_bcast(phase,         p_io)
    CALL p_bcast(driver_call,   p_io)
    
    DO jm=1,100 ! mz_dk_20120203
      DO jc=1,7
        CALL p_bcast(EMIS_CASK(jm,jc),p_io)
      END DO
    END DO
    
    l_io=.FALSE.
    IF(p_parallel_io .AND. loutput) l_io=.TRUE.
!--- Initialize core:
    CALL gmxe_initialize_core(status)
    IF (status /= 0) CALL error_bi("error in initialisation", substr)

    IF (l_aerchem) THEN
      DO jm=1,nmod
        IF ( (JM >= lmode) .AND. (JM <= umode) ) AERCHEM(JM) = .TRUE.
      END DO
      CALL init_aerchem
    END IF

    IF (L_soa) THEN
      IF (L_GASPREC(7,1) .OR. L_GASPREC(8,1)) L_GASPREC(6,1) =.TRUE.
      DO jm=1,nmod
        IF ( (JM >= lmode_soa) .AND. (JM <= umode_soa) ) LSOA(JM) = .TRUE.
      END DO
      call set_soa_params
    ENDIF

    CALL gmxe_initialize_species(p_pe, p_io)

    SELECT CASE (neqm)
    CASE(1)
       CAll init_species_eqsam4clim(1)    
    CASE(2)
       CAll init_species_isorropia(1)
    END SELECT



#ifdef MESSYTENDENCY
    my_handle_td = mtend_get_handle(modstr//'_td')
    my_handle_em = mtend_get_handle(modstr//'_em')
#endif

    IF (p_parallel_io) &
    CALL end_message_bi(modstr,'INITIALIZATION',substr)
    
  END SUBROUTINE gmxe_initialize

!==============================================================================

  SUBROUTINE gmxe_new_tracer

    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer,        ONLY: new_tracer, get_tracer, set_tracer,   &
                                        AIR, ON, OFF, MODAL, AEROSOL,         &
                                        AMOUNTFRACTION, NUMBERDENSITY,        &
                                        I_ADVECT, I_CONVECT,                  &
                                        I_VDIFF,                              &
                                        I_DRYDEP, I_SEDI,                     &
                                        I_SCAV, I_MIX, I_CHARGE,              &
                                        I_AEROSOL_METHOD, I_AEROSOL_MODE,     &
                                        I_AEROSOL_SOL, S_AEROSOL_MODEL,       &
                                        R_MOLARMASS, R_AEROSOL_DENSITY,       &
                                        R_Henry_T0, R_henry_Tdep,             &
                                        R_alpha_T0, R_alpha_Tdep,             &
                                        get_chemprop
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,   ONLY: ti_gp, GPTRSTR
    USE messy_main_mpi_bi,          ONLY: p_parallel_io, p_pe, p_io
    USE MESSY_MAIN_TOOLS,         ONLY: strcrack
    USE MESSY_GMXE_AERCHEM_LIQ
    USE MESSY_GMXE_AERCHEM_KPP,   ONLY: str_field_kpp_l => spc_names
    USE MESSY_GMXE_OC_AGING,      ONLY: num_wsoc
    USE MESSY_GMXE_SOA,           ONLY: ngas_soa => ngas, nsoa, noxi, mw,     &
                                        names

    IMPLICIT NONE

  ! LOCAL
    INTEGER :: status

    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_new_tracer'
    
    INTEGER                     :: i,iz,it,idt,jc,jm, idx1, jt
    INTEGER                     :: ierr, imedium, nquantity, isol, nmedium
    CHARACTER(LEN=10)           :: cunit, cname
    INTEGER                     :: imode, dummy, idtoc
    CHARACTER(LEN=26), POINTER  :: strname(:) => NULL()
    CHARACTER(LEN=26)           :: c_name
    CHARACTER(LEN=2)            :: str_wsoc, str_pa ! mz_dk_20120118
    
    CALL start_message_bi(modstr,'Tracer definition',substr)
    
  
    nquantity = AMOUNTFRACTION
    cunit='mol/mol'
    nmedium=AIR

    DO jm=0,nmod
      nmedium=AIR
      csubname = ''
      IF (jm > 0) THEN
        csubname=TRIM(cmodes(jm))
        nmedium=AEROSOL
      END IF
      isol = 1
      IF (JM > nsoluble) isol = 0
      
      DO jc=1,ngas
        IF (.NOT. td%gas(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%gas(jc)%name)
        CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
          medium=imedium, idx=idt)
        IF(p_parallel_io .AND. ierr == 0) &
          WRITE(*,*) ' Tracer >> ',TRIM(cname)//'_'//TRIM(csubname),' << already defined!',idt
        IF (ierr == 0) CYCLE
        CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
          modstr, subname=csubname, quantity = nquantity,   &
          unit = TRIM(cunit), medium = nmedium, idx = idt)

        CALL tracer_halt(substr,status)
        CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)

        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
        CALL tracer_halt(substr, status)
     END DO

     tracer_init_required = .true.

      DO jc=1,nanions
        IF (jm == 0) CYCLE
        IF (.NOT. td%anion(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%anion(jc)%name)
        
        CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
          medium=imedium, idx=idt)
        IF(p_parallel_io .AND. ierr == 0) &
          WRITE(*,*) ' Tracer >> ',TRIM(cname)//'_'//TRIM(csubname),' << already defined!',idt
        IF (ierr == 0) CYCLE
        CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
          modstr, subname=csubname, quantity = nquantity,   &
          unit = TRIM(cunit), medium = nmedium, idx = idt)
        CALL tracer_halt(substr,status)
        CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)

        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
        CALL tracer_halt(substr, status)
      END DO

      DO jc=1,ncations
        IF (jm == 0) CYCLE
        IF (.NOT. td%cation(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%cation(jc)%name)
        CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
          medium=imedium, idx=idt)
        IF(p_parallel_io .AND. ierr == 0) &
          WRITE(*,*) ' Tracer >> ',TRIM(cname)//'_'//TRIM(csubname),' << already defined!',idt
        IF (ierr == 0) CYCLE

        CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
          modstr, subname=csubname, quantity = nquantity,   &
          unit = TRIM(cunit), medium = nmedium, idx = idt)
        CALL tracer_halt(substr,status)

        CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)

        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
        CALL tracer_halt(substr, status)
      END DO
 
      DO jc=1,nsolutes
! op_pj_20181029+ bug-fix to avoid creation of H2O tracer
!!$     IF (jm > nsoluble) CYCLE
        IF ((jm > nsoluble).OR.(jm==0)) CYCLE
! op_pj_20181029-
        IF (.NOT. td%solute(jc)%ltreat(jm)) CYCLE
        cname = TRIM(td%solute(jc)%name)
        CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
          medium=imedium, idx=idt)
        IF(p_parallel_io .AND. ierr == 0) &
          WRITE(*,*) ' Tracer >> ',TRIM(cname)//'_'//TRIM(csubname),' << already defined!',idt
        IF (ierr == 0) CYCLE
        CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
          modstr, subname=csubname, quantity = nquantity,   &
          unit = TRIM(cunit), medium = nmedium, idx = idt)
        CALL tracer_halt(substr,status)
        CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)

        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
        CALL tracer_halt(substr, status)
      END DO

      DO jc=1,nbulk
        IF (jm == 0) CYCLE
        IF (.NOT. bulk(jc)%ltreat(jm)) CYCLE
        cname = TRIM(bulk(jc)%name)
        CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
          medium=imedium, idx=idt)
        IF(p_parallel_io .AND. ierr == 0) &
          WRITE(*,*) ' Tracer >> ',TRIM(cname)//'_'//TRIM(csubname),' << already defined!',idt
        IF (ierr == 0) CYCLE
        CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
          modstr, subname=csubname, quantity = nquantity,   &
          unit = TRIM(cunit), medium = nmedium, idx = idt)
        CALL tracer_halt(substr,status)
        CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , &
          bulk(jc)%molmass)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)

        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, R_aerosol_density , &
          bulk(jc)%density*1.e3_dp )
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
        CALL tracer_halt(substr, status)
        CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
        CALL tracer_halt(substr, status)
      END DO

    END DO

    DO jm=1,nmod
      csubname=TRIM(cmodes(jm))
      nmedium=AEROSOL
      isol = 1
      IF (JM > nsoluble) isol = 0
      nquantity = Numberdensity
      cname = "N"
      cunit='1/mol'
      CALL get_tracer(ierr, GPTRSTR, TRIM(cname), subname=csubname, &
        medium=imedium, idx=idt)
      IF(p_parallel_io .AND. ierr == 0) &
        WRITE(*,*) ' Tracer >> ',TRIM(ctracername),' << already defined!',idt
      IF (ierr == 0) CYCLE
      CALL new_tracer(status, GPTRSTR, TRIM(cname),       &
        modstr, subname=csubname, quantity = nquantity,   &
        unit = TRIM(cunit), medium = nmedium, idx = idt)
      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , 1._dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , TRIM(modstr))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, R_aerosol_density , 1._dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
      CALL tracer_halt(substr, status)
      
    END DO



! add tracers for various sub - submodels , e.g. aerchem
  
    IF (L_aerchem) THEN

      nspec_aer = 0              ! counter for additional aerosol phase compounds 
                               ! in paerml array
      do jm = 1, nsoluble
        IF (.NOT.AERCHEM(jm) ) CYCLE
        ! use all modes in which aerchem shall be applied
        ! set subname according to the actual modes
        csubname=TRIM(cmodes(jm))
        imode = jm
        do jt = 1,lspec_liq
          idx1 = kpp_l_idx%liq_spec(jt, liq_idx)
          ! check if a tracer with the same name already exists for each mode 
          ! treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          IF (TRIM(strname(1)) == "Prod") CYCLE
          c_name=TRIM(strname(1))
          !        print*, "new tracer aer: ",  TRIM(c_name), csubname
          CALL get_tracer(ierr, GPTRSTR, TRIM(c_name), &
            subname = csubname, idx = idt)
          IF (ierr /= 0) THEN
          ! Define new tracer
!          print*, TRIM(c_name)
            CALL new_tracer(status, GPTRSTR, TRIM(c_name),  &
              modstr, subname=csubname, quantity = AMOUNTFRACTION,      &
              unit = 'mol/mol', medium = AEROSOL, idx = idt)
          
            CALL tracer_halt(substr,status)
            CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , imode)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , &
              TRIM(modstr))
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , 1)
            CALL tracer_halt(substr, status)
            status = get_chemprop(TRIM(c_name), R_molarmass &
               , KPP_L_IDX%LIQ_ATTR(JT,liq_MW))
            IF (status /=0 )  CALL error_bi (&
                 'species '//TRIM(c_name)//' not part of chemprop'&
                 , substr)
            nspec_aer = nspec_aer + 1
          END IF
        END do
      ENDDO
      
      ! check for gas phase compounds necessary to guarantee mass conservation
      ! for kpp chemistry
      nspec_gas = 0               ! counter for additional gas phase compounds 
      ! in paerml array
      do jt = 1, lspec_gas  
        !      print*, "new tracer gas: ", str_field_kpp_l
        IDX1 = kpp_l_idx%gas_spec(jt,gas_idx)
        CALL get_tracer(ierr, GPTRSTR, TRIM(str_field_kpp_l(idx1)), idx=idt)
        
        IF (ierr /= 0) THEN
          !        print*,  TRIM(str_field_kpp_l(idx1))
          CALL new_tracer(status, GPTRSTR, TRIM(str_field_kpp_l(idx1)), &
            modstr, quantity=AMOUNTFRACTION, unit='mol/mol',     &
            medium = AIR, idx = idt)
          
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_drydep          , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_sedi            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
          CALL tracer_halt(substr, status)
          
          nspec_gas = nspec_gas + 1

        END IF
      END do
    END IF

    IF (L_OC_AGING) THEN
      ! tracers should be initialised via the bulk species already
!!$     DO jm=1,nsoluble
!!$        csubname=TRIM(cmodes(jm))
!!$        CALL get_tracer(ierr, GPTRSTR, "OC", subname=TRIM(csubname), idx=idtoc)
!!$        IF (IERR == 0) THEN
!!$           DO jc = 1,num_wsoc
!!$              IF (jc < 10) then
!!$                 write (str_wsoc,'(A,I1)') "0",jc
!!$              ELSE
!!$                 write (str_wsoc,'(I2)') jc
!!$              END IF
!!$              CALL get_tracer(ierr, GPTRSTR, "WSOC"//str_wsoc, &
!!$                subname=TRIM(csubname), idx=idt)
!!$              IF (IERR /= 0 ) THEN
!!$                 CALL new_tracer(status, GPTRSTR, "WSOC"//str_wsoc,  &
!!$                      modstr, subname=TRIM(csubname), quantity = AMOUNTFRACTION,&
!!$                      unit = 'mol/mol', medium = AEROSOL, idx = idt)             
!!$                 CALL tracer_halt(substr,status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
!!$                 CALL tracer_halt(substr, status)  
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , &
!!$                      tbulk(MS,iOC))
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , &
!!$                      TRIM(modstr))
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, R_aerosol_density , &
!!$                      tbulk(De,iOC)*1.e3_dp ) 
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
!!$                 CALL tracer_halt(substr, status)
!!$                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , 1)
!!$                 CALL tracer_halt(substr, status)
!!$              END IF
!!$           END DO
!!$        END IF
!!$     END DO
    END IF

! mz_dk_20120118+
! TRACER initialization for PASSIVE AEROSOLS
    passive:IF (L_PASSIVE_AER) THEN
      mode_loop: DO jm=pamode1,pamode2
        csubname=TRIM(cmodes(jm)) ! mode name, e.g. ki, etc.
        isol = 1
        IF ( jm > nsoluble ) isol = 0 ! if insoluble tracer
           DO jc = 1,num_pa
              IF (jc < 10) then
                 write (str_pa,'(A,I1)') "0",jc
              ELSE
                 write (str_pa,'(I2)') jc
              END IF
              CALL get_tracer(ierr, GPTRSTR, "PASSAER"//str_pa, &
                subname=TRIM(csubname), idx=idt)
              IF (IERR /= 0 ) THEN
                 CALL new_tracer(status, GPTRSTR, "PASSAER"//str_pa,  &
                      modstr, subname=TRIM(csubname), quantity = AMOUNTFRACTION,&
                      unit = 'mol/mol', medium = AEROSOL, idx = idt)
                 CALL tracer_halt(substr,status)
                 CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
                 CALL tracer_halt(substr, status)  
                 CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_drydep          , ON)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_sedi            , ON)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , 1._dp)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_mode    , jm)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, S_aerosol_model   , &
                      TRIM(modstr))
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, R_aerosol_density , 1000._dp)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_method  , MODAL)
                 CALL tracer_halt(substr, status)
                 CALL set_tracer(status, GPTRSTR, idt, I_aerosol_sol     , isol)
                 CALL tracer_halt(substr, status)
              END IF
           END DO
     END DO mode_loop
  END IF passive
! mz_dk_20120118-                 

  IF (L_SOA) THEN
     ! gaseous species first
     DO jt=1,NSOA ! gaseous SOA compounds
        CALL get_tracer(ierr, GPTRSTR, TRIM(names(jt)), idx=idt)
        IF (ierr /= 0) THEN

           CALL new_tracer(status, GPTRSTR, TRIM(names(jt)),      &
                modstr, quantity=AMOUNTFRACTION, unit='mol/mol',     &
                medium = AIR, idx = idt)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
           CALL tracer_halt(substr, status)
           ! beware drydep for this species is not included
           CALL set_tracer(status, GPTRSTR, idt, I_drydep          , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_sedi            , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , MW(jt))
           CALL tracer_halt(substr, status)
        END IF
     END DO
     DO jt=1,NGAS_SOA ! gaseous SOA precursors
        CALL get_tracer(ierr, GPTRSTR, TRIM(names(NSOA + NOXI + jt)), idx=idt)
        IF (ierr /= 0) THEN

           CALL new_tracer(status, GPTRSTR, TRIM(names(NSOA + NOXI + jt)),      &
                modstr, quantity=AMOUNTFRACTION, unit='mol/mol',     &
                medium = AIR, idx = idt)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_advect          , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_convect         , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_vdiff           , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_scav            , ON)
           CALL tracer_halt(substr, status)
           ! beware drydep for this species is not included
           CALL set_tracer(status, GPTRSTR, idt, I_drydep          , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_sedi            , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, I_mix             , OFF)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, GPTRSTR, idt, R_molarmass       , MW(jt))
           CALL tracer_halt(substr, status)
        END IF
     END DO
     ! aerosol tracers should be added in bulk structure
  END IF

  CALL end_message_bi(modstr,'Tracer definition',substr)

  END SUBROUTINE gmxe_new_tracer

!=============================================================================!

  SUBROUTINE gmxe_init_tracer

!!  USE messy_main_tracer_bi, ONLY: tracer_init

    IMPLICIT NONE

!!  IF (tracer_init_required) CALL tracer_init(modstr)

  END SUBROUTINE gmxe_init_tracer

!=============================================================================

  SUBROUTINE gmxe_init_memory

    USE messy_main_tools,            ONLY: int2str, strcrack, match_wild
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, ngpblks, nproma
    USE messy_main_timer,            ONLY: lstart
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, DC_GP  &
                                         , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                         , DC_BC &
                                         , gp_nseg, gp_start, gp_cnt &
                                         , gp_meml, gp_memu &
                                         , GP_2D_HORIZONTAL
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO  &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL        &
                                           , repr_def_axes

    USE MESSY_GMXE_AERCHEM_LIQ,        ONLY: nprod, diag_aerchem, lspec
    USE MESSY_GMXE_AERCHEM_KPP,        ONLY: str_field_kpp_l => spc_names

    IMPLICIT NONE

  !--- Declare dimensions for the channels: -------------------------------------!

  ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'gmxe_init_memory'
    CHARACTER(LEN=15)                     :: name, wl

    INTEGER                               :: jc, jm, n, i, jt
    
    CHARACTER(LEN=1)                      :: char1
    CHARACTER(LEN=2)                      :: char2

    ! auxiliary pointer for data managment
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  p4 => NULL()
    REAL(dp), DIMENSION(:,:,:),     POINTER ::  p3 => NULL()

    INTEGER                               :: jmod, jsol, status
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem => NULL()
    INTEGER                               :: DIMID_NMODE
    INTEGER                               :: REPR_GMXE_4D_NMOD
    INTEGER                               :: REPR_GMXE_1D
    CHARACTER(LEN=1)                      :: ichar
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CHARACTER(LEN=*), PARAMETER      :: channel_name=modstr//'_gp'
    CHARACTER(LEN=*), PARAMETER      :: aerchem_channel_name=modstr//'_aerchem_gp'

    INTEGER                          :: kprod
    INTEGER                          :: dummy
    CHARACTER(LEN=26), POINTER       :: strname(:) => NULL()

#ifdef MESSYTENDENCY
    CALL mtend_register (my_handle_td,mtend_id_t)
    CALL mtend_register (my_handle_td,mtend_id_q)   
    CALL mtend_register (my_handle_td,mtend_id_xl)
    CALL mtend_register (my_handle_td,mtend_id_xi)
    CALL mtend_register (my_handle_td,mtend_id_tracer)
    !
    CALL mtend_register (my_handle_em,mtend_id_tracer)
#endif

    IF(p_parallel_io) &
    CALL start_message_bi(modstr,'Channel definition',substr)
    
    ! POINTER
    ALLOCATE(anum(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(rwet(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(rdry(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    ALLOCATE(ddry(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
!    ALLOCATE(cph(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
!    ALLOCATE(ccn(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))

    !--- 1) Construct the gmxe channel: --------------------------------------!
    CALL new_dimension(status, DIMID_NMODE, 'GMXE_NMODE', nmod)
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATIONS
    CALL new_representation(status, REPR_GMXE_4D_NMOD, &
         'REPR_GMXE_4D_NMOD'    &
         , rank = 4, link = 'xxxx', dctype = DC_GP               &
         , dimension_ids = (/ &
            _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NMODE) /) &
         , ldimlen       = (/ &
            _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /) &
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
    
    cnt(:,_IN_XYZN_)  = nmod
    memu(:,_IN_XYZN_) = nmod
    
    CALL set_representation_decomp(status, REPR_GMXE_4D_NMOD &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL new_representation(status, REPR_GMXE_1D,   &
         'REPR_GMXE_1D'                               &
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
    
    CALL set_representation_decomp(status, REPR_GMXE_1D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)
    ! mz_pj_20061112-

    CALL new_channel(status, channel_name, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)


    !------------------------------------------------------------- 
    ! GMXe DIAGNOSTIC OUTPUT
    !------------------------------------------------------------- 
 
    IF (p_parallel_io) &
    WRITE(*,*) 'add new channel ', TRIM(channel_name),' ...'

    IF (p_parallel_io) &
    WRITE(*,*) 'add channel objects temp, press, rh, sh, zgh2o, zclh2o, zcih2o, zaclc ...'
    CALL new_channel_object(status, channel_name, 'temp' &
         , p3 = t_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'temp' &
         , 'long_name', c='air temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'temp' &
         , 'units', c='deg. C' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'press' &
         , p3 = p_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'press' &
         , 'long_name', c='air pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'press' &
         , 'units', c='hPa'    )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'airdens' &
         , p3 = arho_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'airdens' &
         , 'long_name', c='air density')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'airdens' &
         , 'units', c='kg m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'rh' &
         , p3 = rh_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'rh' &
         , 'long_name', c='relative humdity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'rh' &
         , 'units', c='%'      )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'sh' &
         , p3 = sh_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'sh' &
         , 'long_name', c='specific humdity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'sh' &
         , 'units', c='kg kg-1 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zsat' &
         , p3 = zsat_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zsat' &
         , 'long_name', c='saturation water mass')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zsat' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zgh2o' &
         , p3 = zgh2o_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zgh2o' &
         , 'long_name', c='water vapor')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zgh2o' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zclh2o' &
         , p3 = zclh2o_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zclh2o' &
         , 'long_name', c='cloud water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zclh2o' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zcih2o' &
         , p3 = zcih2o_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcih2o' &
         , 'long_name', c='cloud ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcih2o' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zaclc' &
         , p3 = zaclc_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclc' &
         , 'long_name', c='cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclc' &
         , 'units', c='%')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zaclcac' &
         , p3 = zaclcac_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclcac' &
         , 'long_name', c='accum. cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclcac' &
         , 'units', c='%')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zaclcov' &
         , p2 = zaclcov_2D, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclcov' &
         , 'long_name', c='total cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaclcov' &
         , 'units', c='%')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, channel_name, 'zaopt' &
         , p3 = zaopt_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zaopt' &
     , 'long_name', c='fog type')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zcdnc' &
         , p3 = zcdnc_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcdnc' &
         , 'long_name', c='cloud droplet number conc.')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcdnc' &
         , 'units', c='1 cm-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zicnc' &
         , p3 = zicnc_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zicnc' &
         , 'long_name', c='ice crystal number conc.')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zicnc' &
         , 'units', c='1 cm-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zcrain' &
         , p3 = zcrain_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcrain' &
         , 'long_name', c='cloud rain water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zcrain' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zprec' &
         , p2 = zprec_2D, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zprec' &
         , 'long_name', c='total precipitation (accum.)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zprec' &
         , 'units', c='kg m-2(air) s-1')
    CALL channel_halt(substr, status)

!!$    CALL new_channel_object(status, channel_name, 'e5aclc' &
!!$         , p3 = e5aclc_3D)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclc' &
!!$         , 'long_name', c='echam5 cloud cover')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclc' &
!!$         , 'units', c='%')
!!$    CALL channel_halt(substr, status)
!!$    
!!$    CALL new_channel_object(status, channel_name, 'e5aclcac' &
!!$         , p3 = e5aclcac_3D)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclcac' &
!!$         , 'long_name', c='echam5 accum. cloud cover')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclcac' &
!!$         , 'units', c='%')
!!$    CALL channel_halt(substr, status)
!!$    
!!$    CALL new_channel_object(status, channel_name, 'e5aclcov' &
!!$         , p2 = e5aclcov_2D, reprid=GP_2D_HORIZONTAL)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclcov' &
!!$         , 'long_name', c='echam5 total cloud cover')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5aclcov' &
!!$         , 'units', c='%')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, channel_name, 'e5clh2o' &
!!$         , p3 = e5clh2o_3D)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5clh2o' &
!!$         , 'long_name', c='echam5 cloud water')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5clh2o' &
!!$         , 'units', c='g m-3 (air)')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, channel_name, 'e5cih2o' &
!!$         , p3 = e5cih2o_3D)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5cih2o' &
!!$         , 'long_name', c='echam5 cloud ice')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5cih2o' &
!!$         , 'units', c='g m-3 (air)')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, channel_name, 'e5prec' &
!!$         , p2 = e5prec_2D, reprid=GP_2D_HORIZONTAL)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5prec' &
!!$         , 'long_name', c='echam5 total rain+snow')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5prec' &
!!$         , 'units', c='kg m-2(air) s-1')
!!$    CALL channel_halt(substr, status)
!!$
!!$    CALL new_channel_object(status, channel_name, 'e5cdnc' &
!!$         , p3 = e5cdnc_3D)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5cdnc' &
!!$         , 'long_name', c='echam5 cloud droplet number conc.')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, channel_name, 'e5cdnc' &
!!$         , 'units', c='1 cm-3')
!!$    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'val_strat' &
         , p3 = val_strat_3D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'val_strat' &
         , 'long_name', c='vertical flag (0) for PSC/stratosphere region')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'val_strat' &
         , 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zvap' &
         , p3 = zvap)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zvap' &
         , 'long_name', c='water vapour content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zvap' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zvap2' &
         , p3 = zvap2)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zvap2' &
         , 'long_name', c='water vapour content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zvap2' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zprod' &
         , p3 = zprod)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zprod' &
         , 'long_name', c='water vapour content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zprod' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, channel_name, 'zminc' &
         , p3 = zminc)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zminc' &
         , 'long_name', c='water vapour content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, channel_name, 'zminc' &
         , 'units', c='g m-3 (air)')
    CALL channel_halt(substr, status)

    ! aerosol composition diagnostics (OUTPUT)

    IF(nlowermode < 1 .OR. nlowermode > nmod) THEN
       IF (p_parallel_io) &
       WRITE(*,*) 'Caution: nlowermode = ',nlowermode
       nlowermode = 1 
       IF (p_parallel_io) &
       WRITE(*,*) 'Now    : nlowermode = ',nlowermode
    END IF
    IF(nuppermode < 1 .OR. nuppermode > nmod) THEN
       IF (p_parallel_io) &
       WRITE(*,*) 'Caution: nuppermode = ',nuppermode
       nuppermode = 3
       IF (p_parallel_io) &
       WRITE(*,*) 'Now    : nuppermode = ',nuppermode
    END IF

    IF(nlowermode < 1 .OR. nlowermode > nmod) THEN
       IF (p_parallel_io) &
       WRITE(*,*) 'nlowermode = ',nlowermode
       CALL error_bi('error: nlowermode out of range', substr)
    END IF
    IF(nuppermode < 1 .OR. nuppermode > nmod) THEN
       IF (p_parallel_io) &
       WRITE(*,*) 'nuppermode = ',nuppermode
       CALL error_bi('error: nuppermode out of range', substr)
    END IF

    ALLOCATE(diagaer(0:naerodiag,nlowermode:nuppermode))
    call init_paerosol
    IF(p_parallel_io) &
    WRITE(*,*) 'add channel objects DIAGAER_M',nlowermode, '-',nuppermode
    DO jm = nlowermode, nuppermode
      CALL int2str(char1,jm)
      DO jc=1,naerodiag
        CALL int2str(char2, jc)
        WRITE(name,'(A9,2A1,A2)') 'DIAGAER_M',char1,'_',char2
        If(p_parallel_io )                             &
          WRITE(*,*) TRIM(name), ' ', TRIM(cdiagaer(jc,1)),'  ', &
          TRIM(cdiagaer(jc,2))
        CALL new_channel_object(status, channel_name, TRIM(name)  &
          , p3 = diagaer(jc,jm)%ptr)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, TRIM(name)  &
          , 'long_name', c= cdiagaer(jc,1))
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, TRIM(name)  &
          , 'units', c= cdiagaer(jc,2))
        CALL channel_halt(substr, status)
        
        diagaer(jc,jm)%ptr(:,:,:) = zero ! initialization
        If(p_parallel_io) &
          WRITE(*,*) ' ... ', TRIM(name)
      END DO ! naerodiag
    END DO ! jm

    IF (l_aerchem) THEN
      CALL new_channel(status, aerchem_channel_name, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      ALLOCATE(DIAG_AERCHEM(lmode:umode))
      DO jm=lmode,umode
        CALL int2str(char1,jm)
        CALL new_channel_object(status, aerchem_channel_name, 'LWC_M'//char1  &
          , p3 = diag_aerchem(jm)%lwc)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, aerchem_channel_name, 'LWC_M'//char1  &
          , 'long_name', c='aerosol liquid water content for mode '//char1 )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, aerchem_channel_name, 'LWC_M'//char1  &
          , 'units', c='kg/m^3')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, aerchem_channel_name, 'pH_M'//char1  &
          , p3 = diag_aerchem(jm)%ph)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, aerchem_channel_name, 'pH_M'//char1  &
          , 'long_name', c='aerosol pH value for mode '//char1 )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, aerchem_channel_name, 'pH_M'//char1  &
          , 'units', c='-')
        CALL channel_halt(substr, status)
        
        
        IF (nprod > 0) THEN
          ALLOCATE(DIAG_AERCHEM(jm)%PROD_IDX(nprod))
          ALLOCATE(DIAG_AERCHEM(jm)%PROD_CHAR(nprod))
          ALLOCATE(DIAG_AERCHEM(jm)%PROD(nprod))
          DIAG_AERCHEM(jm)%PROD_CHAR(:) = ""
          DIAG_AERCHEM(jm)%PROD_IDX(:)  = 0

          kprod = 0
          DO jc=1,lspec
             IF (ASSOCIATED(strname)) DEALLOCATE (strname)
             NULLIFY(strname)
             !        print*, "in specloop ", jt, str_field_kpp_l(jt)
             IF (MATCH_WILD( '*Prod*', str_field_kpp_l(jc) )) THEN
                kprod = kprod + 1
                CALL strcrack(str_field_kpp_l(jc), '_', strname, dummy)
                !          print*, "correctly found a prod rate ", kprod, strname, strname(3), dummy
                DIAG_AERCHEM(jm)%PROD_CHAR(kprod) = strname(2)
                DIAG_AERCHEM(jm)%PROD_IDX(kprod) = jc
             END IF
          END DO

          DO jc=1,nprod
            char2 = DIAG_AERCHEM(jm)%PROD_CHAR(jc)

            CALL new_channel_object(status, aerchem_channel_name, 'PROD_R'//char2//'_M'//char1  &
              , p3 = diag_aerchem(jm)%prod(jc)%ptr_3d)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, aerchem_channel_name, 'PROD_R'//char2//'_M'//char1  &
              , 'long_name', c='reaction rate from rate '//char2 )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, aerchem_channel_name, 'PROD_R'//char2//'_M'//char1 &
              , 'units', c='molecules/(cm^3 s)')
            CALL channel_halt(substr, status)
            
          END DO
        END IF
      END DO
    END IF
         
    ! aerosol number

        ! define channel element
        mem => anum(:,:,:,:,1)

        CALL new_channel_object(status, channel_name, 'anumber' &
             , p4 = anum_str, reprid=REPR_GMXE_4D_NMOD, mem=mem)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'anumber' &
             , 'long_name', c='aerosol number' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'anumber' &
             , 'units', c='1/cm3' )
        CALL channel_halt(substr, status)

        DO jm=1,nmod
           CALL int2str(char1, jm)
           ! aerosol number (OUTPUT)
           WRITE(name,'(A9,A1)') 'AERNUMB_M',char1
           IF (p_parallel_io) &
           WRITE(*,*) ' ... ',TRIM(name)
           mem => anum(_RI_XYZN_(:,:,:,jm),:)
           CALL new_channel_object(status, channel_name, TRIM(name), mem = mem)
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name) &
                , 'long_name', c='aerosol number' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name)  &
                , 'units', c='1/cm3' )
           CALL channel_halt(substr, status)
           ! initialization ...
           IF(lnumber) &
           anum(_RI_XYZN_(:,:,:,jm),:) = 0._dp ! [1/cm3]
        END DO

    ! aerosol wet radius

        ! define channel element
        mem => rwet(:,:,:,:,1)
        CALL new_channel_object(status, channel_name, 'wetradius' &
             , p4 = rwet_str, reprid=REPR_GMXE_4D_NMOD, mem=mem)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'wetradius' &
             , 'long_name', c='aerosol wet radius' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'wetradius' &
             , 'units', c='m' )
        CALL channel_halt(substr, status)

        DO jm=1,nmod
           CALL int2str(char1, jm)
           ! ambient aerosol radius (OUTPUT)
           WRITE(name,'(A9,A1)') 'RWETAER_M',char1
           IF (p_parallel_io) &
           WRITE(*,*) ' ... ',TRIM(name)
           mem => rwet(_RI_XYZN_(:,:,:,jm),:)
           CALL new_channel_object(status, channel_name, TRIM(name), mem = mem)
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name) &
                , 'long_name', c='aerosol wet radius' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name)  &
                , 'units', c='m' )
           CALL channel_halt(substr, status)
           ! initialization ...
           rwet(_RI_XYZN_(:,:,:,jm),:) = 0._dp ! [m]
        END DO

    ! aerosol dry radius

        IF (p_parallel_io) &
        WRITE(*,*) ' add channel element:'
        ! define channel element
        mem => rdry(:,:,:,:,1)
        CALL new_channel_object(status, channel_name, 'dryradius' &
             , p4 = rdry_str, reprid=REPR_GMXE_4D_NMOD, mem=mem)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'dryradius' &
             , 'long_name', c='aerosol dry radius' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'dryradius' &
             , 'units', c='m' )
        CALL channel_halt(substr, status)
        DO jm=1,nmod
           CALL int2str(char1, jm)
           ! dry aerosol radius (OUTPUT)
           WRITE(name,'(A9,A1)') 'RDRYAER_M',char1
           IF (p_parallel_io) &
           WRITE(*,*) ' ... ',TRIM(name)
           mem => rdry(_RI_XYZN_(:,:,:,jm),:)
           CALL new_channel_object(status, channel_name, TRIM(name), mem = mem)
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name) &
                , 'long_name', c='aerosol dry radius' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name)  &
                , 'units', c='m' )
           CALL channel_halt(substr, status)
           ! initialization ...
           rdry(_RI_XYZN_(:,:,:,jm),:) = 0._dp ! [m]
        END DO

    ! aerosol dry density

    ! check for channel
        ! define channel element
        IF (p_parallel_io) &
        WRITE(*,*) ' add channel element:'
        mem => ddry(:,:,:,:,1)
        CALL new_channel_object(status, channel_name, 'densaer' &
             , p4 = ddry_str, reprid=REPR_GMXE_4D_NMOD, mem=mem)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'densaer' &
             , 'long_name', c='aerosol dry density' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'densaer' &
             , 'units', c='kg/m3' )
        CALL channel_halt(substr, status)
        DO jm=1,nmod
           CALL int2str(char1, jm)
           ! ambient aerosol density (OUTPUT)
           WRITE(name,'(A9,A1)') 'DDRYAER_M',char1
           IF (p_parallel_io) &
           WRITE(*,*) ' ... ',TRIM(name)
           mem => ddry(_RI_XYZN_(:,:,:,jm),:)
           CALL new_channel_object(status, channel_name, TRIM(name), mem = mem)
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name) &
                , 'long_name', c='aerosol dry density' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, channel_name, TRIM(name)  &
                , 'units', c='kg/m3' )
           CALL channel_halt(substr, status)
           ! initialization ...
           ddry(_RI_XYZN_(:,:,:,jm),:) = 0._dp ! [kg m-3]
        END DO

    ! standard deviation of the aerosol modes

    ! check for channel
        IF (p_parallel_io) THEN
        WRITE(*,*) ' add channel element:'
        WRITE(*,*) ' ... sigma'
        END IF
        CALL new_channel_object(status, channel_name, 'sigma' &
             , p1 = sigma_str, reprid=REPR_GMXE_1D)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'sigma' &
             , 'long_name', c=' standard deviation of aerosol modes' )
        CALL channel_halt(substr, status)

        ! initialization ...
        sigma_str(1:nmod) = sigma(1:nmod)
        
        CALL new_channel_object(status, channel_name, 'crdiv_mid' &
             , p1 = crdiv_str, reprid=REPR_GMXE_1D)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, channel_name, 'crdiv_mid' &
             , 'long_name', c='fixed modal mean size' )
        CALL channel_halt(substr, status)

        ! initialization ...
        crdiv_str(1:nmod) = crdiv_mid(1:nmod)

    IF(p_parallel_io) &
    CALL end_message_bi(modstr,'channel definition',substr)

END SUBROUTINE gmxe_init_memory

!==============================================================================

  SUBROUTINE gmxe_global_start

    IMPLICIT NONE
    
    RETURN

  END SUBROUTINE gmxe_global_start

!==============================================================================

  SUBROUTINE gmxe_radiation

    ! This subroutine is called from messy_radiation (messy_main_control_e5)
    ! and calls gmxe_driver to optionally update various cloud properties
    ! (CDNC, ICNC, CLWC, CIWC, CLC) that are need before the radiation 
    ! calculation.
    IF(TRIM(ADJUSTL(driver_call)) /= 'radiation') RETURN
    CALL gmxe_driver
  
  END SUBROUTINE gmxe_radiation

!==============================================================================

  SUBROUTINE gmxe_physc

  ! Alternatively the driver can be called from physc to be comparable with the 
  ! M7 setup.

    IF(TRIM(ADJUSTL(driver_call)) /= 'physc') RETURN

    CALL gmxe_driver

  END SUBROUTINE gmxe_physc

!==============================================================================



  SUBROUTINE gmxe_driver


#ifndef MESSYTENDENCY
    USE messy_main_data_bi,       ONLY: qm1, tm1, xlm1, xim1, &
                                        qte_3d, tte_3d, xlte_3d, xite_3d
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif 

    USE messy_main_grid_def_mem_bi, ONLY:       nlev       &
                                       ,nproma, kproma     &
                                       ,nrow => ngpblks    &
                                       ,jrow               

#if defined(ECHAM5) || defined(CESM1)
    USE messy_main_data_bi,      ONLY: vervel_3d
#endif

    USE messy_main_timer,    ONLY: delta_time,time_step_len &
                                  ,nstep=>current_time_step
    USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer_mem_bi,  ONLY: ntrac => ntrac_gp
    USE messy_main_mpi_bi,         ONLY: p_io, p_bcast, &
                                         p_parallel_io, p_pe
    USE messy_gmxe_aerchem_liq,    ONLY: jval_h2o2, jval_h2o2_2d, &
                                         jval_o3,   jval_o3_2d,   &
                                         jval_no3,  jval_no3_2d,  &
                                         diag_aerchem, nprod
!
    IMPLICIT NONE

    !--- Parameter list:
    
    ! I/O
    
    ! Local variables:
    
    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_driver'
    !
    LOGICAL,  SAVE              :: entered = .FALSE.
    LOGICAL                     :: loabort = .FALSE.
    
    INTEGER                     :: it, jl, jk, jc, jt, jm, i, j
    
    REAL(dp)                    :: zmass_pre,  zmass_post
    REAL(dp), DIMENSION(kproma,nlev) :: zmass_pre2, zmass_post2
    
    REAL(dp), DIMENSION(kproma,nlev) :: zrhoa,  zpress, ztemp,  zrhum,  zshum, &
      zaopt,  zgh2o,  zaclc,  zclh2o, zcih2o,&
      zcdnc,  zicnc,  zsat,   zair,          &
      zcrain, zdpress
    REAL(dp)                           :: zconv
    REAL(dp), DIMENSION(kproma,nlev+1) :: zpressi
    
    REAL(dp), DIMENSION(kproma,nlev,0:nmod,0:naerodiag) :: zaerosol
    
    REAL(dp), DIMENSION(kproma,nlev,nmod)       :: zaernl, zcph, zccn
    REAL(dp), DIMENSION(kproma,nlev,nmod)       :: zrdry, zrwet, zddry
    REAL(dp), DIMENSION(kproma,nlev,0:naertot)  :: zaerml
    REAL(dp), DIMENSION(kproma,nlev,0:naertot)  :: zxtm
    
    REAL(dp), DIMENSION(kproma,nlev)      :: pscregion, pscregions
    REAL(dp), DIMENSION(kproma)           :: tp_i, zhelp
    REAL(dp)                              :: zdtime, zepsec, zxsec, ztmst
    
    INTEGER                               :: idx1, idx2, idx3, idt
    REAL(dp)                              :: val_coup
    

! local tendencies
    REAL(dp), DIMENSION(kproma,nlev)      :: ztte, zqte, zxlte, zxite
    REAL(dp), DIMENSION(kproma,nlev,ntrac):: zxtte
! start values
    REAL(dp), DIMENSION(kproma,nlev)       :: ztp1_0, zqp1_0, zxlp1_0, zxip1_0
    REAL(dp), DIMENSION(kproma,nlev,ntrac) :: zxtp1_0


    IF(p_parallel_io .AND. .NOT. entered ) &
      CALL start_message_bi(modstr,' prepare aerosol composition and dynamics calculation ',substr)
    
    !--- Initialisations: -----------------------------------------------

    l_io=.FALSE.
    IF(p_parallel_io .AND. loutput .AND. .NOT. entered) l_io=.TRUE.
    
    zconv             = zero
    zair        (:,:) = zero
    ztemp       (:,:) = zero
    zpress      (:,:) = zero
    zpressi     (:,:) = zero
    zdpress     (:,:) = zero
    zrhoa       (:,:) = zero
    zrhum       (:,:) = zero
    zshum       (:,:) = zero
    zsat        (:,:) = zero
    zgh2o       (:,:) = zero
    zclh2o      (:,:) = zero
    zcih2o      (:,:) = zero
    zaclc       (:,:) = zero
    zaopt       (:,:) = zero
    zcrain      (:,:) = zero
    zcdnc       (:,:) = zero
    zicnc       (:,:) = zero
    zcph      (:,:,:) = NaN
    zccn      (:,:,:) = zero
    zaerml    (:,:,:) = zero
    zaernl    (:,:,:) = zero
    zxtm      (:,:,:) = zero
    zaerosol(:,:,:,:) = ZERO

    zqte(:,:)    = 0._dp
    zxlte(:,:)   = 0._dp
    zxite(:,:)   = 0._dp
    ztte(:,:)    = 0._dp
    zxtte(:,:,:) = 0._dp

    zrdry = 0._dp
    zrwet = 0._dp
    zddry = 0._dp
    !--- Security parameters
    
    zepsec = FACm12
    zxsec  = FACp0-zepsec
    
    !--- Computational constants
    
    zdtime = delta_time
    ztmst  = time_step_len
    
    !--- GRID SPACE
    
    zmass_pre         = zero
    zmass_post        = zero
    zmass_pre2(:,:)   = zero
    zmass_post2(:,:)  = zero
    
    !--- Ambient properties: ----------------------------------
    
    !--- Temperature [K]:
#ifndef MESSYTENDENCY    
    ztp1_0(1:kproma,1:nlev)       = tm1  (_RI_XYZ__(1:kproma,jrow,1:nlev))   +       &
                                     tte_3d (_RI_XYZ__(1:kproma,jrow,1:nlev)) * time_step_len
    !--- Specific humidity [kg kg-1 (air)]:

    zqp1_0(1:kproma,1:nlev)   = &
         MAX(qm1(_RI_XYZ__(1:kproma,jrow,1:nlev)) + &
         qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * time_step_len, zero)
    ! liquid water [kg kg-1 (air)]:
    zxlp1_0(1:kproma,1:nlev)   =&
         MAX(xlm1(_RI_XYZ__(1:kproma,jrow,1:nlev)) + &
         xlte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * time_step_len, zero)
    ! ice water [kg kg-1 (air)]:
    zxip1_0(1:kproma,1:nlev)   = &
         MAX(xim1(_RI_XYZ__(1:kproma,jrow,1:nlev)) + &
         xite_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) * time_step_len, zero)
    ! tracer mixing ratios [mol mol-1]
    !um_mb_20171103+
    DO jk=1,nlev
        DO jt=1,ntrac
            zxtp1_0(1:kproma,jk,jt) = pxtm1(_RI_X_ZN_(1:kproma,jk,jt)) &
                                     + pxtte(_RI_X_ZN_(1:kproma,jk,jt)) * time_step_len
        ENDDO
    ENDDO
    !um_mb_20171103-

#else
    call mtend_get_start_l (mtend_id_t,  v0 = ztp1_0)
    call mtend_get_start_l (mtend_id_q,  v0 = zqp1_0)
    call mtend_get_start_l (mtend_id_xl, v0 = zxlp1_0)
    call mtend_get_start_l (mtend_id_xi, v0 = zxip1_0)
    ! ub_ak_20190613+
    !call mtend_get_start_l (mtend_id_tracer, v0t = zxtp1_0)
    DO jt = 1, ntrac
       call mtend_get_start_l (jt, v0 = zxtp1_0(:,:,jt))
    END DO
    ! ub_ak_20190613-
#endif
    zshum(1:kproma,1:nlev) = zqp1_0(1:kproma,1:nlev)
    ztemp(1:kproma,1:nlev) = ztp1_0(1:kproma,1:nlev)
    

    DO jm = 1, nmod
      zaerosol(1:kproma,1:nlev,jm,iTT) = ztemp(1:kproma,1:nlev)
    END DO
    
    !--- Pressure [Pa]:
    
    zpress  (1:kproma,1:nlev)      = press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))     ! [Pa]
    zpressi (1:kproma,1:nlev+1)    = pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev+1))  ! [Pa]

    DO jk=1,nlev
      zdpress(1:kproma,jk)         = pressi_3d(_RI_XYZ__(1:kproma,jrow,jk+1)) - &
                                     pressi_3d(_RI_XYZ__(1:kproma,jrow,jk)) 
    END DO


    !--- Air density [kg (air) m-3 (air)]:

    zrhoa  (1:kproma,1:nlev)  = grmass(_RI_XYZ__(1:kproma,jrow,1:nlev)) / &
                                grvol(_RI_XYZ__(1:kproma,jrow,1:nlev))
    zair   (1:kproma,1:nlev)  = zrhoa (1:kproma,1:nlev)         * (avo*FACm6)

    !--- Relative humidity [0-1]:

    zrhum  (1:kproma,1:nlev)             = rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))
    DO jm = 1, nmod
      zaerosol(1:kproma,1:nlev,jm,iRH)     = zrhum  (1:kproma,1:nlev)
    END DO
    zrhum  (1:kproma,1:nlev)             = zrhum  (1:kproma,1:nlev)        * FACm2
    

    !--- Saturation water mass [g m-3 (air)]:

    zsat (1:kproma,1:nlev)         = zsat_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))
    
    !--- Water vapor [kg kg-1 (air)] ->  [g m-3 (air)]:
    !    (assuming that model specific humidity represents only water vapor and no aqueous phase)
    zgh2o(1:kproma,1:nlev)         = MAX(0._dp, zshum(1:kproma,1:nlev) &
                                   * zrhoa(1:kproma,1:nlev)          * FACp3)
    
    
    !--- Cloud liquid water [g m-3 (air)]:
    zclh2o(1:kproma,1:nlev)         = MAX(zclh2o_3D (_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)
    
    !--- Cloud ice water [g m-3 (air)]:
    zcih2o(1:kproma,1:nlev)        = MAX(zcih2o_3D (_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)
    
    !--- Cloud rain water [g/m3 (air)] <- [kg/m**2s]:   
    zcrain(1:kproma,1:nlev)        = MAX(zcrain_3D(_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)


    !--- Cloud cover [%]:
    zaclc(1:kproma,1:nlev)         = MAX(zaclc_3D(_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)

    !--- Fog type; 0=no fog, 0-50=brown fog, 50-100=fog [%]:
    zaopt(1:kproma,1:nlev)         = MAX(zaopt_3D(_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)

    !--- Cloud droplet number concentration [1 cm-3 (air)]:
    zcdnc(1:kproma,1:nlev)         = MAX(zcdnc_3D(_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)

    !--- Cloud ice crystal number concentration [1 cm-3 (air)]:
    zicnc(1:kproma,1:nlev)         = MAX(zicnc_3D(_RI_XYZ__(1:kproma,jrow,1:nlev)), zero)

    ! Tropopause-mask
    tp_i(1:kproma) = ZERO
    IF (LSTRAT .AND. ASSOCIATED(tp_i0)) tp_i(1:kproma) = tp_i0(1:kproma,jrow)
    
    val_strat_3D(_RI_XYZ__(:,jrow,:)) = FACp0
    ! Don't update tendencies for stratosphere, if val_strat_3D = zero
    DO jk=1,nlev
      DO jl=1,kproma
        IF(jk < NINT(tp_i(jl))) val_strat_3D(_RI_XYZ__(jl,jrow,jk)) = zero
      END DO
    END DO
    
    ! PSC-mask  modified cb-01-2010
    pscregion(1:kproma,:)  = ZERO
    pscregions(1:kproma,:) = ZERO
    IF (LPSC .AND. ASSOCIATED(flt_pscreg0)) pscregion(1:kproma,:) =  &
      flt_pscreg0(_RI_XYZ__(1:kproma,jrow,:))
    IF (LPSC .AND. ASSOCIATED(flt_pscreg1)) pscregions(1:kproma,:) = &
      flt_pscreg1(_RI_XYZ__(1:kproma,jrow,:))
    ! Don't update tendencies for PSC region, if val_strat_3D = zero
    DO jk=1,nlev
      DO jl=1,kproma
        IF(pscregions(jl,jk) > 1.01_dp) val_strat_3D(_RI_XYZ__(jl,jrow,jk)) = zero
      END DO
    END DO

    !--- Assign and convert tracers:
    
    ! Trace gases + aerosols [mol mol-1] -> [umol m-3 (air)]

    ! Assign tracer values and conversion factors (if neccessary)

    ! molecules/cm^3  ----> micro-mol/m^3 
    zconv    = 1._dp / (avo * facm6 * facm6)
  
    DO jt=1,spec_number
      DO jm=0,nmod
        idx1 = species(jt)%tracidx(jm)
        idx2 = species(jt)%aermlidx(jm)
        IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
        zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2)
        IF(lrhoa(idx2)) &
          zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2) &
                                         * zrhoa(1:kproma,1:nlev)
        idt = idx1

        zxtm(1:kproma,1:nlev,idx2) = MAX(0._dp, zxtp1_0(1:kproma,1:nlev,idt))

        zaerml(1:kproma,1:nlev,idx2) = zxtm(1:kproma,1:nlev,idx2)       &  
                                     ! tracer mass       [mol  mol-1]
                                     * zconvert(1:kproma,1:nlev,idx2)   & 
                                     ! conversion factor [umol m-3 (air)]
                                     * val_strat_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))
                                     ! stratosphere/PSC mask

        IF(.NOT. lnumber) THEN
          idx2 = species(jt)%aernlidx(jm)
          IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
          zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2)
          IF(lrhoa(idx2)) &
            zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2) &
                                           * zrhoa(1:kproma,1:nlev)
          idt = idx1

          zxtm(1:kproma,1:nlev,idx2) = MAX(0._dp, zxtp1_0(1:kproma,1:nlev,idt))

          zaernl(1:kproma,1:nlev,jm) = zxtm(1:kproma,1:nlev,idx2)        &
                                     ! tracer number       [#  mol-1]
                                     * zconvert(1:kproma,1:nlev,idx2)    &
                                     ! conversion factor [umol m-3 (air)]
                                     * val_strat_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))
                                     ! stratosphere/PSC mask
          zaerml(1:kproma,1:nlev,idx2) = 0._dp
        END IF
      END DO
    END DO


    !--- Assign water vapor [umol m-3 (air)]:
    !    [umol m-3 (air)]             =      [g m-3 (air)]     / [g mol-1]  -> [g umol-1]
    zaerml(1:kproma,1:nlev,nwh2o(0)) = zgh2o(1:kproma,1:nlev) / (mwh2o * FACm6)

    !--- Assign aerosol quantities from channels:
    IF(nstep > 0) THEN
      !--- Aerosol Number [1 cm-3]:
      IF (lnumber) THEN
        DO jm=1,nmod
          zaernl  (1:kproma,1:nlev,jm)        = anum_str(_RI_XYZN_(1:kproma,jrow,1:nlev,jm))
        END DO
      END IF
    END IF

    !--- Sum total mass of all compounds for mass diagnostics:
    IF (l_aerchem) THEN
      jval_h2o2_2d => jval_h2o2(_RI_XYZ__(:,jrow,:))
      jval_no3_2d  => jval_no3(_RI_XYZ__(:,jrow,:))
      jval_o3_2d   => jval_o3(_RI_XYZ__(:,jrow,:))

      DO jm=lmode,umode
        diag_aerchem(jm)%lwc_2d => diag_aerchem(jm)%lwc(_RI_XYZ__(:,jrow,:))
        diag_aerchem(jm)%ph_2d  => diag_aerchem(jm)%ph(_RI_XYZ__(:,jrow,:))
        diag_aerchem(jm)%lwc_2d(:,:) = 0._dp
        diag_aerchem(jm)%pH_2d(:,:)  = 0._dp   
        DO jt=1,nprod
          diag_aerchem(jm)%prod(jt)%ptr_2d => diag_aerchem(jm)%prod(jt)%ptr_3d(_RI_XYZ__(:,jrow,:))
          diag_aerchem(jm)%prod(jt)%ptr_2d = 0._dp
        END DO
      END DO
    END IF
    
    IF(lmass_diag) CALL sum_mass(zmass_pre)
    IF(lmass_diag) CALL sum_mass2(zmass_pre2)

    vervel => vervel_3d(_RI_XYZ__(:,jrow,:))

    !--- Call GMXe main -----------------------------------------------------

!CDIR NOIEXPAND
    IF(lgmxe .AND. nstep > 0) &  ! skip first time step to avoid numerical garbage
      CALL gmxe_main(kproma,         nlev,           time_step_len,               &
     ! ECHAM indices
                 zrhum,          ztemp,          zpress,          zsat,           &
     !   "   thermodynamics
                 zaclc,          zaopt,          zcih2o,          zclh2o,         &
     !   "   aerosol-cloud properties
                 zcdnc,          zicnc,          zcph,            zcrain,         &
     !   "   aerosol-cloud properties
                 zccn,           zaerml,         zaernl,          zaerosol,       &
     !   "   aerosol-cloud properties
                 zconv,     zvap(_RI_XYZ__(:,jrow,:)), zvap2(_RI_XYZ__(:,jrow,:)), zprod(_RI_XYZ__(:,jrow,:)),  &
     !   "   aerosol-cloud properties
                 zminc(_RI_XYZ__(:,jrow,:)),zrhoa*FACp3,   zpressi,         vervel,         &
                 zrdry,           zrwet,         zddry )


    !--- Reconvert masses and numbers and assign tracer tendencies:

    !--- Assign water vapor [g m-3 (air)]:
    !    [g m-3 (air)]     =           [umol m-3 (air)]        * [g mol-1]  -> [g umol-1]
    zgh2o(1:kproma,1:nlev) = zaerml(1:kproma,1:nlev,nwh2o(0)) * (mwh2o * FACm6)
    IF (p_parallel_io .AND. .NOT. entered ) &
      CALL end_message_bi(modstr,' prepare aerosol composition and dynamics calculation ',substr)
    entered = .TRUE.
    ! Trace gases + aerosols [umol m-3 (air)] -> [mol mol-1]
    
    DO jt=1,spec_number
      DO jm=0,nmod
        idx1 = species(jt)%tracidx(jm)
        idx2 = species(jt)%aermlidx(jm)
        val_coup = species(jt)%zcoup(jm)
        IF (IDX1 == IDT_H2O) CYCLE
        IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
        idt = idx1
        IF (idx2 /= species(jt)%aernlidx(jm) ) &
             zxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =  (zaerml(1:kproma,1:nlev,idx2)           &
                                      /  zconvert(1:kproma,1:nlev,idx2)         &
                                      -  zxtm(1:kproma,1:nlev,idx2) )           &
                                      /  time_step_len                          &
                                      *  val_strat_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))        &
                                      *  val_coup

        IF(.NOT. lnumber) THEN
          idx2 = species(jt)%aernlidx(jm)
          IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
          idt = idx1
          zxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))  =  (zaernl(1:kproma,1:nlev,jm)       &
                                    /   zconvert(1:kproma,1:nlev,idx2)   &
                                    -   zxtm(1:kproma,1:nlev,idx2) )     &
                                    /   time_step_len                    &
                                    *   val_strat_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))

        END IF
      END DO
    END DO

    !--- Update specific humidity [kg kg-1 (air)]:

    IF (lshum) THEN
      !   [kg kg-1 (air)]            =      [g m-3 (air)]      /  [kg m-3 (air)]  /   [kg (air) m-3 (air)]
      zshum(1:kproma,1:nlev) = zgh2o(1:kproma,1:nlev) * FACm3 / &
                               zrhoa(1:kproma,1:nlev)

      zqte(1:kproma,1:nlev) = (zshum(1:kproma,1:nlev) - zqp1_0(1:kproma,1:nlev)) &
                            / time_step_len 

    END IF
  
  !--- Cloud liquid water [kg kg-1 (air)]:

    IF (lclwc) THEN
      !   [kg kg-1 (air)]            =       [g m-3 (air)]      /  [kg m-3 (air)] /   [kg (air) m-3 (air)]
       
       zxlte(1:kproma,1:nlev) = (zclh2o(1:kproma,1:nlev) * FACm3            &
                              /  zrhoa(1:kproma,1:nlev) -                   &
                                 zxlp1_0(1:kproma,1:nlev))                  &
                              /  time_step_len

    END IF

    !--- Cloud ice water [kg kg-1 (air)]:

    IF (lciwc) THEN
      !   [kg kg-1 (air)]            =      [g m-3 (air)]      /  [kg m-3 (air)]  /   [kg (air) m-3 (air)]
       zxite(1:kproma,1:nlev) = (zcih2o(1:kproma,1:nlev) * FACm3            &
                              /  zrhoa(1:kproma,1:nlev) -                   &
                                 zxip1_0(1:kproma,1:nlev))                  &
                              /  time_step_len

    END IF

#ifndef MESSYTENDENCY
    qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))  = &
         qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))  + zqte(1:kproma,1:nlev)
    xlte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
         xlte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) + zxlte(1:kproma,1:nlev)
    xite_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) = &
         xite_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) + zxite(1:kproma,1:nlev)
    DO jk=1,nlev
        DO jt=1,ntrac
            pxtte(_RI_X_ZN_(1:kproma,jk,jt)) = pxtte(_RI_X_ZN_(1:kproma,jk,jt)) + &
                                zxtte(_RI_X_ZN_(1:kproma,jk,jt))  !!um_mb_20171110 
        ENDDO
    ENDDO
#else
    CALL mtend_add_l (my_handle_td, mtend_id_q, px = zqte)
    CALL mtend_add_l (my_handle_td, mtend_id_xl, px = zxlte)
    CALL mtend_add_l (my_handle_td, mtend_id_xi, px = zxite)
    ! ub_ak_20190613+
    !CALL mtend_add_l (my_handle_td, mtend_id_tracer, pxt = zxtte)
    DO jt = 1, ntrac
       CALL mtend_add_l (my_handle_td, jt, px = zxtte(:,:,jt))
    END DO
    ! ub_ak_20190613-
#endif

    !--- Assign aerosol quantities to channels:

    !--- Aerosol diagnostics:

    DO jm=nlowermode, nuppermode
      DO jc = 1, naerodiag
        diagaer(jc,jm)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) = zaerosol(1:kproma,1:nlev,jm,jc)
      END DO
    END DO

#ifdef _INFO
    IF(loutput .AND. p_parallel_io .AND. jrow == nrow/2) THEN
      DO jm=1, nmod
        WRITE(*,*) ' GMXe: Max val before update - after update - aerosol mode = ',  cmodes(jm)
        IF(jm >= nlowermode .AND. jm <= nuppermode) THEN
          
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'T  [K]        = ',MAXVAL(diagaer (iTT, jm)%ptr(_RI_XYZ__(:,jrow,:))),  &
            ' ',MAXVAL(zaerosol(:,:,jm,iTT        ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'RH [0-1]      = ',MAXVAL(diagaer (iRH, jm)%ptr(_RI_XYZ__(:,jrow,:))),  &
            ' ',MAXVAL(zaerosol(:,:,jm,iRH        ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'sPM[umol m-3] = ',MAXVAL(diagaer (isPM, jm)%ptr(_RI_XYZ__(:,jrow,:))), &
            ' ',MAXVAL(zaerosol(:,:,jm,isPM       ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'aPM[umol m-3] = ',MAXVAL(diagaer (iaPM, jm)%ptr(_RI_XYZ__(:,jrow,:))), &
            ' ',MAXVAL(zaerosol(:,:,jm,iaPM       ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'PMs  [ug m-3] = ',MAXVAL(diagaer (iPMs, jm)%ptr(_RI_XYZ__(:,jrow,:))), &
            ' ',MAXVAL(zaerosol(:,:,jm,iPMs       ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'PMt  [ug m-3] = ',MAXVAL(diagaer (iPMt, jm)%ptr(_RI_XYZ__(:,jrow,:))), &
            ' ',MAXVAL(zaerosol(:,:,jm,iPMt       ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'WH2O [ug m-3] = ',MAXVAL(diagaer (iWH2O, jm)%ptr(_RI_XYZ__(:,jrow,:))),&
            ' ',MAXVAL(zaerosol(:,:,jm,iWH2O      ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'GF   [------] = ',MAXVAL(diagaer (iGF, jm)%ptr(_RI_XYZ__(:,jrow,:))),  &
            ' ',MAXVAL(zaerosol(:,:,jm,iGF        ))      ,jm
          WRITE(*,'(A,E12.3,A,E12.3,I4)') 'pH   [------] = ',MAXVAL(diagaer (iPH, jm)%ptr(_RI_XYZ__(:,jrow,:))),  &
            ' ',MAXVAL(zaerosol(:,:,jm,iPH        ))      ,jm
        END IF

        WRITE(*,'(A,E12.3,A,E12.3,I4)') 'anum [N cm-3] = ',MAXVAL(anum    (_RI_XYZN_(:,jrow,1:nlev,jm),   1)),       &
          ' ',MAXVAL(zaernl  (:,:,jm          ))       ,jm
        WRITE(*,'(A,E12.3,A,E12.3,I4)') 'ddry [g cm-3] = ',MAXVAL(zddry    (:,:,jm))*FACm3,jm
        WRITE(*,'(A,E12.3,A,E12.3,I4)') 'rwet     [um] = ',MAXVAL(zrwet    (:,:,jm))*FACp6,jm
      END DO
    END IF
#endif
 
    DO jm=1,nmod
      !--- Aerosol Number [1 cm-3]:
      anum    (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) = zaernl  (1:kproma,1:nlev,jm)
      !--- Dry Count Mean Radius [m]:
      rdry    (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) = zrdry(1:kproma,1:nlev,jm) * FACm2
      !--- Wet Count Mean Radius [m]:
      rwet    (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) = zrwet(1:kproma,1:nlev,jm) * Facm2
      !--- Wet Mean Mode Density [kg/m3]:
      ddry    (_RI_XYZN_(1:kproma,jrow,1:nlev,jm),1) = zddry(1:kproma,1:nlev,jm) * Facp3
      
      !--- Standard deviation for the modes:
      sigma_str(jm) = sigma(jm)
    END DO
      
    !--- Temperature [K]:
    t_3D    (_RI_XYZ__(1:kproma,jrow,1:nlev))         = ztemp   (1:kproma,1:nlev)  - 273.15_dp    ! [deg. C]

    !--- Pressure [hPa]:
    p_3D    (_RI_XYZ__(1:kproma,jrow,1:nlev))         = zpress  (1:kproma,1:nlev)   * FACm2       ! [hPa]

    !--- Air density [kg (air) m-3 (air)]:
    arho_3D (_RI_XYZ__(1:kproma,jrow,1:nlev))         = zrhoa   (1:kproma,1:nlev)                 ! [kg (air) m-3 (air)]

    !--- Relative humidity [%]:
    rh_3D   (_RI_XYZ__(1:kproma,jrow,1:nlev))         = zrhum   (1:kproma,1:nlev) * FACp2         ! [%]

    !--- Specific humidity [kg kg-1 (air)]:
    sh_3D   (_RI_XYZ__(1:kproma,jrow,1:nlev))         = zshum   (1:kproma,1:nlev)                 ! [kg kg-1 (air)]

    !--- Fog type; 0=no fog, 0-50=brown fog, 50-100=fog [%]:
    zaopt_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))         = zaopt   (1:kproma,1:nlev)                 ! [%]

    !--- Saturation water mass [g m-3 (air)]:
    zsat_3d (_RI_XYZ__(1:kproma,jrow,1:nlev))         = zsat    (1:kproma,1:nlev)                 ! [g m-3 (air)]

    !--- Cloud liquid water [g m-3 (air)]:
    zclh2o_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))        = zclh2o(1:kproma,1:nlev)                   ! [g m-3 (air)]

    !--- Cloud ice water water [g m-3 (air)]:
    zcih2o_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))        = zcih2o(1:kproma,1:nlev)                   ! [g m-3 (air)]

    !--- Water vapor [g m-3 (air)]:
    zgh2o_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))         = zgh2o (1:kproma,1:nlev)                   ! [g m-3 (air)]

    !--- Cloud cover [%]:
    IF (lclc) &
     !              [0-1]                =             [%]
      aclc    (_RI_XYZ__(1:kproma,jrow,1:nlev))        =   MAX(zaclc(1:kproma,1:nlev)*FACm2, zero)
    ! Output diagnostic value in any case

    ! Instantanous cloud cover [%]:
    zaclc_3D  (_RI_XYZ__(1:kproma,jrow,1:nlev))      =       zaclc(1:kproma,1:nlev)               ! [%]
    ! Accumulated  cloud cover [%]:
    zaclcac_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))      =  zaclcac_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))    &
                                       +  zaclc(1:kproma,1:nlev)  * zdtime

    ! Total cloud cover [%]:
    DO jl = 1,kproma
      zhelp(jl) = FACp2-zaclc(jl,1)
    END DO
    DO jk = 2,nlev
      DO jl = 1,kproma
        zhelp(jl) = zhelp(jl)*(FACp2-MAX(zaclc(jl,jk),  zaclc(jl,jk-1))) &
          /(FACp2-MIN(zaclc(jl,jk-1),zxsec))
      END DO
    END DO
    DO jl = 1,kproma
      zhelp(jl)  = FACp2-zhelp(jl)
      zaclcov_2D(jl,jrow) = zaclcov_2D(jl,jrow)+zdtime*zhelp(jl)
    END DO

    !--- Cloud rain water [g m-3 (air)]:
    zcrain_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))       = MAX(zcrain(1:kproma,1:nlev) , zero)

    !--- Total precipitation (accumulated) [kg m-2(air) s-1]:

    zhelp(:) = zero
    DO jk = 1,nlev
      DO jl = 1,kproma
      ! g_water/m^3 * 1e-3 * rho[kg_air/m^3] -> kg_water/kg_air
      ! XXX   * deltap / g [kg_air/m^2] -> kg_water/m^2
      ! XXX  / ztmst                    -> kg_water/(m^2*s)
        zhelp(jl) = zhelp(jl) + zcrain(jl,jk) * FACm3 * &
                                zrhoa(jl,jk) * zdpress(jl,jk)/(g*ztmst)
      END DO
    END DO
    zprec_2D(1:kproma,jrow) = zprec_2D(1:kproma,jrow) + zhelp(1:kproma) * zdtime

    !--- Cloud droplet number concentration [1 m-3 (air)]:

    IF (lcdnc) &
     !            [1 m-3 (air)]        =             [1 cm-3 (air)]
      acdnc   (_RI_XYZ__(1:kproma,jrow,1:nlev))      =  MAX(zcdnc(1:kproma,1:nlev), zero) * FACm6
  ! Output diagnostic value in any case
    zcdnc_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))        =      zcdnc(1:kproma,1:nlev)               ! [1 cm-3 (air)]

    !--- Cloud ice crystal number concentration [1 cm-3 (air)]:
!    IF (licnc) &
    !            [1 m-3(air) ]          =            [1 cm-3 (air)]
!      aicnc   (_RI_XYZ__(1:kproma,jrow,1:nlev))         = MAX(zicnc(1:kproma,1:nlev), zero)  * FACm6
    ! Output diagnostic value in any case
    zicnc_3D(_RI_XYZ__(1:kproma,jrow,1:nlev))         =     zicnc(1:kproma,1:nlev)                ! [1 cm-3 (air)]

    !--- Perform mass conservation check: 
    
    IF(lmass_diag) THEN

      !--- Sum total mass of all compounds for mass diagnostics:

      CALL sum_mass(zmass_post)
      CALL sum_mass2(zmass_post2)
      
      !--- Perform mass conservation check:

      IF( ABS(zmass_pre-zmass_post) > FACm1*ABS(MAX(zmass_pre,zmass_post))) THEN

        IF (p_parallel_io) &
          WRITE(*,*) 'GMXe ','violates the mass conservation : Error > 10%'
        
        IF (MAX(zmass_pre,zmass_post) > zero) THEN
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            (zmass_post-zmass_pre)/MAX(zmass_pre,zmass_post)*FACp2,"%"
        ELSE
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            zmass_post-zmass_pre
        END IF
        
        loabort=.TRUE.
        
      ELSE IF( ABS(zmass_pre-zmass_post) > FACm2*ABS(MAX(zmass_pre,zmass_post))) THEN
        
        IF (p_parallel_io) &
          WRITE(*,*) 'GMXe ','violates the mass conservation : Error > 1% < 10%'
        
        IF (MAX(zmass_pre,zmass_post) >zero) THEN
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            (zmass_post-zmass_pre)/MAX(zmass_pre,zmass_post)*FACp2,"%"
        ELSE
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            zmass_post-zmass_pre
        END IF
        
      ELSE IF( ABS(zmass_pre-zmass_post) > FACm3*ABS(MAX(zmass_pre,zmass_post))) THEN
        
        IF (p_parallel_io) &
          WRITE(*,*) 'GMXe ','violates the mass conservation : Error > 0.1% < 1%'
        
        IF (MAX(zmass_pre,zmass_post) >zero) THEN
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            (zmass_post-zmass_pre)/MAX(zmass_pre,zmass_post)*FACp2,"%"
        ELSE
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            zmass_post-zmass_pre
        END IF
        
        loabort=.TRUE.
      ELSE
        
        IF (ABS(zmass_pre-zmass_post) > (100._dp * SPACING(zmass_post)) ) THEN
          !           IF (ABS(zmass_pre-zmass_post) > FACm3) THEN
          IF (p_parallel_io) &
            WRITE(*,*) 'GMXe ','violates the mass conservation : Error < 0.1%'
          IF (p_parallel_io) &
            WRITE(*,*) '              change in mass : ', &
            (zmass_post-zmass_pre)/MAX(zmass_pre,zmass_post)*FACp2,"%"
          IF (p_parallel_io) &
            WRITE(*,*)'post =',zmass_post,' pre=',zmass_pre,' denom=',MAX(zmass_pre,zmass_post)*FACp2
          
        END IF
        
      END IF  !zmass_pre-zmass_post 
      
      IF (loabort)THEN
        IF (p_parallel_io) THEN
          WRITE(*,*) 'GMXe',TRIM(modver),' violates mass conservation (in gmxe_driver)'
        END IF
      END IF
      DO jk=1,nlev
        DO jl=1,kproma
          IF (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) >                &
            (10._dp * SPACING(zmass_post2(jl,jk))) ) THEN
            print*, "mass changes boxwise; box: ",jl,jk," Error: (%)", &
              (ABS(zmass_pre2(jl,jk)-zmass_post2(jl,jk)) /               &
              zmass_post2(jl,jk)*100._dp),                               &
              zmass_pre2(jl,jk), zmass_post2(jl,jk)
          END IF
        END DO
      END DO
    END IF

  

  !------------------------------------------------------------------------------

  !--- Auxiliary routine:

  CONTAINS

  !------------------------------------------------------------------------------

    SUBROUTINE sum_mass(pmasssum)

      IMPLICIT NONE

      INTEGER     :: jl, jk, jm, jc, it
      REAL(dp)    :: pmasssum, zscale

! this one tests the conservation of molecules

!!$      pmasssum = ZERO
!!$      DO jm = 0, nmod
!!$        DO jc = 1, mcomp(jm)
!!$          it = icomp(jc,jm)
!!$          zscale = FACp6
!!$          IF(it == nwh2o(jm)) zscale = FACm4 ! scale down water vapor
!!$          DO jk = 1,nlev
!!$            DO jl = 1,kproma
!!$              pmasssum = pmasssum + zaerml(jl,jk,it) * zscale
!!$            END DO
!!$          END DO
!!$        END DO
!!$      END DO

! this one tests the conservation of masses

      pmasssum = ZERO
      DO jm = 0, nmod
        DO jc = 1, spec_number
          it = species(jc)%aermlidx(jm)
          IF (it == species(jc)%aernlidx(jm)) CYCLE
          zscale = 1._dp
          DO jk = 1,nlev
            DO jl = 1,kproma
              pmasssum = pmasssum + zaerml(jl,jk,it) &
                       * species(jc)%molmass * zscale
              
            END DO
          END DO
        END DO
      END DO
      
    END SUBROUTINE sum_mass
    SUBROUTINE sum_mass2(pmasssum)

      IMPLICIT NONE

      INTEGER     :: jl, jk, jm, jc, it
      REAL(dp)    :: pmasssum(kproma,nlev), zscale

! this one tests the conservation of molecules


      pmasssum = ZERO
      DO jm = 0, nmod
!!$            DO jc = 1, mcomp(jm)
!!$              it = icomp(jc,jm)
!!$              zscale = 1._dp
!!$              pmasssum(jl,jk) = pmasssum(jl,jk) + zaerml(jl,jk,it) * &
!!$                                molarmass(jc,jm) * zscale
        DO jc = 1, spec_number
          it = species(jc)%aermlidx(jm)
          IF (it == species(jc)%aernlidx(jm)) CYCLE
          zscale = 1._dp
          DO jl=1,kproma
            DO jk=1,nlev
              pmasssum(jl,jk) = pmasssum(jl,jk) + zaerml(jl,jk,it) * &
                                species(jc)%molmass * zscale
            END DO
          END DO
        END DO
      END DO
      
    END SUBROUTINE sum_mass2


!------------------------------------------------------------------------------

  END SUBROUTINE gmxe_driver

!=============================================================================!

  SUBROUTINE gmxe_free_memory

    IMPLICIT NONE


    IF (ASSOCIATED(cph))      DEALLOCATE(cph)
    IF (ASSOCIATED(ccn))      DEALLOCATE(ccn)
    IF (ASSOCIATED(anum))     DEALLOCATE(anum)
    IF (ASSOCIATED(rwet))     DEALLOCATE(rwet)
    IF (ASSOCIATED(rdry))     DEALLOCATE(rdry)
    IF (ASSOCIATED(ddry))     DEALLOCATE(ddry)
    IF (ALLOCATED(diagaer))   DEALLOCATE(diagaer)
    IF (ALLOCATED(xconvert))  DEALLOCATE(xconvert)
    IF (ALLOCATED(zconvert))  DEALLOCATE(zconvert)
    IF (ALLOCATED(lrhoa))     DEALLOCATE(lrhoa)
    IF (ALLOCATED(cdiagaer))   DEALLOCATE(cdiagaer)
    
  END SUBROUTINE gmxe_free_memory

!=============================================================================!

  SUBROUTINE gmxe_read_nml_cpl(status, iou)

    ! GMXe MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to online tracers

    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, &
                                        read_nml_close
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    
    INTRINSIC :: TRIM
    
    IF(p_parallel_io) &
      CALL start_message_bi(modstr,'Reading coupling namelist',substr)
    
    status = 1
    
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    !mz_kp_20080112 +
    !       Added choice of call routine location
    IF (p_parallel_io) THEN  
      IF(TRIM(ADJUSTL(driver_call)) /= 'radiation' .AND. &
        TRIM(ADJUSTL(driver_call)) /= 'physc')    THEN
        CALL error_bi(&
          'Choice of driver call (in gmxe.nml) not possible', substr)  
      ENDIF
    ENDIF
    !mz_kp_20080112 - 
    
    CALL read_nml_close(substr, iou, modstr)
    
    status = 0  ! no ERROR
    
    IF(p_parallel_io) &
      CALL end_message_bi(modstr,'Reading coupling namelist',substr)
    
  END SUBROUTINE gmxe_read_nml_cpl

!==============================================================================!

  SUBROUTINE gmxe_init_coupling

    USE messy_main_data_bi,          ONLY: basemod => modstr
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nproma, kproma, nlevp1
#if defined(ECHAM5) || defined(CESM1)           
    USE messy_main_grid_def_mem_bi,  ONLY:nvclev, vct
#endif
#if defined(COSMO) || defined(MESSYDWARF)
    USE messy_main_grid_def_bi,      ONLY:h_a => hyai, h_b=>hybi
#endif    
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, &
                                           p_io, p_bcast
    USE messy_main_tools,            ONLY: find_next_free_unit, strcrack
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info,&
                                           set_channel_object_restreq
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_constants_mem,    ONLY: M_air, STRLEN_MEDIUM
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer,           ONLY: get_tracer, AEROSOL, R_molarmass, &
                                           I_AEROSOL_MODE,                   &
                                           R_MOLARMASS, R_AEROSOL_DENSITY,   &
                                           R_Henry_T0, R_henry_Tdep,         &
                                           R_alpha_T0, R_alpha_Tdep,         &
                                           I_CHARGE, get_chemprop

    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer_mem_bi,    ONLY: ti_gp, GPTRSTR
    USE messy_gmxe_aerchem_liq,      ONLY: jval_h2o2, jval_no3, jval_o3 &
                                         , hs0, dht, alpha0, alpha_T    &
                                         , kpp_l_idx, gas_idx, gas_mw   &
                                         , liq_dens, liq_idx, liq_mw    &
                                         , lspec_liq, lspec_gas
    USE MESSY_GMXE_AERCHEM_KPP,      ONLY: str_field_kpp_l => spc_names
    USE messy_gmxe_kappa,            ONLY: initialise_salts
    USE messy_gmxe_isorropia2,       ONLY: init_species_isorropia
    USE messy_gmxe_eqsam4clim,       ONLY: init_species_eqsam4clim
 
    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'gmxe_init_coupling'

    INTEGER  :: it, jc, jt, jm, klev, ierr, dummy_idt, jk
    INTEGER  :: counter, counter2, dummy, status
    INTEGER  :: idx1, idx2, idx, idt
    INTEGER  :: icharge
    REAL(dp) :: hypi(nlevp1), sfpress, altitude
#if defined(ECHAM5) || defined(CESM1)
    REAL(dp) :: h_a(nvclev), h_b(nvclev) 
#endif 
    CHARACTER(LEN=10)                        :: cname
    CHARACTER(LEN=2)                         :: csubname = ''  ! tracer subname
    CHARACTER(LEN=26), POINTER, DIMENSION(:) :: strname => NULL()
    IF(p_parallel_io) &
      CALL start_message_bi(modstr,'INIT COUPLING',substr) 

    ! get tracer properties from CHEMPROP + optionally modidified by tracer.nml)
    DO jt=1,ngas
       cname = TRIM(td%gas(jt)%name)
       CALL get_tracer(status, GPTRSTR, TRIM(cname), idx=idt)  
       IF (status == 0) THEN
          CALL get_tracer(status,GPTRSTR, idt, R_molarmass, td%gas(jt)%molmass)
          CALL tracer_halt('GAS: molar mass not available '//substr, status)
          CALL get_tracer(status,GPTRSTR, idt, R_aerosol_density &
               , td%gas(jt)%density)
          CALL tracer_halt('GAS: density not available '//substr, status)
       ELSE
          status = get_chemprop(TRIM(cname), R_molarmass, td%gas(jt)%molmass)
          IF (status /= 0) CALL error_bi (&
               'species '//TRIM(cname)//' not part of chemprop',substr) 
          status = get_chemprop(TRIM(cname), R_aerosol_density, td%gas(jt)%density)
       END IF
       td%gas(jt)%density = td%gas(jt)%density * 1.e-3
       IF ( trim(cname) == "H2SO4") td%gas_sulph_idx = jt
    END DO

    DO jm=1,nmod
       csubname=TRIM(cmodes(jm))
       DO jt=1,nanions
          cname = TRIM(td%anion(jt)%name)
          CALL get_tracer(status, GPTRSTR, TRIM(cname), subname=csubname &
               , idx=idt)
          IF (status == 0) THEN
             CALL get_tracer(status,GPTRSTR, idt, R_molarmass &
                  , td%anion(jt)%molmass)
             CALL tracer_halt('ANIONS: molmass not available '//substr,status)
             CALL get_tracer(status,GPTRSTR, idt, R_aerosol_density &
                  , td%anion(jt)%density)
             CALL tracer_halt('ANIONS: density not available '//substr,status)
             CALL get_tracer(status,GPTRSTR, idt, I_CHARGE,icharge)
             CALL tracer_halt('ANIONS: charge not available '//substr,status)
          ELSE
             status = get_chemprop(TRIM(cname), R_molarmass, td%anion(jt)%molmass)
             IF (status /= 0) CALL error_bi (&
                  'species '//TRIM(cname)//' not part of chemprop',substr) 
             status = get_chemprop(TRIM(cname), R_aerosol_density &
                  , td%anion(jt)%density)
             status = get_chemprop(TRIM(cname), I_CHARGE, icharge)
          END IF
          td%anion(jt)%charge  = REAL(icharge,dp)
          td%anion(jt)%density = td%anion(jt)%density * 1.e-3
       END DO


       DO jt=1,ncations
          cname = TRIM(td%cation(jt)%name)
          CALL get_tracer(status, GPTRSTR, TRIM(cname), subname=csubname &
               , idx=idt)
          IF (status == 0) THEN
             CALL get_tracer(status,GPTRSTR, idt, R_molarmass &
                  , td%cation(jt)%molmass)
             CALL tracer_halt('CATIONS: molmass not available '//substr,status)
             CALL get_tracer(status,GPTRSTR, idt, R_aerosol_density &
                  , td%cation(jt)%density)
             CALL tracer_halt('CATIONS: density not available '//substr,status)
             CALL get_tracer(status,GPTRSTR, idt, I_CHARGE, icharge)
             CALL tracer_halt('CATIONS: charge not available '//substr,status)
          ELSE
             status = get_chemprop(TRIM(cname), R_molarmass, td%cation(jt)%molmass)
             IF (status /= 0) CALL error_bi (&
                  'species '//TRIM(cname)//' not part of chemprop', substr) 
             status = get_chemprop(TRIM(cname), R_aerosol_density &
                  , td%cation(jt)%density)
             status = get_chemprop(TRIM(cname), I_CHARGE, icharge)
          END IF
          td%cation(jt)%charge  = REAL(icharge,dp)
          td%cation(jt)%density = td%cation(jt)%density * 1.e-3
          IF ( trim(td%cation(jt)%name) == "Hp") THEN
             td%hp_idx = jt
          END IF
       END DO

       DO jt=1,nsolutes
          cname = TRIM(td%solute(jt)%name)
          CALL get_tracer(status, GPTRSTR, TRIM(cname), subname=csubname &
               , idx=idt)
          IF (status == 0) THEN
             CALL get_tracer(status,GPTRSTR, idt, R_molarmass &
                  , td%solute(jt)%molmass)
             CALL tracer_halt('SOLUTES: molmass not available '//substr,status)
             CALL get_tracer(status,GPTRSTR, idt, R_aerosol_density &
                  , td%solute(jt)%density)
             CALL tracer_halt('SOLUTES: density not available '//substr,status)
          ELSE
             status = get_chemprop(TRIM(cname), R_molarmass, td%solute(jt)%molmass)
             IF (status /= 0) CALL error_bi (&
                  'species '//TRIM(cname)//' not part of chemprop', substr) 
             status = get_chemprop(TRIM(cname), R_aerosol_density &
                  , td%solute(jt)%density)
          END IF
          td%solute(jt)%density = td%solute(jt)%density * 1.e-3
       END DO
       
    END DO

    hs0(:)     = 0._dp
    dht(:)     = 0._dp
    alpha0(:)  = 0._dp
    alpha_T(:) = 0._dp

    DO jt=1,lspec_gas
       IDX = kpp_l_idx%gas_spec(jt,gas_idx)
       cname = str_field_kpp_l(idx)
       CALL get_tracer(ierr, GPTRSTR, TRIM(cname), idx=idt)
       IF (status == 0) THEN
          CALL get_tracer(status,GPTRSTR, idt, R_molarmass &
               , KPP_L_IDX%GAS_ATTR(JT,GAS_MW))
          CALL tracer_halt('CATIONS: molmass not available '//substr,status)
          CALL get_tracer(status,GPTRSTR, idt, R_HENRY_T0,HS0(idx))
          CALL tracer_halt('AERCHEM GAS: Henry_T0 not available '//substr,status)
          CALL get_tracer(status,GPTRSTR, idt, R_HENRY_Tdep,DHT(idx))
          CALL tracer_halt('AERCHEM GAS: Henry_Tdep not available '//substr,status)
          CALL get_tracer(status,GPTRSTR, idt, R_ALPHA_T0,alpha0(idx))
          CALL tracer_halt('AERCHEM GAS: alpha_T0 not available '//substr,status)
          CALL get_tracer(status,GPTRSTR, idt, R_alpha_Tdep, alpha_T(idx))
          CALL tracer_halt('AERCHEM GAS: Alpha_Tdep not available '//substr,status)
        ELSE
           status = get_chemprop( TRIM(cname), R_molarmass &
                , KPP_L_IDX%GAS_ATTR(JT,GAS_MW))
           IF (status /= 0) CALL error_bi (&
                'species '//TRIM(cname)//' not part of chemprop', substr) 
           status = get_chemprop(TRIM(cname), R_Henry_T0,   HS0(idx))
           status = get_chemprop(TRIM(cname), R_Henry_Tdep, DHT(idx))
           status = get_chemprop(TRIM(cname), R_ALpha_T0,   alpha0(idx))
           status = get_chemprop(TRIM(cname), R_alpha_Tdep, alpha_T(idx))
        END IF
     ENDDO
    do jm = 1, nsoluble
       IF (.NOT.AERCHEM(jm) ) CYCLE
       ! use all modes in which aerchem shall be applied
       ! set subname according to the actual modes
       csubname=TRIM(cmodes(jm))
       DO jt=1,lspec_liq
          IDX = kpp_l_idx%liq_spec(jt,liq_idx)
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname)     
          call strcrack(STR_FIELD_KPP_L(idx),'_', strname, dummy)
          cname=TRIM(strname(1))
          CALL get_tracer(status, GPTRSTR, TRIM(cname),  subname = csubname &
               , idx = idt)
          IF (status == 0) THEN
             CALL get_tracer(status,GPTRSTR, idt, R_molarmass &
                  , KPP_L_IDX%LIQ_ATTR(JT,LIQ_MW))
             CALL tracer_halt('CATIONS: molmass not available '//substr, status)
             CALL get_tracer(status,GPTRSTR, idt, R_aerosol_density &
                  , KPP_L_IDX%LIQ_ATTR(JT,LIQ_DENS))
             CALL tracer_halt('CATIONS: density not available '//substr,status)
          ELSE
             status = get_chemprop(TRIM(cname), R_molarmass &
                  , KPP_L_IDX%LIQ_ATTR(JT,LIQ_MW))
             IF (status /= 0) CALL error_bi (&
                  'species '//TRIM(strname(1))//' not part of chemprop'&
                  , substr) 
             status = get_chemprop(TRIM(cname), R_aerosol_density &
                  , KPP_L_IDX%LIQ_ATTR(JT,LIQ_DENS))
          END IF
          KPP_L_IDX%LIQ_ATTR(JT,LIQ_DENS) =  KPP_L_IDX%LIQ_ATTR(JT,LIQ_DENS) * 1.e-3
       END DO
    ENDDO

    !--- Initialisations: -----------------------------------------------
    klev = nlev

    ! get h2o tracer idt
    CALL info_bi('Looking for H2O tracer and check its idt:', substr)
    CALL get_tracer(ierr, GPTRSTR, 'H2O', subname='', idx=idt_H2O)
    IF (idt_H2O == 0) THEN
      CALL info_bi('H2O tracer not available.', substr)
    ELSE
      CALL info_bi('H2O tracer found.', substr)
    ENDIF
    
    ! get so42m_ns tracer idt
    CALL info_bi('Looking for sulphate tracer in the nucleation mode:', substr)
    CALL get_tracer(ierr, GPTRSTR, 'SO4mm', subname='ns', idx=dummy_idt)
    IF (dummy_idt == 0) THEN
      CALL info_bi('Nucleation mode sulphate tracer not available.', substr)
      CALL info_bi('Therefore, no sulphate aerosol nucleation allowed!', substr)
      CALL info_bi('!!! WARNING !!! LNUCL is set to FALSE !!!', substr)
      lnucl =.FALSE.
    ELSE
      CALL info_bi('Nucleation mode sulphate tracer found.', substr)
    ENDIF

    ! Setting up species structure and corresponding arrays 
    ! (conversion factors)

    CALL info_bi('Setting up species structure.', substr)
!mz_se_20170313 moved the IF(L_ORACLE) block to here. Was originally just before the end of this subroutine
! the channels here are needed in INIT_CPL_SPECIES_E5. Changed IF(status == 1) to IF(status /= 0)
!mz_sb_20160211+
    IF (L_ORACLE) THEN
      CALL get_channel_object(status, 'oracle','nmode', p0=nomode)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NfPOA', p0=NfPOA)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NbbPOA', p0=NbbPOA)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NfSOAsv', p0=NfSOAsv)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NbbSOAsv', p0=NbbSOAsv)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NfSOAiv', p0=NfSOAiv)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NbbSOAiv', p0=NbbSOAiv)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','NSOAv', p0=NSOAv)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','tmode', p1=tomode)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
      CALL get_channel_object(status, 'oracle','kMOM', p1=kMOM)
      IF (status /= 0) &
           CALL error_bi('channel object from ORACLE not found', substr)
    ENDIF
!mz_sb_20160211-


    CALL INIT_CPL_SPECIES_E5
    CALL initialise_salts(nproma,nlev,nmod)
    SELECT CASE (neqm)
    CASE(1)
       CALL init_species_eqsam4clim(2)
    CASE(2)
       CALL init_species_isorropia(2)
    END SELECT
    CALL info_bi('Species structure finished.', substr)
    
    ALLOCATE(zconvert (nproma,nlev,0:naertot)); zconvert(:,:,:) = FACp0
    ALLOCATE(xconvert (nproma,nlev,0:naertot)); xconvert(:,:,:) = FACp0
    ALLOCATE(lrhoa    (0:naertot))            ; lrhoa       (:) = .TRUE.  
    
    DO jt=1,spec_number
      DO jm = 0,nmod
        
        idx1 = species(jt)%tracidx(jm)
        IF (idx1 == 0) cycle
        idx2 = species(jt)%aermlidx(jm)
        lrhoa(idx2) =.true.
        
        !--- Assign conversion factors depending on the unit of the tracer:
        
        
        SELECT CASE (ti_gp(idx1)%tp%ident%unit)
          
          !--- Particle mass:
        CASE ('kg/kg')
          ! [kg kg^-1(air)]   => [umol m^-3 (air)]
          xconvert(1:nproma,1:klev,idx2) = FACp9    & 
            / ti_gp(idx1)%tp%meta%cask_r(R_molarmass)
          
        CASE ('mol/mol')
          !  [mol mol^-1 (air)] => [umol g^-1 (air)]
          xconvert(1:nproma,1:klev,idx2) = FACp9 / M_air 
          
          !--- Particle number:
          
        CASE ('1/kg')
          !  [# kg^-1 (air)]  => [# mg^- (air)]
          xconvert(1:nproma,1:klev,idx2) = FACm6 
          
        CASE ('1/mol')
          !  [# mol^-1 (air)] => [# cm^-3 (air)]
          xconvert(1:nproma,1:klev,idx2) = FACm3 / M_air
          
        CASE ('1/cm3')
          !  [# cm^-3 (air)] => [# cm^-3 (air)]
          lrhoa(idx2)    = .false.
          xconvert(1:nproma,1:klev,idx2) = FACp0
          
        END SELECT
      END DO
    END DO


    IF (p_parallel_io) THEN
      WRITE(*,*) 'tropopause channel: ', TRIM(Tropopchannel)
      WRITE(*,*) 'channel element containing tropopause index: ', &
        TRIM(TropopIndex)
    END IF
    
    ierr_tropop = 1
    CALL get_channel_object(ierr_tropop, TRIM(Tropopchannel), TRIM(TropopIndex) &
      , p2=tp_i0)
    IF (p_parallel_io .AND. ierr_tropop /= 0) &
      CALL warning_bi('channel object for tropopause index not found ', substr)
    
    IF (p_parallel_io) THEN
      WRITE(*,*) 'PSC region channel: ', TRIM(Pscchannel)
      WRITE(*,*) 'channel element containing PSC region index: ', TRIM(Pscreg)
    END IF
    
    CALL get_channel_info(ierr_psc, TRIM(Pscchannel))
    
    IF (ierr_psc == 0) THEN
      CALL get_channel_object(ierr, TRIM(Pscchannel), TRIM(phase), p3=flt_pscreg1)
      IF(p_parallel_io .AND. ierr /= 0) &
        CALL warning_bi('channel object for PSC phase index not found', substr)
      
      CALL get_channel_object(ierr, TRIM(Pscchannel), TRIM(Pscreg), p3=flt_pscreg0)
      IF(p_parallel_io .AND. ierr /= 0) &
        CALL warning_bi('channel object for PSC region index not found', substr)
    ELSE
      IF (p_parallel_io) &
        CALL warning_bi('PSC CHANNEL NOT AVAILABLE, USING DEFAULT VALUES', substr)
    END IF ! ierr_psc
    
    IF (p_parallel_io .AND. ierr_tropop /= 0 .AND. ierr_psc /= 0) THEN
      CALL warning_bi(&
        'Caution. NO information on tropopause height and PSC region available!', &
        substr)
      CALL warning_bi('Caution. GMXe will calculate stratospheric aerosols.', substr)
    END IF ! ierr_tropop .and. ierr_psc


! This part deals with the emission fluxes
    
    emis: IF (l_calc_emis) THEN
      CALL gmxe_emis_init_e5        
    END IF emis
    
    IF (l_aerchem) THEN
      CALL info_bi('... getting photolysis rates in liquid phase !', substr)
      
      CALL get_channel_object(status, TRIM(photol(1)), 'J_H2O2', p3=jval_h2o2)
      IF (status /= 0) THEN
        CALL info_bi('J_H2O2 channel object not found !', substr)
        CALL info_bi('no aqueous phase photolysis of H2O2', substr)
        CALL info_bi('There will be no photolysis in aqueous phase!', substr)
      ELSE
        CALL info_bi( &
          'WARNING !!! Rerun activated for photolysis rate of H2O2', substr)
        CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_H2O2' )
        CALL channel_halt(substr, status)
      END IF
      
      CALL get_channel_object(status, TRIM(photol(1)), 'J_NO2', p3=jval_no3)
      IF (status /= 0) THEN
        CALL info_bi('J_NO2 channel object not found !', substr)
        CALL info_bi('no aqueous phase photolysis of NO3', substr)
        CALL info_bi('There will be no photolysis in aqueous phase!', substr)
      ELSE
        CALL info_bi( &
          'WARNING !!! Rerun activated for photolysis rate of NO3', substr)
        CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_NO2' )
        CALL channel_halt(substr, status)
      END IF
      
      CALL get_channel_object(status, TRIM(photol(1)), 'J_O1D', p3=jval_o3)
      IF (status /= 0) THEN
        CALL info_bi('J_O1D channel object not found !', substr)
        CALL info_bi('no aqueous phase photolysis of O3', substr)
        CALL info_bi('There will be no photolysis in aqueous phase!', substr)
      ELSE
        CALL info_bi( &
          'WARNING !!! Rerun activated for photolysis rate of O3', substr)
        CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_O1D' )
        CALL channel_halt(substr, status)
      END IF
      
    END IF
    
    ! determine altitude (level) for stratospheric evaporation
    altitude = 25000._dp
#if !defined(COSMO) && !defined(MESSYDWARF)
       DO jk=1,nvclev
         h_a(jk) = vct(jk)
         h_b(jk) = vct(jk+nvclev)
       ENDDO
#endif
    sfpress     = 101325._dp   ! reference pressure of 1000 hPa
    do jk=1,nlev+1
      hypi(jk)      = h_a(jk) + h_b(jk) * sfpress
    enddo
    k300 = 1
    do jk=1,nlev
      if (hypi(jk) < altitude .and. hypi(jk+1) >= altitude) then
        k300 = jk
        exit
      endif
    end do

    
    IF(p_parallel_io) &
         CALL info_bi('Coupling to main data',substr) 
    CALL get_channel_object(status, 'grid_def','grvol', p3=grvol)
    IF (status == 1) &
         CALL error_bi('channel object for grvol not found', substr)
    CALL get_channel_object(status, 'grid_def','grmass', p3=grmass)
    IF (status == 1) &
         CALL error_bi('channel object for grmass not found', substr)

    CALL get_channel_object(status, basemod,'aclc', p3=aclc)
    IF (status == 1) &
         CALL error_bi('channel object for aclc not found', substr)
    CALL get_channel_object(status, basemod,'acdnc', p3=acdnc)
    IF (status == 1) &
         CALL error_bi('channel object for acdnc not found', substr)
    
    CALL get_channel_object(status, basemod,'rhum', p3=rhum_3d)
    IF (status == 1) &
         CALL error_bi('channel object for rhum not found', substr)
    CALL get_channel_object(status, basemod,'press', p3=press_3d)
    IF (status == 1) &
         CALL error_bi('channel object for press not found', substr)
    CALL get_channel_object(status, basemod,'pressi', p3=pressi_3d)
    IF (status == 1) &
         CALL error_bi('channel object for pressi not found', substr)
    CALL get_channel_object(status, 'grid_def','deltaz', p3=deltaz)
    IF (status == 1) &
         CALL error_bi('channel object for deltaz not found', substr)

#if defined(COSMO) || defined(MESSYDWARF)
   CALL get_channel_object(status, basemod,'vervel', p3=vervel_3D)
    IF (status == 1) &
         CALL error_bi('channel object for vervel not found', substr)
#endif

    IF(p_parallel_io) &
      CALL end_message_bi(modstr,'INIT COUPLING',substr) 
    
  END SUBROUTINE gmxe_init_coupling

!==============================================================================!

  SUBROUTINE INIT_CPL_SPECIES_E5

    USE messy_main_tracer_mem_bi,      ONLY: ti_gp, ntrac => ntrac_gp
    USE messy_main_tracer,             ONLY: AEROSOL, AIR, i_aerosol_mode, &
                                             R_molarmass, R_aerosol_density, &
                                             CLOUD, s_aerosol_model
    USE MESSY_MAIN_TOOLS,              ONLY: strcrack
    USE messy_main_mpi_bi,             ONLY: p_io, p_bcast, p_parallel_io, p_pe
    USE MESSY_GMXE_AERCHEM_LIQ
    USE MESSY_GMXE_AERCHEM_KPP,        ONLY: str_field_kpp_l => spc_names
    USE MESSY_GMXE_OC_AGING,           ONLY: num_wsoc, kappa, init_oc_aging, &
                                             spec_idx_wsoc
    USE MESSY_GMXE_SOA,                ONLY: nsoa, noxi, ngas_soa => ngas, names


    INTEGER :: jm ,jc, jt, jk !mz_sb_20161020: jk added
    INTEGER :: counter
    LOGICAL :: FOUND=.FALSE.

    CHARACTER(LEN=STRLEN_MEDIUM) :: names_array(NMAX_SPEC)
    CHARACTER(LEN=STRLEN_MEDIUM) :: cname

    INTEGER                      :: dummy, idx1
    CHARACTER(LEN=26), POINTER   :: strname(:) => NULL()
    CHARACTER(LEN=2)             :: str_wsoc
    !mz_sb_20161020+
    CHARACTER(LEN=2)             :: str_soa
    CHARACTER(LEN=2)             :: str_mom
    !mz_sb_20161020-

    CHARACTER(LEN=2)             :: str_pa ! mz_dk_20120120
    INTEGER                      :: spec_idx_pa(1:num_pa) ! mz_dk_20120120
!mz_sb_20160211+
!    INTEGER                      :: spec_idx_fPOA(NfPOA),     spec_idx_bbPOA(NbbPOA),      &
!                                    spec_idx_fSOAsv(NfSOAsv), spec_idx_bbSOAsv(NbbSOAsv),  &
!                                    spec_idx_fSOAiv(NfSOAiv), spec_idx_bbSOAiv(NbbSOAiv),  &
!                                    spec_idx_SOAv(NSOAv,maxMOM) 
!mz_sb_20160211-

    cname          = ""
    names_array(:) = ""

    counter = 0
      ! thermodynamic / inorganic partitioning compounds
    DO jc = 1, ngas
      cname = TRIM(td%gas(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,nanions
      cname = TRIM(td%anion(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,ncations
      cname = TRIM(td%cation(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO

    DO jc=1,nsolutes
      cname = TRIM(td%solute(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,nbulk
      cname = TRIM(bulk(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
! This one is for numbers per mode
    IF (.NOT. lnumber) THEN
      DO jm=1,nmod
        cname = 'N'
        DO jt=1,counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
            found=.true.
        END DO
        IF (.NOT.FOUND) THEN
          counter = counter + 1
          names_array(counter) = cname
        END IF
      END DO
    END IF
      
    IF (l_aerchem) THEN
      ! additional species from aerchem (aerosol phase)
      DO jm=1,nsoluble
        IF (.NOT.AERCHEM(jm)) CYCLE
        
        DO jc = 1, lspec_liq
          idx1 = kpp_l_idx%liq_spec(jc, liq_idx)
          ! check if a tracer with the same name already exists for each mode 
          ! treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          IF (TRIM(strname(1)) == 'Prod') CYCLE
          cname=TRIM(strname(1))
          found=.false.
          DO jt=1,counter
            IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
          END DO
          IF (.NOT.FOUND) THEN
            counter = counter + 1
            names_array(counter) = cname
          END IF
        END DO
      END DO
      
      ! additional species from aerchem (gas phase)
      jm = 0
      
      DO jc = 1, lspec_gas
        idx1 = kpp_l_idx%gas_spec(jc, gas_idx)
        ! check if a tracer with the same name already exists for each mode 
        ! treated by aerchem
        if (associated(strname)) DEALLOCATE (strname)
        NULLIFY(strname) 
        CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
        cname=TRIM(strname(1))
        found=.false.
        DO jt=1,counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.FOUND) THEN
          counter = counter + 1
          names_array(counter) = cname
        END IF
      END DO
      
    END IF ! aerchem

!mz_sb_20160211+
! additional species from oracle submodel (aerosol phase) 
    IF (L_ORACLE) THEN
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NfPOA)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "fPOA"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NbbPOA)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "bbPOA"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NfSOAsv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "fSOAsv"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NbbSOAsv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "bbSOAsv"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NfSOAiv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "fSOAiv"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NbbSOAiv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             found=.false.
             cname = "bbSOAiv"//str_soa
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
          END DO
       END DO
       DO jm = NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc = 1,NINT(NSOAv)
           DO jk = 1,NINT(kMOM(jc))
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             IF (jk < 10) then
                write (str_mom,'(A,I1)') "0",jk
             ELSE
                write (str_mom,'(I2)') jk
             END IF
             found=.false.
             cname = "SOAv"//str_soa//str_mom
             DO jt=1,counter
                IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
             END DO
             IF (.NOT.FOUND) THEN
                counter = counter + 1
                names_array(counter) = cname
             END IF
            END DO
          END DO
       END DO
    ENDIF
!mz_sb_20160211-

    IF (L_OC_AGING) THEN
      ! get gas phase OH concentration into the species array
      found = .false.
      cname = "OH"
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
       
    ENDIF ! oc_aging
    
    ! mz_dk_20120119+
    ! PASSIVE AEROSOLS
    ! set names_array and counter
    IF (L_PASSIVE_AER) THEN
      DO jm = 1, nmod
        DO jc = 1, num_pa
          IF (jc < 10) then
            write (str_pa,'(A,I1)') "0",jc
          ELSE
            write (str_pa,'(I2)') jc
          END IF
          found=.false.
          cname = "PASSAER"//str_pa
          DO jt = 1, counter
            IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
          END DO
          IF (.NOT.found) THEN
            counter = counter+1
            names_array(counter) = cname
          END IF
        END DO
      END DO
    END IF  ! pass_aer
    ! mz_dk_20120119-

    IF (L_SOA) THEN
      DO jc=1,ngas_soa
        found=.false.
        cname = TRIM(names(NSOA+NOXI+jc))
        DO jt = 1, counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.found) THEN
          counter = counter+1
          names_array(counter) = cname
        END IF
      END DO
      DO jc=1,noxi
        found=.false.
        cname = TRIM(names(NSOA+jc))
        DO jt = 1, counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.found) THEN
          counter = counter+1
          names_array(counter) = cname
        END IF
      END DO
       ! aerosol species should be already included in the bulk structure
!!$       DO jm=1,nmod
!!$          IF (.NOT.LSOA(jm)) CYCLE
!!$          DO jc=1,nsoa
!!$             cname = TRIM(names(jc))
!!$             found=.false.
!!$             DO jt = 1, counter
!!$                IF ( TRIM(names_array(jt)) == TRIM(cname)) found=.true.
!!$             END DO
!!$             IF (.NOT.found) THEN
!!$                counter = counter+1
!!$                names_array(counter) = cname
!!$             END IF
!!$          END DO
!!$       END DO
    END IF ! SOA

    ! allocate the species index field
    spec_number = counter
    ALLOCATE(SPECIES(spec_number))
    DO jt=1,spec_number
      ALLOCATE(SPECIES(jt)%tracidx(0:nmod))
      ALLOCATE(SPECIES(jt)%aermlidx(0:nmod))
      ALLOCATE(SPECIES(jt)%aernlidx(0:nmod))
      ALLOCATE(SPECIES(jt)%kppidx(0:nmod))
      ALLOCATE(SPECIES(jt)%soaidx(0:nmod))
      ALLOCATE(SPECIES(jt)%zcoup(0:nmod))
      ALLOCATE(SPECIES(jt)%npassive(0:nmod)) ! mz_dk_20120120 PASSIVE TRACER or NOT?
      SPECIES(jt)%name = names_array(jt)
      DO jm=0,nmod
        SPECIES(jt)%tracidx(jm)  = 0
        SPECIES(jt)%aermlidx(jm) = 0
        SPECIES(jt)%aernlidx(jm) = 0
        SPECIES(jt)%kppidx(jm)   = 0
        SPECIES(jt)%soaidx(jm)   = 0
        SPECIES(jt)%zcoup(jm)    = 1.0_dp
        SPECIES(jt)%molmass      = 0.0_dp
        SPECIES(jt)%density      = 0.0_dp
        SPECIES(jt)%kappa        = 0.0_dp
        SPECIES(jt)%charge       = 0.0_dp
        SPECIES(jt)%npassive(jm) = .false. ! mz_dk_20120120
      END DO
    END DO


    ! determine tracer indices
    DO jt=1,ntrac
      IF (ti_gp(jt)%tp%ident%medium == AEROSOL) THEN
        IF (ti_gp(jt)%tp%meta%cask_s(S_aerosol_model) == 'gmxe') THEN
          jm = ti_gp(jt)%tp%meta%cask_i(I_aerosol_mode)
        ELSE
          CYCLE
        ENDIF
      ENDIF

      IF (ti_gp(jt)%tp%ident%medium == AIR) &
        jm = 0

      IF (ti_gp(jt)%tp%ident%medium == CLOUD) CYCLE
        
      DO jc=1,spec_number
        IF ( TRIM(ti_gp(jt)%tp%ident%basename) == TRIM(species(jc)%name) ) THEN
          SPECIES(jc)%tracidx(jm) = jt
          SPECIES(jc)%molmass     = ti_gp(jt)%tp%meta%cask_r(R_molarmass)
          IF (jm > 0) &
            SPECIES(jc)%density   = ti_gp(jt)%tp%meta%cask_r(R_aerosol_density) / 1.e3_dp
        END IF
      ENDDO
    ENDDO

    ! determine paerml indices and coupling switches

    counter = 0
    DO jm = 0, nmod
      DO jc=1, ngas
        IF (td%gas(jc)%ltreat(jm)) THEN
          cname = TRIM(td%gas(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == TRIM(species(jt)%name)) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%gas(jc)%aerml_idx(jm) = counter
              IF (lgas) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, nanions
        IF (jm == 0) CYCLE
        IF (td%anion(jc)%ltreat(jm)) THEN
          cname = TRIM(td%anion(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%anion(jc)%aerml_idx(jm) = counter
              species(jt)%charge = td%anion(jc)%charge
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, ncations
        IF (jm == 0) CYCLE
        IF (td%cation(jc)%ltreat(jm)) THEN
          cname = TRIM(td%cation(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%cation(jc)%aerml_idx(jm) = counter
              species(jt)%charge = td%cation(jc)%charge
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
       DO jc=1, nsolutes
        IF (td%solute(jc)%ltreat(jm)) THEN
          cname = TRIM(td%solute(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%solute(jc)%aerml_idx(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, nbulk
        IF (jm == 0) CYCLE
        IF (bulk(jc)%ltreat(jm)) THEN
          cname = TRIM(bulk(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              bulk(jc)%aerml_idx(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
              species(jt)%kappa = bulk(jc)%kappa
            END IF
          END DO
        END IF
      END DO
      IF (.NOT. lnumber) THEN
        IF (jm /= 0) THEN
          cname='N'
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              species(jt)%aernlidx(jm) = counter
              nnum(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END IF

    END DO
    
    counter = 0
    do jt=1,spec_number
      DO jm=0,nmod
        counter = MAX(counter, SPECIES(jt)%aermlidx(jm))
      ENDDO
    ENDDO
    
    IF ( l_aerchem ) THEN

      ! aerchem liquid compounds
      
      DO jm= 1,nsoluble
        IF (.NOT. AERCHEM(jm)) CYCLE
        DO jc = 1,lspec_liq
          idx1 = kpp_l_idx%liq_spec(jc, liq_idx)
        ! check if a compound with the same name already exists for each mode 
        ! within aerml array treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          cname=TRIM(strname(1))
          DO jt=1,spec_number
            IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
              SPECIES(jt)%kppidx(jm)   = idx1
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              SPECIES(jt)%aermlidx(jm) = counter
            ENDIF
          ENDDO
        ENDDO

      ENDDO   ! aerchem modes

        ! aerchem gas compounds

      jm = 0
      DO jc = 1,lspec_gas
        idx1 = kpp_l_idx%gas_spec(jc, gas_idx)
        ! check if a compound with the same name already exists for each mode 
        ! within aerml array treated by aerchem
        if (associated(strname)) DEALLOCATE (strname)
        NULLIFY(strname) 
        CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
        cname=TRIM(strname(1))
        DO jt=1,spec_number
          IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
            SPECIES(jt)%kppidx(jm)   = idx1
            IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
            counter = counter + 1
            SPECIES(jt)%aermlidx(jm) = counter
          ENDIF
        ENDDO
      ENDDO
      
    ENDIF

!mz_sb_20160211+
    IF (L_ORACLE) THEN
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NfPOA)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "fPOA"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NbbPOA)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "bbPOA"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NfSOAsv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "fSOAsv"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NbbSOAsv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "bbSOAsv"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NfSOAiv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "fSOAiv"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NbbSOAiv)
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             cname = "bbSOAiv"//str_soa
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
          END DO
       END DO
       DO jm=NINT(tomode(1)),NINT(tomode(NINT(nomode)))
          DO jc=1,NINT(NSOAv)
           DO jk=1,NINT(kMOM(jc))
             IF (jc < 10) then
                write (str_soa,'(A,I1)') "0",jc
             ELSE
                write (str_soa,'(I2)') jc
             END IF
             IF (jk < 10) then
                write (str_mom,'(A,I1)') "0",jk
             ELSE
                write (str_mom,'(I2)') jk
             END IF
             cname = "SOAv"//str_soa//str_mom
             DO jt=1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                ENDIF
             END DO
           END DO
          END DO
       END DO
    END IF ! L_ORACLE
!mz_sb_20160211-

    IF (L_OC_AGING) THEN
       ! WSOC compounds treated already in the bulk structure
       ! gas phase OH
       jm = 0
       cname="OH"
       DO jt=1,spec_number
          IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
             IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
             IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
             counter = counter + 1
             SPECIES(jt)%aermlidx(jm) = counter
             SPECIES(JT)%zcoup(jm) = 0._dp
          ENDIF
       END DO
    END IF
    
    ! mz_dk_20120119+
    IF (L_PASSIVE_AER) THEN
       ! passive aerosol tracers
       ! set aerosol mass index
       DO jm = 1, nmod
          DO jc = 1, num_pa
             IF (jc < 10) then
                write (str_pa,'(A,I1)') "0",jc
             ELSE
                write (str_pa,'(I2)') jc
             END IF
             cname = "PASSAER"//str_pa
             DO jt = 1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                END IF
             END DO
          END DO
       END DO
    END IF
    ! mz_dk_20120119-

    IF (L_SOA) THEN
       ! gas phase first
       jm = 0
       DO jc=1,NSOA+NOXI+NGAS_SOA
          cname = TRIM(names(jc))
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%npassive(jm) = .false.
                IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                counter = counter + 1
                SPECIES(jt)%aermlidx(jm) = counter
                species(jt)%SOAIDX(jm) = jc
             END IF
          END DO
       END DO
       ! aerosols treated already in the bulk structure
       ! aermlidx is therefore already set, but SOAIDX 
       ! must still be determined
       DO jm = 1,nmod
         IF (.NOT.LSOA(jm)) CYCLE
         DO jc=1,NSOA
           cname = TRIM(names(jc))
           DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) &
               species(jt)%SOAIDX(jm) = jc
           END DO
         END DO
       END DO
         
!!$       DO jt=1,spec_number
!!$          DO jm=0,nmod
!!$             IF (species(jt)%SOAIDX(jm) == 0 ) CYCLE
!!$             IF (p_pe /= p_io) CYCLE
!!$             print*, "names, soa", species(jt)%name, jt, species(jt)%aermlidx(jm), &
!!$                  species(jt)%SOAIDX(jm), jm, species(jt)%tracidx(jm), &
!!$                  species(jt)%npassive(jm)
!!$          END DO
!!$       END DO
    ENDIF ! L_SOA

    DO jt=1,spec_number
      SELECT CASE (TRIM(species(jt)%name))
      CASE("SS")
!        spec_idx_ss = jt
      CASE("OC")
        spec_idx_oc = jt
      CASE("DU")
!        spec_idx_du = jt     
      CASE("BC")
!        spec_idx_bc = jt       
      CASE("Hp")   
        spec_idx_hp = jt
      CASE("OHm")   
        spec_idx_ohm = jt
      CASE("HSO4m")  
        spec_idx_hso4m = jt
      CASE("H2SO4")
        spec_idx_h2so4 = jt
      CASE("SO4mm")
        spec_idx_so4mm = jt
      CASE("NH3")
        spec_idx_nh3 = jt
      CASE("NH4p")
        spec_idx_nh4p = jt
        
      CASE("H2O")
        spec_idx_h2o = jt
        species(jt)%molmass = mwh2o
      CASE("OH")
        spec_idx_oh = jt
      CASE("Number", "N")
        number_idx = jt
      END SELECT
    END DO

!mz_sb_20160211+
    IF (L_ORACLE) THEN
       DO jc=1,NINT(NfPOA)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "fPOA"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_fPOA(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NbbPOA)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "bbPOA"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_bbPOA(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NfSOAsv)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "fSOAsv"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_fSOAsv(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NbbSOAsv)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "bbSOAsv"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_bbSOAsv(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NfSOAiv)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "fSOAiv"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_fSOAiv(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NbbSOAiv)
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          cname = "bbSOAiv"//str_soa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_bbSOAiv(jc) = jt
             END IF
          END DO
       END DO
       DO jc=1,NINT(NSOAv)
        DO jk=1,NINT(kMOM(jc))
          IF (jc < 10) then
             write (str_soa,'(A,I1)') "0",jc
          ELSE
             write (str_soa,'(I2)') jc
          END IF
          IF (jk < 10) then
             write (str_mom,'(A,I1)') "0",jk
          ELSE
             write (str_mom,'(I2)') jk
          END IF
          cname = "SOAv"//str_soa//str_mom
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0.14_dp
                !spec_idx_SOAv(jc,jk) = jt
             END IF
          END DO
        END DO
       END DO
    END IF ! L_ORACLE
!mz_sb_20160211-

    IF (L_OC_AGING) THEN
       jt = spec_idx_oc
       species(jt)%kappa = 0.01_dp
       DO jc=1,nbulk
         IF ( TRIM(species(jt)%name) == TRIM(bulk(jc)%name) ) &
           bulk(jc)%kappa = 0.01_dp
       END DO
       spec_idx_wsoc(0) = spec_idx_oc
       CALL init_oc_aging

       DO jc=1,num_wsoc
         IF (jc < 10) then
           write (str_wsoc,'(A,I1)') "0",jc
         ELSE
           write (str_wsoc,'(I2)') jc
         END IF
         cname = "WSOC"//str_wsoc
         DO jt=1,spec_number
           IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
             species(jt)%kappa = kappa(jc)
             spec_idx_wsoc(jc) = jt
           END IF
         END DO
         DO jt=1,nbulk
           IF (TRIM(bulk(jt)%name) == TRIM(cname) ) &
             bulk(jt)%kappa = kappa(jc)
         END DO
       END DO
    END IF
    
    ! mz_dk_2012019+
    IF (L_PASSIVE_AER) THEN
       ! set kappa, index, and npassive flag
       DO jc=1,num_pa
          IF (jc < 10) then
             write (str_pa,'(A,I1)') "0",jc
          ELSE
             write (str_pa,'(I2)') jc
          END IF
          cname = "PASSAER"//str_pa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0._dp
                spec_idx_pa(jc) = jt
                species(jt)%npassive(1:nmod)= .true.
                !print*,species(jt)%name,jt,species(jt)%npassive(1:nmod)
             END IF
          END DO
       END DO
    END IF
    ! mz_dk_2012019-
    naertot = counter


    DO jc=1,spec_number
      DO jm=0,nmod
        jt = species(jc)%aermlidx(jm)
        SELECT CASE (TRIM(species(jc)%name))
        CASE ("H2O")
          nwh2o(jm) = jt
        CASE ("H2SO4")
          nh2so4(jm) = jt
        CASE ("HSO4m")
          nhso4m(jm) = jt
        CASE ("SO4mm")  
          nso42m(jm) = jt
        CASE ("Hp")  
          nhp(jm) = jt
        CASE ("OHm")  
          nohm(jm) = jt
        END SELECT
      END DO
    END DO
    
!    counter = 0
!    DO jt=1,spec_number
!      DO jm=0,nmod
!        IF (SPECIES(jt)%tracidx(jm) /= 0) counter = counter + 1
!      END DO
!    END DO
!
!    print*, "end counter: ",counter, naertot

    IF (p_parallel_io) then 
      do jt=1,spec_number
        print*, "name: ",TRIM(species(jt)%name)
        DO jm=0,nmod
          print*, "mode: ", jm, "tracer_index: ", SPECIES(jt)%tracidx(jm),  &
                              "aerml_index: ",  SPECIES(jt)%aermlidx(jm), &
                              "kpp_index: ",  SPECIES(jt)%kppidx(jm),     &
                              "density: ",  SPECIES(jt)%density,          &
                              "molarmass: ",  SPECIES(jt)%molmass,        &
                              "kappa: ",  SPECIES(jt)%kappa,              &
                              "npassive: ", species(jt)%npassive(jm)
        ENDDO
      ENDDO
    END IF



  END SUBROUTINE INIT_CPL_SPECIES_E5

!==============================================================================

!==============================================================================!

  SUBROUTINE gmxe_vdiff

#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte
#endif
    USE messy_main_tracer_mem_bi, ONLY: ti_gp
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, avo => N_A
    USE messy_main_grid_def_bi,   ONLY: philat_2d
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, jrow, kproma

#ifdef ECHAM5
    USE messy_main_data_bi,       ONLY: pxtems
#endif
    USE messy_main_timer,         ONLY: time_step_len &
                                      , nstep=>current_time_step, lstart
    USE messy_main_mpi_bi,        ONLY: p_io, p_bcast, p_parallel_io
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='gmxe_vdiff'
    INTEGER                     :: it,jc,jm,jt,jk,idt,jmod,jl, mlev, ktop
    REAL(dp), POINTER           :: zxtems(:,:)
    REAL(dp)                    :: diameter3
    REAL(dp)                    :: convMtoN
    REAL(dp)                    :: convMolectoKg
    REAL(dp)                    :: zdp    (nproma,nlev)
    REAL(dp)                    :: zmass  (nproma,nlev)
    REAL(dp)                    :: znumber(nproma,nlev)
    REAL(dp)                    :: fac_N  (nproma,nlev)
    REAL(dp)                    :: conv_unit(nproma,nlev)
    REAL(dp)                    :: zdz    (nproma,nlev)
    REAL(dp)                    :: flux(nproma,nlev)
#ifdef MESSYTENDENCY
    REAL(dp)                    :: zxtte(nproma,nlev) ! op_pj_20181115
#endif

    IF (.not. l_calc_emis) RETURN

#if defined (ECHAM5)
   zxtems => pxtems(:,1,:,jrow)
#endif

    ! default unit for channel defintion of aerosol N [1/cm3]
    fac_N(1:kproma,1:nlev) = ZERO
    DO jk=1,nlev
      fac_N(1:kproma,jk) = grmass(_RI_XYZ__(1:kproma,jrow,jk)) &
                         / grvol(_RI_XYZ__(1:kproma,jrow,jk)) * FACm6
    END DO
    DO jk=1,nlev
      zdp(1:kproma,jk) = pressi_3d(_RI_XYZ__(1:kproma,jrow,jk+1)) - &
                         pressi_3d(_RI_XYZ__(1:kproma,jrow,jk))
    END DO
    zdz(:,:) = 0.0_dp

    zdz(1:kproma,2:nlev) = deltaz(_RI_XYZ__(1:kproma,jrow,2:nlev))
    ! 
    IF(.NOT. lnumber) THEN
      DO jm = 1,nmod
        jt = species(number_idx)%tracidx(jm)
        fac_N(1:kproma,1:nlev) = FACp0
        SELECT CASE (TRIM(ti_gp(jt)%tp%ident%unit))
        CASE ('1/cm3')
          fac_N(1:kproma,1:nlev) = grmass(_RI_XYZ__(1:kproma,jrow,1:nlev)) &
                                 / grvol(_RI_XYZ__(1:kproma,jrow,1:nlev)) * FACm6
        CASE ('1/mol')
          fac_N(1:kproma,1:nlev) = M_air * FACm3
        CASE ('1/kg')
          fac_N(1:kproma,1:nlev) = FACp0
        CASE DEFAULT 
          IF (p_parallel_io) &
            CALL error_bi('no valid tracer unit for aerosol number!',&
            substr)
        END SELECT
      END DO
    END IF

! Loop over emission fluxes

    DO jc = 1,num_fluxes
      flux(:,:)  = 0._dp
      
      jm = emis_flux_array(jc)%mode
      convMtoN = 1._dp
      conv_unit(1:kproma,1:nlev) = 1._dp
      ! tracer numbers
      idt = species(number_idx)%tracidx(jm)
      conv_unit(1:kproma,1:nlev) = conv_unit(1:kproma,1:nlev) &
                                 * fac_N(1:kproma,1:nlev)
      IF (.NOT. ASSOCIATED(emis_flux_array(jc)%nflux) ) THEN
        diameter3 = emis_flux_array(jc)%diameter * emis_flux_array(jc)%diameter &
                  * emis_flux_array(jc)%diameter

        convMtoN  = 6._dp / pi &
                  / (emis_flux_array(jc)%density * 1.e3_dp) &
                  / (diameter3 * sigma_exp_ln(jm))
          
        SELECT CASE (TRIM(emis_flux_array(jc)%unit))
        CASE("kg/(m^2 s)")
          !conv_unit(1:kproma,2:nlev) = conv_unit(1:kproma,2:nlev)  &
          !                           * 1._dp
        CASE("kg/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = conv_unit(1:kproma,2:nlev)  &
                                     * zdz(1:kproma,2:nlev)
        CASE("molecules/(m^2 s)")
          conv_unit(1:kproma,2:nlev) = emis_flux_array(jc)%specs(1)%molarmass &
                                     * 1.e-3_dp * conv_unit(1:kproma,2:nlev)  &
                                     / avo 
        CASE("molecules/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = emis_flux_array(jc)%specs(1)%molarmass &
                                     * 1.e-3_dp * conv_unit(1:kproma,2:nlev)  &
                                     / avo * zdz(1:kproma,2:nlev)
        END SELECT
        
        IF (emis_flux_array(jc)%dim == 3) THEN         
          IF (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,_IZ_XYZ__)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%flux(_RI_XYZ__(1:kproma,jrow,1:mlev))
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%flux_2D(1:kproma,jrow)
        ENDIF
      ELSE
        IF ( emis_flux_array(jc)%dim == 3) THEN
          IF  (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,_IZ_XYZ__)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%nflux(_RI_XYZ__(1:kproma,jrow,1:mlev))
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%nflux_2D(1:kproma,jrow)
        ENDIF
      ENDIF
      
!!$      print*, "Flux before: ",jc, emis_flux_array(jc)%name, minval(flux(:,:)), MAXVAL(flux(:,:)), &
!!$        emis_flux_array(jc)%fac_num_emis, emis_flux_array(jc)%total_frac, jm, &
!!$        MINVAL(pxtte(:,:,idt)), MAXVAL(pxtte(:,:,idt))
      IF ( (l_tendency) .OR. &
        (emis_flux_array(jc)%dim == 3) ) THEN
#ifndef MESSYTENDENCY
        IF (emis_flux_array(jc)%NxD) THEN
          DO ktop=1,mlev
            DO jl=1,kproma
              jk= NINT(emis_flux_array(jc)%VIND(_RI_XYZ__(jl,jrow,ktop)))
              pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))      &
                + conv_unit(jl,jk) * convMtoN            &
                * flux(jl,ktop)                          &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
            END DO
          END DO
        ELSE
          DO jk=1,nlev
            DO jl=1,kproma
              pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))      &    
                + conv_unit(jl,jk) * convMtoN            &
                * flux(jl,jk)                            &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
            END DO
          END DO
        END IF
#else
        zxtte(:,:) = 0.0_dp
        IF (emis_flux_array(jc)%NxD) THEN
          DO ktop=1,mlev
            DO jl=1,kproma
              jk= NINT(emis_flux_array(jc)%VIND(_RI_XYZ__(jl,jrow,ktop)))
              zxtte(jl,jk) = &
                  conv_unit(jl,jk) * convMtoN            &
                * flux(jl,ktop)                          &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
            END DO
          END DO
        ELSE
          DO jk=1,nlev
            DO jl=1,kproma
              zxtte(jl,jk) = &
                 conv_unit(jl,jk) * convMtoN            &
                * flux(jl,jk)                            &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
            END DO
          END DO
        END IF
        CALL mtend_add_l(my_handle_em, idt, px = zxtte)
#endif
      ELSE
#if defined (ECHAM5)
        DO jl=1,kproma
          zxtems(jl,idt) = zxtems(jl,idt)              &
            + conv_unit(jl,nlev) * convMtoN            &
            * flux(jl,nlev)                            &
            * emis_flux_array(jc)%fac_num_emis         &
            * emis_flux_array(jc)%total_frac
        END DO
#endif
      ENDIF
!!$      print*, "Flux after: ",jc, emis_flux_array(jc)%name, minval(flux(:,:)), MAXVAL(flux(:,:)), &
!!$        emis_flux_array(jc)%fac_num_emis, emis_flux_array(jc)%total_frac, jm, &
!!$        MINVAL(pxtte(:,:,idt)), MAXVAL(pxtte(:,:,idt)), "add term: ", &
!!$        MINVAL( conv_unit(:,:) *convmton * flux(:,:) *  emis_flux_array(jc)%fac_num_emis       &
!!$                * emis_flux_array(jc)%total_frac / zdp(:,:) * g ), &
!!$        MAXVAL( conv_unit(:,:) *convmton * flux(:,:) *  emis_flux_array(jc)%fac_num_emis       &
!!$                * emis_flux_array(jc)%total_frac / zdp(:,:) * g )


      ! tracer mass
      convMtoN = 1._dp
      
      DO jt = 1, emis_flux_array(jc)%num_spec_emis
        flux(:,:)  = 0._dp
        conv_unit(:,:) = 1._dp
        idt = emis_flux_array(jc)%specs(jt)%trac_idx
        IF (idt == 0) CYCLE
        SELECT CASE (TRIM(emis_flux_array(jc)%unit))
        CASE("kg/(m^2 s)")
          conv_unit(1:kproma,1:nlev) = &
            M_air / emis_flux_array(jc)%specs(jt)%molarmass
        CASE("kg/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = &
            M_air / emis_flux_array(jc)%specs(jt)%molarmass * &
            zdz(1:kproma,2:nlev)
        CASE("molecules/(m^2 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp
        CASE("molecules/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp * &
                                       zdz(1:kproma,2:nlev)
        END SELECT

        IF (emis_flux_array(jc)%dim == 3) THEN
          IF (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,_IZ_XYZ__)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%flux(_RI_XYZ__(1:kproma,jrow,1:mlev))
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%flux_2D(1:kproma,jrow)
        ENDIF
        
        IF ( (l_tendency) .OR. &
             (emis_flux_array(jc)%dim == 3) ) THEN

#ifndef MESSYTENDENCY
          IF (emis_flux_array(jc)%NxD) THEN
            DO ktop=1,mlev
              DO jl=1,kproma
                jk= NINT(emis_flux_array(jc)%VIND(_RI_XYZ__(jl,jrow,ktop)))
                pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))      &
                  + conv_unit(jl,jk) * convMtoN            &
                  * flux(jl,ktop)                          &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ELSE
            DO jk=2,nlev
              DO jl=1,kproma
                pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))      &    
                  + conv_unit(jl,jk) * convMtoN            &
                  * flux(jl,jk)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ENDIF
#else
          zxtte(:,:) = 0.0_dp
          IF (emis_flux_array(jc)%NxD) THEN
            DO ktop=1,mlev
              DO jl=1,kproma
                jk= NINT(emis_flux_array(jc)%VIND(_RI_XYZ__(jl,jrow,ktop)))
                zxtte(jl, jk) = &
                    conv_unit(jl,jk) * convMtoN            &
                  * flux(jl,ktop)                          &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ELSE
            DO jk=2,nlev
              DO jl=1,kproma
                zxtte(jl, jk) = &
                    conv_unit(jl,jk) * convMtoN            &
                  * flux(jl,jk)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO            
          ENDIF
          CALL mtend_add_l(my_handle_em, idt, px = zxtte)
#endif
        ELSE
#if defined (ECHAM5)
          DO jl=1,kproma
            zxtems(jl,idt) = zxtems(jl,idt)              &
              + conv_unit(jl,nlev) * convMtoN            &
              * flux(jl,nlev)                            &
              * emis_flux_array(jc)%specs(jt)%frac       &
              * emis_flux_array(jc)%total_frac
          END DO
#endif
        ENDIF
      END DO
    END DO

  END SUBROUTINE gmxe_vdiff

!-------------------------------------------------------------------------------
  SUBROUTINE gmxe_emis_init_e5


    USE MESSY_MAIN_TRACER,              ONLY: r_molarmass, i_aerosol_mode, &
                                              numberdensity
    USE MESSY_MAIN_TRACER_MEM_BI,       ONLY: ti_gp, ntrac => ntrac_gp
    USE messy_main_tools,               ONLY: strcrack, str2num

    USE messy_main_channel,             ONLY: get_channel_object,  &
                                              get_channel_object_info
    USE messy_main_channel_error_bi,    ONLY: channel_halt
    USE messy_main_channel_bi,          ONLY: GP_3D_MID, &
                                              GP_3D_1LEV, GP_2D_HORIZONTAL
    USE messy_main_mpi_bi,              ONLY: p_parallel_io

    IMPLICIT NONE
    INTEGER :: jm, jc, jt, dummy, counter, status, id_repr
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring2(:) => NULL()

    LOGICAL                                   :: found
    CHARACTER(LEN=*), PARAMETER               :: substr='gmxe_emis_e5'
    REAL(dp)                                  :: val

    DO jt=1,100
      emis_flux_list(jt)%name         = ""
      emis_flux_list(jt)%flux_name    = ""
      emis_flux_list(jt)%nflux_name   = ""
      emis_flux_list(jt)%channel_name = ""
      emis_flux_list(jt)%density      = 0._dp
      emis_flux_list(jt)%diameter     = 0._dp
      emis_flux_list(jt)%mode         = 0
      emis_flux_list(jt)%scal_fac     = 1._dp
      emis_flux_list(jt)%total_frac   = 1._dp
      emis_flux_list(jt)%fac_num_emis = 1._dp
    END DO
! seasalt

    emis_flux_list(1)%name          = "seasalt_mass_ks"
    emis_flux_list(1)%density       = 2.170_dp
    emis_flux_list(1)%diameter      = 2._dp * 0.035e-6_dp
    emis_flux_list(1)%unit          = "kg/(m^2 s)"
    emis_flux_list(1)%mode          = det_mode(nmod,"ks","KS")

    emis_flux_list(2)%name          = "monahan_seasalt_mass_as"
    emis_flux_list(2)%density       = 2.170_dp
    emis_flux_list(2)%diameter      = 2._dp * 0.258E-6_dp
    emis_flux_list(2)%unit          = "kg/(m^2 s)"
    emis_flux_list(2)%mode          = det_mode(nmod,"as","AS")
    emis_flux_list(2)%fac_num_emis  = 1._dp / sigma_exp_ln(emis_flux_list(2)%mode)

    emis_flux_list(3)%name          = "monahan_seasalt_mass_cs"
    emis_flux_list(3)%density       = 2.170_dp
    emis_flux_list(3)%diameter      = 2._dp * 1.65e-6_dp
    emis_flux_list(3)%unit          = "kg/(m^2 s)"
    emis_flux_list(3)%mode          = det_mode(nmod,"cs","CS")
    emis_flux_list(3)%fac_num_emis  = 1._dp / sigma_exp_ln(emis_flux_list(3)%mode)

    emis_flux_list(4)%name          = "seasalt_mass_as"
    emis_flux_list(4)%density       = 2.170_dp
    emis_flux_list(4)%diameter      = 2._dp * 0.258E-6_dp
    !emis_flux_list(4)%diameter      = 2._dp * 0.156e-6_dp !mz_vk_20170317
    emis_flux_list(4)%unit          = "kg/(m^2 s)"
    emis_flux_list(4)%mode          = det_mode(nmod,"as","AS")

    emis_flux_list(5)%name          = "seasalt_mass_cs"
    emis_flux_list(5)%density       = 2.170_dp
!    emis_flux_list(5)%diameter      = 2._dp * 1.65e-6_dp suitable for CB setup
!    emis_flux_list(5)%diameter      = 2._dp * 0.75e-6_dp
    emis_flux_list(5)%diameter      = 2._dp * 1.65e-6_dp
    !emis_flux_list(5)%diameter      = 2._dp * 0.85e-6_dp !mz_vk_20170317
    emis_flux_list(5)%unit          = "kg/(m^2 s)"
    emis_flux_list(5)%mode          = det_mode(nmod,"cs","CS")

    emis_flux_list(6)%name          = "seasalt_mass_cs_as"
    emis_flux_list(6)%density       = 2.170_dp
    emis_flux_list(6)%diameter      = 2._dp * 0.45e-6_dp
    emis_flux_list(6)%unit          = "kg/(m^2 s)"
    emis_flux_list(6)%mode          = det_mode(nmod,"as","AS")

! organic carbon
    emis_flux_list(10)%name         = "oc_mass_soa_ks"
    emis_flux_list(10)%density      = 2.0_dp
    emis_flux_list(10)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(10)%unit         = "kg/(m^2 s)"
!  soa should condense on existing particles
    emis_flux_list(10)%fac_num_emis = 0.1_dp
!    emis_flux_list(10)%fac_num_emis = 1._dp
    emis_flux_list(10)%mode         = det_mode(nmod,"ks","KS")

    emis_flux_list(11)%name         = "oc_mass_ff_ks"
    emis_flux_list(11)%density      = 2.0_dp
    emis_flux_list(11)%diameter     = 2._dp * 0.03e-6_dp ! aerocom
!    emis_flux_list(11)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(11)%unit         = "kg/(m^3 s)"
    emis_flux_list(11)%mode         = det_mode(nmod,"ks","KS")

    emis_flux_list(12)%name         = "oc_mass_bb_ks"
    emis_flux_list(12)%density      = 2.0_dp
!    emis_flux_list(12)%diameter     = 2._dp * 0.035e-6_dp !test_ht_121205
    emis_flux_list(12)%diameter     = 2._dp * 0.075e-6_dp ! aerocom
!    emis_flux_list(12)%diameter     = 2._dp * 0.09e-6_dp 
    emis_flux_list(12)%unit         = "kg/(m^3 s)"
    emis_flux_list(12)%mode         = det_mode(nmod,"ks","KS")

    emis_flux_list(13)%name         = "oc_mass_soa_ki"
    emis_flux_list(13)%density      = 2.0_dp
    emis_flux_list(13)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(13)%unit         = "kg/(m^2 s)"
!    emis_flux_list(13)%fac_num_emis = 0._dp
    emis_flux_list(13)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(14)%name         = "oc_mass_ff_ki"
    emis_flux_list(14)%density      = 2.0_dp
    emis_flux_list(14)%diameter     = 2._dp * 0.03e-6_dp ! aerocom
!    emis_flux_list(14)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(14)%unit         = "kg/(m^3 s)"
    emis_flux_list(14)%mode         = det_mode(nmod,"ki","KI")
    
    emis_flux_list(15)%name         = "oc_mass_bb_ki"
    emis_flux_list(15)%density      = 2.0_dp
!    emis_flux_list(15)%diameter     = 2._dp * 0.035e-6_dp  ! test_ht_121205
    emis_flux_list(15)%diameter     = 2._dp * 0.075e-6_dp ! aerocom
!    emis_flux_list(15)%diameter     = 2._dp * 0.09e-6_dp
    emis_flux_list(15)%unit         = "kg/(m^3 s)"
    emis_flux_list(15)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(16)%name         = "oc_mass_ks"
    emis_flux_list(16)%density      = 2.0_dp
    emis_flux_list(16)%diameter     = 2._dp * 0.258E-7_dp
    emis_flux_list(16)%mode         = det_mode(nmod,"ks","KS")    
    emis_flux_list(16)%unit         = "kg/(m^2 s)"
    emis_flux_list(16)%fac_num_emis = 1. / (emis_flux_list(16)%density * 1.e3_dp &
                                     * sigma_exp_ln(emis_flux_list(16)%mode))

    emis_flux_list(18)%name         = "oc_mass_ki"
    emis_flux_list(18)%density      = 2.0_dp
    emis_flux_list(18)%diameter     = 2._dp * 0.258E-7_dp
    emis_flux_list(18)%mode         = det_mode(nmod,"ki","KI")
    emis_flux_list(18)%unit         = "kg/(m^2 s)"
    emis_flux_list(18)%fac_num_emis = 1. / (emis_flux_list(18)%density * 1.e3_dp &
                                    * sigma_exp_ln(emis_flux_list(18)%mode))
    
    emis_flux_list(19)%name         = "oc_mass_road_ki"
    emis_flux_list(19)%density      = 2.0_dp
    emis_flux_list(19)%diameter     = 2._dp * 0.03 ! aerocom
    emis_flux_list(19)%unit         = "kg/(m^2 s)"
    emis_flux_list(19)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(20)%name         = "oc_ss_mass_as"
    emis_flux_list(20)%density      = 2.0_dp
    emis_flux_list(20)%diameter     = 2._dp * 0.258E-6_dp
    emis_flux_list(20)%unit         = "kg/(m^2 s)"
    emis_flux_list(20)%fac_num_emis = 0._dp ! oc part of seasalt particles
!    emis_flux_list(20)%fac_num_emis = 1._dp  ! oc separated from sea salt particles
    emis_flux_list(20)%mode         = det_mode(nmod,"as","AS")


! black carbon / soot

    emis_flux_list(30)%name         = "bc_mass_ff_ki"
    emis_flux_list(30)%density      = 2.0_dp
    emis_flux_list(30)%diameter     = 2._dp * 0.03e-6_dp !aerocom
!    emis_flux_list(30)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(30)%unit         = "kg/(m^3 s)"
    emis_flux_list(30)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(31)%name         = "bc_mass_bb_ki"
    emis_flux_list(31)%density      = 2.0_dp
    emis_flux_list(31)%diameter     = 2._dp * 0.075e-6_dp ! aerocom
!    emis_flux_list(31)%diameter     = 2._dp * 0.09e-6_dp ! aerocom
    emis_flux_list(31)%unit         = "kg/(m^3 s)"
    emis_flux_list(31)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(32)%name         = "bc_mass_ki"
    emis_flux_list(32)%density      = 2.0_dp
    emis_flux_list(32)%diameter     = 2._dp * 0.258E-7_dp
    emis_flux_list(32)%mode         = det_mode(nmod,"ki","KI")
    emis_flux_list(32)%unit         = "kg/(m^2 s)"
    emis_flux_list(32)%fac_num_emis = 1. / (emis_flux_list(32)%density * 1.e3_dp &
                                    * sigma_exp_ln(emis_flux_list(32)%mode))

    emis_flux_list(33)%name         = "bc_mass_air_ki"
    emis_flux_list(33)%density      = 2.0_dp
    emis_flux_list(33)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(33)%unit         = "kg/(m^3 s)"
    emis_flux_list(33)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(34)%name         = "bc_mass_road_ki"
    emis_flux_list(34)%density      = 2.0_dp
    emis_flux_list(34)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(34)%mode         = det_mode(nmod,"ki","KI")
    emis_flux_list(34)%unit         = "kg/(m^2 s)"

!mz_ap_20190507+
! WARNING: 
! for emissions imported and vertically distributed by OFFEMIS
! the [X]/(m^3s) units should be used
! for emissions imported by IMPORT and are 3D 
! the [X]/(m^3s) units should be used
! for emissions on a surface (example: BIOBURN)
! the [X]/(m^2s) units should be used

    emis_flux_list(35)%name         = "bc_mlc_2d_ki"
    emis_flux_list(35)%density      = 2.0_dp
    emis_flux_list(35)%diameter     = 2._dp * 0.03e-6_dp !aerocom
!    emis_flux_list(35)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(35)%unit         = "molecules/(m^2 s)"
    emis_flux_list(35)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(36)%name         = "bc_mlc_3d_ki"
    emis_flux_list(36)%density      = 2.0_dp
    emis_flux_list(36)%diameter     = 2._dp * 0.03e-6_dp !aerocom
!    emis_flux_list(36)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(36)%unit         = "molecules/(m^3 s)"
    emis_flux_list(36)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(37)%name         = "bc_mass_2d_ki"
    emis_flux_list(37)%density      = 2.0_dp
    emis_flux_list(37)%diameter     = 2._dp * 0.03e-6_dp !aerocom
!    emis_flux_list(37)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(37)%unit         = "kg/(m^2 s)"
    emis_flux_list(37)%mode         = det_mode(nmod,"ki","KI")

    emis_flux_list(38)%name         = "bc_mass_3d_ki"
    emis_flux_list(38)%density      = 2.0_dp
    emis_flux_list(38)%diameter     = 2._dp * 0.03e-6_dp !aerocom
!    emis_flux_list(38)%diameter     = 2._dp * 0.04e-6_dp
    emis_flux_list(38)%unit         = "kg/(m^3 s)"
    emis_flux_list(38)%mode         = det_mode(nmod,"ki","KI")

!mz_ap_20190507-

! dust
    
    emis_flux_list(40)%name         = "dust_mass_ci"
    emis_flux_list(40)%density      = 2.650_dp
!    emis_flux_list(40)%diameter     = 2.0_dp * 1.35e-6_dp !suitable for CB setup
!    emis_flux_list(40)%diameter     = 2.0_dp * 1.e-6_dp
    emis_flux_list(40)%diameter     = 2.0_dp * 0.65e-6_dp
    emis_flux_list(40)%unit         = "kg/(m^2 s)"
    emis_flux_list(40)%mode         = det_mode(nmod,"ci","CI")

    emis_flux_list(41)%name         = "dust_mass_ai"
    emis_flux_list(41)%density      = 2.650_dp
    emis_flux_list(41)%diameter     = 2.0_dp * 0.16e-6_dp
!    emis_flux_list(41)%diameter     = 2.0_dp * 0.35e-6_dp
    emis_flux_list(41)%unit         = "kg/(m^2 s)"
    emis_flux_list(41)%mode         = det_mode(nmod,"ai","AI")

!mz_ap_20171019+
    emis_flux_list(42)%name         = "du_nap_mass_ci"
    emis_flux_list(42)%density      = 2.650_dp
    emis_flux_list(42)%diameter     = 2.0_dp * 1.35e-6_dp
!!    emis_flux_list(42)%diameter     = 2.0_dp * 0.65e-6_dp
    emis_flux_list(42)%mode         = det_mode(nmod,"ci","CI") 
    emis_flux_list(42)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(43)%name         = "du_nap_mass_ai"
    emis_flux_list(43)%density      = 2.650_dp
    emis_flux_list(43)%diameter     = 2.0_dp * 0.16e-6_dp
    emis_flux_list(43)%mode         = det_mode(nmod,"ai","AI") 
    emis_flux_list(43)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(44)%name         = "du_kp_mass_ci"
    emis_flux_list(44)%density      = 2.650_dp
    emis_flux_list(44)%diameter     = 2.0_dp * 1.35e-6_dp
!!     emis_flux_list(44)%diameter     = 2.0_dp * 0.65e-6_dp
    emis_flux_list(44)%mode         = det_mode(nmod,"ci","CI") 
    emis_flux_list(44)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(45)%name         = "du_kp_mass_ai"
    emis_flux_list(45)%density      = 2.650_dp
    emis_flux_list(45)%diameter     = 2.0_dp * 0.16e-6_dp
    emis_flux_list(45)%mode         = det_mode(nmod,"ai","AI") 
    emis_flux_list(45)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(46)%name         = "du_capp_mass_ci"
    emis_flux_list(46)%density      = 2.650_dp
    emis_flux_list(46)%diameter     = 2.0_dp * 1.35e-6_dp
!!    emis_flux_list(46)%diameter     = 2.0_dp * 0.65e-6_dp
    emis_flux_list(46)%mode         = det_mode(nmod,"ci","CI") 
    emis_flux_list(46)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(47)%name         = "du_capp_mass_ai"
    emis_flux_list(47)%density      = 2.650_dp
    emis_flux_list(47)%diameter     = 2.0_dp * 0.16e-6_dp
    emis_flux_list(47)%mode         = det_mode(nmod,"ai","AI") 
    emis_flux_list(47)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(48)%name         = "du_mgpp_mass_ci"
    emis_flux_list(48)%density      = 2.650_dp
    emis_flux_list(48)%diameter     = 2.0_dp * 1.35e-6_dp
!!    emis_flux_list(48)%diameter     = 2.0_dp * 0.65e-6_dp
    emis_flux_list(48)%mode         = det_mode(nmod,"ci","CI") 
    emis_flux_list(48)%unit         = "kg/(m^2 s)"
 
    emis_flux_list(49)%name         = "du_mgpp_mass_ai"
    emis_flux_list(49)%density      = 2.650_dp
    emis_flux_list(49)%diameter     = 2.0_dp * 0.16e-6_dp
    emis_flux_list(49)%mode         = det_mode(nmod,"ai","AI") 
    emis_flux_list(49)%unit         = "kg/(m^2 s)"
!mz_ap_20171019-

! SO2 -> SO4

    emis_flux_list(50)%name         = "so2_mass_ks"
    emis_flux_list(50)%density      = species(spec_idx_SO4mm)%density
    emis_flux_list(50)%diameter     = 2.0_dp * 3.0e-8_dp
    emis_flux_list(50)%unit         = "molecules/(m^3 s)"
    emis_flux_list(50)%mode         = det_mode(nmod,"ks","KS")

    emis_flux_list(51)%name         = "so2_mass_as"
    emis_flux_list(51)%density      = species(spec_idx_SO4mm)%density
    emis_flux_list(51)%diameter     = 2.0_dp * 3.0e-7_dp
    emis_flux_list(51)%unit         = "molecules/(m^3 s)"
    emis_flux_list(51)%mode         = det_mode(nmod,"as","AS")

    emis_flux_list(52)%name         = "so2_mass_2d_ks"
    emis_flux_list(52)%density      = species(spec_idx_SO4mm)%density
    emis_flux_list(52)%diameter     = 2.0_dp * 3.0e-8_dp
    emis_flux_list(52)%unit         = "molecules/(m^2 s)"
    emis_flux_list(52)%mode         = det_mode(nmod,"ks","KS")

    emis_flux_list(53)%name         = "so2_mass_2d_as"
    emis_flux_list(53)%density      = species(spec_idx_SO4mm)%density
    emis_flux_list(53)%diameter     = 2.0_dp * 3.0e-7_dp
    emis_flux_list(53)%unit         = "molecules/(m^2 s)"
    emis_flux_list(53)%mode         = det_mode(nmod,"as","AS")

    num_fluxes = 0
!    DO jm=1,50 
    DO jm=1,100 ! mz_dk_20120203
      found = .false.
      IF ( ADJUSTL(TRIM(EMIS_CASK(jm,1))) == "") CYCLE
      DO jt=1,100
        IF (TRIM(emis_flux_list(jt)%name) == TRIM(EMIS_CASK(jm,1)) ) THEN
          num_fluxes = num_fluxes + 1
          found = .TRUE.
        ENDIF
      END DO
      IF (.NOT. FOUND) CALL warning_bi("EMIS_CASK named "//&
        TRIM(EMIS_CASK(jm,1))//&
        " not found in the list of fluxes and is therefore ignored!", substr)
    END DO

    counter = 0
    ALLOCATE(emis_flux_array(num_fluxes))
!    DO jm=1,50 
    DO jm=1,100 ! mz_dk_20120203
      IF ( TRIM(EMIS_CASK(jm,1)) == "") CYCLE
      DO jt=1,100
        IF (TRIM(emis_flux_list(jt)%name) == TRIM(EMIS_CASK(jm,1)) ) THEN
          counter = counter + 1
          emis_flux_array(counter)%name         = emis_flux_list(jt)%name
          emis_flux_array(counter)%channel_name = TRIM(EMIS_CASK(jm,3))
          emis_flux_array(counter)%flux_name    = TRIM(EMIS_CASK(jm,4))
          emis_flux_array(counter)%nflux_name   = TRIM(EMIS_CASK(jm,5))
          emis_flux_array(counter)%density      = emis_flux_list(jt)%density
          emis_flux_array(counter)%diameter     = emis_flux_list(jt)%diameter
          emis_flux_array(counter)%mode         = emis_flux_list(jt)%mode
          emis_flux_array(counter)%unit         = emis_flux_list(jt)%unit
          emis_flux_array(counter)%fac_num_emis = emis_flux_list(jt)%fac_num_emis
          IF ( TRIM(EMIS_CASK(jm,2)) == "") THEN
            emis_flux_array(counter)%total_frac = emis_flux_list(jt)%total_frac
          ELSE
            call str2num(TRIM(EMIS_CASK(jm,2)), val)
!            print*, "str2num 1", jm, val, TRIM(EMIS_CASK(jm,2))
            emis_flux_array(counter)%total_frac = val
          END IF

          call strcrack(EMIS_CASK(jm,6), ';', outstring,  &
            emis_flux_array(counter)%num_spec_emis)
          call strcrack(EMIS_CASK(jm,7), ';', outstring2, dummy)
    

          ALLOCATE(emis_flux_array(counter)%specs(emis_flux_array(counter)%num_spec_emis))
          
          DO jc = 1, emis_flux_array(counter)%num_spec_emis
            val = 0._dp
            emis_flux_array(counter)%specs(jc)%name = TRIM(outstring(jc))
            call str2num(TRIM(outstring2(jc)), val)
!            print*, "str2num 2", jm, val, TRIM(outstring2(jc))
            emis_flux_array(counter)%specs(jc)%frac = val
          ENDDO

        END IF
      END DO
    END DO

    DO jm = 1,num_fluxes
      DO jc = 1,emis_flux_array(jm)%num_spec_emis
         emis_flux_array(jm)%specs(jc)%trac_idx = 0
         emis_flux_array(jm)%specs(jc)%molarmass = 1._dp
         emis_flux_array(jm)%specs(jc)%mode = 0
         emis_flux_array(jm)%specs(jc)%l_numb = .FALSE.
        DO jt = 1,ntrac
          IF (TRIM(emis_flux_array(jm)%specs(jc)%name) == &
            ti_gp(jt)%tp%ident%fullname ) THEN
            emis_flux_array(jm)%specs(jc)%trac_idx = jt

            emis_flux_array(jm)%specs(jc)%molarmass = &
              ti_gp(jt)%tp%meta%cask_r(R_molarmass)
            emis_flux_array(jm)%specs(jc)%mode = &
              ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE)
            IF (ti_gp(jt)%tp%ident%quantity == numberdensity) &
              emis_flux_array(jm)%specs(jc)%l_numb = .TRUE.
          END IF
        END DO
      END DO
      IF ((emis_flux_array(jm)%mode == 0)    .AND. &
          (emis_flux_array(jm)%num_spec_emis > 0)) &
          emis_flux_array(jm)%mode = emis_flux_array(jm)%specs(1)%mode
    END DO

    DO jm = 1,num_fluxes
      CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
        TRIM(emis_flux_array(jm)%flux_name), p3=emis_flux_array(jm)%flux)
      IF (status /= 0) &
        CALL error_bi(&
        'requested object name '//TRIM(emis_flux_array(jm)%flux_name)//&
        ' for element '//TRIM(emis_flux_array(jm)%name)//&
        ' not found in channel '//TRIM(emis_flux_array(jm)%channel_name), substr)

      CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
        TRIM(emis_flux_array(jm)%nflux_name), p3=emis_flux_array(jm)%nflux)

      IF (status /= 0) THEN
        CALL info_bi( &
        'requested object name for number flux for element '// &
        TRIM(emis_flux_array(jm)%name)//&
        ' not found; calculating number from mass!', substr)
      ENDIF

      CALL get_channel_object_info(status,                     &
        TRIM(emis_flux_array(jm)%channel_name),                &
        TRIM(emis_flux_array(jm)%flux_name), reprid = id_repr )

      emis_flux_array(jm)%dim = 0
      emis_flux_array(jm)%NxD = .FALSE.
      IF (id_repr == GP_2D_HORIZONTAL) emis_flux_array(jm)%dim      = 2
      IF (id_repr == GP_2D_HORIZONTAL) emis_flux_array(jm)%dim_orig = 2
      IF (id_repr == GP_3D_MID)        emis_flux_array(jm)%dim      = 3
      IF (id_repr == GP_3D_MID)        emis_flux_array(jm)%dim_orig = 3
      IF (id_repr == GP_3D_1LEV)       emis_flux_array(jm)%dim_orig = 3
      IF (id_repr == GP_3D_1LEV)       emis_flux_array(jm)%dim      = 2

      IF (emis_flux_array(jm)%dim_orig == 2) THEN
        emis_flux_array(jm)%flux_2D  => emis_flux_array(jm)%flux(:,:,1)
        IF (ASSOCIATED(emis_flux_array(jm)%nflux)) &
          emis_flux_array(jm)%nflux_2D => emis_flux_array(jm)%nflux(:,:,1)
      ELSE IF (emis_flux_array(jm)%dim_orig == 3) THEN
        IF (id_repr == GP_3D_1LEV) THEN
          emis_flux_array(jm)%dim = 2
          emis_flux_array(jm)%flux_2D  => emis_flux_array(jm)%flux(_RI_XYZ__(:,:,1))
          IF (ASSOCIATED(emis_flux_array(jm)%nflux)) &
            emis_flux_array(jm)%nflux_2D => emis_flux_array(jm)%nflux(_RI_XYZ__(:,:,1))
        ENDIF
      END IF
      
      ! NxD emissions 
      IF ( emis_flux_array(jm)%dim == 0 ) THEN
        CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
          TRIM(emis_flux_array(jm)%flux_name)//'_vind', &
          p3=emis_flux_array(jm)%vind)
         IF (status /= 0) &
           CALL error_bi(&
           'requested object index array for element '//TRIM(emis_flux_array(jm)%name)//&
           ' not found in channel '//TRIM(emis_flux_array(jm)%channel_name), substr)
        emis_flux_array(jm)%dim = 3
        emis_flux_array(jm)%NxD = .TRUE.
      END IF

    END DO


!!$    IF (p_parallel_io) THEN
!!$      DO jm=1,num_fluxes
!!$        print*, "flux: ", emis_flux_array(jm)%name, jm
!!$        print*, "mode: ", emis_flux_array(jm)%mode
!!$        print*, "density: ", emis_flux_array(jm)%density
!!$        print*, "sigma: ", sigma(emis_flux_array(jm)%mode), sigma_exp_ln(emis_flux_array(jm)%mode)
!!$      END DO
!!$    END IF
      


  END SUBROUTINE gmxe_emis_init_e5

!---------------------------------------------------------------------------------

  INTEGER FUNCTION det_mode(nmod,str1,str2)

    INTEGER, INTENT(IN) :: nmod
    CHARACTER(LEN=2), INTENT(IN) :: str1 
    CHARACTER(LEN=2), INTENT(IN) :: str2
    
    INTEGER :: jt

    DO jt=1,nmod
      IF ((TRIM(cmodes(jt)) == str1) .OR.  (TRIM(cmodes(jt)) == str2)) THEN
        det_mode = jt
        EXIT
      END IF
    END DO
  END FUNCTION det_mode

!---------------------------------------------------------------------------------
END MODULE messy_gmxe_si

