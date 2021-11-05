#include "messy_main_ppd_bi.inc"

  MODULE messy_oracle_si
  !
  ! DESCRIPTION
  ! -----------
  ! oracle INTERFACE LAYER FOR ECHAM5/MESSY
  !
  ! AUTHOR
  ! ------
  ! Alexandra Tsimpidi, MAX Planck Institute for Chemistry, Mainz, Germany
  ! Vlassis Karydis,    MAX Planck Institute for Chemistry, Mainz, Germany
  ! questions/suggestions: a.tsimpidi@mpic.de
  !
  ! LAST MODIFICATIONS - 
  !*****************************************************************************

  USE messy_oracle
  USE messy_main_tracer,        ONLY: t_ident
  USE messy_main_constants_mem, ONLY: DP,STRLEN_MEDIUM, STRLEN_ULONG 
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  ! POINTER TO STREAM ELEMENTS
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  ! USE SPECIAL NCREGRID EVENT TRIGGER
  USE messy_main_blather_bi,    ONLY: info_bi, warning_bi, error_bi
  IMPLICIT NONE
  PRIVATE
 
  ! SUBROUTINES
  PUBLIC :: oracle_initialize              ! initialization
  PUBLIC :: oracle_new_tracer              ! define tracers
  PUBLIC :: oracle_init_memory             ! allocate memory
  PUBLIC :: oracle_init_tracer             ! initialize tracers
  PUBLIC :: oracle_init_coupling           ! coupling to echam5
  PUBLIC :: oracle_vdiff                   ! distribute online emissions
  PUBLIC :: oracle_physc                   ! calls oracle core layer
  PUBLIC :: oracle_radiation               ! 
  PUBLIC :: oracle_free_memory             ! deallocate memory 

  INTRINSIC ABS, ASSOCIATED, ALLOCATED, MAX, TRIM

  INTEGER, PUBLIC,  SAVE :: npre !,nmode,tmode(3)
  INTEGER, PUBLIC,  SAVE :: trac_id_gp(1) = 0  !location of gas phase speies in tracer array
!idt

   INTEGER :: idt_LfPOG01  = 0 ! created here as non reactive (no reaction in MECCA!)
   INTEGER :: idt_LbbPOG01 = 0 ! created here as non reactive (no reaction in MECCA!)  
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_fPOA
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_fPOG
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_bbPOA
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_bbPOG
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_fSOAsv
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_fSOGsv
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_bbSOAsv
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_bbSOGsv
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_fSOAiv
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_fSOGiv
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: idt_bbSOAiv
   INTEGER, DIMENSION(:),   ALLOCATABLE :: idt_bbSOGiv
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: idt_SOAv
   INTEGER, DIMENSION(:,:),   ALLOCATABLE :: idt_SOGv
!carbon emission mass fluxes insoluble part in kg/kg
!  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_insol => NULL()
!carbon emission mass fluxes soluble part in kg/kg
!  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_sol => NULL()

! carbon emissions
  CHARACTER (LEN=STRLEN_MEDIUM) :: Cemis_channel  = ''
! organic carbon
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_soa_sol   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_ff_sol    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_bb_sol    = ''
   !--- CPL namelist fields:
!mz_ap_20170912+
  LOGICAL, SAVE:: l_tendency      = .FALSE.
!mz_ap_20170912-
 
 !SCALAR CHANNEL FOR GMXE COUPLING------------------------------------------------------------------------
  REAL(DP), POINTER :: nmode_out   => NULL()! number of modes in ORACLE
  REAL(DP), POINTER :: NfPOA_out   => NULL()  
  REAL(DP), POINTER :: NbbPOA_out  => NULL()  
  REAL(DP), POINTER :: NfSOAsv_out => NULL() 
  REAL(DP), POINTER :: NbbSOAsv_out=> NULL()    
  REAL(DP), POINTER :: NfSOAiv_out => NULL()  
  REAL(DP), POINTER :: NbbSOAiv_out=> NULL()   
  REAL(DP), POINTER :: NSOAv_out   => NULL()     
  REAL(dp), DIMENSION(:), POINTER :: tmode_out => NULL() ! type of modes in ORACLE 
  REAL(dp), DIMENSION(:), POINTER :: kMOM_out  => NULL() ! number of MOM species in each volatility bin
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
   INTEGER, SAVE                               :: num_fluxes,idt_N(7) !N_APT!

   REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dryradius   => NULL() !gmxelink
   REAL(dp), POINTER, DIMENSION(:) :: sigma   => NULL() !gmxelink : standard deviation of aerosol mode 
   CHARACTER(LEN=26), POINTER     :: MOM_n_01(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_02(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_03(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_04(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_05(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_06(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_07(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_08(:) => NULL() !names of MOM species - MOMlink
   CHARACTER(LEN=26), POINTER     :: MOM_n_09(:) => NULL() !names of MOM species - MOMlink
   INTEGER, PUBLIC                :: maxMOM=30   !max number of MOM species for each volatility bin- MOMlink
   INTEGER, PUBLIC                :: kMOM(9)     !number of MOM species in each volatility bin- MOMlink

   ! EMIS_CASK = Character array which contains for each of maximum 500 emission 
   !              fields the necessary information for emission assignment:
   !          1 = name of the emission object (for identification)
   !          2 = total scaling factor for incoming flux
   !          3 = channel / channel name of emission flux
   !          4 = channel object name of mass emission flux
   !          5 = corresponding name emission flux (if exists)
   !          6 = list of tracers which should get a value from this emission
   !              ";" separated list of tracers (fullname, CASE SENSITIVE)
   !          7 = list of scaling factors for each tracer
   !              ";" separated list of REAL values
   CHARACTER(LEN=STRLEN_ULONG), DIMENSION(50,7), SAVE :: EMIS_CASK


  NAMELIST /CPL/ Cemis_channel,  EMIS_CASK,  l_tendency

CONTAINS
!==============================================================================

  SUBROUTINE oracle_initialize
!
   USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi
   USE messy_main_tools,        ONLY: find_next_free_unit
   USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, p_pe
 
   IMPLICIT NONE
   ! LOCAL
   CHARACTER(LEN=*), PARAMETER     :: substr = 'oracle_initialize' ! name of subroutine
   INTEGER                         :: iou    ! I/O unit
   INTEGER                         :: status ! error status
   INTEGER                         :: j,jm,jc,k
  
   IF (p_parallel_io) &
   CALL start_message_bi(modstr,'INITIALIZATION',substr)
 
   !--- Read namelist and control variables:
  
   ! INITIALIZE MAIN-CTRL
   IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      ! *** CALL CORE ROUTINE:
      CALL oracle_read_nml_ctrl(status, iou)
      IF (status /= 0)  CALL error_bi('error oracle_read_nml_ctrl', substr)
   END IF

   !--- Read CPL namelist
   IF (p_parallel_io) THEN
      EMIS_CASK(:,:) =''
      iou = find_next_free_unit(100,200)
      CALL oracle_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi("error in coupling namelist", substr)
   END IF



   !--- Broadcast over processors:
   CALL p_bcast (NfPOA,       p_io)
   CALL p_bcast (NbbPOA,       p_io)
   CALL p_bcast (NfSOAsv,      p_io)
   CALL p_bcast (NbbSOAsv,      p_io)
   CALL p_bcast (NfSOAiv,      p_io)
   CALL p_bcast (NbbSOAiv,      p_io)
   CALL p_bcast (NSOAv,       p_io)
   CALL p_bcast (NSOAP,        p_io)
   CALL p_bcast (aermod,       p_io)
   CALL p_bcast (nmode,        p_io)
   CALL p_bcast (tmode,        p_io)
   CALL p_bcast (Cemis_channel,p_io)
   CALL p_bcast (l_tendency,   p_io)

   DO jm=1,50
    DO jc=1,7
     CALL p_bcast(EMIS_CASK(jm,jc),p_io)
    END DO
   END DO


!   call flush (200+p_pe)
    DO j=1,NSOAv
     DO k=1,MAX_SPECIES
      CALL p_bcast (SOGv_mw(j,k), p_io)
     END DO
    END DO

   DO j=1,NSOAP
    CALL p_bcast (mwsoap(j), p_io)
    CALL p_bcast (csat(j), p_io)
    CALL p_bcast (cstemp(j), p_io)
    CALL p_bcast (deltah(j), p_io)
    CALL p_bcast (flagsoap(j), p_io)
!   call flush (200+p_pe)

   END DO

    CALL p_bcast(SOGv01, p_io)
    CALL p_bcast(SOGv02, p_io)
    CALL p_bcast(SOGv03, p_io)
    CALL p_bcast(SOGv04, p_io)
    CALL p_bcast(SOGv05, p_io)
    CALL p_bcast(SOGv06, p_io)
    CALL p_bcast(SOGv07, p_io)
    CALL p_bcast(SOGv08, p_io)
    CALL p_bcast(SOGv09, p_io)

   !--- Initialize core:

!   CALL oracle_initialize_core
 
   IF (p_parallel_io) &
   CALL end_message_bi(modstr,'INITIALIZATION',substr)
  
  END SUBROUTINE oracle_initialize
 
!***************************************************************************

  SUBROUTINE oracle_read_nml_cpl(status, iou)
 
! oracle MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
! read namelist for 'coupling' to channel containing organic emissions
! Authors: Alexandra Tsimpidi, MPIC, 2013
!          Vlassis Karydis,    MPIC, 2013
 
! MESSy
   USE messy_main_tools,          ONLY: read_nml_open, read_nml_check &
                                       , read_nml_close
   USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
   USE messy_main_mpi_bi,        ONLY: p_parallel_io

  IMPLICIT NONE


! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    INTRINSIC TRIM

   IF(p_parallel_io) &
     CALL start_message_bi(modstr,'Reading coupling namelist',substr)
 
    status = 1
 
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)

    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

     IF(p_parallel_io) &
     CALL end_message_bi(modstr,'Reading coupling namelist',substr)

  END SUBROUTINE oracle_read_nml_cpl
 
!*****************************************************************************



!*****************************************************************************

  SUBROUTINE oracle_new_tracer

   USE messy_main_blather_bi,      ONLY: start_message_bi, end_message_bi
   USE messy_main_tracer_tools_bi, ONLY: tracer_halt 
   USE messy_main_tracer,          ONLY: new_tracer, get_tracer, set_tracer,   &
                                         AIR, ON, OFF, MODAL, AEROSOL,         &
                                         AMOUNTFRACTION, NUMBERDENSITY,        &
                                         I_ADVECT, I_CONVECT,                  &
                                         I_VDIFF,                              &
                                         I_DRYDEP, I_SEDI,                     &
                                         I_SCAV, I_MIX,                        &
                                         I_AEROSOL_METHOD, I_AEROSOL_MODE,     &
                                         I_AEROSOL_SOL, S_AEROSOL_MODEL,       &
                                         R_MOLARMASS, R_AEROSOL_DENSITY,       &
                                         R_pss  ,     R_dryreac_sf
   USE messy_main_tracer_mem_bi,   ONLY: ti_gp, GPTRSTR
   USE messy_main_mpi_bi,          ONLY: p_parallel_io!,p_io, p_bcast, p_pe
   USE MESSY_MAIN_TOOLS,           ONLY: strcrack

   IMPLICIT NONE
! LOCAL
   INTEGER :: i, j, k, status
 
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_new_tracer'
   CHARACTER(LEN=2)            :: str_soa 
   CHARACTER(LEN=2)            :: str_mom 
   CHARACTER(LEN=2)            :: str_mode(3)
   INTEGER                     :: ierr, imedium, idt 

     IF(p_parallel_io) &
   CALL start_message_bi(modstr, 'TRACER DEFINITION', substr)

   do i=1,NSOAv
     IF (i.eq.1) CALL strcrack(SOGv01, ';', MOM_n_01, kMOM(i))
     IF (i.eq.2) CALL strcrack(SOGv02, ';', MOM_n_02, kMOM(i))
     IF (i.eq.3) CALL strcrack(SOGv03, ';', MOM_n_03, kMOM(i))
     IF (i.eq.4) CALL strcrack(SOGv04, ';', MOM_n_04, kMOM(i))
     IF (i.eq.5) CALL strcrack(SOGv05, ';', MOM_n_05, kMOM(i))
     IF (i.eq.6) CALL strcrack(SOGv06, ';', MOM_n_06, kMOM(i))
     IF (i.eq.7) CALL strcrack(SOGv07, ';', MOM_n_07, kMOM(i))
     IF (i.eq.8) CALL strcrack(SOGv08, ';', MOM_n_08, kMOM(i))
     IF (i.eq.9) CALL strcrack(SOGv09, ';', MOM_n_09, kMOM(i))
   end do

   maxMOM=MAXVAL(kMOM(1:NSOAv))
!   CALL p_bcast (kMOM,   p_io)
!   CALL p_bcast (maxMOM, p_io)

   ALLOCATE (idt_fPOA(NfPOA,nmode))
   ALLOCATE (idt_bbPOA(NbbPOA,nmode))
   ALLOCATE (idt_fSOAsv(NfSOAsv,nmode))
   ALLOCATE (idt_bbSOAsv(NbbSOAsv,nmode))
   ALLOCATE (idt_fSOAiv(NfSOAiv,nmode))
   ALLOCATE (idt_bbSOAiv(NbbSOAiv,nmode))
   ALLOCATE (idt_fPOG(NfPOA))
   ALLOCATE (idt_bbPOG(NbbPOA))
   ALLOCATE (idt_fSOGsv(NfSOAsv))
   ALLOCATE (idt_bbSOGsv(NbbSOAsv))
   ALLOCATE (idt_fSOGiv(NfSOAiv))
   ALLOCATE (idt_bbSOGiv(NbbSOAiv))
   ALLOCATE (idt_SOAv(NSOAv,nmode,maxMOM))
   ALLOCATE (idt_SOGv(NSOAv,maxMOM))
!mz_ap_20150311+
! Add tracers non reactive  
      CALL new_tracer(status, GPTRSTR, "LfPOG01",  &
           modstr, quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AIR, idx = idt_LfPOG01)
      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, R_molarmass       , 250._dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, R_pss             , 1.0e5_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LfPOG01, R_dryreac_sf      , 0._dp)
      CALL tracer_halt(substr, status)

      CALL new_tracer(status, GPTRSTR, "LbbPOG01",  &
           modstr, quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AIR, idx = idt_LbbPOG01)
      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_molarmass       , 250._dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_pss             , 1.0e5_dp)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR,  idt_LbbPOG01, R_dryreac_sf      , 0._dp)
      CALL tracer_halt(substr, status)
!mz_ap_20150311-

  str_mode(1)= 'ks'
  str_mode(2)= 'as'
  str_mode(3)= 'cs'
   DO j=1,nmode
  npre = 0
! Add aerosols for fPOA species 
    DO i = 1,NfPOA
    idt_fPOA(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "fPOA"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_fPOA(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), R_molarmass       , mwsoap(i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fPOA(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
     npre=NfPOA
! Add aerosols for bbPOA species 
    DO i = 1,NbbPOA
    idt_bbPOA(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "bbPOA"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_bbPOA(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbPOA(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
    npre=npre+NbbPOA
! Add aerosols for fSOAsv species 
    DO i = 1,NfSOAsv
    idt_fSOAsv(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "fSOAsv"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_fSOAsv(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAsv(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
    npre=npre+NfSOAsv
! Add aerosols for bbSOAsv species 
    DO i = 1,NbbSOAsv
    idt_bbSOAsv(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "bbSOAsv"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_bbSOAsv(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAsv(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
    npre=npre+NbbSOAsv
! Add aerosols for fSOAiv species 
    DO i = 1,NfSOAiv
    idt_fSOAiv(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "fSOAiv"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_fSOAiv(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_fSOAiv(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
      npre=npre+NfSOAiv
! Add aerosols for bbSOAiv species 
    DO i = 1,NbbSOAiv
    idt_bbSOAiv(i,j)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF
      CALL new_tracer(status, GPTRSTR, "bbSOAiv"//str_soa,  &
           modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
           unit = 'mol/mol', medium = AEROSOL, idx = idt_bbSOAiv(i,j))

      CALL tracer_halt(substr,status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_advect          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_convect         , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_vdiff           , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_scav            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_drydep          , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_sedi            , ON)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_mix             , OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), R_molarmass       , mwsoap(npre+i))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_aerosol_mode    , tmode(j))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), S_aerosol_model   , TRIM(aermod))
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), R_aerosol_density , 1.0E+03_dp )
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_aerosol_method  , MODAL)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_bbSOAiv(i,j), I_aerosol_sol     , 1)
      CALL tracer_halt(substr, status)
    END DO
    npre=npre+NbbSOAiv 
! Add aerosols for SOAv species
    DO i = 1,NSOAv
      DO k= 1,kMOM(i)
      idt_SOAv(i,j,k)= 0
       IF (i < 10) then
        write (str_soa,'(A,I1)') "0",i
       ELSE
        write (str_soa,'(I2)') i
       END IF

       IF (k < 10) then
        write (str_mom,'(A,I1)') "0",k
       ELSE
        write (str_mom,'(I2)') k
       END IF


       CALL new_tracer(status, GPTRSTR, "SOAv"//str_soa//str_mom,  &
            modstr, subname=str_mode(tmode(j)-1),quantity = AMOUNTFRACTION,&
            unit = 'mol/mol', medium = AEROSOL, idx = idt_SOAv(i,j,k))
       CALL tracer_halt(substr,status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_advect          , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_convect         , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_vdiff           , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_scav            , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_drydep          , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_sedi            , ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_mix             , OFF)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), R_molarmass       , SOGv_mw(i,k))
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_aerosol_mode    , tmode(j))
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), S_aerosol_model   , TRIM(aermod))
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), R_aerosol_density , 1.0E+03_dp )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_aerosol_method  , MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_SOAv(i,j,k), I_aerosol_sol     , 1)
       CALL tracer_halt(substr, status)

      END DO
    END DO
   END DO
     IF(p_parallel_io) &
    CALL end_message_bi(modstr, 'TRACER DEFINITION', substr)

  END SUBROUTINE oracle_new_tracer
!***************************************************************************

  SUBROUTINE oracle_init_memory

!   ! oracle MODULE ROUTINE (ECHAM-5 INTERFACE)
!   !
!   ! define oracle specific channel(s) and allocate memory for
!   ! global fields
!   !
!   ! Authors: Alexandra Tsimpidi, MPIC, 2013
!   !          Vlassis Karydis,    MPIC, 2013

    USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp     
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, DC_BC
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_repr,     ONLY: new_representation, AUTO &
                                         , set_representation_decomp &
                                         , IRANK, PIOTYPE_COL
    USE messy_main_channel_dimensions, ONLY: new_dimension
! 
!   ! ECHAM5/MESSy
   USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
   USE messy_main_mpi_bi,     ONLY: p_parallel_io
!
   IMPLICIT NONE

   INTEGER :: jt,i,k

!   ! LOCAL
   CHARACTER(LEN=2)            :: str_soa 
   CHARACTER(LEN=2)            :: str_mom 
   CHARACTER(LEN=1)            :: str_mod 
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_init_memory'
   CHARACTER(LEN=26)           :: MOM_name
   INTEGER                     :: status ! error status
!  !new representation for channel output
   INTEGER                               :: DIMID_TMODE
   INTEGER                               :: DIMID_KMOM
   INTEGER                               :: REPR_ORACLE_1D
   INTEGER                               :: REPR_ORACLE_KMOM_1D
   ! PARALLEL DECOMPOSITION
   INTEGER                          :: nseg = 0
   INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
   INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
   INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
   INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

     IF(p_parallel_io) &
   CALL start_message_bi(modstr, 'Get gasses from MECCA', substr)
!

!!!!!!!!!!!!!!!!!!!!!!!!!N_APT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO jt = 1,ntrac_gp

       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "ns") idt_N(1) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "ks") idt_N(2) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "as") idt_N(3) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "cs") idt_N(4) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "ki") idt_N(5) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "ai") idt_N(6) = jt
       IF(ti_gp(jt)%tp%ident%basename == "N".AND.ti_gp(jt)%tp%ident%subname == "ci") idt_N(7) = jt
    END DO

    DO i=1,7
      IF  (idt_N(i) == 0) THEN
       write (str_mod,'(I1)') i
       CALL error_bi("N_"//str_mod//" not present ", substr)
      ENDIF
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!N_APT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CHECK FOR CGs from MECCA

    DO i = 1,NfPOA
    idt_fPOG (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LfPOG"//str_soa) THEN
       idt_fPOG(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "fPOG"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF  (idt_fPOG(i) == 0) THEN
    CALL error_bi('fPOG not present in chemical mechanism!', substr)
   ENDIF
   END DO 

    DO i = 1,NbbPOA
    idt_bbPOG (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LbbPOG"//str_soa) THEN
       idt_bbPOG(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "bbPOG"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_bbPOG(i) == 0) THEN
    CALL error_bi('bbPOG not present in chemical mechanism!', substr)
   ENDIF
   END DO 
 
    DO i = 1,NfSOAsv
    idt_fSOGsv (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LfSOGsv"//str_soa) THEN
       idt_fSOGsv(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "fSOGsv"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_fSOGsv(i) == 0) THEN
    CALL error_bi('fSOGsv not present in chemical mechanism!', substr)
   ENDIF
   END DO 
    
    DO i = 1,NbbSOAsv
    idt_bbSOGsv (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LbbSOGsv"//str_soa) THEN
       idt_bbSOGsv(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "bbSOGsv"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_bbSOGsv(i) == 0) THEN
    CALL error_bi('bbSOGsv not present in chemical mechanism!', substr)
   ENDIF
   END DO 
    
    DO i = 1,NfSOAiv
    idt_fSOGiv (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LfSOGiv"//str_soa) THEN
       idt_fSOGiv(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "fSOGiv"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_fSOGiv(i) == 0) THEN
    CALL error_bi('fSOGiv not present in chemical mechanism!', substr)
   ENDIF
   END DO 
    
    DO i = 1,NbbSOAiv
    idt_bbSOGiv (i)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

    DO jt = 1,ntrac_gp
     IF(ti_gp(jt)%tp%ident%basename == "LbbSOGiv"//str_soa) THEN
       idt_bbSOGiv(i) = jt
       IF (p_parallel_io) THEN
           WRITE(*,*) "bbSOGiv"//str_soa, 'exist! '
       END IF
     ENDIF
    END DO
   IF (idt_bbSOGiv(i) == 0) THEN
    CALL error_bi('bbSOGiv not present in chemical mechanism!', substr)
   ENDIF
   END DO
 



    DO i = 1,NSOAv
     DO k=1,kMOM(i)
      idt_SOGv (i,k)= 0
      IF (i < 10) then
       write (str_soa,'(A,I1)') "0",i
      ELSE
       write (str_soa,'(I2)') i
      END IF

      IF (k < 10) then
       write (str_mom,'(A,I1)') "0",k
      ELSE
       write (str_mom,'(I2)') k
      END IF

      IF (i.eq.1) MOM_name = MOM_n_01(k)
      IF (i.eq.2) MOM_name = MOM_n_02(k)
      IF (i.eq.3) MOM_name = MOM_n_03(k)
      IF (i.eq.4) MOM_name = MOM_n_04(k)
      IF (i.eq.5) MOM_name = MOM_n_05(k)
      IF (i.eq.6) MOM_name = MOM_n_06(k)
      IF (i.eq.7) MOM_name = MOM_n_07(k)
      IF (i.eq.8) MOM_name = MOM_n_08(k)
      IF (i.eq.9) MOM_name = MOM_n_09(k)
   
      DO jt = 1,ntrac_gp
       IF(ti_gp(jt)%tp%ident%basename == MOM_name) THEN
!       IF(ti_gp(jt)%tp%ident%basename == "LSOGv"//str_soa) THEN
         idt_SOGv(i,k) = jt
         IF (p_parallel_io) THEN
           WRITE(*,*) "SOGv"//str_soa//str_mom, 'is', MOM_name 
         END IF
       ENDIF
      END DO
      IF (idt_SOGv(i,k) == 0) THEN
       write(*,*) i,k,MOM_name
       CALL error_bi('SOGv not present in chemical mechanism!', substr)
      ENDIF
     END DO
    END DO

    !------------------------------------------------------------- 
    ! ORACLE CHANNEL NEEDED FOR GMXE COUPLING
    !------------------------------------------------------------- 
    IF (p_parallel_io) &
    WRITE(*,*) 'add new channel ', TRIM(modstr),' ...'

    CALL new_channel(status, modstr, reprid=SCALAR, lrestreq =.FALSE.)
    CALL channel_halt(substr, status)

    ! nmode
    CALL new_channel_object(status, modstr, 'nmode' &
         , p0 = nmode_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'nmode' &
         , 'long_name', c='number of modes used for organic aerosol')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'nmode' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    nmode_out=nmode
    ! NfPOA
    CALL new_channel_object(status, modstr, 'NfPOA' &
         , p0 = NfPOA_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfPOA' &
         , 'long_name', c='number of Primary OA species from Fossil Fuel emission')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfPOA' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NfPOA_out=NfPOA
    ! NbbPOA
    CALL new_channel_object(status, modstr, 'NbbPOA' &
         , p0 = NbbPOA_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbPOA' &
         , 'long_name', c='number of Primary OA species from Biomass Burning emission')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbPOA' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NbbPOA_out=NbbPOA
    ! NfSOAsv
    CALL new_channel_object(status, modstr, 'NfSOAsv' &
         , p0 = NfSOAsv_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfSOAsv' &
         , 'long_name', c='number of Secondary OA species from the oxidation of fPOA with csat <  1000')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfSOAsv' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NfSOAsv_out=NfSOAsv
    ! NbbSOAsv
    CALL new_channel_object(status, modstr, 'NbbSOAsv' &
         , p0 = NbbSOAsv_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbSOAsv' &
         , 'long_name', c='number of Secondary OA species from the oxidation of bbPOA with csat <  10')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbSOAsv' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NbbSOAsv_out=NbbSOAsv
    ! NfSOAiv
    CALL new_channel_object(status, modstr, 'NfSOAiv' &
         , p0 = NfSOAiv_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfSOAiv' &
         , 'long_name', c='number of Secondary OA species from the oxidation of fPOA with csat >= 100')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NfSOAiv' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NfSOAiv_out=NfSOAiv
    ! NbbSOAiv
    CALL new_channel_object(status, modstr, 'NbbSOAiv' &
         , p0 = NbbSOAiv_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbSOAiv' &
         , 'long_name', c='number of Secondary OA species from the oxidation of bbPOA with csat >= 1000')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NbbSOAiv' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NbbSOAiv_out=NbbSOAiv
    ! NSOAv
    CALL new_channel_object(status, modstr, 'NSOAv' &
         , p0 = NSOAv_out)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NSOAv' &
         , 'long_name', c='number of Secondary OA species from the oxidation of traditional biogenic VOC')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NSOAv' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    NSOAv_out=NSOAv

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_dimension(status, DIMID_TMODE, 'ORACLE_NMODE', nmode)
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_ORACLE_1D,   &
         'REPR_ORACLE_1D'                               &
         , rank = 1, link = 'x---', dctype = DC_BC  &
         , dimension_ids = (/ DIMID_TMODE /)        &
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
    cnt(:,1)   = nmode
    meml(:,1)  = 1
    memu(:,1)  = nmode
    
    CALL set_representation_decomp(status, REPR_ORACLE_1D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! tmode (1-dimensional!)
    CALL new_channel_object(status, modstr, 'tmode' &
         , p1 = tmode_out, reprid=REPR_ORACLE_1D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tmode' &
    , 'long_name', c='type of modes in ORACLE')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tmode' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    tmode_out(1:nmode)=tmode(1:nmode)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_dimension(status, DIMID_KMOM, 'ORACLE_NMODE_KMOM', size(kMOM))
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_ORACLE_KMOM_1D,   &
         'REPR_ORACLE_KMOM_1D'                               &
         , rank = 1, link = 'x---', dctype = DC_BC  &
         , dimension_ids = (/ DIMID_KMOM /)        &
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
    cnt(:,1)   = size(kMOM)
    meml(:,1)  = 1
    memu(:,1)  = size(kMOM)
    
    CALL set_representation_decomp(status, REPR_ORACLE_KMOM_1D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! k-MOM (1-dimensional!)
    CALL new_channel_object(status, modstr, 'kMOM' &
         , p1 = kMOM_out, reprid=REPR_ORACLE_KMOM_1D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'kMOM' &
         , 'long_name', c='number of MOM species in each volatility bin- MOMlink ')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'kMOM' &
         , 'units', c='-' )
    CALL channel_halt(substr, status)
    kMOM_out(1:SIZE(KMOM))=kMOM(1:SIZE(KMOM))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE oracle_init_memory
!***************************************************************************

!***************************************************************************
  SUBROUTINE oracle_init_coupling

     USE messy_main_channel,          ONLY: get_channel_object, get_channel_info
     USE messy_main_channel_error_bi, ONLY: channel_halt
     USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
     USE messy_main_mpi_bi,           ONLY: p_parallel_io
     USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, nproma,nlev

   IMPLICIT NONE
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_init_coupling'
   INTEGER  :: ierr
   INTEGER  :: status ! error status         !gmxelink


     IF(p_parallel_io) &
     CALL start_message_bi(modstr,'INIT COUPLING',substr)

!Link to other aerosol model
     if(TRIM(aermod).ne."oracle") then
        write(6,*) TRIM(aermod//'_gp')
         CALL get_channel_object(status, TRIM(aermod//'_gp'),&
              'dryradius', p4=dryradius)
         IF (status /= 0) & 
              CALL error_bi('dry radius channel object not found !',substr)
         CALL get_channel_object(status, TRIM(aermod//'_gp'),&
              'sigma', p1=sigma)
         IF (status /= 0) & 
              CALL error_bi('sigma channel object not found !',substr)
     else
! Oracle uses its own distribution
         ALLOCATE(dryradius(nproma,nlev,7,ngpblks))
         ALLOCATE(sigma(7))
         dryradius(1:nproma,1:nlev,1,1:ngpblks) = 0.158E-8_dp  ! nucleation [m]
         dryradius(1:nproma,1:nlev,2,1:ngpblks) = 0.258E-7_dp  ! aitken     [m]
         dryradius(1:nproma,1:nlev,3,1:ngpblks) = 0.258E-6_dp  ! accum      [m]
         dryradius(1:nproma,1:nlev,4,1:ngpblks) = 0.258E-5_dp  ! coarse     [m]
         dryradius(1:nproma,1:nlev,5,1:ngpblks) = 0.258E-7_dp  ! aitken     [m]
         dryradius(1:nproma,1:nlev,6,1:ngpblks) = 0.258E-6_dp  ! accum      [m]
         dryradius(1:nproma,1:nlev,7,1:ngpblks) = 0.258E-5_dp  ! coarse     [m]
         sigma(1:7) = (/ 1.59_dp,1.59_dp,1.59_dp,2.0_dp,1.59_dp,1.59_dp,2.0_dp/)
     endif
!link

! This part deals with the emission fluxes
   CALL oracle_emis_init_si

 
    IF(p_parallel_io) &
    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE oracle_init_coupling

!***************************************************************************

   SUBROUTINE oracle_emis_init_si


    USE MESSY_MAIN_TRACER,              ONLY: r_molarmass, i_aerosol_mode, &
                                              numberdensity
    USE MESSY_MAIN_TRACER_MEM_BI,       ONLY: ti_gp, ntrac => ntrac_gp
    USE messy_main_tools,               ONLY: strcrack, str2num

    USE messy_main_channel,             ONLY: get_channel_object,  &
                                              get_channel_object_info
    USE messy_main_channel_error_bi,    ONLY: channel_halt
    USE messy_main_channel_bi,          ONLY: GP_3D_MID, &
                                              GP_3D_1LEV, GP_2D_HORIZONTAL
    USE messy_main_mpi_bi,              ONLY: p_io, p_bcast, p_parallel_io

    IMPLICIT NONE
    INTEGER :: jm, jc, jt, dummy, counter, status, id_repr
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER     :: outstring2(:) => NULL()

    LOGICAL                                   :: found
    CHARACTER(LEN=*), PARAMETER               :: substr='oracle_emis_si'
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
! organic carbon
    emis_flux_list(1)%name         = "oc_mass_ff_ks"
    emis_flux_list(1)%density      = 2.0_dp
    emis_flux_list(1)%diameter     = 2._dp * 0.03e-6_dp
    emis_flux_list(1)%mode         = 2
    emis_flux_list(1)%unit         = "kg/(m^3 s)"
!    emis_flux_list(1)%unit         = "kg/(m^2 s)"


    emis_flux_list(2)%name         = "oc_mass_ff_as"
    emis_flux_list(2)%density      = 2.0_dp
    emis_flux_list(2)%diameter     = 2._dp * 0.3e-6_dp 
    emis_flux_list(2)%mode         = 3
    emis_flux_list(2)%unit         = "kg/(m^3 s)"
!    emis_flux_list(2)%unit         = "kg/(m^2 s)"

    emis_flux_list(3)%name         = "oc_mass_bb_ks"
    emis_flux_list(3)%density      = 2.0_dp
    emis_flux_list(3)%diameter     = 2._dp * 0.075e-6_dp
    emis_flux_list(3)%mode         = 2
    emis_flux_list(3)%unit         = "kg/(m^2 s)"
!    emis_flux_list(3)%unit         = "kg/(m^3 s)"

    emis_flux_list(4)%name         = "oc_mass_bb_as"
    emis_flux_list(4)%density      = 2.0_dp
    emis_flux_list(4)%diameter     = 2._dp * 0.75e-6_dp
    emis_flux_list(4)%mode         = 3
    emis_flux_list(4)%unit         = "kg/(m^2 s)"
!    emis_flux_list(4)%unit         = "kg/(m^3 s)"

    emis_flux_list(5)%name         = "oc_mass_ks"
    emis_flux_list(5)%density      = 2.0_dp
    emis_flux_list(5)%diameter     = 2._dp * 0.258E-7_dp
    emis_flux_list(5)%mode         = 2
    emis_flux_list(5)%unit         = "kg/(m^2 s)"

    emis_flux_list(6)%name         = "oc_mass_as"
    emis_flux_list(6)%density      = 2.0_dp
    emis_flux_list(6)%diameter     = 2._dp * 0.258E-6_dp
    emis_flux_list(6)%mode         = 3
    emis_flux_list(6)%unit         = "kg/(m^2 s)"

    emis_flux_list(7)%name         = "oc_ss_mass_as"
    emis_flux_list(7)%density      = 2.0_dp
    emis_flux_list(7)%diameter     = 2._dp * 0.258E-6_dp
    emis_flux_list(7)%mode         = 3
    emis_flux_list(7)%unit         = "kg/(m^2 s)"
 
    num_fluxes = 0
     DO jm=1,50
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
    DO jm=1,50
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
          emis_flux_array(counter)%total_frac   = emis_flux_list(jt)%total_frac
           ELSE
            call str2num(TRIM(EMIS_CASK(jm,2)), val)
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
    END DO

    DO jm = 1,num_fluxes
      CALL get_channel_object(status,TRIM(emis_flux_array(jm)%channel_name), &
        TRIM(emis_flux_array(jm)%flux_name), p3=emis_flux_array(jm)%flux)
      IF (status /= 0) &
       CALL error_bi(&
       'requested object name for element '//TRIM(emis_flux_array(jm)%name)//&
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
          emis_flux_array(jm)%flux_2D  =>  emis_flux_array(jm)%flux(_RI_XYZ__(:,:,1))
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
      IF (p_parallel_io) THEN
          write(*,*)  'flux',emis_flux_array(jm)%flux_name,' of dim ',emis_flux_array(jm)%dim_orig
          write(*,*)  'assigned to ', emis_flux_array(jm)%dim
          write(*,*)  'id_repr= ', id_repr, '"GP_3D_1LEV=',GP_3D_1LEV
      ENDIF


    END DO

  END SUBROUTINE oracle_emis_init_si


!***************************************************************************

  SUBROUTINE oracle_init_tracer
  
   IMPLICIT NONE

  END SUBROUTINE oracle_init_tracer

!***************************************************************************
   SUBROUTINE oracle_global_start
 
     IMPLICIT NONE
 
   END SUBROUTINE oracle_global_start
!***************************************************************************
 
  SUBROUTINE oracle_vdiff
! oracle MODULE ROUTINE (ECHAM-5 INTERFACE)
!
! distribute online organic emissions
!
! Authors: Alexandra Tsimpidi, MPIC, 2013
!          Vlassis Karydis,    MPIC, 2013
 
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, ti_gp
    USE messy_main_constants_mem, ONLY: g,M_air,STRLEN_MEDIUM, pi, avo => N_A
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, jrow, kproma 
    USE messy_main_grid_def_bi,     ONLY: philat_2d
! mz_rj_20140918
#if defined(ECHAM5)
    USE messy_main_data_bi,       ONLY: pxtems
#endif
    USE messy_main_data_bi,       ONLY: tm1_3d, tte_3d   &
                                      , qte_3d, qm1_3d, press_3d &
                                      , pressi_3d
    USE messy_main_grid_def_bi,   ONLY: grmass, grvol, deltaz
    USE messy_main_timer,         ONLY: time_step_len &
                                      , nstep=>current_time_step, lstart
    USE messy_main_mpi_bi,        ONLY: p_io, p_bcast, p_parallel_io
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi

   IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='oracle_vdiff'
    INTEGER                     :: jc,jt,jk,idx,jl, mlev, ji, idt, jm  !N_APT!
    REAL(dp)                    :: diameter3 !N_APT!
    REAL(dp)                    :: convMtoN !N_APT!
#if defined (ECHAM5)
    REAL(dp), POINTER           :: zxtems(:,:)
#endif
    REAL(dp)                    :: zdp    (nproma,nlev)
    REAL(dp)                    :: fac_N  (nproma,nlev) !N_APT!
    REAL(dp)                    :: conv_unit(nproma,nlev)
    REAL(dp)                    :: zdz    (nproma,nlev)
    REAL(dp)                    :: flux(nproma,nlev)
    REAL(dp)                    :: sigma_exp_ln 
 
 
!      IF(p_parallel_io) & !N_APT!
!      CALL start_message_bi(modstr,'EMIS COUPLING',substr)!N_APT!

! mz_rj_20140918+
#if defined (ECHAM5)
    zxtems => pxtems(:,1,:,jrow)
#endif
! mz_rj_20140918-
     DO jk=1,nlev
       zdp(1:kproma,jk) = pressi_3d(_RI_XYZ__(1:kproma,jrow,jk+1)) - &
                          pressi_3d(_RI_XYZ__(1:kproma,jrow,jk))
     END DO
     zdz(:,:) = 0.0_dp
     zdz(1:kproma,2:nlev) = deltaz(_RI_XYZ__(1:kproma,jrow,2:nlev))

      DO jm = 1,7
        jt = idt_N(jm)
        fac_N(1:kproma,1:nlev) = 1._dp
        SELECT CASE (TRIM(ti_gp(jt)%tp%ident%unit))
        CASE ('1/cm3')
          fac_N(1:kproma,1:nlev) = grmass(1:kproma,1:nlev,jrow) &
                                 / grvol(1:kproma,1:nlev,jrow) * 1.e-6_dp
        CASE ('1/mol')
          fac_N(1:kproma,1:nlev) = M_air * 0.001_dp
        CASE ('1/kg')
          fac_N(1:kproma,1:nlev) = 1._dp
        CASE DEFAULT 
          IF (p_parallel_io) &
            CALL error_bi('no valid tracer unit for aerosol number!',&
            substr)
        END SELECT
      END DO
!!!!!!!!!!!!!!!!!!!!!!!N_APT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 ! Loop over emission fluxes

    DO jc = 1,num_fluxes
      flux(:,:)  = 0._dp

      jm = emis_flux_array(jc)%mode !N_APT!
      convMtoN = 1._dp !N_APT!
      conv_unit(1:kproma,1:nlev) = 1._dp

!#ifndef VERTICO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!N_APT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! tracer numbers
      idx = idt_N(jm)
      conv_unit(1:kproma,1:nlev) = conv_unit(1:kproma,1:nlev) &
                                 * fac_N(1:kproma,1:nlev)
      IF (.NOT. ASSOCIATED(emis_flux_array(jc)%nflux) ) THEN
        diameter3 = emis_flux_array(jc)%diameter * emis_flux_array(jc)%diameter &
                  * emis_flux_array(jc)%diameter
        sigma_exp_ln = EXP(4.5*LOG(sigma(jm))*LOG(sigma(jm)))
        convMtoN  = 6._dp / pi &
                  / (emis_flux_array(jc)%density * 1.e3_dp) &
                  / (diameter3 * sigma_exp_ln)

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
            mlev = SIZE(emis_flux_array(jc)%VIND,2)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%flux(1:kproma,1:mlev,jrow)
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%flux_2D(1:kproma,jrow)
        ENDIF
      ELSE
        IF ( emis_flux_array(jc)%dim == 3) THEN
          IF  (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,2)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%nflux(1:kproma,1:mlev,jrow)
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%nflux_2D(1:kproma,jrow)
        ENDIF
      ENDIF

!mz_ap_20170912+
      IF ( (l_tendency) .OR. &
        (emis_flux_array(jc)%dim == 3) ) THEN
!mz_ap_20170912-
        IF (emis_flux_array(jc)%NxD) THEN
          DO ji=1,mlev
            DO jl=1,kproma
              jk= NINT(emis_flux_array(jc)%VIND(jl,ji,jrow))
              pxtte(_RI_X_ZN_(jl,jk,idx)) = pxtte(_RI_X_ZN_(jl,jk,idx))   &
                + conv_unit(jl,jk) * convMtoN            &
                * flux(jl,ji)                            &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
!            IF (p_parallel_io) write(6,*) 'N_dbg_NxD', lnumber,conv_unit(jl,jk), convMtoN, emis_flux_array(jc)%fac_num_emis, emis_flux_array(jc)%total_frac
!            IF (p_parallel_io) write(6,*) 'N_dbg_NxD', idx, pxtte(jl,jk,idx) 
             END DO
          END DO
        ELSE
          DO jk=1,nlev
            DO jl=1,kproma
              pxtte(_RI_X_ZN_(jl,jk,idx)) = pxtte(_RI_X_ZN_(jl,jk,idx))        &
                + conv_unit(jl,jk) * convMtoN            &
                * flux(jl,jk)                            &
                * emis_flux_array(jc)%fac_num_emis       &
                * emis_flux_array(jc)%total_frac         &
                / zdp(jl,jk) * g
!            IF (p_parallel_io) write(6,*) 'N_dbg_3D', lnumber,conv_unit(jl,jk), convMtoN, emis_flux_array(jc)%fac_num_emis, emis_flux_array(jc)%total_frac
!            IF (p_parallel_io) write(6,*) 'N_dbg_3D', idx, pxtte(jl,jk,idx) 
            END DO
          END DO
        END IF
#if defined (ECHAM5)
      ELSE
        DO jl=1,kproma
          zxtems(jl,idx) = zxtems(jl,idx)              &
            + conv_unit(jl,nlev) * convMtoN            &
            * flux(jl,nlev)                            &
            * emis_flux_array(jc)%fac_num_emis         &
            * emis_flux_array(jc)%total_frac
!            IF (p_parallel_io) write(6,*) 'N_dbg_2D', lnumber,conv_unit(jl,nlev), convMtoN, emis_flux_array(jc)%fac_num_emis, emis_flux_array(jc)%total_frac
!            IF (p_parallel_io) write(6,*) 'N_dbg_2D', idx, zxtems(jl,idx) 
        END DO
#endif
      ENDIF
!#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!N_APT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! tracer mass 
      convMtoN = 1._dp
      DO jt = 1, emis_flux_array(jc)%num_spec_emis
        flux(:,:)  = 0._dp
        conv_unit(:,:) = 1._dp
        idx = emis_flux_array(jc)%specs(jt)%trac_idx
        IF (idx == 0) CYCLE
        SELECT CASE (TRIM(emis_flux_array(jc)%unit))
        CASE("kg/(m^2 s)")
          conv_unit(1:kproma,1:nlev) = &
            M_air / emis_flux_array(jc)%specs(jt)%molarmass
!mz_ap_20140921+
         CASE("kg/(m^3 s)")
           conv_unit(1:kproma,2:nlev) = &
             M_air / emis_flux_array(jc)%specs(jt)%molarmass*zdz(1:kproma,2:nlev) 
!mz_ap_20140921-
        CASE("molecules/(m^2 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp
        CASE("molecules/(m^3 s)")
          conv_unit(1:kproma,2:nlev) = M_air / avo / 1.e3_dp * &
                                       zdz(1:kproma,2:nlev)
        END SELECT

        IF (emis_flux_array(jc)%dim == 3) THEN
          IF (emis_flux_array(jc)%NxD) THEN
            mlev = SIZE(emis_flux_array(jc)%VIND,2)
          ELSE
            mlev = nlev
          ENDIF
          flux(1:kproma,1:mlev) = emis_flux_array(jc)%flux(_RI_XYZ__(1:kproma,jrow,1:mlev))
        ELSE IF (emis_flux_array(jc)%dim == 2) THEN
          flux(1:kproma,nlev)   = emis_flux_array(jc)%flux_2D(1:kproma,jrow)
        ENDIF
 
          idt = idx
!mz_ap_20170912+
      IF ( (l_tendency) .OR. &
        (emis_flux_array(jc)%dim == 3) ) THEN
!mz_ap_20170912-
          IF (emis_flux_array(jc)%NxD) THEN
            DO ji=1,mlev
              DO jl=1,kproma
                jk= NINT(emis_flux_array(jc)%VIND(jl,ji,jrow))
                pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))        &
                  + conv_unit(jl,jk)                       &
                  * flux(jl,ji)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ELSE
! mz_rj_20140410+
            DO jk=2,nlev
! mz_rj_20140410-
              DO jl=1,kproma
                pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))        &
                  + conv_unit(jl,jk)                       &
                  * flux(jl,jk)                            &
                  * emis_flux_array(jc)%specs(jt)%frac     &
                  * emis_flux_array(jc)%total_frac         &
                  / zdp(jl,jk) * g
              END DO
            END DO
          ENDIF
       ELSE
#if defined (ECHAM5)
          DO jl=1,kproma
            zxtems(jl,idx) = zxtems(jl,idx)              &
              + conv_unit(jl,nlev) * convMtoN            &
              * flux(jl,nlev)                            &
              * emis_flux_array(jc)%specs(jt)%frac       &
              * emis_flux_array(jc)%total_frac
          END DO
#endif
       ENDIF
      END DO
    END DO

!       IF(p_parallel_io) & !N_APT!                                                                                       
!       CALL end_message_bi(modstr,'EMIS COUPLING',substr)!N_APT!

 END SUBROUTINE oracle_vdiff
!***************************************************************************

!***************************************************************************
  SUBROUTINE oracle_physc 
   USE messy_main_data_bi,        ONLY:           tm1       &    ! dry air temperature (at time step -1)  [K]
                                        ,         tte_3d    &    ! dry air temperature (tendency)         [K/s]
                                        ,         press_3d  &    ! air pressure                           [Pa]
                                        ,         pressi_3d     ! air pressure (interface)               [Pa]
   USE messy_main_grid_def_bi,     ONLY:grvol, grmass          ! grid volume [m3] and mass [kg]
   USE messy_main_grid_def_mem_bi, ONLY:klev => nlev       &
                                        ,nproma, kproma     &
                                        ,nrow => ngpblks    &
                                        ,jrow             
   USE messy_main_constants_mem,  ONLY:   M_air,dp
   USE messy_main_timer,          ONLY: delta_time,time_step_len &
                                      , nstep=>current_time_step
   USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi
   USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
   USE messy_main_mpi_bi,         ONLY: p_io, p_bcast, &
                                        p_parallel_io, p_pe

   IMPLICIT NONE
   CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_physc'
   real(dp) caer(NSOAP), cgas(NSOAP), csatT(NSOAP),qaer(NSOAP,nmode),qgas(NSOAP)
   real(dp) maer(NSOAv,nmode,maxMOM),mgas(NSOAv,maxMOM),mtot(NSOAv),mfrac(NSOAv,maxMOM)
   real(dp) tempk,xconvert(kproma,klev),idt_tes(4)    !gmxelink
   REAL(dp), DIMENSION(nmode) :: rsec     ! mode radius - gmxelink
   REAL(dp), DIMENSION(kproma,klev) :: ztemp,zpress
   REAL(dp), DIMENSION(kproma,klev,nmode) :: zdryrad              !gmxelink
   integer i,j,jk,jl,k 
   integer idt,nlev
   LOGICAL,  SAVE              :: entered = .FALSE.

      IF(p_parallel_io.AND. .NOT. entered) &
     CALL start_message_bi(modstr,'CALL oracle core layer',substr)

  nlev = klev ! mz_rj_20140909
  ztemp       (:,:) = 0.0_dp
  zpress      (:,:) = 0.0_dp

   !--- Assign conversion factor to convert mol mol-1 to umol m-3:
   xconvert(1:kproma,1:klev)=(1.E+9_dp/M_air)* &                        !  [mol mol-1] => [umol g-1]
          (grmass(_RI_XYZ__(1:kproma,jrow,1:nlev))/grvol(_RI_XYZ__(1:kproma,jrow,1:nlev))) !  Air density [g m-3]: [umol g-1] => [umol m-3]
          

   !--- Ambient properties: ----------------------------------
   !--- Temperature [K]:
 
   ztemp  (1:kproma,1:klev)       = tm1  (_RI_XYZ__(1:kproma,jrow,1:nlev))   +       &
                                    tte_3d(_RI_XYZ__(1:kproma,jrow,1:klev)) * time_step_len
   !--- Pressure [Pa]:
 
   zpress  (1:kproma,1:klev)      = press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))     ! [Pa]

   !--- Dry radius [m]:

!mz_ap_20170914+
! change momentum in log-normal distribution
! fomr number distribution to surface distribution:
! radius change accordingly : eq 8.50 Senfield and Pandis 2nd edition
! OLD:
!   zdryrad (1:kproma,1:klev,1:nmode) = dryradius(1:kproma,1:klev,tmode(1):tmode(nmode),jrow) ! [m] gmxelink
! NEW
! log = ln ;  log10 =log base 10
   do j=1,nmode
      !zdryrad (1:kproma,1:klev,j) = 0.5 *(exp(log(2*dryradius(1:kproma,1:klev,tmode(j),jrow))+ &
      !                                           2*((log(sigma(tmode(j))))**2))) ! [m] gmxelink
      ! this is equivalent to the above one, just nicer and faster 
      zdryrad (1:kproma,1:klev,j) = dryradius(1:kproma,1:klev,tmode(j),jrow)*exp(2*((log(sigma(tmode(j))))**2))
   enddo                                              

   do jk = 1,klev
    do jl = 1,kproma

     tempk = ztemp  (jl,jk)
     caer=0.0_dp
     qaer=0.0_dp
     maer=0.0_dp
     mtot=0.0_dp
     cgas=0.0_dp
     qgas=0.0_dp
     mgas=0.0_dp
     csatT=0.0_dp
     mfrac=0.0_dp
     do j=1,nmode
      rsec(j) = zdryrad(jl,jk,j)     !gmxelink
      npre=0
      do i = 1,NfPOA      
      idt = idt_fPOA(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3 
                 * mwsoap(i)                                      !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      
      idt = idt_fPOG(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3
                 * mwsoap(i)                                      !convert umol m-3 to ug m-3
      end do

      npre=npre+NfPOA 
      do i = 1,NbbPOA      
      idt = idt_bbPOA(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3 
                 * mwsoap(npre+i)                                 !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt = idt_bbPOG(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                              & ! conversion factor umol m-3
                 * mwsoap(npre+i)                                 !convert umol m-3 to ug m-3
      end do
      
      npre=npre+NbbPOA
      do i = 1,NfSOAsv      
      idt = idt_fSOAsv(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3 
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt = idt_fSOGsv(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3
      end do
     
      npre=npre+NfSOAsv
      do i = 1,NbbSOAsv      
      idt =idt_bbSOAsv(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3 
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt = idt_bbSOGsv(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3
      end do
      
      npre=npre+NbbSOAsv
      do i = 1,NfSOAiv      
      idt = idt_fSOAiv(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3 
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt = idt_fSOGiv(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3
      end do
     
      npre=npre+NfSOAiv
      do i = 1,NbbSOAiv      
      idt = idt_bbSOAiv(i,j)
      qaer(npre+i,j) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))    &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3 
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + qaer(npre+i,j)      

      idt =idt_bbSOGiv(i)
      if(j.eq.1) cgas(npre+i) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                 + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len ))  & ! tracer mass mol mol-1
                 * xconvert(jl,jk)                               & ! conversion factor umol m-3
                 * mwsoap(npre+i)                                  !convert umol m-3 to ug m-3
      end do
      
      npre=npre+NbbSOAiv
      do i = 1,NSOAv
       do k= 1,kMOM(i)      
      idt = idt_SOAv(i,j,k)
      maer(i,j,k) =  MAX(0._dp,(pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                  + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len))   & ! tracer mass mol mol-1
                  * xconvert(jl,jk)                                & ! conversion factor umol m-3
                  * SOGv_mw(i,k)                                   !convert umol m-3 to ug m-3

      caer(npre+i) = caer(npre+i) + maer(i,j,k)      
      qaer(npre+i,j) = qaer(npre+i,j) + maer(i,j,k)     
      mtot(i) = mtot(i) + maer(i,j,k)     

      idt = idt_SOGv(i,k)
      if(j.eq.1) then
        mgas(i,k) = MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))               &
                            + pxtte(_RI_X_ZN_(jl,jk,idt)) * time_step_len )) & ! tracer mass mol mol-1
                            * xconvert(jl,jk)                              & ! conversion factor umol m-3
                            * SOGv_mw(i,k)                                 !convert umol m-3 to ug m-3

        cgas(npre+i) = cgas(npre+i) + mgas(i,k)
        mtot(i) = mtot(i) + mgas(i,k)     
      end if

       end do
      end do
     end do

     do i = 1,NSOAv
       mwsoap(npre+i)=0.d0
        do k= 1,kMOM(i)  
     mfrac(i,k)= MAX(0._dp,((mgas(i,k)+sum(maer(i,:,k)))/max(1.e-50_dp,mtot(i))))
        mwsoap(npre+i)=mwsoap(npre+i)+mfrac(i,k)*SOGv_mw(i,k)
        end do
        mwsoap(npre+i)=MAX(1._dp,mwsoap(npre+i))
     end do

    do i=1,NSOAP
     qgas(i)=cgas(i)
    end do
!!!!!!!!!!!!!!!!!!!!!CALL gas/aerosol partitioning subroutine!!!!!!!!!!!!!!!!!!!!!!!!
     CALL oracle_soap(caer,cgas,tempk,p_pe,jl,jk,csatT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!!!!!!!!!!!!!CALL mode distribution subroutine!!!!!!!!!!!!!!!!!!!!!!!!
     if(nmode.ne.1) then
       CALL oracle_mode(caer,cgas,qaer,qgas,csatT,rsec,p_pe,jl,jk) !gmxelink
     else
       do i=1,NSOAP
        qaer(i,1)=caer(i)
       end do
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do j=1,nmode
      npre=0

      do i = 1,NfPOA      
       idt = idt_fPOA(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              &
                           + (qaer(i,j)/xconvert(jl,jk)/mwsoap(i)           &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_fPOG(i)                    
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))                & 
                           + (cgas(i)/xconvert(jl,jk)/mwsoap(i)           &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))       &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))    &
                           / time_step_len 
       end if
      end do

      npre=npre+NfPOA
      do i = 1,NbbPOA      
       idt = idt_bbPOA(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i) &
!                           + pxtte(jl,jk,idt_BBPOA(i,j))*time_step_len)))   &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_bbPOG(i)                    
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))                & 
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i) &
!                           + pxtte(jl,jk,idt_BBCG(i))*time_step_len)))    &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))       &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))    &
                           / time_step_len 
       end if
      end do

      npre=npre+NbbPOA
      do i = 1,NfSOAsv      
       idt = idt_fSOAsv(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))            &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i) &
!                           + pxtte(jl,jk,idt_FFsSOA(i,j))*time_step_len))) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len))) &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_fSOGsv(i)                    
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              & 
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i) &
!                           + pxtte(jl,jk,idt_FFsCG(i))*time_step_len)))   &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 
       end if
      end do

      npre=npre+NfSOAsv
      do i = 1,NbbSOAsv      
       idt = idt_bbSOAsv(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))            &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))  &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_bbSOGsv(i)                     
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              & 
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 
       end if
      end do

      npre=npre+NbbSOAsv
      do i = 1,NfSOAiv      
       idt = idt_fSOAiv(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))            &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_fSOGiv(i)                      
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              & 
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 
       end if
      end do

      npre=npre+NfSOAiv
      do i = 1,NbbSOAiv      
       idt = idt_bbSOAiv(i,j)
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))            &
                           + (qaer(npre+i,j)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))     &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))  &
                           / time_step_len 

       if(j.eq.1) then
        idt = idt_bbSOGiv(i)                    
        pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))              & 
                           + (cgas(npre+i)/xconvert(jl,jk)/mwsoap(npre+i) &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))      &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))   &
                           / time_step_len 
       end if
      end do

      npre=npre+NbbSOAiv
      do i = 1,NSOAv      
        do k=1,kMOM(i)
         idt = idt_SOAv(i,j,k)
         pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))                &
                           + (mfrac(i,k)*qaer(npre+i,j)/xconvert(jl,jk)/SOGv_mw(i,k)     &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))         &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))      &
                           / time_step_len 


      if(j.eq.1) then
       idt = idt_SOGv(i,k)                    
       pxtte(_RI_X_ZN_(jl,jk,idt)) = pxtte(_RI_X_ZN_(jl,jk,idt))                  &
                           + (mfrac(i,k)*cgas(npre+i)/xconvert(jl,jk)/SOGv_mw(i,k)     &
                           - MAX(0._dp, (pxtm1(_RI_X_ZN_(jl,jk,idt))          &
                           + pxtte(_RI_X_ZN_(jl,jk,idt))*time_step_len)))       &
                           / time_step_len
        end if
        end do
      end do
     end do
    end do
   end do
    
      IF(p_parallel_io.AND. .NOT. entered) &
     CALL end_message_bi(modstr,'CALL oracle core layer',substr)
   entered = .TRUE.

  END SUBROUTINE oracle_physc 
!***************************************************************************
   SUBROUTINE oracle_radiation

    IMPLICIT NONE

   END SUBROUTINE oracle_radiation
!***************************************************************************
   SUBROUTINE oracle_free_memory

    IMPLICIT NONE

   IF (ALLOCATED (idt_fPOA)) DEALLOCATE (idt_fPOA)
   IF (ALLOCATED (idt_bbPOA)) DEALLOCATE (idt_bbPOA)
   IF (ALLOCATED (idt_fSOAsv))DEALLOCATE (idt_fSOAsv)
   IF (ALLOCATED (idt_bbSOAsv))DEALLOCATE (idt_bbSOAsv)
   IF (ALLOCATED (idt_fSOAiv))DEALLOCATE (idt_fSOAiv)
   IF (ALLOCATED (idt_bbSOAiv))DEALLOCATE (idt_bbSOAiv)
   IF (ALLOCATED (idt_fPOG))  DEALLOCATE (idt_fPOG)
   IF (ALLOCATED (idt_bbPOG))  DEALLOCATE (idt_bbPOG)
   IF (ALLOCATED (idt_fSOGsv)) DEALLOCATE (idt_fSOGsv)
   IF (ALLOCATED (idt_bbSOGsv)) DEALLOCATE (idt_bbSOGsv)
   IF (ALLOCATED (idt_fSOGiv)) DEALLOCATE (idt_fSOGiv)
   IF (ALLOCATED (idt_bbSOGiv)) DEALLOCATE (idt_bbSOGiv)
   IF (ALLOCATED (idt_SOAv))DEALLOCATE (idt_SOAv)
   IF (ALLOCATED (idt_SOGv)) DEALLOCATE (idt_SOGv)
   IF(TRIM(aermod).eq."oracle") then
     IF (ASSOCIATED(dryradius))THEN
        DEALLOCATE(dryradius)
        NULLIFY(dryradius)
     ENDIF
     IF (ASSOCIATED(sigma)) THEN 
        DEALLOCATE(sigma)
        NULLIFY(sigma)
     ENDIF
   ENDIF

   END SUBROUTINE oracle_free_memory
!***************************************************************************
!***************************************************************************
  END MODULE messy_oracle_si
!***************************************************************************

