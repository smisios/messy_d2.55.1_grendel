#include "messy_main_ppd_bi.inc"

MODULE messy_bioburn_si

  !----------------------------------------------------------------------------
  !  bioburn : Model of Biomass Burning Emissions
  !
  !  AUTHOR:  Davic Cabrera, MPICH, SEP 2013
  !
  !  MESSy Interface  for the Submodel kernel
  !----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                      !
  !     Scientific algorithm                                             !
  !                                                                      !
  !             Emission = [EF][Firetype][DM]                            !
  !               where [EF]         = emission factor (gr/kg)           !
  !                     [Firetype]   = type of dyre based on locations   !
  !                     [DM]         = Dry matter burned                 !
  !                                                                      !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY,PTR_1D_ARRAY
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_timer_bi,      ONLY: timer_event_init
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_LAST, time_event
  USE messy_bioburn

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC NULL

  !final emissions
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  :: emis_flux => NULL()
  ! number of supported emission factors
  INTEGER, PARAMETER :: IEF=25
  ! number of actual emission factors (= fire types)
  INTEGER            :: AEF = 0

  TYPE T_BB
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_name     = '' ! name of tracer
     ! unit (gr/m2/s or mlc/m2/s)
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_unit     = ''
     REAL(DP)                     :: MW               ! Molar weight
     REAL(DP)                     :: GS               ! Global scale
     REAL(DP), DIMENSION(IEF)     :: EF               ! emission factors
  END TYPE T_BB

  ! max number of tracers / emissions pairs
  INTEGER, PUBLIC, PARAMETER               :: NMAXNTRAC_bioburn =100
  TYPE(T_BB), DIMENSION(NMAXNTRAC_bioburn) :: BB_input
  ! actual number of tracers
  INTEGER, PUBLIC                          :: nbioburn = 0
  !actual number of vertical levels
  INTEGER, PUBLIC                          :: vlev = 0

  ! dry matter burned
  REAL(dp), DIMENSION(:,:,:), POINTER  :: input_dry_mass   => NULL()
  ! fire type classification
  REAL(dp), DIMENSION(:,:,:), POINTER  :: input_fire_type  => NULL()

  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: fire_type(2) = ''
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: dry_mass(2)  = ''

  ! Vertical fraction after Dentener et al. (2006)
  REAL(dp), DIMENSION(:,:,:), POINTER   :: input_vert_frac  => NULL()
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: vert_frac(2)  = ''

  ! PUBLIC INTERFACE ROUTINES
  PUBLIC :: bioburn_initialize
  PUBLIC :: bioburn_init_memory
  PUBLIC :: bioburn_init_coupling
  PUBLIC :: bioburn_global_end
  PUBLIC :: bioburn_vdiff
  PUBLIC :: bioburn_free_memory

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE bioburn_initialize

    ! bioburn MODULE ROUTINE (SUBMODEL INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! DEFINE NCREGRID EVENT TRIGGERs
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: jt
    INTEGER                     :: jt2

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL bioburn_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast( l_verbose,      p_io)

    !  INIT

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL bioburn_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    ! CPL NAMELIST
    CALL p_bcast(fire_type(1),  p_io)
    CALL p_bcast(fire_type(2),  p_io)
    CALL p_bcast(dry_mass(1),   p_io)
    CALL p_bcast(dry_mass(2),   p_io)
    CALL p_bcast(vert_frac(1),   p_io)
    CALL p_bcast(vert_frac(2),   p_io)

    ! INITIALISE VALUES
    DO jt=1, NMAXNTRAC_bioburn
       BB_input(jt)%tr_name  = ''
       BB_input(jt)%tr_unit  = ''
       BB_input(jt)%MW       = 0.0_dp
       BB_input(jt)%GS       = 0.0_dp
       DO jt2=1, IEF
          BB_input(jt)%EF(jt2) = 0.0_dp
       END DO
    END DO

    ! READ BB namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL bioburn_read_nml_bb(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       nbioburn = 0
       DO jt=1, NMAXNTRAC_bioburn
          IF (TRIM(BB_input(jt)%tr_name) == '') CYCLE
          nbioburn = nbioburn + 1
          !
          BB_input(nbioburn)%tr_name = TRIM(BB_input(jt)%tr_name)
          BB_input(nbioburn)%tr_unit = TRIM(BB_input(jt)%tr_unit)
          BB_input(nbioburn)%MW      = BB_input(jt)%MW
          BB_input(nbioburn)%GS      = BB_input(jt)%GS
          DO jt2=1, IEF
             BB_input(nbioburn)%EF(jt2)   = BB_input(jt)%EF(jt2)
          END DO
       END DO
       write (*,*) " Total number of tracers in bioburn_nml = " , nbioburn
    END IF

    ! BB namelist
    CALL p_bcast(nbioburn, p_io)
    DO jt=1, nbioburn
       CALL p_bcast(BB_input(jt)%tr_name, p_io)
       CALL p_bcast(BB_input(jt)%tr_unit, p_io)
       CALL p_bcast(BB_input(jt)%MW,      p_io)
       CALL p_bcast(BB_input(jt)%GS,      p_io)
       DO jt2=1, IEF
          CALL p_bcast(BB_input(jt)%EF(jt2),     p_io)
       END DO
    END DO

    DO jt=1, nbioburn
       IF (BB_input(jt)%tr_unit == 'mlc/m2/s' .AND. &
            BB_input(jt)%MW == 0.0_dp) THEN
          CALL error_bi('Units are mlc/m2/s --> molar mass is not defined!!' &
               , substr)
       END IF
    END DO

  END SUBROUTINE bioburn_initialize

  ! ---------------------------------------------------------------------------

  SUBROUTINE bioburn_init_memory

    IMPLICIT NONE
    INTRINSIC TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'bioburn_init_memory'
    INTEGER                          :: status

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

    ! allocate emis_flux with nbioburn
    ALLOCATE(emis_flux(nbioburn))

    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

  END SUBROUTINE bioburn_init_memory

  ! ---------------------------------------------------------------------------

  SUBROUTINE bioburn_init_coupling

    USE messy_main_channel_repr,  ONLY: new_representation            &
                                      , get_representation            &
                                      , AUTO, IRANK                   &
                                      , repr_def_axes
    USE messy_main_channel,          ONLY: get_channel_object         &
                                         , get_channel_info           &
                                         , get_channel_object_info    &
                                         , new_channel                &
                                         , new_channel_object         &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_tracer,           ONLY: R_molarmass
    USE messy_main_constants_mem,    ONLY: iouerr

    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER(LEN=*), PARAMETER :: substr = 'bioburn_init_coupling'
    INTEGER :: status
    INTEGER :: jt
    INTEGER :: reprid
    INTEGER :: REPR_BB_3D
    INTEGER :: bioburn

    CALL start_message_bi(modstr, 'COUPLING', substr)

    CALL info_bi('IMPORT DATA: ', substr)
    call get_channel_info(status, fire_type(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
          WRITE(iouerr,*) ' no fire type at surface available  '
          WRITE(iouerr,*) ' channel not present ! (check namelist)'
          CALL error_bi('fire type not present!', substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=fire_type(1) &
            , oname=fire_type(2), reprid=reprid )
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
             WRITE(iouerr,*) ' no fire type at surface available '
             WRITE(iouerr,*) ' object not present ! (check namelist)'
             CALL error_bi('surface fire type not present!', substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) ' Fire type at surface available from channel ', &
                  fire_type(2)
          END IF
          CALL get_channel_object(status, cname=fire_type(1) &
               , oname=fire_type(2) &
               , p3=input_fire_type )
          CALL channel_halt(substr, status)
       END IF
    END IF
    ! set actual number of emission factors (= fire types)
    AEF = SIZE(input_fire_type, DIM=_IN_XY_N_)

    call get_channel_info(status, dry_mass(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
          WRITE(iouerr,*) ' no dry matter burned available  '
          WRITE(iouerr,*) ' channel not present ! (check namelist)'
          CALL error_bi('no dry matter burned) not present!', substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=dry_mass(1) &
            , oname=dry_mass(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
             WRITE(iouerr,*) ' no dry matter burned available  '
             WRITE(iouerr,*) ' object not present ! (check namelist)'
             CALL error_bi('no dry matter burned) not present!', substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) ' dry matter burned available from channel ', &
                  dry_mass(2)
          END IF
          CALL get_channel_object(status, cname=dry_mass(1) &
               , oname=dry_mass(2) &
               , p3=input_dry_mass )
          CALL channel_halt(substr, status)
       END IF
    END IF

   CALL get_channel_info(status, vert_frac(1))
   IF (status /= 0) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) ' no vertical fraction data available  '
          WRITE(*,*) ' channel not present ! (check namelist)'
          CALL error_bi('no vertical fraction data) not present!', substr)
       END IF
    ELSE
       CALL get_channel_object_info(status &
            , cname=vert_frac(1), oname=vert_frac(2), reprid=REPR_BB_3D)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) ' no vertical fraction data available  '
             WRITE(*,*) ' object not present ! (check namelist)'
             CALL error_bi('no vertical fraction data) not present!', substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) ' vertical fraction data available from channel ' &
                  , vert_frac(2)
          END IF
          CALL get_channel_object(status &
               , cname=vert_frac(1), oname=vert_frac(2) &
               , p3=input_vert_frac )
          CALL channel_halt(substr, status)
       END IF
    END IF

    ! create new channel with the representation ID from get_channel_object_info

    CALL new_channel(status,modstr//'_gp', reprid=REPR_BB_3D)
    CALL channel_halt(substr, status)

    DO jt=1,nbioburn
       CALL new_channel_object(status, modstr//'_gp'             &
            , TRIM(BB_input(jt)%tr_name)//'_flux'                &
            , p3=emis_flux(jt)%ptr)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp'                   &
            , TRIM(BB_input(jt)%tr_name)//'_flux'                 &
            , 'longname', c='flux for '//TRIM(BB_input(jt)%tr_name))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp'                   &
            , TRIM(BB_input(jt)%tr_name)//'_flux'                 &
            , 'units', c=TRIM(BB_input(jt)%tr_unit))
       CALL channel_halt(substr, status)
    END DO

    ! read the vertical dimension size of input_vert_frac
    vlev = SIZE(input_vert_frac, DIM=_IZ_XYZ__ )

    CALL end_message_bi(modstr, 'COUPLING', substr)

  END SUBROUTINE bioburn_init_coupling

  ! ---------------------------------------------------------------------------

  SUBROUTINE bioburn_vdiff

    USE messy_main_blather_bi,      ONLY: error_bi
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma
    USE messy_main_constants_mem,   ONLY: g, M_air, N_A

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_vdiff'
    INTEGER              :: jt2
    INTEGER              :: jt, jk ! index for level loop

    IF (l_verbose) THEN
       DO jt=1,nbioburn
          write(*,*) ' nbioburn', nbioburn
          write(*,*) ' size of ef matrix ef_g', shape(emis_flux(jt)%ptr)
          write(*,*) ' size of input firetype size', shape(input_fire_type)
          write(*,*) ' size of input dry matter', shape(input_dry_mass)
          write(*,*) ' size of input vertical fraction', shape(input_vert_frac)
          write(*,*) ' size of jk, jt', jk,jt
          write(*,*) ' size of nbioburn', shape(nbioburn)
          write(*,*) ' kproma, jrow',  kproma,jrow
       ENDDO
    ENDIF

    DO jt=1, nbioburn
       !
       emis_flux(jt)%ptr(_RI_XYZ__(1:kproma,jrow,1:vlev)) = 0.0_dp
       !
       DO jk=1, vlev
          !
          DO jt2=1, AEF
             !
             ! now sum up the emissions from the individual fire types
             !
             emis_flux(jt)%ptr(_RI_XYZ__(1:kproma,jrow,jk)) =           &
                  emis_flux(jt)%ptr(_RI_XYZ__(1:kproma,jrow,jk)) +      &
                  input_vert_frac(_RI_XYZ__(1:kproma,jrow,jk)) *        &
                  input_dry_mass(_RI_XYZ__(1:kproma,jrow,1)) *          &
                  BB_input(jt)%GS *                                     &
                  BB_input(jt)%EF(jt2) *                                &
                  input_fire_type(_RI_XY_N_(1:kproma,jrow,jt2))
          ENDDO
          !
       END DO

       SELECT CASE(TRIM(BB_input(jt)%tr_unit))
            !
         CASE ("gr/m2/s")
            ! do nothing, calculations are fine
         CASE ("mlc/m2/s")
            ! convert units            
            emis_flux(jt)%ptr(_RI_XYZ__(1:kproma,jrow,1:vlev)) =        &
                 emis_flux(jt)%ptr(_RI_XYZ__(1:kproma,jrow,1:vlev))     &
                 * (N_A/BB_input(jt)%MW)                                
         CASE DEFAULT
            !
            CALL error_bi('Emission units not recognised.'//&
                 &' Please correct bioburn.nml', substr)
       END SELECT
       !
    ENDDO

  END SUBROUTINE bioburn_vdiff

  !---------------------------------------------------------------------------

  SUBROUTINE bioburn_global_end

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_global_end'

  END SUBROUTINE bioburn_global_end

  ! ---------------------------------------------------------------------------

  SUBROUTINE bioburn_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(emis_flux)) DEALLOCATE(emis_flux)
    NULLIFY(emis_flux)

  END SUBROUTINE bioburn_free_memory

  ! ===========================================================================
  ! PRIVATE ROUTINES
  ! ===========================================================================

  SUBROUTINE bioburn_read_nml_cpl(status, iou)

    ! read namelist for 'coupling'
    !
    ! Author: Andrea Pozzer, MPICH, Feb 2011

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ fire_type, dry_mass, vert_frac

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE bioburn_read_nml_cpl

  ! ===========================================================================

  SUBROUTINE bioburn_read_nml_bb(status, iou)

    ! read namelist for 'coupling'
    !
    ! Author: Andrea Pozzer, MPICH, Feb 2011

    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    !
    NAMELIST /BB/  BB_input

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='bioburn_read_nml_bb'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt
    INTEGER                     :: jt2

    status = 1

    CALL read_nml_open(lex, substr, iou, 'BB', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=BB, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'BB', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

    WRITE(*,*) '.........................................................'
    WRITE(*,*) '             BIOBURN  TRACER DEFINED                       '
    WRITE(*,*) '.........................................................'

    DO jt=1, NMAXNTRAC_bioburn
       IF (TRIM(BB_input(jt)%tr_name) == '') CYCLE
       WRITE(*,*) '  TRACER NO.          ', jt
       WRITE(*,*) '  NAME              = ', TRIM(BB_input(jt)%tr_name)
       WRITE(*,*) '  UNIT              = ', TRIM(BB_input(jt)%tr_unit)
       WRITE(*,*) '  MW                = ', BB_input(jt)%MW
       WRITE(*,*) '  GS                = ', BB_input(jt)%GS
       DO jt2=1, IEF
          IF (BB_input(jt)%EF(jt2) >= 0.0_dp ) THEN
             WRITE(*,*) '  EF ',jt2,' = ', BB_input(jt)%EF(jt2)
          ELSE
             WRITE(*,*) 'ERROR (EF < 0.0):  EF ',jt2,' = ', BB_input(jt)%EF(jt2)
             status = 1
             RETURN
          END IF
       END DO
       WRITE(*,*) '.........................................................'

    END DO
  END SUBROUTINE bioburn_read_nml_bb

  ! ===========================================================================

END MODULE messy_bioburn_si
