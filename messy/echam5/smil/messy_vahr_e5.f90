! ***********************************************************************
MODULE messy_vahr_e5
! ***********************************************************************

  ! MESSy-SMIL FOR SUBMODEL VAHR
  !
  ! VAHR - Volcanic Aerosol Heating Rates 
  ! Description: Inclusion of heating rates from volcanic 
  ! aerosols calculated by Georgiy Stenchikov  
  !
  ! Author: Andreas Baumgaertner, MPICH, Aug 2007

  ! TODO:
  ! use tropopause (channel object) as minimum altitude 

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_vahr

  IMPLICIT NONE
  PRIVATE

  ! GLOBAL PARMETERS
  TYPE IO_VAHR
     ! name of source-channel
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel     = ''
     ! name of source-object
     CHARACTER(LEN=STRLEN_OBJECT)  :: object    = 'VAHR_VOLC_hrates'
     REAL(DP)                      :: coeff  =   0.0_DP ! scale coeff (fraction)
     REAL(DP)                      :: latmin = -90.0_DP
     REAL(DP)                      :: latmax =  90.0_DP
     INTEGER                       :: kmin   =     0
     INTEGER                       :: kmax   =     0
  END TYPE IO_VAHR

  TYPE(IO_VAHR), SAVE :: VAHR

  ! pointer to imported heating rate channel object
  REAL(dp), POINTER, DIMENSION(:,:,:) :: VOLC_hrates  => NULL()

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: vahr_initialize
  PUBLIC :: vahr_init_coupling
  PUBLIC :: vahr_physc
  !PRIVATE :: vahr_read_nml_cpl

CONTAINS

! ========================================================================
  SUBROUTINE vahr_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'vahr_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL vahr_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    CALL p_bcast(VAHR%channel, p_io)
    CALL p_bcast(VAHR%object,  p_io)
    CALL p_bcast(VAHR%coeff,   p_io)
    CALL p_bcast(VAHR%latmin,  p_io)
    CALL p_bcast(VAHR%latmax,  p_io)
    CALL p_bcast(VAHR%kmin,    p_io)
    CALL p_bcast(VAHR%kmax,    p_io)

  END SUBROUTINE vahr_initialize
! ========================================================================

! ========================================================================
  SUBROUTINE vahr_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object
    
    ! LOCAL
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'vahr_init_coupling'

    ! get channel object for vahr
    CALL get_channel_object(status, TRIM(VAHR%channel), TRIM(VAHR%object) &
         , p3=VOLC_hrates)
    CALL channel_halt(substr, status)

    END SUBROUTINE vahr_init_coupling
! ========================================================================

! ========================================================================
  SUBROUTINE vahr_physc

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi,  ONLY: nlev, jrow, kproma 
  USE messy_main_grid_def_bi,      ONLY:philat_2d
  USE messy_main_data_bi,          ONLY:tte_3d                   ! temp tendency
  IMPLICIT NONE

  ! LOCAL
  REAL(DP) :: dayspersecond=1./(3600.*24.)
  INTEGER  :: jp, jk
  INTEGER  :: kmin, kmax, klev
  INTEGER  :: kmin_tmp, kmax_tmp
  LOGICAL  :: llat, llev

  level_loop: DO jk=1, nlev  ! LOOP OVER LEVELS

     vector_loop: DO jp=1, kproma

        llat = (philat_2d(jp, jrow) >= VAHR%latmin) .AND. &
             (philat_2d(jp, jrow) <= VAHR%latmax)

        kmin = VAHR%kmin
        kmax = VAHR%kmax       
        kmin_tmp = MAX(MIN(kmin, kmax),1)
        kmax_tmp = MIN(MAX(kmin, kmax),nlev)
        kmin = kmin_tmp
        kmax = kmax_tmp

        llev = (jk >= kmin) .AND. (jk <= kmax)

        ! ADD TO HEATING RATE
        IF (llat.AND.llev) THEN
           tte_3d(jp,jk,jrow)= tte_3d(jp,jk,jrow)&
                + VOLC_hrates(jp,jk,jrow)*dayspersecond
        END IF

     END DO vector_loop

  END DO level_loop        ! LOOP OVER LEVELS


  END SUBROUTINE vahr_physc
! ========================================================================

! ========================================================================
  SUBROUTINE vahr_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ VAHR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='vahr_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE vahr_read_nml_cpl
! ========================================================================

! ***********************************************************************
END MODULE messy_vahr_e5
! ***********************************************************************

