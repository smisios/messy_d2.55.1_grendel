#include "messy_main_ppd_bi.inc"

MODULE messy_qbo_si

#if defined(ECHAM5) || defined(BLANK) || defined(CESM1)

  ! REFERENCES:
  !   see messy_qbo.f90
  !
  ! AUTHORS:
  !   Marco Giorgetta, MPI for Meteorology, Hamburg, October 1999
  !     - original code for ECHAM4
  !   Maarten van Aalst, MPI for Chemistry, Mainz, 2003/2004
  !     - implementation in ECHAM5 (as preliminary MESSy submodel)
  !   Patrick Joeckel, MPI for Chemistry, Mainz, September 2004
  !     - revision for MESSy 0.9
  !

  ! MESSy
  USE messy_qbo
  USE messy_main_channel,    ONLY: t_chaobj_cpl
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &
                                      mtend_id_u
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: ABS, ALLOCATED, NULL

  ! GLOBAL PARAMETERS
  ! --- VERTICAL ---
  ! index of first level of ECHAM grid in range of r_p
  INTEGER :: ifirst
  ! index of last  level of ECHAM grid in range of r_p
  INTEGER :: ilast
  ! index of next following level of r_p with respect to a given ECHAM5 level
  INTEGER,  DIMENSION(:), ALLOCATABLE :: iindex
  ! ECHAM5 full level grid for p0=101325 Pa
  ! (above 100 hPa: level pressure (almost) independent of location)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: zpf0

  ! POINTERS FOR DIAGNOSTIC CHANNEL OBJECTS
  REAL(DP), DIMENSION(:,:,:), POINTER :: uqb ! QBO field used for nudging
  REAL(DP), DIMENSION(:,:,:), POINTER :: dun ! nudging tendency of zonal wind
  REAL(DP), DIMENSION(:,:,:), POINTER :: anu ! nudging amplitude

  ! FLAG-ARRAY FOR SUPRESSING QBO-NUDGING,
  ! IN CASE WIND IS > 99.0 (missing data)
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: sflag => NULL()

  ! CPL NAMELIST
  TYPE(t_chaobj_cpl) :: c_nudge_data ! op_pj_20100827

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! PUBLIC INTERFACE ROUTINES
  PUBLIC :: qbo_initialize
  PUBLIC :: qbo_init_memory
  PUBLIC :: qbo_init_coupling
  PUBLIC :: qbo_physc
  PUBLIC :: qbo_free_memory
  !PRIVATE :: qbo_read_nml_cpl

CONTAINS

! -----------------------------------------------------------------------
  SUBROUTINE qbo_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,   ONLY: error_bi
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'qbo_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr,'INITIALIZATION', substr)

    ! INITIALIZE CTRL-NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL qbo_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF

    CALL p_bcast(r_lat1, p_io)
    CALL p_bcast(r_lat2, p_io)
    CALL p_bcast(r_nudg0, p_io)
    CALL p_bcast(r_nweight, p_io)
    CALL p_bcast(r_hwidth, p_io)

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL qbo_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CPL NAMELIST', substr)
    END IF
    CALL p_bcast(c_nudge_data%cha, p_io)
    CALL p_bcast(c_nudge_data%obj, p_io)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif
    CALL end_message_bi(modstr,'INITIALIZATION', substr)

  END SUBROUTINE qbo_initialize
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE qbo_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID 
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'qbo_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr,'INITIALIZE MEMORY', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'uqb', p3=uqb)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uqb' &
         , 'long_name', c='QBO field used for nudging')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uqb', 'units', c='m/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dun', p3=dun)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dun' &
         , 'long_name', c='nudging tendency of zonal wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dun', 'units', c='m/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'anu', p3=anu)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'anu' &
         , 'long_name', c='nudging amplitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'anu', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'INITIALIZE MEMORY', substr)

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_u)
#endif
  END SUBROUTINE qbo_init_memory
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE qbo_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object &
                                         , get_channel_object_dimvar
    USE messy_main_tools,            ONLY: PTR_1D_ARRAY
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, apzero
    USE messy_main_grid_def_bi,      ONLY: ceta
    USE messy_main_mpi_bi,           ONLY: p_parallel_io

    IMPLICIT NONE

    INTRINSIC :: TRIM, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'qbo_init_coupling'
    INTEGER :: status
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER          :: dvs => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER :: units => NULL()
    INTEGER :: i
    INTEGER,  DIMENSION(nlev)   :: jindex
    INTEGER                     :: jk

    CALL start_message_bi(modstr,'INITIALIZE COUPLING', substr)

    CALL get_channel_object(status &
         , TRIM(c_nudge_data%cha), TRIM(c_nudge_data%obj), p1=qbodd)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status &
         , TRIM(c_nudge_data%cha), TRIM(c_nudge_data%obj)//'_flg', p1=sflag)
    CALL channel_halt(substr, status)

    nqbolev = SIZE(qbodd)
    CALL qbo_add_mem

    CALL get_channel_object_dimvar(status &
         , TRIM(c_nudge_data%cha), TRIM(c_nudge_data%obj) &
         , dvs, units)
    CALL channel_halt(substr, status)

    DO i=1, SIZE(dvs)
       IF (p_parallel_io) &
            WRITE(*,*) 'DIMVAR ',i,' [',TRIM(units(i)),']: ',dvs(i)%ptr
    END DO
    r_p => dvs(1)%ptr

    ! INITIALIZE VERTICAL GRID-MAPPING
    ! ECHAM5 full level grid for p0=101325 Pa
    ALLOCATE(zpf0(nlev))
    zpf0(:) = ceta(:) * apzero
    !
    CALL qbo_interpolate_index(nqbolev, nlev, r_p, zpf0  &
         , ifirst, ilast, jindex)
    !
    ALLOCATE(iindex(nlev))
    iindex(:) = 0
    DO jk=ifirst,ilast
       iindex(jk)=jindex(jk)
    ENDDO

    CALL end_message_bi(modstr,'INITIALIZE COUPLING', substr)

  END SUBROUTINE qbo_init_coupling
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE qbo_physc

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: um1, vom_3d 
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
    USE messy_main_grid_def_bi,     ONLY: philat_2d
    USE messy_main_timer,           ONLY: ztmst=>time_step_len 

    IMPLICIT NONE

    ! LOCAL
    INTEGER  :: iind0, iind1       ! vertical indices
    REAL(DP) :: lati               ! current latitude
    INTEGER  :: jk, jp

    level_loop: DO jk=1, nlev

       ! SKIP IF NOT IN QBO VERTICAL INDEX RANGE
       IF ((jk < ifirst) .OR. (jk > ilast)) CYCLE

       ! set level indices for vertical interpolation
       iind0=iindex(jk)-1
       iind1=iindex(jk)

       vector_loop: DO jp=1, kproma

          ! SKIP IF NOT IN QBO LATITUDE RANGE
          lati = ABS(philat_2d(jp,jrow)) 
          IF (ABS(lati) > r_lat2 ) CYCLE

          ! CALCULATE QBO NUDGING PARAMETERS
          ! compute QBO at latitude lati and level
          CALL qbo_1(uqb(_RI_XYZ__(jp,jrow,jk)), iind0, iind1, zpf0(jk), lati)

          ! compute nudging strength at latitude lati and level
          CALL qbo_2(anu(_RI_XYZ__(jp,jrow,jk)), iind0, iind1, zpf0(jk),lati)

          ! compute tendendy du/dt due to the QBO nudging
          CALL qbo_3(dun(_RI_XYZ__(jp,jrow,jk)), uqb(_RI_XYZ__(jp,jrow,jk)) &
               , anu(_RI_XYZ__(jp,jrow,jk)) &
               , um1(_RI_XYZ__(jp,jrow,jk)), vom_3d(_RI_XYZ__(jp,jrow,jk))  &
               , ztmst)

          ! ZERO NUDGING TENDENCY, IN CASE OF MISSING DATA
          dun(_RI_XYZ__(jp,jrow,jk)) = dun(_RI_XYZ__(jp,jrow,jk)) * sflag(iind0) * sflag(iind1)

          ! APPLY QBO NUDGING TENDENCY
#ifndef MESSYTENDENCY
          vom_3d(_RI_XYZ__(jp,jrow,jk)) = &
               vom_3d(_RI_XYZ__(jp,jrow,jk)) + dun(_RI_XYZ__(jp,jrow,jk))
#endif

       END DO vector_loop

    END DO level_loop

#ifdef MESSYTENDENCY
    CALL mtend_add_l(my_handle, mtend_id_u, px=dun(_RI_XYZ__(1:kproma,jrow,:)))
#endif

  END SUBROUTINE qbo_physc
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE qbo_free_memory

    IMPLICIT NONE

    ! FREE MEMORY (CORE)
    CALL qbo_clean

    ! FREE MEMORY (INTERFACE)
    IF (ALLOCATED(zpf0))   DEALLOCATE(zpf0)
    IF (ALLOCATED(iindex)) DEALLOCATE(iindex)

  END SUBROUTINE qbo_free_memory
! -----------------------------------------------------------------------

! ************************************************************************
! PRIVATE ROUTINES
! ************************************************************************

! -----------------------------------------------------------------------
  SUBROUTINE qbo_read_nml_cpl(status, iou)
   
    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ c_nudge_data

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='qbo_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE qbo_read_nml_cpl
! -----------------------------------------------------------------------

#else

IMPLICIT NONE

#endif

! ************************************************************************
END MODULE messy_qbo_si
! ************************************************************************ 





