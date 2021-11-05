! **************************************************************************
#include "messy_main_ppd_bi.inc"
MODULE messy_visop_si
! **************************************************************************

#if defined (ICON) || defined (COSMO)

#ifdef _MESSY_OMP
#ifdef ICON
#include "omp_definitions.inc"
#endif
#endif

  ! VISOP
  ! visible satellite image forward operator
  ! 
  ! 2014 Leonhard Scheck

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_visop

#ifdef ICON
  ! needs to be a PARAMETER below, since ALLOCATABLE is not allowed for
  ! NAMELIST variables ...
  USE messy_main_bmluse_bi,      ONLY: max_dom
#endif
  
  USE messy_main_tools, ONLY: PTR_2D_ARRAY, PTR_3D_ARRAY

  IMPLICIT NONE
  PRIVATE
  SAVE

#ifndef ICON
  INTEGER, PARAMETER :: max_dom = 1
#endif
  
  INTEGER, PARAMETER :: n2d = 10
  TYPE(PTR_2D_ARRAY), DIMENSION(n2d) :: v2d
  ! VIS006 SEVIRI channel
  INTEGER, PARAMETER :: id2_tau06w = 1  ! column integrated water optical depth
  INTEGER, PARAMETER :: id2_tau06i = 2  ! column integrated ice   optical depth
  ! VIS008 SEVIRI channel
  INTEGER, PARAMETER :: id2_tau08w = 3  ! column integrated water optical depth
  INTEGER, PARAMETER :: id2_tau08i = 4  ! column integrated ice   optical depth
  INTEGER, PARAMETER :: id2_reffw = 5   ! representative effective radius for water droplets [um]
  INTEGER, PARAMETER :: id2_reffi = 6   ! representative effective radius for ice particels  [um]
  INTEGER, PARAMETER :: id2_cwtop = 7   ! water cloud top  [m]
  INTEGER, PARAMETER :: id2_cwbase = 8  !   "      "  base [m]
  INTEGER, PARAMETER :: id2_citop = 9   ! ice   cloud top  [m]
  INTEGER, PARAMETER :: id2_cibase = 10 !  "       "  base [m]


  INTEGER, PARAMETER :: n3d = 6
  TYPE(PTR_3D_ARRAY), DIMENSION(n3d) :: v3d
  INTEGER, PARAMETER :: id3_temp = 1  ! temperature                  [K]
  INTEGER, PARAMETER :: id3_rho =  2  ! density                      [kg/m3]
  INTEGER, PARAMETER :: id3_qc =   3  ! specific cloud water content [kg/kg]
  INTEGER, PARAMETER :: id3_qi =   4  ! specific cloud ice   content [kg/kg]
  INTEGER, PARAMETER :: id3_qni =  5  ! cloud ice number concentration [1/kg]
  INTEGER, PARAMETER :: id3_z =    6  ! geometric height             [m]


  CHARACTER(LEN=30), DIMENSION(3,n2d) :: chaobj_int = RESHAPE( (/ &
       'tau06w                        ','                              ',&
       'column int. water opt. depth  ', &
       'tau06i                        ','                              ',&
       'column int. ice opt. depth    ', &
       'tau08w                        ','                              ',&
       'column int. water opt. depth  ', &
       'tau08i                        ','                              ',&
       'column int. ice opt. depth    ', &
       'reffw                         ','microns                       ',&
       'eff. radius for water droplets', &
       'reffi                         ','microns                       ',&
       'eff. radius for ice particles ', &
       'cwtop                         ','m                             ',&
       'water cloud top               ', &
       'cwbase                        ','m                             ',&
       'water cloud base              ', &
       'citop                         ','m                             ',&
       'ice cloud top                 ', &
       'cibase                        ','m                             ',&
       'ice cloud base                ' &
       /), SHAPE=(/3,n2d/) )

#if defined (ICON)
  CHARACTER(LEN=16), DIMENSION(2,n3d) :: chaobj_ext = RESHAPE( (/ &
       'nh_state_diag   ', 'temp            ', &
       'nh_state_prog   ', 'rho             ', &
       'nh_state_prog   ', 'qc              ', &
       'nh_state_prog   ', 'qi              ', &
       'nh_state_prog   ', 'qni             ', &
       'nh_state_metrics', 'z_ifc           '  &
       /), SHAPE=(/2,n3d/) )
#endif

#if defined (COSMO)
  CHARACTER(LEN=16), DIMENSION(2,n3d) :: chaobj_ext = RESHAPE( (/ &
       'COSMO           ', 'tm1             ', &
       'COSMO           ', 'rho_air_dry     ', &
       'tracer_gp       ', 'QC              ', &
       'tracer_gp       ', 'QI              ', &
       '                ', 'QNI             ', &
       'COSMO_ORI       ', 'HHL             '  &
       /), SHAPE=(/2,n3d/) )
#endif

  INTEGER, DIMENSION(max_dom)        :: patches = -1
  LOGICAL, DIMENSION(:), ALLOCATABLE :: lcalc
  LOGICAL                            :: qni_avail=.TRUE.

  PUBLIC :: visop_initialize
  PUBLIC :: visop_init_memory
  PUBLIC :: visop_init_coupling
  PUBLIC :: visop_local_end
  PUBLIC :: visop_free_memory
  PUBLIC :: visop_set_domain

CONTAINS

  ! -------------------------------------------------------------------------
  SUBROUTINE visop_initialize
    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_main_channel_bi, ONLY: n_dom

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'visop_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! ### initialize submodel / read namelist...
    CALL start_message_bi(modstr,'INITIALIZE',substr)  ! log-output

    !! INITIALIZE MAIN-CTRL
    !IF (p_parallel_io) THEN
    !   iou = find_next_free_unit(100,200)
    !   CALL visop_read_nml_ctrl(status, iou)
    !   IF (status /= 0) CALL error_bi('error in namelist CTRL ',substr)
    !END IF

    ! initialize CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL visop_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    DO i=1, max_dom
       CALL p_bcast(patches(i), p_io)
    END DO

    ALLOCATE(lcalc(n_dom), stat=status)
    lcalc = .FALSE.
    IF (patches(1) == -1) THEN
       lcalc = .TRUE.
    ELSE
       DO i = 1, max_dom
          IF (patches(i) > n_dom) CYCLE
          IF (patches(i) < 1) CYCLE
          lcalc(patches(i)) = .TRUE.                   
       END DO
    END IF

    CALL end_message_bi(modstr,'INITIALIZE',substr)  ! log-output


  END SUBROUTINE visop_initialize
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE visop_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------


    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL
    USE messy_main_channel_mem,       ONLY: dom_current
    USE messy_main_channel,           ONLY: new_channel, new_channel_object  &
                                          , new_attribute
    USE messy_main_channel_mem,       ONLY: dom_unbound

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'visop_init_memory'
    INTEGER                       :: status
    INTEGER                       :: i

    IF (.NOT. (lcalc(dom_current) .OR. (dom_current == dom_unbound))) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! define new channel
    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    DO i=1, n2d
       CALL new_channel_object(status, modstr, &
            TRIM(chaobj_int(1,i)), p2=v2d(i)%ptr, reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, &
            TRIM(chaobj_int(1,i)), 'units', c=TRIM(chaobj_int(2,i)))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, &
            TRIM(chaobj_int(1,i)), 'long_name', c=TRIM(chaobj_int(3,i)))
       CALL channel_halt(substr, status)
    END DO

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
 
  END SUBROUTINE visop_init_memory
  ! -------------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  SUBROUTINE visop_init_coupling

    IMPLICIT NONE

    CALL visop_set_domain

  END SUBROUTINE visop_init_coupling
  ! -------------------------------------------------------------------------

  ! ====================================================================
  SUBROUTINE visop_local_end

    USE messy_main_grid_def_mem_bi,  ONLY: jrow, nlev, startidx, endidx
    USE messy_main_channel_mem,      ONLY: dom_current
    USE messy_main_channel,          ONLY: get_channel_output
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_mem,      ONLY: dom_unbound

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'visop_local_end'
    INTEGER :: status
    INTEGER :: i
    LOGICAL :: lout, lstat

    IF (.NOT. (lcalc(dom_current) .OR. (dom_current == dom_unbound))) RETURN

    CALL get_channel_output(status, modstr, lout, lstat)
    CALL channel_halt(substr, status)

    ! CALCULATE ONLY, IF OUTPUT IS TRIGGERED
    IF (.NOT. (lout .OR. lstat)) RETURN

    IF (qni_avail) THEN
#ifdef _MESSY_OMP
#ifdef ICON
!$OMP PARALLEL
!$OMP DO PRIVATE(i) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(i) SCHEDULE(dynamic,1)
#endif
#endif

    DO i=startidx,endidx
       CALL OPTPROP_COLUMN_QNI( nlev &
            , v3d(id3_z)%ptr(_RI_XYZ__(i,jrow,:))    &
            , v3d(id3_rho)%ptr(_RI_XYZ__(i,jrow,:))  &
            , v3d(id3_qc)%ptr(_RI_XYZ__(i,jrow,:))   &
            , v3d(id3_qi)%ptr(_RI_XYZ__(i,jrow,:))   &
            , v3d(id3_qni)%ptr(_RI_XYZ__(i,jrow,:))  &
            , v2d(id2_tau06w)%ptr(i,jrow),   v2d(id2_tau06i)%ptr(i,jrow)  &
            , v2d(id2_tau08w)%ptr(i,jrow),   v2d(id2_tau08i)%ptr(i,jrow)  &
            , v2d(id2_reffw)%ptr(i,jrow),    v2d(id2_reffi)%ptr(i,jrow)   &
            , v2d(id2_cwtop)%ptr(i,jrow),    v2d(id2_cwbase)%ptr(i,jrow)  &
            , v2d(id2_citop)%ptr(i,jrow),    v2d(id2_cibase)%ptr(i,jrow)   )
    END DO

#ifdef _MESSY_OMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

    ELSE

#ifdef _MESSY_OMP
#ifdef ICON
!$OMP PARALLEL
!$OMP DO PRIVATE(i) ICON_OMP_DEFAULT_SCHEDULE
#else
!$OMP PARALLEL
!$OMP DO PRIVATE(i) SCHEDULE(dynamic,1)
#endif
#endif

    DO i=startidx,endidx
       CALL OPTPROP_COLUMN( nlev &
            , v3d(id3_z)%ptr(_RI_XYZ__(i,jrow,:))    &
            , v3d(id3_temp)%ptr(_RI_XYZ__(i,jrow,:)) &
            , v3d(id3_rho)%ptr(_RI_XYZ__(i,jrow,:))  &
            , v3d(id3_qc)%ptr(_RI_XYZ__(i,jrow,:))   &
            , v3d(id3_qi)%ptr(_RI_XYZ__(i,jrow,:))   &
            , v2d(id2_tau06w)%ptr(i,jrow),   v2d(id2_tau06i)%ptr(i,jrow)  &
            , v2d(id2_tau08w)%ptr(i,jrow),   v2d(id2_tau08i)%ptr(i,jrow)  &
            , v2d(id2_reffw)%ptr(i,jrow),    v2d(id2_reffi)%ptr(i,jrow)   &
            , v2d(id2_cwtop)%ptr(i,jrow),    v2d(id2_cwbase)%ptr(i,jrow)  &
            , v2d(id2_citop)%ptr(i,jrow),    v2d(id2_cibase)%ptr(i,jrow)   )
    END DO

#ifdef _MESSY_OMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif
    END IF

  END SUBROUTINE visop_local_end
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE visop_free_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to deallocate the memory, which has
    ! been "manually" allocated in visop_init_memory.
    ! Note: channel object memory must not be deallocated! This is
    !       performed centrally.
    ! ------------------------------------------------------------------

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'visop_free_memory'

    IF (ALLOCATED(lcalc)) DEALLOCATE(lcalc)

  END SUBROUTINE visop_free_memory
  ! ====================================================================

  ! -------------------------------------------------------------------------
  SUBROUTINE visop_set_domain

    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_mem,      ONLY: dom_unbound, dom_current

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER   :: substr = 'visop_set_domain'
    INTEGER                       :: status
    INTEGER                       :: i


    IF (.NOT. (lcalc(dom_current) .OR. (dom_current == dom_unbound))) RETURN

    ! INTERNAL OBJECTS
    DO i=1, n2d
       CALL get_channel_object(status, modstr &
            , TRIM(chaobj_int(1,i)), p2=v2d(i)%ptr )
       CALL channel_halt(substr, status)       
    END DO

    ! EXTERNAL OBJECTS
    DO i=1, n3d
       CALL get_channel_object(status, TRIM(chaobj_ext(1,i)) &
            , TRIM(chaobj_ext(2,i)), p3=v3d(i)%ptr )
       IF (i == id3_qni) THEN
          qni_avail = (status == 0)
       ELSE
          CALL channel_halt(substr, status)
       END IF
    END DO

  END SUBROUTINE visop_set_domain
  ! -------------------------------------------------------------------------

  !==========================================================================
  ! PRIVATE ROUTINES
  !==========================================================================

  !--------------------------------------------------------------------------
  SUBROUTINE visop_read_nml_cpl(status, iou)

    ! VISOP MODULE ROUTINE (PRIVATE)
    !
    ! read namelist for 'coupling'
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_visop,      ONLY : modstr

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'visop_read_nml_cpl'

    NAMELIST /CPL/ patches

    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not. lex) RETURN ! namelist file (<modstr>.nml) not available

    READ (iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! ...done, without error


  END SUBROUTINE visop_read_nml_cpl
  !------------------------------------------------------------------------

#endif

! **************************************************************************
END MODULE messy_visop_si
! **************************************************************************
