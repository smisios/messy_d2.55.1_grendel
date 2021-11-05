#include "messy_main_ppd_bi.inc"

! **************************************************************************
MODULE messy_vaxtra_si
! **************************************************************************

  ! MODULE FOR VERTICAL AXES TRANSFORMATION
  !
  ! Authors: Patrick Joeckel, DLR, December 2016
  !

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi
  USE messy_main_channel_bi,     ONLY: n_dom
  USE messy_main_channel,        ONLY: STRLEN_CHANNEL, STRLEN_OBJECT, REPR_UNDEF
  USE messy_main_constants_mem,  ONLY: STRLEN_MEDIUM, FLAGGED_BAD
  USE messy_vaxtra

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! ----- VERTICAL AXES ----------------------------------------------------
  ! MAX NUMBER OF VERTICAL AXES IN NAMELIST
  INTEGER, PARAMETER :: NMAXVAX    = 50
  ! MAX DIMENSION LENGTH OF AXES
  INTEGER, PARAMETER :: NMAXVAXLEN = 100

  TYPE IO_VAXIS
     CHARACTER(LEN=STRLEN_MEDIUM)   :: axname = ''
     CHARACTER(LEN=STRLEN_CHANNEL)  :: channel  = ''
     CHARACTER(LEN=STRLEN_OBJECT)   :: object = ''
     REAl(DP)                       :: scal = 1.0_dp
     CHARACTER(LEN=STRLEN_MEDIUM)   :: unit = ''
     LOGICAL                        :: topdown = .false.
     LOGICAL                        :: llog = .false.
     INTEGER                        :: dimlen
     REAL(DP),DIMENSION(NMAXVAXLEN) :: axis = 0.0_dp
  END type IO_VAXIS

  TYPE VAXIS
     TYPE(IO_VAXIS) :: io
     REAL(DP), DIMENSION(:,:,:), POINTER :: ptr => NULL() ! channel object
     ! representation ID ...
     ! ... of object on which trafo is based
     INTEGER                             :: reprid_obj = REPR_UNDEF 
     ! ... of new representation with new vertical axis
     INTEGER                             :: reprid     = REPR_UNDEF
     INTEGER                             :: domain_idx = 0
  END type VAXIS

  TYPE(IO_VAXIS), DIMENSION(NMAXVAX)     :: VAX
  INTEGER                                :: NVAX = 0
  TYPE(VAXIS), DIMENSION(:), ALLOCATABLE :: XVAX
  ! ----- END VERTICAL AXES ------------------------------------------------

  ! ----- TRANSFORMED OBJECTS ----------------------------------------------
  INTEGER, PARAMETER :: NMAXTRA = 300
  INTEGER, PARAMETER :: IPO_NONE    =  0  ! no vertical interpolation
  INTEGER, PARAMETER :: IPO_MID2INT =  1  ! GP_3D_MID -> GP_3D_INT
  INTEGER, PARAMETER :: IPO_INT2MID = -1  ! GP_3D_INT -> GP_3D_MID

  TYPE IO_TRAOBJ
     CHARACTER(LEN=STRLEN_OBJECT)   :: name = ''  ! new object name
     CHARACTER(LEN=STRLEN_CHANNEL)  :: cha  = ''  ! channel name
     CHARACTER(LEN=STRLEN_OBJECT)   :: obj  = ''  ! object name
     CHARACTER(LEN=STRLEN_MEDIUM)   :: axn  = ''  ! axis name
     REAL(DP)                       :: missing = FLAGGED_BAD
  END type IO_TRAOBJ

  TYPE TRAOBJ
     TYPE(IO_TRAOBJ)                     :: io
     REAL(DP), DIMENSION(:,:,:), POINTER :: inp => NULL() ! input object
     ! representation ID
     INTEGER                             :: reprid = REPR_UNDEF  
     INTEGER                             :: ipo  = IPO_NONE
     INTEGER                             :: axid = 0      ! axis ID
     INTEGER                             :: domain_idx = 0
     REAL(DP), DIMENSION(:,:,:), POINTER :: out => NULL() ! ouput object
  END type TRAOBJ

  TYPE(IO_TRAOBJ), DIMENSION(NMAXTRA)        :: TRA
  INTEGER                                    :: NTRA = 0
  TYPE(TRAOBJ),    DIMENSION(:), ALLOCATABLE :: XTRA
  ! ----- END TRANSFORMED OBJECTS ------------------------------------------

  PUBLIC :: vaxtra_initialize
  PUBLIC :: vaxtra_init_memory
  PUBLIC :: vaxtra_init_coupling
  PUBLIC :: vaxtra_global_end
  PUBLIC :: vaxtra_free_memory
  !PRIVATE :: vaxtra_read_nml_cpl

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_initialize

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit, quick_sort &
                                   , domains_from_string

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'vaxtra_initialize'
    INTEGER                        :: iou    ! I/O unit
    INTEGER                        :: status ! error status
    INTEGER                        :: i, j
    LOGICAL                        :: lfound
    INTEGER, DIMENSION(NMAXVAXLEN) :: order = 0
    CHARACTER(LEN=32)              :: varname = '  '
    INTEGER, DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                        :: nd, num

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL vaxtra_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    ! ----- VERTICAL AXES --------------------------------
    IF (p_parallel_io) THEN
       ! GET NUMBER OF VERTICAL AXES
       NVAX = 1
       DO i=1, NMAXVAX
          IF (TRIM(ADJUSTL(VAX(i)%axname)) == '') CYCLE
          IF (TRIM(ADJUSTL(VAX(i)%channel)) == '') CYCLE
          IF (TRIM(ADJUSTL(VAX(i)%object)) == '') CYCLE
          ! analyze name for domains
          CALL domains_from_string(status,VAX(i)%axname, n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing (VAX)' &
               , substr)
          NVAX = NVAX + num
       END DO
       NVAX = NVAX - 1
    END IF
    CALL p_bcast(NVAX, p_io)
    ALLOCATE(XVAX(NVAX))

    IF (p_parallel_io) THEN
       NVAX = 1
       ! CHECK AND COPY DATA
       DO i=1, NMAXVAX
          IF (TRIM(ADJUSTL(VAX(i)%axname)) == '') CYCLE

          CALL domains_from_string(status, VAX(i)%axname, n_dom, num &
                ,varname, dnums=domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing (VAX)' &
               , substr)
          
          ! XVAX(NVAX)%io%axname = TRIM(ADJUSTL(VAX(i)%axname))
          domain_loop1: DO nd = 1, SIZE(domnum)
             XVAX(NVAX)%io%axname  = TRIM(varname)
             XVAX(NVAX)%domain_idx = domnum(nd)
             WRITE(*,*) 'VERTICAL AXIS: ''',TRIM(XVAX(NVAX)%io%axname),''''
             !
             IF (TRIM(ADJUSTL(VAX(i)%channel)) == '') THEN
                WRITE(*,*) ' ... empty channel ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '      CHANNEL: ',TRIM(ADJUSTL(VAX(i)%channel))
                XVAX(NVAX)%io%channel = TRIM(ADJUSTL(VAX(i)%channel))
             END IF
             !
             IF (TRIM(ADJUSTL(VAX(i)%object)) == '') THEN
                WRITE(*,*) ' ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '      OBJECT : ',TRIM(ADJUSTL(VAX(i)%object))
                XVAX(NVAX)%io%object = TRIM(ADJUSTL(VAX(i)%object))
             END IF
             !
             XVAX(NVAX)%io%scal    = VAX(i)%scal
             WRITE(*,*) '      SCALING: ', XVAX(NVAX)%io%scal
             !
             XVAX(NVAX)%io%unit    = TRIM(ADJUSTL(VAX(i)%unit))
             WRITE(*,*) '      UNIT   : ',TRIM(VAX(i)%unit)
             !
             XVAX(NVAX)%io%topdown = VAX(i)%topdown
             WRITE(*,*) '      TOPDOWN: ', XVAX(NVAX)%io%topdown
             !
             XVAX(NVAX)%io%llog = VAX(i)%llog
             WRITE(*,*) '  LOGARITHMIC: ', XVAX(NVAX)%io%llog
             !
             XVAX(NVAX)%io%dimlen  = VAX(i)%dimlen
             WRITE(*,*) '      LENGTH : ', XVAX(NVAX)%io%dimlen
             !
             XVAX(NVAX)%io%axis(1:XVAX(NVAX)%io%dimlen) = &
                  VAX(i)%axis(1:VAX(i)%dimlen)
             CALL quick_sort(XVAX(NVAX)%io%axis(1:XVAX(NVAX)%io%dimlen) &
                  , order(1:XVAX(NVAX)%io%dimlen))
             WRITE(*,*) '      VALUES : ', &
                  XVAX(NVAX)%io%axis(1:XVAX(NVAX)%io%dimlen)

             ! NEXT AXIS
             NVAX = NVAX + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop1
          DEALLOCATE(domnum)
          NULLIFY(domnum)
       END DO
       NVAX = NVAX - 1 
    END IF
    CALL p_bcast(NVAX, p_io)
    IF (NVAX /= SIZE(XVAX)) &
         CALL error_bi('error determing number of vertical axes' , substr)
    
    ! BROADCAST RESULTS
    DO i=1, NVAX
       CALL p_bcast(XVAX(i)%io%axname,  p_io)
       CALL p_bcast(XVAX(i)%io%channel, p_io)
       CALL p_bcast(XVAX(i)%io%object,  p_io)
       CALL p_bcast(XVAX(i)%io%scal,    p_io)
       CALL p_bcast(XVAX(i)%io%unit,    p_io)
       CALL p_bcast(XVAX(i)%io%topdown, p_io)
       CALL p_bcast(XVAX(i)%io%llog,    p_io)
       CALL p_bcast(XVAX(i)%io%dimlen,  p_io)
       CALL p_bcast(XVAX(i)%io%axis(:), p_io)
       CALL p_bcast(XVAX(i)%domain_idx, p_io)
    END DO
    ! ----- END VERTICAL AXES --------------------------------

    ! ----- TRANSFORMED OBJECTS ------------------------------
    IF (p_parallel_io) THEN
       ! GET NUMBER OF OBJECTS
       NTRA = 1
       DO i=1, NMAXTRA
          IF (TRIM(ADJUSTL(TRA(i)%name)) == '') CYCLE
          IF (TRIM(ADJUSTL(TRA(i)%cha)) == '') CYCLE
          IF (TRIM(ADJUSTL(TRA(i)%obj)) == '') CYCLE
          IF (TRIM(ADJUSTL(TRA(i)%axn)) == '') CYCLE
          ! analyze name for domains
          CALL domains_from_string(status,TRA(i)%name,n_dom,num)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)

          NTRA = NTRA + num
       END DO
       NTRA = NTRA - 1
    END IF
    CALL p_bcast(NTRA, p_io)
    ALLOCATE(XTRA(NTRA))
    
    IF (p_parallel_io) THEN
       NTRA = 1
       ! CHECK AND COPY DATA
       DO i=1, NMAXTRA
          IF (TRIM(ADJUSTL(TRA(i)%name)) == '') CYCLE

          CALL domains_from_string(status,TRA(i)%name, n_dom, num &
                ,varname, dnums=domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          
          domain_loop2: DO nd = 1, SIZE(domnum)
             XTRA(NTRA)%io%name    = TRIM(varname)
             XTRA(NTRA)%domain_idx = domnum(nd)
             WRITE(*,*) ' OBJECT       : ''',TRIM(XTRA(NTRA)%io%name),''''
             !
             IF (TRIM(ADJUSTL(TRA(i)%cha)) == '') THEN
                WRITE(*,*) ' ... empty input channel ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) ' INPUT CHANNEL: ',TRIM(TRA(i)%cha)
                XTRA(NTRA)%io%cha = TRIM(ADJUSTL(TRA(i)%cha))
             END IF
             !
             IF (TRIM(ADJUSTL(TRA(i)%obj)) == '') THEN
                WRITE(*,*) ' ... empty input object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) ' INPUT OBJECT : ',TRIM(TRA(i)%obj)
                XTRA(NTRA)%io%obj = TRIM(ADJUSTL(TRA(i)%obj))
             END IF
             !
             IF (TRIM(ADJUSTL(TRA(i)%axn)) == '') THEN
                WRITE(*,*) ' ... empty axis name ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) ' AXIS NAME    : ',TRIM(TRA(i)%axn)
                XTRA(NTRA)%io%axn = TRIM(ADJUSTL(TRA(i)%axn))
             END IF
             !
             XTRA(NTRA)%io%missing = TRA(i)%missing
             WRITE(*,*) ' MISSING VALUE: ',TRA(i)%missing
             !
             ! NEXT OBJECT
             NTRA = NTRA + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop2
          DEALLOCATE(domnum)
          NULLIFY(domnum)
       END DO
       NTRA = NTRA - 1
    END IF
    CALL p_bcast(NTRA, p_io)
    IF (NTRA /= SIZE(XTRA)) &
         CALL error_bi('error determing number of objects' &
         , substr)

    ! BROADCAST RESULTS
    DO i=1, NTRA
       CALL p_bcast(XTRA(i)%io%name,    p_io)
       CALL p_bcast(XTRA(i)%io%cha,     p_io)
       CALL p_bcast(XTRA(i)%io%obj,     p_io)
       CALL p_bcast(XTRA(i)%io%axn,     p_io)
       CALL p_bcast(XTRA(i)%io%missing, p_io)
       CALL p_bcast(XTRA(i)%domain_idx, p_io)
    END DO
    ! ----- END TRANSFORMED OBJECTS --------------------------

    ! ----- MAP AXIS ID --------------------------------------
    DO i=1, NTRA
       lfound = .FALSE.
       DO j=1, NVAX
          IF (XVAX(j)%domain_idx /= XTRA(i)%domain_idx) CYCLE
          IF (TRIM(XVAX(j)%io%axname) == TRIM(XTRA(i)%io%axn)) THEN
             XTRA(i)%axid = j
             lfound = .TRUE.
             EXIT ! exit inner do loop
          END IF
       END DO
       IF (.NOT. lfound) THEN
          CALL error_bi('AXIS '//TRIM(XTRA(i)%io%axn)//' REQUESTED FOR '//&
               &' OBJECT '//TRIM(XTRA(i)%io%name)//' NOT DEFINED',substr)
       END IF
    END DO
    ! --------------------------------------------------------

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NVAX,' VERTICAL AXES DEFINED !'
       WRITE(*,*) ' ---> ',NTRA,' OBJECTS TO BE TRANSFORMED !'
    END IF

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE vaxtra_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_init_memory

    USE messy_main_channel_error_bi,   ONLY: channel_halt
    USE messy_main_channel_bi,         ONLY: DIMID_LAT, DIMID_LON, gp_nseg &
                                           , DC_GP                         &
                                           , gp_start, gp_cnt, gp_memu     &
                                           , gp_meml
    USE messy_main_grid_def_mem_bi,    ONLY: nproma, ngpblks
    USE messy_main_channel_dimensions, ONLY: new_dimension &
                                           , add_dimension_variable        &
                                           , add_dimension_variable_att
    USE messy_main_channel_repr,       ONLY: new_representation            &
                                           , set_representation_decomp     &
                                           , AUTO, IRANK, PIOTYPE_COL      &
                                           , repr_def_axes
    USE messy_main_blather_bi,        ONLY: info_bi
    USE messy_main_mpi_bi,            ONLY: p_parallel_io
    USE messy_main_channel_mem,       ONLY: dom_curid

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'vaxtra_init_memory'
    INTEGER                     :: status
    INTEGER                     :: i   
    INTEGER                     :: dimid
    !
    ! FOR PARALLEL DECOMPOSITION
    INTEGER                          :: nseg    = 0  ! max no. of segments
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()


    CALL start_message_bi(modstr, 'REPRESENTATION DEFINITIONS', substr)

    ! for parallel I/O
    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    DO i=1, NVAX
       IF (XVAX(i)%domain_idx /= dom_curid) CYCLE
       
       ! define new dimension
       CALL info_bi(' DIMENSION     : '//TRIM(XVAX(i)%io%axname),substr)
       CALL new_dimension(status, dimid &
            , TRIM(XVAX(i)%io%axname), XVAX(i)%io%dimlen)
       CALL channel_halt(substr, status)

       ! add dimension variable (of same name!)
       CALL info_bi(' DIM.VARIABLE  : '//TRIM(XVAX(i)%io%axname),substr)
       CALL add_dimension_variable(status, dimid, TRIM(XVAX(i)%io%axname) &
            , XVAX(i)%io%axis(1:XVAX(i)%io%dimlen))
       CALL channel_halt(substr, status)

       ! add dimension variable attributes
       CALL info_bi(' UNIT          : '//TRIM(XVAX(i)%io%unit),substr)
       CALL add_dimension_variable_att(status, dimid, TRIM(XVAX(i)%io%axname) &
            , 'units', c=TRIM(XVAX(i)%io%unit))
       CALL channel_halt(substr, status)
       !
       IF (XVAX(i)%io%topdown) THEN
          CALL info_bi(' ATTRIBUTE     : positive=down',substr)
          CALL add_dimension_variable_att(status, dimid &
               , TRIM(XVAX(i)%io%axname) &
               , 'positive', c='down')
          CALL channel_halt(substr, status)
       ENDIF
       !
       CALL info_bi(' ATTRIBUTE     : interpolation',substr)
       IF (XVAX(i)%io%llog) THEN
          CALL add_dimension_variable_att(status, dimid &
               , TRIM(XVAX(i)%io%axname) &
               , 'interpolation', c='logarithmic')
       ELSE
          CALL add_dimension_variable_att(status, dimid &
               , TRIM(XVAX(i)%io%axname) &
               , 'interpolation', c='linear')
       END IF
       CALL channel_halt(substr, status)


       CALL info_bi(' REPRESENTATION: GP_3D_VAX_'//&
            &TRIM(XVAX(i)%io%axname),substr)
       CALL new_representation(status, XVAX(i)%reprid            &
            , 'GP_3D_VAX_'//TRIM(XVAX(i)%io%axname)              &
            , rank = 3, link = 'xxx-', dctype = DC_GP            &
            , dimension_ids = (/_RI_XYZ__(DIMID_LON, DIMID_LAT, dimid) /)    &
            , ldimlen       = (/ _RI_XYZ__(nproma , ngpblks, AUTO)    /)     &
            , output_order  = (/ _IX_XYZ__ , _IY_XYZ__ , _IZ_XYZ__ /) &
            , axis = repr_def_axes(_RI_XYZ__('X','Y','Z'),'-')               &
           )
       CALL channel_halt(substr, status)

       start(:,:)= gp_start(:,:)
       meml(:,:) = gp_meml(:,:)

       cnt(:,:)  = gp_cnt(:,:)
       memu(:,:) = gp_memu(:,:)

       start(:,4) = 1
       cnt(:,4) = 1
       meml(:,4) = 1
       memu(:,4) = 1

       start(:,_IY_XYZ__) = gp_start(:,_IY_XYZN_)
       cnt(:,_IY_XYZ__)   = gp_cnt(:,_IY_XYZN_)
       meml(:,_IY_XYZ__)  = gp_meml(:,_IY_XYZN_)
       memu(:,_IY_XYZ__)  = gp_memu(:,_IY_XYZN_)

       cnt(:,_IZ_XYZ__)  = XVAX(i)%io%dimlen
       memu(:,_IZ_XYZ__) = XVAX(i)%io%dimlen

       CALL set_representation_decomp(status, XVAX(i)%reprid &
            , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       CALL channel_halt(substr, status)

       IF (p_parallel_io) &
       WRITE(*,*) '------------------------------------------------------'
    END DO

    !----- DEALLOCATE------------------------------------------
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)
    ! ---------------------------------------------------------

    CALL end_message_bi(modstr, 'REPRESENTATION DEFINITIONS', substr)

  END SUBROUTINE vaxtra_init_memory
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_init_coupling

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, warning_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_3D_INT 
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , get_channel_object_info &
                                         , new_channel, new_attribute &
                                         , new_channel_object &
                                         , copy_channel_object_atts
    USE messy_main_channel_mem,       ONLY: dom_curid
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'vaxtra_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: i
    INTEGER                     :: reprid
    INTEGER                     :: axis_id

    CALL start_message_bi(modstr, 'COUPLING TO EXTERNAL OBJECTS', substr)

    ! FOR AXES
    DO i=1, NVAX
       IF (XVAX(i)%domain_idx /= dom_curid) CYCLE
       
       IF (p_parallel_io) THEN
          WRITE(*,*) ' CHANNEL / OBJECT ' &
               , TRIM(XVAX(i)%io%channel),' / ',TRIM(XVAX(i)%io%object) &
               , ' FOR AXIS ', TRIM(XVAX(i)%io%axname)
       ENDIF
       CALL get_channel_object(status &
            , TRIM(XVAX(i)%io%channel), TRIM(XVAX(i)%io%object) &
            , p3=XVAX(i)%ptr)
       CALL channel_halt(substr, status)
       CALL get_channel_object_info(status &
            , TRIM(XVAX(i)%io%channel), TRIM(XVAX(i)%io%object) &
            , reprid = reprid)
       CALL channel_halt(substr, status)
       IF ((reprid /= GP_3D_MID) .AND. (reprid /= GP_3D_INT)) THEN
          CALL error_bi('REPRESENTATION (AXIS) MUST BE GP_3D_MID/INT' ,substr)
       ENDIF
       XVAX(i)%reprid_obj = reprid
    END DO

    ! .................................................................

    ! OBJECTS TO BE TRANSFORMED
    DO i=1, NTRA
       IF (XTRA(i)%domain_idx /= dom_curid) CYCLE
       
       IF (p_parallel_io) THEN
          WRITE(*,*) ' CHANNEL / OBJECT ' &
               , TRIM(XTRA(i)%io%cha),' / ',TRIM(XTRA(i)%io%obj) &
               , ' FOR NEW OBJECT ', TRIM(XTRA(i)%io%name)
       ENDIF
       CALL get_channel_object(status &
            , TRIM(XTRA(i)%io%cha), TRIM(XTRA(i)%io%obj) &
            , p3=XTRA(i)%inp)
       !CALL channel_halt(substr, status)
       IF (status /= 0) THEN
          CALL warning_bi(&
               '  '//TRIM(XTRA(i)%io%cha)//' - '//TRIM(XTRA(i)%io%obj)//&
               &' not found ... skipping', substr)
          NULLIFY(XTRA(i)%inp)
          CYCLE
       ELSE
          CALL get_channel_object_info(status &
               , TRIM(XTRA(i)%io%cha), TRIM(XTRA(i)%io%obj) &
               , reprid = reprid)
          CALL channel_halt(substr, status)
          XTRA(i)%reprid = reprid
       END IF
       !
       ! CHECK CONFORMITY OF REPRESENTATIONS
       axis_id = XTRA(i)%axid
       !
       IF ((XTRA(i)%reprid /= GP_3D_MID) .AND. &
            (XTRA(i)%reprid /= GP_3D_INT)) THEN
          CALL error_bi('REPRESENTATION (OBJECT) MUST BE GP_3D_MID/INT' ,substr)
       ENDIF
       ! XVAX(axis_id)%reprid_obj is EITHER GP_3D_MID OR GP_3D_INT
       ! (see above)
       IF (XTRA(i)%reprid /= XVAX(axis_id)%reprid_obj) THEN
          IF (XTRA(i)%reprid == GP_3D_MID) THEN
             XTRA(i)%ipo = IPO_MID2INT
          ENDIF
          IF (XTRA(i)%reprid == GP_3D_INT) THEN
             XTRA(i)%ipo = IPO_INT2MID
          END IF
       ELSE
          XTRA(i)%ipo = IPO_NONE
       ENDIF
    END DO

    CALL end_message_bi(modstr, 'COUPLING TO EXTERNAL OBJECTS', substr)

    ! .................................................................

    CALL start_message_bi(modstr, 'CHANNEL OBJECT DEFINITIONS', substr)

    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    DO i=1, NTRA
       IF (XTRA(i)%domain_idx /= dom_curid) CYCLE
       
       ! SKIP, OF INPUT IS NOT PRESENT
       IF (.NOT. ASSOCIATED(XTRA(i)%inp)) CYCLE
       !
       axis_id = XTRA(i)%axid
       CALL new_channel_object(status, modstr, TRIM(XTRA(i)%io%name) &
            , p3=XTRA(i)%out, reprid = XVAX(axis_id)%reprid)
       CALL channel_halt(substr, status)
       !
       ! copy all attributes from original object to new object
       CALL copy_channel_object_atts(status &
            , TRIM(XTRA(i)%io%cha), TRIM(XTRA(i)%io%obj)  &
            , modstr, TRIM(XTRA(i)%io%name) )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr, TRIM(XTRA(i)%io%name) &
            , 'missing_value', r=XTRA(i)%io%missing, loverwrite=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr, TRIM(XTRA(i)%io%name) &
            , '_FillValue', r=XTRA(i)%io%missing, loverwrite=.TRUE.)
       CALL channel_halt(substr, status)
       !
    END DO

    CALL end_message_bi(modstr, 'CHANNEL OBJECT DEFINITIONS', substr)

  END SUBROUTINE vaxtra_init_coupling
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_global_end

    USE messy_main_grid_def_mem_bi,ONLY: nproma, npromz, ngpblks, nlev
    USE messy_main_tools,       ONLY: iso2ind, ind2val
    USE messy_main_channel_mem, ONLY: dom_curid
    
    IMPLICIT NONE
    INTRINSIC :: MINVAL, MAXVAL, LOG

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'vaxtra_global_end'
    INTEGER                             :: i, j, jk, jrow, kproma, jp, mlev
    REAL(DP), DIMENSION(nproma)         :: iso, f
    INTEGER,  DIMENSION(nproma)         :: kindex
    REAL(DP), DIMENSION(nproma, nlev+1) :: col, acol

    ! ASSUME THAT MANY OBJECTS ARE TRANSFORMED TO THE SAME VERTICAL AXIS
    ! -> CALCULATE LEVEL INDICES AND WEIGHTS (IN AXIS FIELD) ONLY 
    !    ONCE PER AXIS AND LEVEL

    axis_loop: DO i=1, NVAX
       IF (XVAX(i)%domain_idx /= dom_curid) CYCLE
       
       taxis_level: DO jk=1, XVAX(i)%io%dimlen

          ! SCALE AXIS VALUES (!)
          iso(:) = XVAX(i)%io%axis(jk) / XVAX(i)%io%scal
          IF (XVAX(i)%io%llog) THEN
             iso(:) = LOG(iso(:))
          END IF
             
          row_loop: DO jrow=1, ngpblks
#ifndef CESM1
             IF (jrow == ngpblks) THEN
                kproma = npromz
             ELSE
                kproma = nproma
             END IF
#else
             kproma = npromz(jrow)
#endif

             mlev = SIZE(XVAX(i)%ptr, _IZ_XYZ__)  ! nlev or nlev+1
             acol(1:kproma,1:mlev) = XVAX(i)%ptr(_RI_XYZ__(1:kproma,jrow,:))
             IF (XVAX(i)%io%llog) THEN
                acol(1:kproma,1:mlev) = LOG(acol(1:kproma,1:mlev))
             END IF

             CALL iso2ind(kproma, acol(1:kproma,1:mlev)              &
                  , iso(1:kproma) , kindex(1:kproma), f(1:kproma))

             object_loop: DO j=1, NTRA
                IF (XTRA(j)%domain_idx /= dom_curid) CYCLE
                
                ! CHECK, IF AXIS MATCHES ...
                IF (XTRA(j)%axid /= i) CYCLE   ! NO
                !
                ! CHECK, IF INPUT CHANNEL OBJECT IS PRESENT
                IF (.NOT. ASSOCIATED(XTRA(j)%inp)) CYCLE ! NO
                !
                ! INITIALIZE WITH MISSING VALUE
                ! (:,jk,jrow) or (:,jrow,jk)
                XTRA(j)%out(_RI_XYZ__(:,jrow,jk)) = XTRA(j)%io%missing
                !
                mlev = SIZE(XTRA(j)%inp, _IZ_XYZ__)
                !
                ! PER-INTERPOLATION
                SELECT CASE(XTRA(j)%ipo)
                CASE(IPO_NONE)
                   col(1:kproma,1:mlev) = XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,:))
                CASE(IPO_INT2MID)
                   col(1:kproma,1:nlev) = 0.5_dp * ( &
                        XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,1:nlev)) + &
                        XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) )
                CASE(IPO_MID2INT)
                   ! 1,jrow
                   col(1:kproma,1)      = XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,1))
                    ! nlev,jrow
                   col(1:kproma,nlev+1) = XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,nlev))
                   col(1:kproma,2:nlev) = 0.5_dp * ( &
                        XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,1:nlev-1)) + &
                        XTRA(j)%inp(_RI_XYZ__(1:kproma,jrow,2:nlev)) )
                   ! reset number of available vertical levels
                   mlev = nlev + 1 ! ju_ak_20191030
                END SELECT
                !
                ! INTERPOLATION
                CALL ind2val(kproma, XTRA(j)%out(_RI_XYZ__(1:kproma,jrow,jk)) &
                     , col(1:kproma,1:mlev) &
                     , kindex(1:kproma), f(1:kproma) )
                ! 
                ! RESET TO MISSING VALUES BEYOND BOUNDARIES
                mlev = SIZE(XVAX(i)%ptr, _IZ_XYZ__)  ! nlev or nlev+1
                DO jp=1, kproma
                   IF (iso(jp) < MINVAL(acol(jp,1:mlev))) THEN
                      XTRA(j)%out(_RI_XYZ__(jp,jrow,jk)) = XTRA(j)%io%missing
                   END IF
                   IF (iso(jp) > MAXVAL(acol(jp,1:mlev))) THEN
                      XTRA(j)%out(_RI_XYZ__(jp,jrow,jk)) = XTRA(j)%io%missing
                   END IF
                END DO
                !
             END DO object_loop

          END DO row_loop
          
       END DO taxis_level

    END DO axis_loop

  END SUBROUTINE vaxtra_global_end
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_free_memory

    IMPLICIT NONE

    INTEGER :: i

    DO i=1, NVAX
       NULLIFY(XVAX(i)%ptr)
    END DO
    DEALLOCATE(XVAX)

    DO i=1, NTRA
       NULLIFY(XTRA(i)%inp)
       NULLIFY(XTRA(i)%out)
    END DO
    DEALLOCATE(XTRA)
    
  END SUBROUTINE vaxtra_free_memory
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE vaxtra_read_nml_cpl(status, iou)

    ! vaxtra MODULE ROUTINE (INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to basemodel
    !
    ! Author: Patrick Joeckel, DLR, Dec 2016

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'vaxtra_read_nml_cpl'

    NAMELIST /CPL/ VAX, TRA

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE vaxtra_read_nml_cpl
! ----------------------------------------------------------------------
 
! **************************************************************************
END MODULE messy_vaxtra_si
! **************************************************************************
