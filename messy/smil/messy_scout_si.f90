#include "messy_main_ppd_bi.inc"

! ***********************************************************************
MODULE MESSY_SCOUT_SI
! ***********************************************************************

  ! Selectable Column OUTput
  ! MODULE FOR HF-OUTPUT OF TRACERS/CHANNEL OBJECTS (VERTICAL COLUMNS)
  ! AT SELECTED GEOGRAPHIC LOCATIONS
  !
  ! INTERFACE FOR ECHAM5 (MESSy/SMIL)
  !
  ! Author: Patrick Joeckel, MPICH, December 2003
  !

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_main_constants_mem, ONLY: STRLEN_XLONG, STRLEN_VLONG
  USE messy_scout

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TRIM, ASSOCIATED, NULL

  ! MAX. NUMBER OF LOCATIONS IN NAMELIST
  INTEGER, PARAMETER :: NMAXLOC        = 500

  TYPE IO_SCOUT
     CHARACTER(LEN=STRLEN_VLONG)  :: name = ''
     REAL(DP)                     :: lat  = -999.0_DP
     REAL(DP)                     :: lon  = -999.0_DP
     CHARACTER(LEN=STRLEN_XLONG)  :: str  = ''
  END TYPE IO_SCOUT

  TYPE LOCATION
     TYPE(IO_SCOUT)        :: loc
     INTEGER               :: n = 0              ! no. of requests/objects
     CHARACTER(LEN=STRLEN_CHANNEL), DIMENSION(:) &
          , POINTER        :: channel => NULL()  ! channel name
     CHARACTER(LEN=STRLEN_OBJECT),  DIMENSION(:) &
          , POINTER        :: object => NULL()   ! object name
     TYPE(PTR_3D_ARRAY),           DIMENSION(:) &
          , POINTER        :: field  => NULL()   ! ptr to channel object data
     INTEGER, DIMENSION(:)                       &
          , POINTER        :: reprid =>NULL()    ! representation ID
     TYPE(PTR_1D_ARRAY),           DIMENSION(:) &
          , POINTER        :: column => NULL()   ! column data
     LOGICAL               :: ldo = .FALSE.      ! location in domain ?
     INTEGER               :: pe = 0             ! on which CPU
     INTEGER               :: jp = 0             ! which column
     INTEGER               :: jrow = 0           ! which row
     INTEGER               :: domain_idx = 0
  END TYPE LOCATION

  TYPE(IO_SCOUT), DIMENSION(NMAXLOC), SAVE :: LOC
  INTEGER                                  :: NLOC
  TYPE(LOCATION), DIMENSION(:), POINTER, SAVE :: XLOC => NULL()

  ! NEW REPRESENTATIONS
  INTEGER, SAVE                     :: SCT_GP_1D_COLUMN_MID
  INTEGER, SAVE                     :: SCT_GP_1D_COLUMN_INT

  PUBLIC :: scout_initialize
  PUBLIC :: scout_init_memory
  PUBLIC :: scout_init_coupling
  PUBLIC :: scout_write_output
  PUBLIC :: scout_free_memory

  !PRIVATE :: scout_read_nml_cpl

CONTAINS

! ---------------------------------------------------------------------
  SUBROUTINE scout_initialize

    ! scout MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_channel_bi, ONLY: n_dom
    USE messy_main_tools,      ONLY: str2chob, find_next_free_unit &
                                   , domains_from_string

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, LEN_TRIM
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'scout_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j
    !
    CHARACTER(LEN=32)                        :: varname = '  '
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                                  :: nd, num

    ALLOCATE(XLOC(NMAXLOC))

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL scout_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('',substr)
    END IF

    CALL start_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       ! GET NUMBER OF LOCATIONS
       NLOC = 1
       ! COPY DATA AND PARSE STR
       DO i=1, NMAXLOC
          IF (TRIM(LOC(i)%name) == '') CYCLE
          IF ((LOC(i)%lon < -180.0_dp) .OR. (LOC(i)%lon > 360.0_dp)) CYCLE
          CALL domains_from_string(status, LOC(i)%name, n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          NLOC = NLOC + num
       END DO
       NLOC = NLOC -1
    END IF
    CALL p_bcast(NLOC, p_io)
    ALLOCATE(XLOC(NLOC))

    IF (p_parallel_io) THEN
       ! GET NUMBER OF LOCATIONS
       NLOC = 1
       ! COPY DATA AND PARSE STR
       location_loop: DO i=1, NMAXLOC
          IF (TRIM(LOC(i)%name) == '') CYCLE

          CALL domains_from_string(status,LOC(i)%name,n_dom,num &
                ,varname,dnums=domnum )
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          IF (LEN_TRIM(ADJUSTL(varname)) > 9) &
               CALL error_bi('location name too long (max. 5(+4) characters)' &
               , substr)
          
          domain_loop: DO nd = 1, SIZE(domnum)
             XLOC(NLOC)%loc%name   = TRIM(ADJUSTL(varname))
             XLOC(NLOC)%domain_idx = domnum(nd)

             WRITE(*,*) 'LOCATION ''',TRIM(XLOC(NLOC)%loc%name),''' ...'
             !
             IF ((LOC(i)%lat < -90.0_dp) .OR. (LOC(i)%lat > 90.0_dp)) THEN
                WRITE(*,*) '... LATITUDE out of range [-90,90] : ' &
                     ,LOC(i)%lat,''' ... skipping'
                CYCLE
             ELSE
                WRITE(*,*) '    LATITUDE : ',LOC(i)%lat
             END IF
             XLOC(NLOC)%loc%lat  = LOC(i)%lat
             !
             IF ((LOC(i)%lon < -180.0_dp) .OR. (LOC(i)%lon > 360.0_dp)) THEN
                WRITE(*,*) '... LONGITUDE out of range [-180, 360] : ' &
                     ,LOC(i)%lon,''' ... skipping'
                CYCLE
             ELSE
                WRITE(*,*) '    LONGITUDE: ',LOC(i)%lon
             END IF
             XLOC(NLOC)%loc%lon  = LOC(i)%lon
             !
             XLOC(NLOC)%loc%str  = LOC(i)%str
             ! SET PRELIMINARY %n
             CALL str2chob(status, XLOC(NLOC)%loc%str, XLOC(NLOC)%n &
                  , XLOC(NLOC)%channel, XLOC(NLOC)%object)
             IF (status /= 0) THEN
                WRITE(*,*) '... ERROR IN STRING ... skipping'
                CYCLE
             END IF
             !
             IF (XLOC(NLOC)%n == 0) THEN
                WRITE(*,*) '... EMPTY OUTPUT LIST ... skipping'
                CYCLE
             ELSE
                WRITE(*,*) '    REQUESTS : ',XLOC(NLOC)%n
             END IF
             !
             WRITE(*,'(1x,a16,1x,a32)') 'CHANNEL','OBJECT(S)'
             WRITE(*,'(1x,a16,1x,a32)') '-------','---------'
             DO j=1, XLOC(NLOC)%n 
                WRITE(*,'(1x,a16,1x,a32)') TRIM(XLOC(NLOC)%channel(j)) &
                     , TRIM(XLOC(NLOC)%object(j))
             END DO
             !
             ! NEXT LOCATION
             NLOC = NLOC + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop 
          DEALLOCATE(domnum); NULLIFY(domnum)
       END DO location_loop
       NLOC = NLOC - 1 
       IF (NLOC > SIZE(XLOC)) CALL error_bi('error in namelist parsing', substr)
    END IF

    ! BROADCAST RESULTS
    DO i=1, NLOC
       CALL p_bcast(XLOC(i)%loc%name, p_io)
       CALL p_bcast(XLOC(i)%loc%lat,  p_io)
       CALL p_bcast(XLOC(i)%loc%lon,  p_io)
       CALL p_bcast(XLOC(i)%loc%str,  p_io)
       CALL p_bcast(XLOC(i)%domain_idx, p_io)
       CALL p_bcast(XLOC(i)%n,       p_io)
       ! SPACE ON NON-IO CPUs
       ! SET PRELIMINARY %n
       IF (.NOT.p_parallel_io) THEN
          ALLOCATE(XLOC(i)%channel(XLOC(i)%n))
          ALLOCATE(XLOC(i)%object(XLOC(i)%n))
       END IF
       DO j=1, XLOC(i)%n
          CALL p_bcast(XLOC(i)%channel(j),  p_io)
          CALL p_bcast(XLOC(i)%object(j), p_io)
       END DO
       !
       ! column, field, reprid, ldo, pe, jp, jrow 
       ! -> set in scout_init_coupling
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NLOC,' LOCATION(S) FOR OUTPUT INITIALIZED !'
    END IF

    CALL end_message_bi(modstr,'OUTPUT INITIALISATION ',substr)

  END SUBROUTINE scout_initialize
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE scout_init_memory

    ! SCOUT MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DIMID_LEV, DIMID_ILEV, DC_BC
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel_repr,     ONLY: new_representation, AUTO &
                                         , set_representation_decomp &
                                         , IRANK, PIOTYPE_COL

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'scout_init_memory'
    INTEGER                      :: status
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    ! -------------------

    ! NEW REPRESENTATIONS
    CALL new_representation(status, SCT_GP_1D_COLUMN_MID, 'SCT_GP_1D_COL_MID' &
         , rank = 1, link = 'x---', dctype = DC_BC                   &
         , dimension_ids = (/ DIMID_LEV /) &
         , ldimlen       = (/ AUTO  /)     &
         , axis = 'Z---'                   &
         )
    CALL channel_halt(substr, status)

    start(:,1) = 1
    cnt(:,1)   = nlev
    meml(:,1)  = 1
    memu(:,1)  = nlev
    
    CALL set_representation_decomp(status, SCT_GP_1D_COLUMN_MID &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -------------------

    CALL new_representation(status, SCT_GP_1D_COLUMN_INT, 'SCT_GP_1D_COL_INT' &
         , rank = 1, link = 'x---', dctype = DC_BC                       &
         , dimension_ids = (/ DIMID_ILEV /) &
         , ldimlen       = (/ AUTO  /)     &
         , axis = 'Z---'                   &
         )
    CALL channel_halt(substr, status)

    start(:,1) = 1
    cnt(:,1)   = nlev+1
    meml(:,1)  = 1
    memu(:,1)  = nlev+1
    
    CALL set_representation_decomp(status, SCT_GP_1D_COLUMN_INT &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    ! -------------------

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL end_message_bi(modstr, 'NEW REPRESENTATIONS', substr)

  END SUBROUTINE scout_init_memory
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE scout_init_coupling

    ! SCOUT MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! define specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, warning_bi, info_bi
    USE messy_main_transform_bi,     ONLY: locate_in_decomp
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_GP, DC_BC                 &
                                         , SCALAR, GP_3D_MID, GP_3D_INT &
                                         , GP_3D_1LEV, GP_2D_HORIZONTAL &
                                         , DIMID_LON, DIMID_LAT
    USE messy_main_channel,          ONLY: new_channel, new_channel_object   &
                                         , new_attribute, get_channel_object &
                                         , get_attribute, get_channel_info   &
                                         , get_channel_object_info
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_channel_repr,     ONLY: get_representation_id, IRANK, AUTO &
                                         , PIOTYPE_COL, get_representation_info&
                                         , new_representation &
                                         , set_representation_decomp
    USE messy_main_channel_dimensions, ONLY: get_dimension_info
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG, STRLEN_MEDIUM
    USE messy_main_tools,              ONLY: match_wild

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'scout_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: reprid
    INTEGER                      :: reprid_new
    INTEGER                      :: i, jr, n, nptr, jo
    CHARACTER(LEN=STRLEN_ULONG)  :: charatt
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER :: ONAMES => NULL()

    INTEGER                      :: il
    INTEGER                      :: rank
    INTEGER, DIMENSION(IRANK)    :: dimids
    INTEGER                      :: dctype 
    INTEGER                      :: aaxlen = -99 ! length of additional axis
    CHARACTER(LEN=STRLEN_MEDIUM) :: dimname  = ''! NAME OF additional DIMENSION
    INTEGER                      :: SCT_GP_1D_COLUMN_AAX
    INTEGER                      :: irx, iry, iraax
    LOGICAL                      :: lstatic
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    IF (p_parallel_io) THEN
       WRITE(*,*)
       WRITE(*,'(4x,a6,1x,a8,1x,a9,1x,a9,1x,4(a7))') &
            'DOMAIN' &
            , 'LOCATION',' LATITUDE','LONGITUDE' &
            , '     PE','    JP','  JROW','STATUS'
       WRITE(*,'(4x,a6,1x,a8,1x,a9,1x,a9,1x,4(a7))')  &
            '------' &
            , '--------','---------','---------' &
            , '-------','------','------','------'
    END IF

    location_loop: DO i=1, NLOC

       IF (XLOC(i)%domain_idx /= dom_curid) CYCLE

       CALL locate_in_decomp(status, XLOC(i)%loc%lon, XLOC(i)%loc%lat &
            , XLOC(i)%pe, XLOC(i)%jp, XLOC(i)%jrow)

       IF (p_parallel_io) THEN
          WRITE(*,'(4x,i6,1x,a8,1x,f9.4,1x,f9.4,1x,4i7)') &
               XLOC(i)%domain_idx  &
               , TRIM(XLOC(i)%loc%name), XLOC(i)%loc%lat,XLOC(i)%loc%lon &
               , XLOC(i)%pe, XLOC(i)%jp, XLOC(i)%jrow , status
       END IF

       IF (status /= 0) THEN
          XLOC(i)%ldo = .FALSE.
          !                                     
          CYCLE  ! -----------------------------|
       ELSE
          XLOC(i)%ldo = .TRUE.
       END IF

       ! OPEN ONE CHANNEL PER LOCATION
       CALL new_channel(status, modstr//'_'//TRIM(XLOC(i)%loc%name))
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//TRIM(XLOC(i)%loc%name) &
            , 'latitude', r=XLOC(i)%loc%lat)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//TRIM(XLOC(i)%loc%name) &
            , 'longitude', r=XLOC(i)%loc%lon)
       CALL channel_halt(substr, status)

       ! CALCULATE REQUIRED NUMBER OF POINTERS
       nptr = 0
       DO jr=1,  XLOC(i)%n  ! LOOP OVER REQUESTS
          CALL get_channel_info(status, TRIM(XLOC(i)%channel(jr)) &
               , ONAMES = ONAMES)
          IF (status /= 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CALL channel_halt(substr, status)
          ELSE
             CALL warning_bi(' ... channel '&
                  &//TRIM(XLOC(i)%channel(jr))//' does not exist.', substr)
             CYCLE
          END IF
          DO jo = 1, SIZE(ONAMES)
             IF ( match_wild(TRIM(XLOC(i)%object(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1
             END IF
          END DO
       END DO
       !
       IF (ASSOCIATED(ONAMES)) THEN
          DEALLOCATE(ONAMES)
          NULLIFY(ONAMES)
       END IF

       ! ALLOCATE SPACE FOR POINTER ARRAYS (FIELD, COLUMN OBJECT)
       ALLOCATE(XLOC(i)%field(nptr))    ! POINTER TO CHANNEL OBJECT
       ALLOCATE(XLOC(i)%reprid(nptr))   ! REPRESENTATION OF CHANNEL OBJECT
       ALLOCATE(XLOC(i)%column(nptr))   ! POINTER TO COLUMN CHANNEL OBJECT

       ! RESET NUMBER OF REQUESTS
       n = XLOC(i)%n    ! SAVE NUMBER OF REQUESTS FOR LOOP BELOW       
       XLOC(i)%n = nptr ! NUMBER OF MATCHING OBJECTS
       ! (SO FAR XLOC(i)%n CONTAINED NUMBER OF REQUESTS)

       ! ADD REQUESTED CHANNEL OBJECTS
       nptr = 0
       request_loop: DO jr=1, n  ! LOOP OVER REQUESTS

          ! GET ALL POTENTIAL OBJECTS
          CALL get_channel_info(status, TRIM(XLOC(i)%channel(jr)) &
               , ONAMES = ONAMES)
          IF (status == 3003) THEN ! CHANNEL (NAME) DOES NOT EXIST
             CYCLE
          ELSE
             CALL channel_halt(substr, status)
          END IF

          object_loop: DO jo = 1, SIZE(ONAMES)

             ! CHECK IF POTENTIAL OBJECT IS MATCHING
             IF ( match_wild(TRIM(XLOC(i)%object(jr)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1      ! NEXT POINTER INDEX
             ELSE
                CYCLE
             END IF
             
             CALL info_bi('          '//TRIM(XLOC(i)%channel(jr))//&
                  &' - '//TRIM(ONAMES(jo)) &
                  , substr)

             CALL get_channel_object(status &
                  , TRIM(XLOC(i)%channel(jr)), TRIM(ONAMES(jo)) &
                  , p3=XLOC(i)%field(nptr)%ptr )
             IF (status /= 0) THEN
                CALL warning_bi( &
                     '            ... channel object not found ... skipping' &
                     , substr)
                NULLIFY(XLOC(i)%column(nptr)%ptr)
                NULLIFY(XLOC(i)%field(nptr)%ptr)
                CYCLE
             END IF

             CALL get_channel_object_info(status &
                  , TRIM(XLOC(i)%channel(jr)), TRIM(ONAMES(jo)) &
                  , reprid=reprid, lstatic=lstatic)
             CALL channel_halt(substr, status)
             !
             XLOC(i)%reprid(nptr) = reprid
             !
             IF (reprid == GP_3D_MID) THEN
                CALL get_representation_id(status, 'SCT_GP_1D_COL_MID' &
                     , reprid=reprid_new)
                CALL channel_halt(substr, status)
             ELSE IF (reprid == GP_3D_INT) THEN
                CALL get_representation_id(status, 'SCT_GP_1D_COL_INT' &
                     , reprid=reprid_new)
                CALL channel_halt(substr, status)
             ELSE IF (reprid == GP_3D_1LEV) THEN
                reprid_new = SCALAR
             ELSE IF (reprid == GP_2D_HORIZONTAL) THEN
                reprid_new = SCALAR
             ELSE
                CALL get_representation_info(status, '' &
                  , ID=reprid, rank=rank, dctype=dctype, dim_ids=dimids)
                CALL channel_halt(substr, status)

                IF (rank == 3 .AND. dctype==DC_GP) THEN
                   irx = -99
                   iry = -99
                   iraax = -99
                   ! search for longitude and latitude rank
                   DO il = 1, IRANK
                      IF (dimids(il) == DIMID_LON) THEN
                         irx = il
                      ELSE IF (dimids(il) == DIMID_LAT) THEN
                         iry = il
                      ELSE IF (dimids(il) > 0) THEN
                         ! index of additional z or n axis
                         iraax = dimids(il)
                      END IF
                   END DO
                   IF (irx > 0 .AND. iry > 0) THEN

                      CALL get_dimension_info(status, iraax &
                           , len=aaxlen, name=dimname)

                      CALL  get_representation_info(status &
                           , 'SCT_GP_1D_COL_'//TRIM(ADJUSTL(dimname)) &
                           , Id = reprid_new)

                      IF (status /= 0) THEN
                         ! NEW REPRESENTATIONS
                         nseg = 1
                         ALLOCATE(start(nseg,IRANK))
                         ALLOCATE(cnt(nseg,IRANK))
                         ALLOCATE(meml(nseg,IRANK))
                         ALLOCATE(memu(nseg,IRANK))
                         CALL new_representation(status, SCT_GP_1D_COLUMN_AAX &
                              , 'SCT_GP_1D_COL_'//TRIM(ADJUSTL(dimname))       &
                              , rank = 1, link = 'x---', dctype = DC_BC       &
                              , dimension_ids = (/ iraax /) &
                              , ldimlen       = (/ AUTO  /)     &
                              , axis = 'Z---'                   &
                              )
                         CALL channel_halt(substr, status)
                         
                         start(:,:)  = 1
                         cnt(:,1)    = aaxlen
                         meml(:,:)   = 1
                         memu(:,1)   = aaxlen
                         cnt(:,2:4)  = 1
                         memu(:,2:4) = 1
                         
                         CALL set_representation_decomp(status &
                              , SCT_GP_1D_COLUMN_AAX &
                              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
                         CALL channel_halt(substr, status)
                         ! -------------------
                         
                         DEALLOCATE(start) ; NULLIFY(start)
                         DEALLOCATE(cnt)   ; NULLIFY(cnt)
                         DEALLOCATE(meml)  ; NULLIFY(meml)
                         DEALLOCATE(memu)  ; NULLIFY(memu)

                         reprid_new =  SCT_GP_1D_COLUMN_AAX
                      END IF
                   ELSE
                      CALL warning_bi( &
                           '            ... representation not supported ... skipping' &
                           , substr)
                      NULLIFY(XLOC(i)%column(nptr)%ptr)
                      NULLIFY(XLOC(i)%field(nptr)%ptr)
                      CYCLE
                   END IF
                       
                ELSE
                   CALL warning_bi( &
                  '            ... representation not supported ... skipping' &
                  , substr)
                   NULLIFY(XLOC(i)%column(nptr)%ptr)
                   NULLIFY(XLOC(i)%field(nptr)%ptr)
                   CYCLE
                END IF
             END IF

             CALL new_channel_object(status &
                  , modstr//'_'//TRIM(XLOC(i)%loc%name) &
                  , TRIM((XLOC(i)%channel(jr)))//'_'//TRIM(ONAMES(jo))  &
                  , p1=XLOC(i)%column(nptr)%ptr &
                  , reprid=reprid_new, lstatic=lstatic)
             CALL channel_halt(substr, status)
             ! COPY ATTRIBUTES
             CALL get_attribute(status &
                  , TRIM(XLOC(i)%channel(jr)), TRIM(ONAMES(jo)) &
                  , 'long_name', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XLOC(i)%loc%name) &
                     , TRIM((XLOC(i)%channel(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'long_name', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF
             CALL get_attribute(status &
                  , TRIM(XLOC(i)%channel(jr)), TRIM(ONAMES(jo)) &
                  , 'units', c=charatt)
             IF (status == 0) THEN
                CALL new_attribute(status &
                     , modstr//'_'//TRIM(XLOC(i)%loc%name) &
                     , TRIM((XLOC(i)%channel(jr)))//'_'//TRIM(ONAMES(jo))  &
                     , 'units', c=TRIM(charatt))
                CALL channel_halt(substr, status)
             END IF

          END DO object_loop

       END DO request_loop
       IF (nptr /= XLOC(i)%n) THEN
          CALL error_bi('something went wrong with pointer counting', substr)
       END IF

       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)

    END DO location_loop

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE scout_init_coupling
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE scout_write_output

    USE messy_main_mpi_bi,      ONLY: p_pe, p_bcast
    USE messy_main_channel_bi,  ONLY: GP_2D_HORIZONTAL
    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'scout_write_output'
    INTEGER                     :: i, j
    INTEGER                     :: p_pe_col, jp, jrow

    location_loop: DO i=1, NLOC
       
       IF (XLOC(i)%domain_idx /= dom_curid) CYCLE

       !                                            ^
       IF (.NOT. XLOC(i)%ldo) CYCLE  ! -------------|

       ! COLUMN IS LOCATED ON THIS PE
       p_pe_col = XLOC(i)%pe
       jp       = XLOC(i)%jp
       jrow     = XLOC(i)%jrow

       object_loop: DO j=1, XLOC(i)%n

          ! CHECK: CHANNEL AND/OR OBJECT NOT AVAIALABLE
          IF (.NOT.ASSOCIATED(XLOC(i)%field(j)%ptr)) CYCLE

          ! TRANSFER DATA ...
          ! ... PART 1: FROM CHANNEL OBJECT TO COLUMN CHANNEL
          !             (ON CPU, WHERE COLUMN IS LOCATED)
          IF (p_pe == p_pe_col) THEN
             IF (XLOC(i)%reprid(j) == GP_2D_HORIZONTAL) THEN
                XLOC(i)%column(j)%ptr(:) = &
                     XLOC(i)%field(j)%ptr(jp, jrow, 1)
             ELSE
                XLOC(i)%column(j)%ptr(:) = &
                     XLOC(i)%field(j)%ptr(_RI_XYZ__(jp,jrow,:))
             END IF
          END IF  ! p_pe == p_pe_col

          ! ... PART 2: FROM CPU, WHERE COLUMN IS LOCATED
          !             TO ALL OTHER CPUs
          CALL p_bcast(XLOC(i)%column(j)%ptr(:), p_pe_col)

       END DO object_loop

    END DO location_loop

  END SUBROUTINE scout_write_output
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE scout_free_memory

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: i

    location_loop: DO i=1, NLOC
       IF (ASSOCIATED(XLOC(i)%channel))  DEALLOCATE(XLOC(i)%channel)
       NULLIFY(XLOC(i)%channel)
       IF (ASSOCIATED(XLOC(i)%object)) DEALLOCATE(XLOC(i)%object)
       NULLIFY(XLOC(i)%object)
       !
       IF (ASSOCIATED(XLOC(i)%field))   DEALLOCATE(XLOC(i)%field)
       NULLIFY(XLOC(i)%field)
       IF (ASSOCIATED(XLOC(i)%reprid))   DEALLOCATE(XLOC(i)%reprid)
       NULLIFY(XLOC(i)%reprid)
       IF (ASSOCIATED(XLOC(i)%column))  DEALLOCATE(XLOC(i)%column)
       NULLIFY(XLOC(i)%column)
    END DO location_loop

    IF (ASSOCIATED(XLOC))  DEALLOCATE(XLOC)
    NULLIFY(XLOC)
    
  END SUBROUTINE scout_free_memory
! ---------------------------------------------------------------------
 
! ----------------------------------------------------------------------
  SUBROUTINE scout_read_nml_cpl(status, iou)

    ! scout MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'scout_read_nml_cpl'

    NAMELIST /CPL/ LOC


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

  END SUBROUTINE scout_read_nml_cpl
! ----------------------------------------------------------------------

! ***********************************************************************
END MODULE MESSY_SCOUT_SI
! ***********************************************************************
