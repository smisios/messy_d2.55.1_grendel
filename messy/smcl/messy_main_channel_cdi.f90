! **********************************************************************
MODULE messy_main_channel_cdi
! **********************************************************************

  ! MESSY DATA TRANSFER AND EXPORT INTERFACE (MEMORY MANAGEMENT)
  !
  ! Author: Bastian Kern, DLR, May 2014
  ! based on netcdf implementation by Patrick Joeckel, MPICH, May 2005
  ! adaption from ICON mo_name_list_output(_init)

  USE messy_main_channel

#ifdef HAVE_CDI

!!$  USE mo_cdi_constants
!!$  USE mo_cdi

  IMPLICIT NONE
  PRIVATE

!  INCLUDE 'cdi.inc'
#include <cdi.inc>

  INTEGER, SAVE :: NCOUT_PREC =  DATATYPE_FLT32 ! default

  PUBLIC :: ch_cdi_init_rst
  PUBLIC :: ch_cdi_init_io
  PUBLIC :: ch_cdi_write_header
  PUBLIC :: ch_cdi_write_time
  PUBLIC :: ch_cdi_write_data
#endif
  PUBLIC :: ch_cdi_finish_io
  !
  ! op_bk_20140506+
  ! read capabilities not implemented yet,
  ! cdi restart file can not be used, use netcdf library restart files (2)
  ! disabled below, so disable here for now
!   PUBLIC :: ch_cdi_read_data
  ! op_bk_20140506-
  !
  !PRIVATE :: cdi_write_attribute_list
  !PRIVATE :: cdi_define_dimvar_list
  !PRIVATE :: cdi_write_dimvar_list
  !PRIVATE :: cdi_check_attributes
  !PRIVATE :: handle_cdi_error

CONTAINS

  ! -------------------------------------------------------------------
  ! PUBLIC SUBROUTINES
  ! -------------------------------------------------------------------

#ifdef HAVE_CDI

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_init_rst(status, fname, att)

    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: fname
    TYPE(t_attribute_list), POINTER     :: att 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ch_cdi_init_rst'
    INTEGER                     :: fileID

    fileID = streamOpenRead(TRIM(fname))
    CALL handle_cdi_error(fileID, status, substr)
    IF (status /= 0) RETURN

    CALL cdi_check_attributes(status, fileID, CDI_GLOBAL, att) 
    IF (status /= 0) RETURN

    CALL streamClose(fileID)

  END SUBROUTINE ch_cdi_init_rst
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_init_io(status, IOMODE, channel, AMODE, att, lclose)

    USE messy_main_channel_attributes, ONLY: t_attribute_list

    IMPLICIT NONE

    INTRINSIC :: TRIM, PRESENT, IOR

    ! I/O
    INTEGER,                      INTENT(OUT) :: status
    INTEGER,                      INTENT(IN)  :: IOMODE
    TYPE(t_channel),              POINTER     :: channel ! INTENT(INOUT)    
    INTEGER,                      INTENT(IN)  :: AMODE
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att     ! INTENT(IN)
    LOGICAL,           INTENT(IN),   OPTIONAL :: lclose

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER               :: substr = 'ch_cdi_init_io'
    INTEGER                                   :: fileID, vlistID
    LOGICAL                                   :: lexist
    INTEGER                                   :: nc_cmode
    LOGICAL, SAVE                             :: lfirst = .TRUE.
    LOGICAL                                   :: llclose

    IF (PRESENT(lclose)) THEN
       llclose = lclose
    ELSE
       llclose = .FALSE.
    END IF

    IF (lfirst) THEN
       SELECT CASE(OUT_PREC(FTYPE_CDI_NC))
       CASE(1)
          NCOUT_PREC = DATATYPE_FLT32
       CASE(2)
          NCOUT_PREC = DATATYPE_FLT64
       CASE DEFAULT
          status = 4014 ! unknown precision flag for cdi
          RETURN
       END SELECT
       lfirst = .false.
    END IF

    SELECT CASE(AMODE)
    CASE(AMODE_READ)
       !
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          !
          status = 3205 ! no input of output files
          RETURN
          !
       CASE(IOMODE_RST)
          !
          IF (I_VERBOSE_LEVEL >= 1) THEN
             WRITE(*,*) substr,': OPENING  FILE: ' &
                  &   , TRIM(channel%int%fname(IOMODE))
          ENDIF
          !
          INQUIRE(file = TRIM(channel%int%fname(IOMODE)), exist = lexist)
          IF (.NOT.lexist) THEN
             IF (channel%int%lrestreq) THEN
                IF (channel%int%lign) THEN
                   ! ignore all missing restart fields
                   IF (I_VERBOSE_LEVEL >= 1) THEN
                      WRITE(*,*) substr,':       WARNING: REQUIRED RESTART ', &
                           &    'FILE NOT PRESENT ... ALL OBJECTS TO BE IGNORED'
                   END IF
                   ! reset: prevent from closing and reading
                   channel%int%fname(IOMODE) = ''
                   status = 0
                ELSE
                   status = 3206 ! restart file required but not present
                END IF
             ELSE
                IF (I_VERBOSE_LEVEL >= 1) THEN
                   WRITE(*,*) substr,':       WARNING: RESTART ', &
                        &    'FILE NOT PRESENT ... NOT REQUIRED'
                END IF
                ! reset: prevent from closing and reading
                channel%int%fname(IOMODE) = ''
                status = 0
             END IF
             RETURN
          END IF
          !
          fileID = streamOpenRead(TRIM(channel%int%fname(IOMODE)))
          CALL handle_cdi_error(fileID, status, substr)
          channel%int%cdi(IOMODE)%fileID = fileID
          !
          ! optional special attributes
          IF (PRESENT(att)) THEN
             CALL cdi_check_attributes(status, fileID, CDI_GLOBAL, att)
             IF (status /= 0) RETURN
          END IF
          CALL cdi_check_attributes(status, fileID, CDI_GLOBAL, GATT)
          IF (status /= 0) RETURN
          CALL cdi_check_attributes(status, fileID, CDI_GLOBAL, channel%att)
          IF (status /= 0) RETURN
          !
          IF (llclose) THEN
             CALL get_attribute(status, channel%att, 'channel_time_slo' &
                  , r=channel%int%tslo)
             IF (status /= 0) RETURN
             CALL streamClose(fileID)
             channel%int%cdi(IOMODE)%fileID = CDI_UNDEFID
          END IF
          ! 
          IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*)
          !
       END SELECT
       !
    CASE(AMODE_WRITE)
       !
       ! close 'old' file first, if still open
       IF (channel%int%cdi(IOMODE)%fileID /= CDI_UNDEFID) THEN
          CALL streamClose(channel%int%cdi(IOMODE)%fileID)
          channel%int%cdi(IOMODE)%fileID = CDI_UNDEFID
       END IF
       !
       ! open 'new' file
       IF (I_VERBOSE_LEVEL >= 2) THEN
          WRITE(*,*) substr,': CREATING FILE: ',TRIM(channel%int%fname(IOMODE))
       END IF
       nc_cmode = FILETYPE_NC2
       fileID = streamOpenWrite(TRIM(channel%int%fname(IOMODE)), nc_cmode)
       CALL handle_cdi_error(fileID, status, substr)
       IF (status /= 0) RETURN
       channel%int%cdi(IOMODE)%fileID = fileID
       !
       ! create a new vlist for channel objects (variables) handling
       ! op_bk_20170116+
       IF (channel%int%cdi(IOMODE)%vlistID == CDI_UNDEFID) THEN
       ! op_bk_20170116-
          vlistID = vlistCreate()
          CALL handle_cdi_error(vlistID, status, substr)
          IF (status /= 0) RETURN
          channel%int%cdi(IOMODE)%vlistID = vlistID
       END IF
       !
       ! op_bk_20170116+
       IF (channel%int%cdi(IOMODE)%taxisID == CDI_UNDEFID) THEN
          channel%int%cdi(IOMODE)%taxisID = taxisCreate(TAXIS_ABSOLUTE)
          CALL vlistDefTaxis(vlistID, channel%int%cdi(IOMODE)%taxisID)
       END IF
       ! op_bk_20170116-
    END SELECT

    status = 0

  END SUBROUTINE ch_cdi_init_io
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_write_header(status, IOMODE, channel, att)

    USE messy_main_channel_attributes,   ONLY: t_attribute_list
    USE mo_name_list_output_init,        ONLY: patch_info
    USE mo_util_uuid,                    ONLY: uuid2char
    USE mo_name_list_output_gridinfo,    ONLY: set_grid_info_netcdf
    USE mo_vertical_coord_table,         ONLY: vct
    USE mo_nonhydrostatic_config,        ONLY: ivctype
    USE mo_lnd_nwp_config,               ONLY: nlev_snow
    USE mo_run_config,                   ONLY: num_lev
    USE mo_impl_constants,               ONLY: zml_soil
    USE mo_name_list_output_gridinfo,    ONLY: GRID_INFO_BCAST

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, PRESENT

    ! I/O
    INTEGER,                INTENT(OUT)       :: status
    INTEGER,                INTENT(IN)        :: IOMODE
    TYPE(t_channel),        POINTER           :: channel    ! INTENT(INOUT)
    TYPE(t_attribute_list), POINTER, OPTIONAL :: att        ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER               :: substr = 'ch_cdi_write_header'
    TYPE(t_channel_object_list), POINTER      :: le
    TYPE(t_channel_object),      POINTER      :: object
    INTEGER                                   :: i
    INTEGER                                   :: i2nd  ! index in secondary data
    INTEGER                                   :: jsnd  ! output data type
    INTEGER                                   :: vlistID, varID
    LOGICAL                                   :: lskip
    INTEGER                                   :: CDINF_PREC
    INTEGER                                   :: idom
    INTEGER                                   :: gridtype
    CHARACTER(len=1)                          :: uuid_string(16)
    INTEGER                                   :: max_cell_connectivity
    INTEGER                                   :: k, nlev, nlevp1, nplev, nzlev, nilev, nzlevp1 &
         &                                     , znlev_soil, idate, itime
    INTEGER                                   :: nsoil_jsbach ! JSBACH number of soil layers
    REAL(dp), ALLOCATABLE                     :: levels_i(:), levels_m(:), p_lonlat(:)
    REAL(dp), ALLOCATABLE                     :: levels(:), lbounds(:), ubounds(:)
    INTEGER                                   :: iret
    ! op_bk_20141106+
    INTEGER                                   :: steptype
    ! op_bk_20141106-

    ! INIT
    vlistID = channel%int%cdi(IOMODE)%vlistID

    ! set precision
    IF (IOMODE == IOMODE_OUT) THEN
       CDINF_PREC = NCOUT_PREC
    ELSE
       CDINF_PREC = DATATYPE_FLT64
    END IF

    ! ---------------------------------------------------------
    ! write global attributes
    ! ---------------------------------------------------------
    ! - common to all channels
    CALL cdi_write_attribute_list(status, vlistID, CDI_GLOBAL, GATT &
         , CDINF_PREC)
    IF (status /= 0) RETURN
    ! - channel specific
    CALL cdi_write_attribute_list(status, vlistID, CDI_GLOBAL &
         , channel%att, CDINF_PREC)
    IF (status /= 0) RETURN
    ! - optional special atributes
    IF (PRESENT(att)) THEN
       CALL cdi_write_attribute_list(status, vlistID, CDI_GLOBAL, att &
            , CDINF_PREC)
       IF (status /= 0) RETURN
    END IF

    ! ---------------------------------------------------------
    ! define time axis
    ! ---------------------------------------------------------
    IF (IOMODE == IOMODE_OUT) THEN
       IF (channel%int%cdi(IOMODE)%taxisID == CDI_UNDEFID) THEN
          channel%int%cdi(IOMODE)%taxisID = taxisCreate(TAXIS_ABSOLUTE)
          CALL vlistDefTaxis(vlistID, channel%int%cdi(IOMODE)%taxisID)
       END IF
    END IF
    
    ! ---------------------------------------------------------
    ! define horizontal and vertical grids
    ! ---------------------------------------------------------
    ! - init
    gridtype = GRID_UNSTRUCTURED

    idom = channel%dom_id
    max_cell_connectivity = patch_info(idom)%max_cell_connectivity

    ! Cells
    channel%int%cdi(IOMODE)%CellGridID = gridCreate(gridtype, patch_info(idom)%cells%n_glb)
    CALL gridDefNvertex(channel%int%cdi(IOMODE)%CellGridID, max_cell_connectivity)
    !
    CALL gridDefXname(channel%int%cdi(IOMODE)%CellGridID, 'clon')
    CALL gridDefXlongname(channel%int%cdi(IOMODE)%CellGridID, 'center longitude')
    CALL gridDefXunits(channel%int%cdi(IOMODE)%CellGridID, 'radian')
    !
    CALL gridDefYname(channel%int%cdi(IOMODE)%CellGridID, 'clat')
    CALL gridDefYlongname(channel%int%cdi(IOMODE)%CellGridID, 'center latitude')
    CALL gridDefYunits(channel%int%cdi(IOMODE)%CellGridID, 'radian')
    !
    CALL uuid2char(patch_info(idom)%grid_uuid, uuid_string)
    CALL gridDefUUID(channel%int%cdi(IOMODE)%CellGridID, uuid_string)
    !
    CALL gridDefNumber(channel%int%cdi(IOMODE)%CellGridID, patch_info(idom)%number_of_grid_used)
    !
    CALL gridDefPosition(channel%int%cdi(IOMODE)%CellGridID, GRID_CELL)

    ! Vertices
    channel%int%cdi(IOMODE)%VertGridID = gridCreate(gridtype, patch_info(idom)%verts%n_glb)
    CALL gridDefNvertex(channel%int%cdi(IOMODE)%VertGridID, 9-max_cell_connectivity)
    !
    CALL gridDefXname(channel%int%cdi(IOMODE)%VertGridID, 'vlon')
    CALL gridDefXlongname(channel%int%cdi(IOMODE)%VertGridID, 'vertex longitude')
    CALL gridDefXunits(channel%int%cdi(IOMODE)%VertGridID, 'radian')
    !
    CALL gridDefYname(channel%int%cdi(IOMODE)%VertGridID, 'vlat')
    CALL gridDefYlongname(channel%int%cdi(IOMODE)%VertGridID, 'vertex latitude')
    CALL gridDefYunits(channel%int%cdi(IOMODE)%VertGridID, 'radian')
    !
    CALL uuid2char(patch_info(idom)%grid_uuid, uuid_string)
    CALL gridDefUUID(channel%int%cdi(IOMODE)%VertGridID, uuid_string)
    !
    CALL gridDefNumber(channel%int%cdi(IOMODE)%VertGridID, patch_info(idom)%number_of_grid_used)
    !
    CALL gridDefPosition(channel%int%cdi(IOMODE)%VertGridID, GRID_VERTEX)

    ! Edges
    channel%int%cdi(IOMODE)%EdgeGridID = gridCreate(gridtype, patch_info(idom)%edges%n_glb)
    CALL gridDefNvertex(channel%int%cdi(IOMODE)%EdgeGridID, 4)
    !
    CALL gridDefXname(channel%int%cdi(IOMODE)%EdgeGridID, 'elon')
    CALL gridDefXlongname(channel%int%cdi(IOMODE)%EdgeGridID, 'edge midpoint longitude')
    CALL gridDefXunits(channel%int%cdi(IOMODE)%EdgeGridID, 'radian')
    !
    CALL gridDefYname(channel%int%cdi(IOMODE)%EdgeGridID, 'elat')
    CALL gridDefYlongname(channel%int%cdi(IOMODE)%EdgeGridID, 'edge midpoint latitude')
    CALL gridDefYunits(channel%int%cdi(IOMODE)%EdgeGridID, 'radian')
    !
    CALL uuid2char(patch_info(idom)%grid_uuid, uuid_string)
    CALL gridDefUUID(channel%int%cdi(IOMODE)%EdgeGridID, uuid_string)
    !
    CALL gridDefNumber(channel%int%cdi(IOMODE)%EdgeGridID, patch_info(idom)%number_of_grid_used)
    !
    CALL gridDefPosition(channel%int%cdi(IOMODE)%EdgeGridID, GRID_EDGE)

    ! additional grid informations
    IF (patch_info(idom)%grid_info_mode == GRID_INFO_BCAST) THEN
       ! op_bk_20140818+
!        CALL set_grid_info_netcdf(channel%int%cdi(IOMODE)%CellGridID, patch_info(idom)%cells%grid_info)
!        CALL set_grid_info_netcdf(channel%int%cdi(IOMODE)%EdgeGridID, patch_info(idom)%edges%grid_info)
!        CALL set_grid_info_netcdf(channel%int%cdi(IOMODE)%VertGridID, patch_info(idom)%verts%grid_info)
       ! op_bk_20140818-
    END IF

    ! 4. add vertical grid descriptions
    !
    channel%int%cdi(IOMODE)%ZaxisID(:) = CDI_UNDEFID ! not all are set

    ! surface level
    channel%int%cdi(IOMODE)%ZaxisID(ZA_surface) = zaxisCreate(ZAXIS_SURFACE, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_surface), levels)
    DEALLOCATE(levels)

    nlev   = num_lev(idom)
    nlevp1 = num_lev(idom) + 1
    znlev_soil = SIZE(zml_soil)

    ! cloud base level
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_cloud_base) = zaxisCreate(ZAXIS_CLOUD_BASE, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_cloud_base), levels)
    DEALLOCATE(levels)

    ! cloud top level
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_cloud_top) = zaxisCreate(ZAXIS_CLOUD_TOP, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_cloud_top), levels)
    DEALLOCATE(levels)

    ! level of 0°C isotherm
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_isotherm_zero) = zaxisCreate(ZAXIS_ISOTHERM_ZERO, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_isotherm_zero), levels)
    DEALLOCATE(levels)

    ! reference_layer
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_reference) = zaxisCreate(ZAXIS_REFERENCE, nlev)
    ALLOCATE(lbounds(nlev), ubounds(nlev), levels(nlev))
    DO k = 1, nlev
       lbounds(k) = REAL(k,dp)
       levels(k)  = REAL(k,dp)
    END DO
    DO k = 2, nlevp1
       ubounds(k-1) = REAL(k,dp)
    END DO
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference), levels )  !necessary for NetCDF
    ! set numberOfVGridUsed
    ! Dependent on the algorithm chosen to generate the vertical grid (ivctype)
    ! op_bk_20140507+
!     CALL zaxisDefNumber(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference), get_numberOfVgridUsed(ivctype) )
    ! op_bk_20140507-
    !
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference),nlevp1)
    !
    ! UUID not yet available - write dummy UUID
    ! op_bk_20140507+
!     CALL zaxisDefUUID(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference), uuidOfVGrid_string ) !uuidOfVGrid
    ! op_bk_20140507-
    DEALLOCATE(lbounds, ubounds, levels)

    ! reference half levels
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half) = zaxisCreate(ZAXIS_REFERENCE, nlevp1)
    ALLOCATE(levels(nlevp1))
    DO k = 1, nlevp1
       levels(k) = REAL(k,dp)
    END DO
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half), levels)
    ! set numberOfVGridUsed
    ! Dependent on the algorithm chosen to generate the vertical grid (ivctype)
    ! op_bk_20140507+
!     CALL zaxisDefNumber(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half), get_numberOfVgridUsed(ivctype) )
    ! op_bk_20140507-
    !
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half),nlevp1)
    !
    ! UUID not yet available - write dummy UUID
    ! op_bk_20140507+
!     CALL zaxisDefUUID(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half), uuidOfVGrid_string ) !uuidOfVGrid
    ! op_bk_20140507-
    DEALLOCATE(levels)

    ! reference (special version for hhl)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl) = zaxisCreate(ZAXIS_REFERENCE, nlevp1)
    ALLOCATE(lbounds(nlevp1), ubounds(nlevp1), levels(nlevp1))
    DO k = 1, nlevp1
       lbounds(k) = REAL(k,dp)
       levels(k)  = REAL(k,dp)
    END DO
    ubounds(1:nlevp1) = 0._dp
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl), levels)  !necessary for NetCDF
    ! set numberOfVGridUsed
    ! Dependent on the algorithm chosen to generate the vertical grid (ivctype)
    ! op_bk_20140507+
!     CALL zaxisDefNumber(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl), get_numberOfVgridUsed(ivctype) )
    ! op_bk_20140507-
    !
    ! Define number of half levels for z-axis 
    CALL zaxisDefNlevRef(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl),nlevp1)
    !
    ! UUID not yet available - write dummy UUID
    ! op_bk_20140507+
!       CALL zaxisDefUUID(channel%int%cdi(IOMODE)%ZaxisID(ZA_reference_half_hhl), uuidOfVGrid_string ) !uuidOfVGrid
    ! op_bk_20140507-
    DEALLOCATE(lbounds, ubounds, levels)

    ! hybrid_layer
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid) = zaxisCreate(ZAXIS_HYBRID, nlev)
    ALLOCATE(lbounds(nlev), ubounds(nlev), levels(nlev))
    DO k = 1, nlev
       lbounds(k) = REAL(k,dp)
       levels(k)  = REAL(k,dp)
    END DO
    DO k = 2, nlevp1
       ubounds(k-1) = REAL(k,dp)
    END DO
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid), levels)   !necessary for NetCDF
    DEALLOCATE(lbounds, ubounds, levels)
    CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid), 2*nlevp1, vct(1:2*nlevp1))

    ! hybrid
    !
    ! Note: "ZAXIS_HYBRID_HALF" is deprecated and will soon be
    ! removed from the CDI (in principle its use should be simply
    ! replaced by ZAXIS_HALF, as long as lbounds and ubounds are set
    ! correctly).
    channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half) = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)
    ALLOCATE(levels(nlevp1))
    DO k = 1, nlevp1
       levels(k) = REAL(k,dp)
    END DO
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half), levels)
    DEALLOCATE(levels)
    CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half), 2*nlevp1, vct(1:2*nlevp1))

    ! hybrid (special version for hhl)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half_hhl) = zaxisCreate(ZAXIS_HYBRID_HALF, nlevp1)
    ALLOCATE(lbounds(nlevp1), ubounds(nlevp1), levels(nlevp1))
    DO k = 1, nlevp1
       lbounds(k) = REAL(k,dp)
       levels(k)  = REAL(k,dp)
    END DO
    ubounds(1:nlevp1) = 0._dp
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half_hhl), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half_hhl), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half_hhl), levels)   !necessary for NetCDF
    DEALLOCATE(lbounds, ubounds, levels)
    CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_hybrid_half_hhl), 2*nlevp1, vct(1:2*nlevp1))

    ! Define axis for output on mean sea level
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_meansea) = zaxisCreate(ZAXIS_MEANSEA, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.0_dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_meansea), levels)
    DEALLOCATE(levels)

    ! Define axes for soil model (DEPTH_BELOW_LAND)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land_p1) = &
         &                   zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil+1)
    ALLOCATE(levels(znlev_soil+1))
    levels(1) = 0._dp
    DO k = 1, znlev_soil
       levels(k+1) = REAL(zml_soil(k)*1000._dp,dp)  ! in mm
    END DO
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land_p1), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land_p1), "mm")
    DEALLOCATE(levels)

    !(depth_below_land_layer)
    !
    ! op_bk_20140507+
!     IF (ALLOCATED(lnd_jsbach_config)) THEN     ! For JSBACH
!        nsoil_jsbach = lnd_jsbach_config(idom)%nsoil
!        channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, nsoil_jsbach)
!        ALLOCATE(levels(nsoil_jsbach), lbounds(nsoil_jsbach), ubounds(nsoil_jsbach))
!        levels = 0._dp
!        levels(1) = REAL(1000._dp * lnd_jsbach_config(idom)%zlev_soil(1), dp)
!        DO k = 2,nsoil_jsbach
!           levels(k) = levels(k-1) + REAL(1000._dp * lnd_jsbach_config(idom)%zlev_soil(k), dp)
!        END DO
!        lbounds(1) = 0._dp  ! surface
!        DO k = 2,nsoil_jsbach
!           lbounds(k) = REAL((levels(k-1) + (levels(k-1) - lbounds(k-1))), dp)
!        END DO
!        DO k = 1,nsoil_jsbach
!           ubounds(k) = REAL((levels(k) + (levels(k) - lbounds(k))), dp)
!        END DO
!     ELSE                                       ! For TERRA
    ! op_bk_20140507-
       channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land) = zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, znlev_soil)
       ALLOCATE(lbounds(znlev_soil), ubounds(znlev_soil), levels(znlev_soil))
       lbounds(1) = 0._dp   ! surface
       DO k = 2, znlev_soil
          lbounds(k)   = REAL((zml_soil(k-1) + (zml_soil(k-1) - lbounds(k-1))),dp)
       ENDDO
       DO k = 1, znlev_soil
          ubounds(k) = REAL((zml_soil(k) + (zml_soil(k) - lbounds(k))),dp)
          levels(k)  = REAL(zml_soil(k)*1000._dp,dp)
       ENDDO
       ubounds(:) = ubounds(:) * 1000._dp        ! in mm
       lbounds(:) = lbounds(:) * 1000._dp        ! in mm
    ! op_bk_20140507+
!     END IF
    ! op_bk_20140507-

    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land), levels)   !necessary for NetCDF
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_below_land), "mm")
    DEALLOCATE(lbounds, ubounds, levels)
    !
    ! Specific soil axis for Runoff_s
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_s) = &
         &                        zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, 1)
    ALLOCATE(levels(1))
    levels(1) = 0._dp  ! in mm
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_s), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_s), "mm")
    DEALLOCATE(levels)
    !
    ! Specific soil axis for Runoff_g
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_g) = &
         &                        zaxisCreate(ZAXIS_DEPTH_BELOW_LAND, 1)
    ALLOCATE(levels(1))
    levels(1) = 0.1_dp * 1000._dp  ! in mm
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_g), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_depth_runoff_g), "mm")
    DEALLOCATE(levels)
    !
    ! SNOW axis (for multi-layer snow model)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_snow_half) = zaxisCreate(ZAXIS_SNOW, nlev_snow+1)
    ALLOCATE(levels(nlev_snow+1))
    DO k = 1, nlev_snow+1
       levels(k) = REAL(k,dp)
    END DO
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_snow_half), levels)
    DEALLOCATE(levels)
    !
    ! SNOW-layer axis (for multi-layer snow model)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_snow) = zaxisCreate(ZAXIS_SNOW, nlev_snow)
    ALLOCATE(levels(nlev_snow), lbounds(nlev_snow), ubounds(nlev_snow))
    DO k = 1, nlev_snow
       lbounds(k) = REAL(k,dp)
       levels(k)  = REAL(k,dp)
    ENDDO
    DO k = 1, nlev_snow
       ubounds(k) = REAL(k+1,dp)
    ENDDO

    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_snow), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_snow), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_snow), levels)   !necessary for NetCDF
    DEALLOCATE(levels, lbounds, ubounds)

    ! Specified height level above ground: 2m
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_height_2m) = zaxisCreate(ZAXIS_HEIGHT, 1)
    ALLOCATE(levels(1))
    levels(1) = 2._dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_height_2m), levels)
    DEALLOCATE(levels)

    ! Specified height level above ground: 10m
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_height_10m) = zaxisCreate(ZAXIS_HEIGHT, 1)
    ALLOCATE(levels(1))
    levels(1) = 10._dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_height_10m), levels)
    DEALLOCATE(levels)

    ! Top of atmosphere
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_toa) = zaxisCreate(ZAXIS_TOA, 1)
    ALLOCATE(levels(1))
    levels(1) = 1._dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_toa), levels)
    DEALLOCATE(levels)

    ! Isobaric surface 800 hPa (layer)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_800) = zaxisCreate(ZAXIS_PRESSURE, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1)= 800._dp   ! hPa
    ubounds(1)= 1013._dp  ! hPa
    levels(1) = 800._dp   ! hPa
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_800), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_800), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_800), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_800), "hPa")
    DEALLOCATE(lbounds, ubounds, levels)

    ! Isobaric surface 400 hPa (layer)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_400) = zaxisCreate(ZAXIS_PRESSURE, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1)= 400._dp   ! hPa
    ubounds(1)= 800._dp   ! hPa
    levels(1) = 400._dp   ! hPa
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_400), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_400), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_400), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_400), "hPa")
    DEALLOCATE(lbounds, ubounds, levels)

    ! Isobaric surface 0 hPa (layer)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_0) = zaxisCreate(ZAXIS_PRESSURE, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1)= 0._dp ! hPa
    ubounds(1)= 400._dp   ! hPa
    levels(1) = 0._dp   ! hPa
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_0), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_0), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_0), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure_0), "hPa")
    DEALLOCATE(lbounds, ubounds, levels)

    ! Specific vertical axis for Lake-model
    ! -------------------------------------
    ! Lake bottom (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom) = zaxisCreate(ZAXIS_LAKE_BOTTOM, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1)= 1._dp
    ubounds(1)= 0._dp
    levels(1) = 1._dp
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom), "m")
    DEALLOCATE(lbounds, ubounds, levels)

    ! Lake bottom half (interface, i.e. only typeOfFirstFixedSurface)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom_half) = zaxisCreate(ZAXIS_LAKE_BOTTOM, 1)
    ALLOCATE(levels(1))
    levels(1) = 0._dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom_half), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_lake_bottom_half), "m")
    DEALLOCATE(levels)

    ! Mixing layer (we define it as a layer in order to be able to re-set
    ! either the first- or secondFixedSurfaces if necessary)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_mix_layer) = zaxisCreate(ZAXIS_MIX_LAYER, 1)
    ALLOCATE(lbounds(1), ubounds(1), levels(1))
    lbounds(1)= 1._dp
    ubounds(1)= 0._dp
    levels(1) = 1._dp
    CALL zaxisDefLbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_mix_layer), lbounds) !necessary for GRIB2
    CALL zaxisDefUbounds(channel%int%cdi(IOMODE)%ZaxisID(ZA_mix_layer), ubounds) !necessary for GRIB2
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_mix_layer), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_mix_layer), "m")
    DEALLOCATE(lbounds, ubounds, levels)

    ! Bottom of sediment layer penetrated by thermal wave (interface, i.e. only typeOfFirstFixedSurface)
    !
    channel%int%cdi(IOMODE)%ZaxisID(ZA_sediment_bottom_tw_half) = zaxisCreate(ZAXIS_SEDIMENT_BOTTOM_TW, 1)
    ALLOCATE(levels(1))
    levels(1) = 0._dp
    CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_sediment_bottom_tw_half), levels)
    CALL zaxisDefUnits(channel%int%cdi(IOMODE)%ZaxisID(ZA_sediment_bottom_tw_half), "m")
    DEALLOCATE(levels)

    ! Define axes for output on p-, i- and z-levels
    !
    ! op_bk_20140507+
!     lwrite_pzlev = (of%name_list%pl_varlist(1) /= ' ')  .OR.  &
!          &         (of%name_list%hl_varlist(1) /= ' ')  .OR.  &
!          &         (of%name_list%il_varlist(1) /= ' ')
!     IF (lwrite_pzlev) THEN

!        ! p-axis
!        !
!        nplev = nh_pzlev_config(idom)%nplev
!        channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure) = zaxisCreate(ZAXIS_PRESSURE, nplev)
!        ALLOCATE(levels(nplev))
!        DO k = 1, nplev
!           levels(k) = REAL(nh_pzlev_config(idom)%plevels(k),dp)
!        END DO
!        CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure), levels)
!        CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_pressure), nplev, levels)
!        DEALLOCATE(levels)
!
!        ! Altitude above mean sea level
!        !
!        nzlev = nh_pzlev_config(idom)%nzlev
!        channel%int%cdi(IOMODE)%ZaxisID(ZA_altitude) = zaxisCreate(ZAXIS_ALTITUDE, nzlev)
!        ALLOCATE(levels(nzlev))
!        DO k = 1, nzlev
!           levels(k) = REAL(nh_pzlev_config(idom)%zlevels(k),dp)
!        END DO
!        CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_altitude), levels)
!        CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_altitude), nzlev, levels)
!        DEALLOCATE(levels)
!
!        ! i-axis (isentropes)
!        !
!        nilev = nh_pzlev_config(idom)%nilev
!        channel%int%cdi(IOMODE)%ZaxisID(ZA_isentropic) = zaxisCreate(ZAXIS_ISENTROPIC, nilev)
!        ALLOCATE(levels(nilev))
!        DO k = 1, nilev
!           levels(k) = REAL(nh_pzlev_config(idom)%ilevels(k),dp)
!        END DO
!        CALL zaxisDefLevels(channel%int%cdi(IOMODE)%ZaxisID(ZA_isentropic), levels)
!        CALL zaxisDefVct(channel%int%cdi(IOMODE)%ZaxisID(ZA_isentropic), nilev, levels)
!        DEALLOCATE(levels)
!     ENDIF
    ! op_bk_20140507-

    ! ---------------------------------------------------------
    ! define variables (including secondary data)
    ! ---------------------------------------------------------
    ! - loop over all channel objects / variables
    le => channel%list
    object_loop1: DO
       IF (.NOT. ASSOCIATED(le)) EXIT
       object => le%this
       ! any output for this object ???
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          lskip = .NOT. (object%int%lout .AND. object%int%cdi(IOMODE)%lout)
       CASE(IOMODE_RST)
          lskip = .NOT. object%int%lrst
       END SELECT
       !
       IF (lskip) THEN
          le => le%next          
          CYCLE
       END IF

       ! define primary and secondary data
       DO jsnd=1, SND_MAXLEN
          ! set horizontal grid
          SELECT CASE(object%repr%hgrid)
          CASE(GRID_UNSTRUCTURED_CELL)
             object%int%cdi(IOMODE)%gridID = channel%int%cdi(IOMODE)%CellGridID
          CASE(GRID_UNSTRUCTURED_EDGE)
             object%int%cdi(IOMODE)%gridID = channel%int%cdi(IOMODE)%EdgeGridID
          CASE(GRID_UNSTRUCTURED_VERT)
             object%int%cdi(IOMODE)%gridID = channel%int%cdi(IOMODE)%VertGridID
          END SELECT
          ! set vertical grid
          object%int%cdi(IOMODE)%zaxisID = channel%int%cdi(IOMODE)%ZaxisID(object%repr%vgrid)

          IF (.NOT. object%int%lexp(jsnd, IOMODE)) CYCLE ! no output
          ! op_bk_20141106+
          SELECT CASE(jsnd)
          CASE(SND_INS)
             steptype = object%steptype
          CASE(SND_AVE)
             steptype = TSTEP_AVG
          CASE(SND_STD)
             steptype = TSTEP_SD
          CASE(SND_MIN)
             steptype = TSTEP_MIN
          CASE(SND_MAX)
             steptype = TSTEP_MAX
          CASE(SND_CNT)
             ! TSTEP_RANGE ??
             steptype = TSTEP_RANGE
          CASE(SND_CAV)
             ! TSTEP_AVG ??
             steptype = TSTEP_AVG
          END SELECT
          ! op_bk_20141106-
          varID =  vlistDefVar(vlistID, object%int%cdi(IOMODE)%gridID    &
               &             , object%int%cdi(IOMODE)%zaxisID, steptype)
          IF (varID == CDI_UNDEFID) THEN
             WRITE(*,*) 'VARID FAILED'
             RETURN
          END IF

          CALL vlistDefVarName(vlistID, varID, TRIM(object%name)//TRIM(SND_TEXT(jsnd,IOMODE)))
          CALL vlistDefVarDatatype(vlistID, varID, CDINF_PREC)

          IF (jsnd == SND_INS) then
             object%int%cdi(IOMODE)%varID = varID        ! primary data
          ELSE
             i2nd = object%int%i2nd(jsnd)
             object%int%cdi(IOMODE)%svarID(i2nd) = varID ! 2ndary data
          END IF

          ! write variable attributes
          CALL cdi_write_attribute_list(status, vlistID, varID, object%att &
               , CDINF_PREC)
          ! write special attributes for secondary data
          IF ((jsnd == SND_CNT) .OR. (jsnd == SND_CAV)) THEN
             SELECT CASE(CDINF_PREC)
             CASE(DATATYPE_FLT64)
                iret = vlistDefAttFlt(vlistID, varID, 'range' &
                     , DATATYPE_FLT64, SIZE(object%io%range), object%io%range)
                CALL handle_cdi_error(iret, status, substr)
                IF (status /= 0) RETURN
             CASE(DATATYPE_FLT32)
                iret = vlistDefAttFlt(vlistID, varID, 'range' &
                     , DATATYPE_FLT32, SIZE(REAL(object%io%range)) &
                     , REAL(object%io%range))
                CALL handle_cdi_error(iret, status, substr)
                IF (status /= 0) RETURN
             END SELECT
          END IF
       END DO
       ! -------------------------------------------------------

       le => le%next
    END DO object_loop1

  END SUBROUTINE ch_cdi_write_header
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_write_time(status, channel)


    USE messy_main_timer,      ONLY: &
    &                                YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
!!$    USE mtime,                 ONLY: setCalendar, PROLEPTIC_GREGORIAN

    IMPLICIT NONE
    
    ! I/O
    INTEGER,           INTENT(OUT) :: status
    TYPE(t_channel),   POINTER     :: channel ! INTENT(INOUT)

    ! local
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_cdi_write_time'
    INTEGER                        :: fileID, vlistID
    INTEGER                        :: idate, itime, iret

    fileID = channel%int%cdi(IOMODE_OUT)%fileID
    IF (fileID == CDI_UNDEFID) THEN
       status = 4011 ! undefined file-id
       RETURN
    END IF

    vlistID = channel%int%cdi(IOMODE_OUT)%vlistID
    IF (vlistID == CDI_UNDEFID) THEN
       status = 4012 ! undefined vlist-id
       RETURN
    END IF

    IF (streamInqVlist(fileID) == CDI_UNDEFID) THEN
       CALL streamDefVlist(fileID, vlistID)
    END IF

    ! set timestep
!!$    CALL setCalendar(PROLEPTIC_GREGORIAN)
    idate = cdiEncodeDate(YEAR, MONTH, DAY)
    itime = cdiEncodeTime(HOUR, MINUTE, SECOND)
    WRITE(*,'(A24,I4,5(A1,I2.2))') "MESSY: time for output: ", YEAR, "-", MONTH, "-", DAY, " ", HOUR, ":", MINUTE, ":", SECOND
    CALL taxisDefVdate(channel%int%cdi(IOMODE_OUT)%taxisID, idate)
    CALL taxisDefVtime(channel%int%cdi(IOMODE_OUT)%taxisID, itime)
    iret = streamDefTimestep(channel%int%cdi(IOMODE_OUT)%fileID, channel%int%ntpfcnt-1)
    CALL handle_cdi_error(iret, status, substr)

    status = 0

  END SUBROUTINE ch_cdi_write_time
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_write_data(status, IOMODE, channel, object &
       , ptr, jsnd, i2nd)

!!$    USE mo_communication,         ONLY: exchange_data, t_comm_gather_pattern
!!$    USE mo_name_list_output_init, ONLY: patch_info
!!$    USE mo_mpi,                   ONLY: my_process_is_mpi_workroot

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    ! I/O
    INTEGER,                INTENT(OUT)   :: status
    INTEGER,                INTENT(IN)    :: IOMODE
    TYPE(t_channel),        POINTER       :: channel ! INTENT(INOUT)
    TYPE(t_channel_object), POINTER       :: object  ! INTENT(INOUT)
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr
    INTEGER,                INTENT(IN)    :: jsnd    ! OUTPUT DATA TYPE
    INTEGER,                INTENT(IN)    :: i2nd    ! index of 2ndary DATA

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'ch_cdi_write_data'
    INTEGER                               :: fileID, varID
    REAL(DP), DIMENSION(:,:), POINTER     :: r_out
!!$    TYPE(t_comm_gather_pattern), POINTER  :: p_pat
    INTEGER                               :: idom, n_points, nlevs
    INTEGER                               :: nrank

    IF (I_VERBOSE_LEVEL >= 2) THEN
       SELECT CASE(IOMODE)
       CASE(IOMODE_OUT)
          IF (I_VERBOSE_LEVEL >= 4) THEN
             WRITE(*,*) '... netCDF-output: '                           &
                  &   , TRIM(channel%int%fname(IOMODE))                 &
                  &   , ' (',channel%int%ntpfcnt,') <- '                &
                  &   , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE)) &
                  &   , ' (',MINVAL(ptr),' - ',MAXVAL(ptr),')'
          ELSE
             WRITE(*,*) '... netCDF-output: '                           &
                  &   , TRIM(channel%int%fname(IOMODE))                 & 
                  &   , ' (',channel%int%ntpfcnt,') <- '                &
                  &   , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))
          ENDIF
          !
       CASE(IOMODE_RST)
          IF (I_VERBOSE_LEVEL >= 4) THEN
             WRITE(*,*) '... netCDF-output: '                           &
                  &   , TRIM(channel%int%fname(IOMODE))                 &
                  &   , ' <- '                                          &
                  &   , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE)) &
                  &   , ' (',MINVAL(ptr),' - ',MAXVAL(ptr),')'
          ELSE
             WRITE(*,*) '... netCDF-output: '                           &
                  &   , TRIM(channel%int%fname(IOMODE))                 &
                  &   , ' <- '                                          &
                  &   , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))
          ENDIF
       END SELECT
    END IF

    fileID = channel%int%cdi(IOMODE)%fileID
    IF (fileID == CDI_UNDEFID) THEN
       status = 4011 ! undefined file-id
       RETURN
    END IF

    SELECT CASE(jsnd)
    CASE(SND_UNDEF)
       status = 3203 ! secondary data type unknown
       RETURN
    CASE(SND_INS)
       varID = object%int%cdi(IOMODE)%varID
    CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
       varID = object%int%cdi(IOMODE)%svarID(i2nd)
    END SELECT

!    idom = channel%int%cdi(IOMODE)%dom_id
!    nrank = object%repr%rank    
!
!    IF (nrank == 2) THEN
!       nlevs = 1
!    ELSE
!       nlevs = object%repr%gdimlen(2)
!    END IF
!
!    SELECT CASE(object%repr%hgrid)
!    CASE(GRID_UNSTRUCTURED_CELL)
!       p_pat => patch_info(idom)%p_pat_c
!       n_points = patch_info(idom)%cells%n_glb
!    CASE(GRID_UNSTRUCTURED_EDGE)
!       p_pat => patch_info(idom)%p_pat_e
!       n_points = patch_info(idom)%edges%n_glb
!    CASE(GRID_UNSTRUCTURED_VERT)
!       p_pat => patch_info(idom)%p_pat_v
!       n_points = patch_info(idom)%verts%n_glb
!    CASE DEFAULT
!       status = 4015
!       RETURN
!    END SELECT
!
!    ALLOCATE(r_out(MERGE(n_points, 0, my_process_is_mpi_workroot()), nlevs))
!    r_out(:,:) = 0._dp
!    r_out(:,:) = ptr(:,:,1,1)
!
!    CALL streamWriteVar(channel%int%cdi(IOMODE)%fileID, object%int%cdi(IOMODE)%varID, r_out, 0)
    ! op_bk_20141106+
!     CALL streamWriteVar(channel%int%cdi(IOMODE)%fileID, object%int%cdi(IOMODE)%varID, ptr(:,:,1,1), 0)
    CALL streamWriteVar(channel%int%cdi(IOMODE)%fileID, varID, ptr(:,:,1,1), 0)
    ! op_bk_20141106-
    IF (L_FLUSH_IOBUFFER) THEN
       CALL streamSync(channel%int%cdi(IOMODE)%fileID)
    END IF

!    DEALLOCATE(r_out)
    status = 0

  END SUBROUTINE ch_cdi_write_data
  ! -------------------------------------------------------------------

#endif

  ! -------------------------------------------------------------------
  SUBROUTINE ch_cdi_finish_io(status, IOMODE, channel, lclose)

!!$    USE mo_name_list_output_types,     ONLY: max_z_axes

    IMPLICIT NONE

#ifdef HAVE_CDI
    INTRINSIC :: TRIM
#endif

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: IOMODE
    TYPE(t_channel),  POINTER     :: channel ! INTENT(INOUT)
    LOGICAL,          INTENT(IN)  :: lclose

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'ch_cdi_finish_io'
#ifdef HAVE_CDI
    INTEGER :: streamid, i
#endif

    IF (channel%int%fname(IOMODE) == '') THEN
       status = 0 ! O.K.: file not opened
       RETURN
    END IF
    
#ifdef HAVE_CDI

    streamid = channel%int%cdi(IOMODE)%fileID
    IF (streamid == CDI_UNDEFID) THEN
       status = 4011 ! undefined file-id
       RETURN
    END IF

    IF (lclose) THEN
       IF (I_VERBOSE_LEVEL >= 2) THEN
          WRITE(*,*) substr,': CLOSING FILE: ',TRIM(channel%int%fname(IOMODE))
       END IF
       CALL streamClose(streamid)
       channel%int%cdi(IOMODE)%fileID = CDI_UNDEFID
       channel%int%fname(IOMODE) = ''
       IF (channel%int%cdi(IOMODE)%vlistID /= CDI_UNDEFID) THEN
          CALL vlistDestroy(channel%int%cdi(IOMODE)%vlistID)
          channel%int%cdi(IOMODE)%vlistID = CDI_UNDEFID
       END IF
       IF (channel%int%cdi(IOMODE)%taxisID /= CDI_UNDEFID) THEN
          CALL taxisDestroy(channel%int%cdi(IOMODE)%taxisID)
          channel%int%cdi(IOMODE)%taxisID = CDI_UNDEFID
       END IF
       DO i=1, SIZE(channel%int%cdi(IOMODE)%zaxisID)
          IF (channel%int%cdi(IOMODE)%zaxisID(i) /= CDI_UNDEFID) THEN
             CALL zaxisDestroy(channel%int%cdi(IOMODE)%zaxisID(i))
             channel%int%cdi(IOMODE)%zaxisID(i) = CDI_UNDEFID
          END IF
       END DO
       IF (channel%int%cdi(IOMODE)%CellGridID /= CDI_UNDEFID) THEN
          CALL gridDestroy(channel%int%cdi(IOMODE)%CellGridID)
          channel%int%cdi(IOMODE)%CellGridID = CDI_UNDEFID
       END IF
       IF (channel%int%cdi(IOMODE)%EdgeGridID /= CDI_UNDEFID) THEN
          CALL gridDestroy(channel%int%cdi(IOMODE)%EdgeGridID)
          channel%int%cdi(IOMODE)%EdgeGridID = CDI_UNDEFID
       END IF
       IF (channel%int%cdi(IOMODE)%VertGridID /= CDI_UNDEFID) THEN
          CALL gridDestroy(channel%int%cdi(IOMODE)%VertGridID)
          channel%int%cdi(IOMODE)%VertGridID = CDI_UNDEFID
       END IF
       status = 0
    ELSE
       status = 0
    END IF
#endif

  END SUBROUTINE ch_cdi_finish_io
  ! -------------------------------------------------------------------

#ifdef HAVE_CDI

  ! op_bk_20140425+
!   ! -------------------------------------------------------------------
!   SUBROUTINE ch_cdi_read_data(status, IOMODE, channel, object &
!        , ptr, jsnd, i2nd)

!     USE messy_main_channel_repr,  ONLY: IRANK

!     IMPLICIT NONE

!     INTRINSIC :: ASSOCIATED, TRIM

!     ! I/O
!     INTEGER,                INTENT(OUT)   :: status
!     INTEGER,                INTENT(IN)    :: IOMODE
!     TYPE(t_channel),        POINTER       :: channel ! INTENT(INOUT)
!     TYPE(t_channel_object), POINTER       :: object  ! INTENT(INOUT)
!     REAL(DP), DIMENSION(:,:,:,:), POINTER :: ptr
!     INTEGER,                INTENT(IN)    :: jsnd    ! OUTPUT DATA TYPE
!     INTEGER,                INTENT(IN)    :: i2nd    ! index of 2ndary DATA

!     ! LOCAL
!     CHARACTER(LEN=*), PARAMETER :: substr = 'ch_cdi_read_data'
!     INTEGER                               :: ncid
!     INTEGER                               :: varID
!     INTEGER                               :: zstat
!     INTEGER, DIMENSION(IRANK)             :: nsz

!     IF (channel%int%fname(IOMODE) == '') THEN
!        status = 0 ! O.K.: FILE NOT OPENED
!        RETURN
!     END IF

!     ncid  = channel%int%netcdf(IOMODE)%fileID
!     !
!     IF (ncid == NC_ID_UNDEF) THEN
!        status = 4001 ! UNDEFINED FILE-ID
!        RETURN
!     END IF

!     IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
!          WRITE(*,*) '... netCDF-input: ' &
!          , TRIM(channel%int%fname(IOMODE)), ' <- ' &
!          , TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE))

!     ! GET VARAIABLE ID
!     zstat = NF90_INQ_VARID( ncid, &
!          TRIM(object%name)//TRIM(SND_TEXT(jsnd, IOMODE)), varID)
!     IF (zstat /= NF90_NOERR) THEN
!        ! netCDF VARIABLE DOES NOT EXIST -> CHECK
!        SELECT CASE(jsnd)
!        CASE(SND_INS)
!           IF (object%lrestreq) THEN
!              IF (object%int%lign) THEN
!                 IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
!                      WRITE(*,*) substr &
!                      ,'    WARNING: REQUIRED RESTART VARIABLE ', &
!                      'NOT PRESENT! ' &
!                      ,'HOWEVER: IGNORE = T FOR THIS OBJECT'
!                 status = 0
!              ELSE
!                 status = 3211 ! RESTART VARIABLE REQUIRED BUT NOT PRESENT
!              END IF
!           ELSE
!              IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
!                   WRITE(*,*) '    WARNING: VARIABLE NOT PRESENT IN RESTART FILE'
!              status = 0
!           END IF
!        CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX
!           IF (I_VERBOSE_LEVEL >= 1) & ! op_pj_20110803
!                WRITE(*,*) '    WARNING: VARIABLE NOT PRESENT IN RESTART FILE'
!           status = 0
!        END SELECT
!        RETURN ! RETURN IF VARIABLE NOT PRESENT
!     END IF
!     !
!     SELECT CASE(jsnd)
!     CASE(SND_UNDEF)
!        status = 3203 ! SECONDARY DATA TYPE UNKNOWN
!        RETURN
!     CASE(SND_INS)
!        object%int%netcdf(IOMODE)%varID = varID
!     CASE DEFAULT ! SND_AVE, SND_STP, SND_MIN, SND_MAX, ...
!        object%int%netcdf(IOMODE)%svarID(i2nd) = varID
!     END SELECT

!     ! CHECK VARIABLE ATTRIBUTES
!     CALL cdi_check_attributes(status, ncid, varID, object%att)
!     IF (status /= 0) RETURN

!     ! SETUP TEMPORARY MEMORY FOR IMPORT
!     nsz(:) = object%repr%shape_out(:)
!     ALLOCATE(ptr(nsz(1),nsz(2),nsz(3),nsz(4)), STAT=zstat)
!     IF (zstat /= 0) THEN
!        status = 1000 ! MEMORY ALLOCATION FAILED
!        RETURN
!     END IF

!     ! IMPORT DATA
!     CALL nf( &
!          NF90_GET_VAR(ncid, varID, ptr) &
!          , status, substr)
  ! op_bk_20140425-
!     IF (status /= 0) RETURN
    
!     status = 0

!   END SUBROUTINE ch_cdi_read_data
!   ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! PRIVATE SUBROUTINES
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE cdi_write_attribute_list(status, vlistID, varID, attlist, prec)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
         , TYPE_STRING, TYPE_INTEGER, TYPE_REAL_DP

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, NULL, TRIM

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    INTEGER,                INTENT(IN)  :: vlistID
    INTEGER,                INTENT(IN)  :: varID
    TYPE(t_attribute_list), POINTER     :: attlist
    INTEGER,                INTENT(IN)  :: prec

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'cdi_write_attribute_list'
    TYPE(t_attribute_list), POINTER :: ai  => NULL()
    TYPE(t_attribute),      POINTER :: att => NULL()
    INTEGER                         :: iret, len

    ! op_bk_20151007+
    INTEGER :: OneInt(1)
    REAL(DP) :: OneReal_dp(1)
    REAL     :: OneReal(1)
    ! op_bk_20151007-

    ai => attlist
    DO
       IF (.NOT. ASSOCIATED(ai)) EXIT

       att => ai%this

       SELECT CASE(ai%this%type)
       CASE(TYPE_STRING)
          len = LEN_TRIM(att%c)
          iret = vlistDefAttTxt(vlistID, varID, TRIM(att%name), len, TRIM(att%c))
          CALL handle_cdi_error(iret, status, substr)
          IF (status /= 0) RETURN
       CASE(TYPE_INTEGER)
          OneInt(1) = att%i
          iret = vlistDefAttInt(vlistID, varID, TRIM(att%name), DATATYPE_INT32, 1, OneInt)
          CALL handle_cdi_error(iret, status, substr)
          IF (status /= 0) RETURN
       CASE(TYPE_REAL_DP)
          SELECT CASE(prec)
          CASE(DATATYPE_FLT64)
             OneReal_dp(1) = att%r
             iret = vlistDefAttFlt(vlistID, varID, TRIM(att%name) &
                  , DATATYPE_FLT64, 1, OneReal_dp)
             CALL handle_cdi_error(iret, status, substr)
             IF (status /= 0) RETURN
          CASE(DATATYPE_FLT32)
             OneReal(1) = REAL(att%r)
             iret = vlistDefAttFlt(vlistID, varID, TRIM(att%name) &
                  , DATATYPE_FLT32, 1, OneReal)
             CALL handle_cdi_error(iret, status, substr)
             IF (status /= 0) RETURN
          END SELECT
       CASE DEFAULT
          status = 806  ! UNKNOWN ATTRIBUTE TYPE
       END SELECT

       IF (status /= 0) RETURN

       ai => ai%next
    END DO    

    status = 0

  END SUBROUTINE cdi_write_attribute_list
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE cdi_check_attributes(status, fileid, varID, attlist)

    USE messy_main_channel_attributes, ONLY: t_attribute_list, t_attribute &
         &                                 , TYPE_INTEGER, TYPE_REAL_DP    &
         &                                 , TYPE_STRING
    USE messy_main_constants_mem,      ONLY: STRLEN_ULONG, STRLEN_MEDIUM

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM, ABS, TINY

    ! I/O
    INTEGER,          INTENT(OUT)   :: status
    INTEGER,          INTENT(IN)    :: fileid
    INTEGER,          INTENT(IN)    :: varID
    TYPE(t_attribute_list), POINTER :: attlist  ! INTENT(IN)

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'cdi_check_attributes'
    TYPE(t_attribute_list), POINTER :: list
    TYPE(t_attribute),      POINTER :: att
    INTEGER                         :: i, vlistID, natts
    INTEGER                         :: zstat
    INTEGER                         :: xtype, len, attnum
    INTEGER                         :: ai(1)
    CHARACTER(LEN=STRLEN_ULONG)     :: ac = ''
    REAL(DP)                        :: ar(1)
    LOGICAL                         :: lok, lfound
    ! op_bk_20140425+
    CHARACTER(LEN=STRLEN_MEDIUM)    :: attname
    ! op_bk_20140425-


    ! vlistID from file
    vlistID = streamInqVlist(fileid)
    ! total number of attributes in file
    CALL handle_cdi_error(vlistInqNatts(vlistID, varID, natts),status,substr)
    IF (status /= 0) RETURN

    list => attlist
    DO
       IF (.NOT. ASSOCIATED(list)) EXIT
       att => list%this
       ! -------------------------
       SELECT CASE(att%iflag)
       CASE(AF_RST_NONE)
          list => list%next
          CYCLE
       CASE(AF_RST_CMP)
          IF (I_VERBOSE_LEVEL >= 1) &
               & WRITE(*,*) '   CHECKING ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE(AF_RST_INP)
          IF (I_VERBOSE_LEVEL >= 1) &
               & WRITE(*,*) '   READING  ATTRIBUTE ''',TRIM(att%name),''' ...'
       CASE DEFAULT
          list => list%next
          CYCLE
       END SELECT

       ! inquire attribute
       lfound = .FALSE.
       DO i = 0, natts-1
          attname = ''
          len = 0
          status = vlistInqAtt(vlistID, varID, i, attname, xtype, len)
          IF (status /= CDI_NOERR) THEN
             IF (TRIM(att%name) == TRIM(attname)) THEN
                lfound = .TRUE.
                EXIT
             END IF
          END IF
          !
       END DO

       IF (.NOT. lfound) THEN
          IF (I_VERBOSE_LEVEL >= 0) &
               & WRITE(*,*) '   ... *** ERROR *** ATTRIBUTE NOT PRESENT'
          status = 3207  ! MISSING ATTRIBUTE IN RESTART FILE             
          RETURN
       END IF

       ! check length
       IF (xtype /= DATATYPE_TXT) THEN
          IF (len /= 1) THEN
             status = 3208 ! restart attribute has non-scalar rank
             RETURN
          END IF
       ELSE
          IF (len > STRLEN_ULONG) THEN
             status = 3209 ! restart character attribute too long
             RETURN
          END IF
       END IF
       
       ! get attribute and set type
       SELECT CASE(xtype)
       CASE(DATATYPE_INT8)
          CALL handle_cdi_error(vlistInqAttInt(vlistID, varID  &
               &              , TRIM(att%name), 1, ai(1)), status, substr)
       CASE(DATATYPE_TXT)
          CALL handle_cdi_error(vlistInqAttTxt(vlistID, varID  &
               &              , TRIM(att%name), len, ac(1:len+1)), status, substr)
       CASE(DATATYPE_INT16)
          CALL handle_cdi_error(vlistInqAttInt(vlistID, varID  &
               &              , TRIM(att%name), 1, ai(1)), status, substr)
       CASE(DATATYPE_INT32)
          CALL handle_cdi_error(vlistInqAttInt(vlistID, varID  &
               &              , TRIM(att%name), 1, ai(1)), status, substr)
       CASE(DATATYPE_FLT32)
          CALL handle_cdi_error(vlistInqAttFlt(vlistID, varID  &
               &              , TRIM(att%name), 1, ar(1)), status, substr)
       CASE(DATATYPE_FLT64)
          CALL handle_cdi_error(vlistInqAttFlt(vlistID, varID  &
               &              , TRIM(att%name), 1, ar(1)), status, substr)
       END SELECT
       IF (status /= 0) RETURN
          
       SELECT CASE(att%iflag)
       CASE(AF_RST_CMP) ! COMPARE
          ! CHECK VALUE
          SELECT CASE(att%type)
          CASE(TYPE_INTEGER)
             lok = (att%i == ai(1))
             IF (.NOT. lok) THEN
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     & WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ai &
                     & ,' /= ', att%i
             ELSE
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     & WRITE(*,*) '      ... OK: ',ai(1)
             END IF
          CASE(TYPE_REAL_DP)
             lok = (ABS(att%r - ar(1)) < TINY(ar(1)))
             IF (.NOT. lok) THEN
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     & WRITE(*,*) '      ... ATTRIBUTE MISMATCH: ',ar(1) &
                     & ,' /= ', att%r
             ELSE
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     & WRITE(*,*) '      ... OK: ',ar(1)
             END IF
          CASE(TYPE_STRING)
             lok = (att%c(1:len) == ac(1:len))
             IF (.NOT. lok) THEN
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     WRITE(*,*) '      ... ATTRIBUTE MISMATCH: '      &
                     & ,ac(1:len), ' /= ',att%c(1:len)
             ELSE
                IF (I_VERBOSE_LEVEL >= 1)                             &
                     & WRITE(*,*) '      ... OK: ',TRIM(ac(1:len))
             END IF
          END SELECT
          !
          IF (.NOT. lok) THEN
             status = 3210  ! RESTART ATTRIBUTE MISMATCH
             RETURN
          END IF
          !
       CASE(AF_RST_INP) ! INPUT
          ! INPUT VALUE
          SELECT CASE(att%type)
          CASE(TYPE_INTEGER)
             att%i = ai(1)
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', att%i
          CASE(TYPE_REAL_DP)
             att%r = ar(1)
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', att%r
          CASE(TYPE_STRING)
             att%c = ac
             IF (I_VERBOSE_LEVEL >= 1) WRITE(*,*) '      ... ', TRIM(att%c)
          END SELECT
          !
       CASE DEFAULT
          !
          !
       END SELECT
       
       ! -------------------------
       list => list%next
    END DO

    CALL vlistDestroy(vlistID)
    ! return
    status = 0

  END SUBROUTINE cdi_check_attributes
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE handle_cdi_error(sin, sout, substr)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER,          INTENT(IN)  :: sin
    INTEGER,          INTENT(OUT) :: sout
    CHARACTER(LEN=*), INTENT(IN)  :: substr

    ! mz_pj_20080807+
    ! LOCAL
    INTEGER :: iou
    LOGICAL :: opened
    ! mz_pj_20080807-

    IF (sin < CDI_NOERR) THEN

       WRITE(*,*) TRIM(substr),': *** CDI ERROR: ',CDI_STRERROR(sin)
       sout = 4010 ! CDI error

       DO iou=100,300
          INQUIRE(unit=iou,opened=opened)
          IF (.NOT.opened) EXIT
       END DO
       OPEN(iou, FILE='ERROR.cdi', STATUS='UNKNOWN')
       WRITE(iou,*) TRIM(substr),': *** CDI ERROR: ',CDI_STRERROR(sin)
       CLOSE(iou)
    ELSE
       sout = 0
    ENDIF
    
  END SUBROUTINE handle_cdi_error
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  FUNCTION CDI_STRERROR(sin)
    INTEGER, INTENT(IN) :: sin
    CHARACTER(LEN=80)   :: CDI_STRERROR

    SELECT CASE(sin)
    CASE (-10)
       CDI_STRERROR = "Operating system error"
       RETURN
    CASE (-20)
       CDI_STRERROR = "Invalid argument"
       RETURN
    CASE (-21)
       CDI_STRERROR = "Unsupported file format"
       RETURN
    CASE (-22)
       CDI_STRERROR = "Library not available"
       RETURN
    CASE (-23)
       CDI_STRERROR = "Unsupported file structure"
       RETURN
    CASE (-24)
       CDI_STRERROR = "Unsupported NetCDF4 structure"
       RETURN
    CASE (-99)
       CDI_STRERROR = "Internal limits exceeded"
       RETURN
    CASE DEFAULT
       CDI_STRERROR = "Unknown error"
       RETURN
    END SELECT
  END FUNCTION CDI_STRERROR
  ! -------------------------------------------------------------------

#endif

! **********************************************************************
END MODULE messy_main_channel_cdi
! **********************************************************************
 
