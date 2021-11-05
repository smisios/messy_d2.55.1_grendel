! ****************************************************************************
MODULE messy_main_grid_trafo_scrp
! ****************************************************************************

  USE messy_main_blather,              ONLY: debug
  USE messy_main_grid_netcdf,          ONLY: t_ncdim
  USE messy_main_constants_mem,        ONLY: dp, sp
  USE messy_main_grid_trafo_scrp_base, ONLY: norm_opt_none     &
                                           , norm_opt_dstarea  &
                                           , norm_opt_frcarea  & 
                                           , map_type_conserv  &
                                           , map_type_bilinear &
                                           , map_type_bicubic  &
                                           , map_type_distwgt 

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: ABS, ALLOCATED, ASSOCIATED, MAX, MAXVAL &
             , MIN, MINVAL, NULL, PRESENT, REAL, SIZE, TRIM

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr='grid_trafo_scrp'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver='0.1'

  ! DEFINE NUMBER OF REMAPPINGS: HERE REMAPPING IS ALWAYS RESTRICTED TO
  ! ONE DIRECTION ONLY, i.e., from source to destination
  INTEGER, PARAMETER, PUBLIC :: num_maps   = 1

  REAL(dp), DIMENSION(2)     :: lon_ranges = (/0._dp, 360._dp/)

  TYPE t_scrip_weights
     ! number of unique address pairs in the remapping == number of entries
     ! in the sparse matrix for the remapping
     INTEGER :: num_links = 0

     ! number of required weights 
     !    - bilinear          num_wgts = 1
     !    - distance-weighted num_wgts = 1
     !    - conservative      num_wgts = 3
     !    - bicubic           num_wgts = 4
     INTEGER                        :: num_wgts  = 0

     ! normally the weights are calculated only at the beginning of the
     ! simulation. If ltimedependent is set to true, the weights are
     ! recalculated each import time step. 
     ! Note: this is much more computing time and memory intensive and is only
     ! required, if time dependent masks (e.g. ice mask) are important for the 
     ! remapping
     LOGICAL                         :: ltimedependent = .FALSE.

     ! THE FOLLOWING DATA SHOULD OPTIMALLY be CHANNEL OBJECTS
     ! (TODO) CHALLENGE: How to associate RESTART elements to SDATA structs?
     ! remap matrix dimensioned by (num_wghts, num_links)
     REAL(dp), DIMENSION(:), POINTER :: weights => NULL()
     REAL(dp), DIMENSION(:), POINTER :: dstfrac => NULL()  ! fraction 
     ! area of destination field
     REAL(dp), DIMENSION(:), POINTER :: dstarea => NULL()
     ! source address of each link (dim: num_links)
     INTEGER,  DIMENSION(:), POINTER :: srcadd  => NULL()
     ! destination address of each link (dim: num_links)
     INTEGER,  DIMENSION(:), POINTER :: dstadd  => NULL()
  END type t_scrip_weights

  TYPE t_scrip_grid
     INTEGER :: ID      = -99
     INTEGER :: size    = 0   ! "Horizontal" size product of lon/lat grid dims 
                              ! i.e. a 1D size horizontal part  of the grid
     INTEGER :: corners = 4   ! number of corners, in our case always 4
     INTEGER :: rank    = 0   ! number of dimensions in "model code"
     
     ! length of dimensions in "model code" (dimensioned by rank)
     TYPE (t_ncdim), DIMENSION(:), POINTER :: dim => NULL()

     LOGICAL,        DIMENSION(:), POINTER :: lmask => NULL() 
                                                    ! dimensioned by grid_size
                                                    ! T for participating points
                                                    ! F for neglected points
                                                       
     ! longitude of grid centers (dim: grid_size ; unit: radians)
     REAL(dp), DIMENSION(:),   POINTER :: center_lon => NULL()
     ! latitude of grid centers (dim: grid_size ; unit: radians)
     REAL(dp), DIMENSION(:),   POINTER :: center_lat => NULL()
     ! longitude of grid corners (dim: (grid_size, grid_corners); unit: radians)
     REAL(dp), DIMENSION(:,:), POINTER :: corner_lon => NULL()
     ! latitude of grid corners (dim: (grid_size, grid_corners); unit: radians)
     REAL(dp), DIMENSION(:,:), POINTER :: corner_lat => NULL()
  END type t_scrip_grid

  TYPE t_scrip_data
     ! SCRIP data ID
     INTEGER            :: id             = 0
     ! remapping type
     INTEGER            :: map_type       = map_type_conserv
     ! normaize option
     INTEGER            :: norm_opt       = norm_opt_none
     LOGICAL            :: luse_sgrd_area = .FALSE.! use source grid area
     LOGICAL            :: luse_dgrd_area = .FALSE.! use destination grid area
     TYPE(t_scrip_weights)          :: wghts
     TYPE(t_scrip_grid), POINTER    :: sgrd => NULL()
     TYPE(t_scrip_grid), POINTER    :: dgrd => NULL()
  END type t_scrip_data

  ! BUILD UP concatenated list of scrip_data
  TYPE t_scrip_data_list
     TYPE(t_scrip_data)               :: this
     TYPE(t_scrip_data_list), POINTER :: next => NULL()
  END type t_scrip_data_list
  PUBLIC :: t_scrip_data_list
  TYPE(t_scrip_data_list), POINTER, PUBLIC :: SCRIPDATALIST => NULL()

  ! number of defined grids
  INTEGER, PUBLIC                          :: NSDATA = 0 

  TYPE t_scrip_grid_list
     TYPE(t_scrip_grid)               :: this
     TYPE(t_scrip_grid_list), POINTER :: next => NULL()
  END type t_scrip_grid_list
  PUBLIC :: t_scrip_grid_list

  TYPE(t_scrip_grid_list), POINTER,  PUBLIC :: SCRIPGRIDLIST => NULL()
  ! number of defined grids
  INTEGER                                   :: NSGRID = 0 

  ! TYPES
  PUBLIC :: T_SCRIP_DATA, T_SCRIP_GRID, T_SCRIP_WEIGHTS

  ! SUBROUTINES
  PUBLIC :: GEOHYB2SCRIPGRID
  PUBLIC :: DEFINE_SCRIPGRID
  PUBLIC :: DEFINE_SCRIPDATA
  PUBLIC :: CALC_SCRIPDATA
  PUBLIC :: CALC_SCRIP_WEIGHTS
  PUBLIC :: APPLY_SCRIP_WEIGHTS
  PUBLIC :: LOCATE_SCRIPDATA

  PUBLIC :: SCRIP_CONTROL
  PUBLIC :: BALANCE_CURVILINEAR_PS
  PUBLIC :: CONSTRUCT_INPUT_SURF_PRESSURE
  PUBLIC :: INTERPOL_GEOHYBGRID_PS

  PUBLIC :: CLEAN_SCRIPGRID_LIST
  PUBLIC :: CLEAN_SCRIPDATA_LIST
  PUBLIC :: CONSTRUCT_VERTICAL_AXIS

  ! PRIVATE
  ! INIT_SCRIPGRID
  ! COPY_SCRIPGRID
  ! INIT_SCRIPDATA
  ! COPY SCRIPDATA

CONTAINS

  ! ---------------------------------------------------------------------------
  SUBROUTINE GEOHYB2SCRIPGRID(status, ggrid, sgrid, psgrid, ID, l_set_ranges &
       , llab)

    USE messy_main_constants_mem, ONLY: dp, pi, dtr &
                                      , DTR_icon
    USE messy_main_grid_trafo,    ONLY: COMPLETE_GEOHYBGRID
    USE messy_main_grid,          ONLY: t_geohybgrid, COPY_GEOHYBGRID      &
                                      , INIT_GEOHYBGRID, RGEMPTY
    USE messy_main_grid_netcdf,   ONLY: GRD_MAXSTRLEN, ELEMENT, NULL_XTYPE &
                                      , DOUBLE_NARRAY, t_narray            &
                                      , COPY_NARRAY, INIT_NARRAY           &
                                      , QDEF_NCVAR

    IMPLICIT NONE
    
    ! I/O
    INTEGER,            INTENT(OUT)             :: status
    TYPE(t_geohybgrid), INTENT(IN)              :: ggrid
    TYPE(t_scrip_grid), OPTIONAL, INTENT(INOUT) :: sgrid
    TYPE(t_scrip_grid), OPTIONAL, POINTER       :: psgrid
    INTEGER,            OPTIONAL, INTENT(OUT)   :: ID
    LOGICAL,            OPTIONAL, INTENT(IN)    :: l_set_ranges
    ! KEEP llab for debugging
    CHARACTER(LEN=9),   OPTIONAL, INTENT(IN)    :: llab
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr='geohyb2scripgrid'
    INTEGER                               :: rank
    INTEGER, DIMENSION(2)                 :: locdims
    INTEGER                               :: gsize
    LOGICAL,  DIMENSION(:),   ALLOCATABLE :: lmask
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: clon
    REAL(dp), DIMENSION(:),   ALLOCATABLE :: clat
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: corlon
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: corlat
    INTEGER, SAVE            :: callnum = 0
    REAL(dp)                 :: fac
    REAL(dp)                 :: dlon, dlat1, dlat2
    CHARACTER(GRD_MAXSTRLEN) :: unit = '' 
    INTEGER                  :: i,j,n, ix, lc ! loop indices
    INTEGER, DIMENSION(:), POINTER :: ivec => NULL()
    ! mids of longitudes etc.defined:
    LOGICAL                  :: l_mids  = .FALSE.
    LOGICAL                  :: l_interfaces  = .FALSE.
    LOGICAL                  :: l_curvilinear = .FALSE.
    LOGICAL                  :: l_unstructured = .FALSE.
    INTEGER, DIMENSION(2)    :: vdim
    TYPE(t_geohybgrid)       :: gi
    TYPE(t_scrip_grid), POINTER :: pgs => NULL()
    TYPE(t_scrip_grid)       :: gs
    INTEGER                  :: SID = -99
    LOGICAL                  :: llsetr
    TYPE(t_narray)           :: llonm
    TYPE(t_narray)           :: llatm
    TYPE(t_narray)           :: lloni
    TYPE(t_narray)           :: llati
    TYPE(t_narray)           :: vdmask
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tcorlon 
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: tcorlat 
    INTEGER                  :: jll, jur, ill, ilr, iul, iur
    INTEGER                  :: corners
    INTEGER                  :: len

    corners = ggrid%corners ! for quadrilateral grids = 4

    IF (PRESENT(l_set_ranges)) THEN
       llsetr= l_set_ranges
    ELSE
       llsetr= .FALSE.
    ENDIF

    unit = '' 
    SID  = -99

    IF (ASSOCIATED(pgs)) THEN
       NULLIFY(pgs)
    ENDIF
    IF (PRESENT(psgrid)) THEN
       NULLIFY(psgrid)
    ENDIF
    CALL INIT_NARRAY(llonm)
    CALL INIT_NARRAY(lloni)
    CALL INIT_NARRAY(llatm)
    CALL INIT_NARRAY(llati)
    
    l_mids        = .FALSE.
    l_interfaces  = .FALSE.
    l_curvilinear = .FALSE.
    l_unstructured = .FALSE.
   
    status = 3999 ! UNKNOWN ERROR

    CALLNUM = CALLNUM + 1
    CALL INIT_GEOHYBGRID(gi)

    CALL COPY_GEOHYBGRID(gi, ggrid)

    CALL COMPLETE_GEOHYBGRID(gi)

    ! initialise grid information as required by SCRIP
    ! determine unit of lon /lat
    ! assume same unit for lon and lat
    !   check if grid mids are defined

    IF (gi%lonm%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%lonm%natts
          IF (TRIM(gi%lonm%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%lonm%att(i)%dat%vc)
                unit(j:j)=gi%lonm%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_mids        = .TRUE.
       ! FORCE gi%lonm%dat%vd to be associated from now on
       CALL COPY_NARRAY(llonm, gi%lonm%dat)
       CALL COPY_NARRAY(llatm, gi%latm%dat)
    ENDIF
    IF (gi%loni%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%loni%natts
          IF (TRIM(gi%loni%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%loni%att(i)%dat%vc)
                unit(j:j)=gi%loni%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_interfaces  = .TRUE.
       CALL COPY_NARRAY(lloni, gi%loni%dat)
       CALL COPY_NARRAY(llati, gi%lati%dat)
    ENDIF
    
    ! just search curvilinear and unstructured  fields if local llon 
    ! and llat not yet defined
    IF (.NOT. l_mids .AND. .NOT. l_interfaces) THEN
    IF  (gi%clonm%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%clonm%natts
          IF (TRIM(gi%clonm%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%clonm%att(i)%dat%vc)
                unit(j:j)=gi%clonm%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_mids        = .TRUE.
       l_curvilinear = .TRUE.
       CALL COPY_NARRAY(llonm, gi%clonm%dat)
       CALL COPY_NARRAY(llatm, gi%clatm%dat)
    END IF
    IF (gi%cloni%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%cloni%natts
          IF (TRIM(gi%cloni%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%cloni%att(i)%dat%vc)
                unit(j:j)=gi%cloni%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_interfaces  = .TRUE.
       l_curvilinear = .TRUE.
       CALL COPY_NARRAY(lloni, gi%cloni%dat)
       CALL COPY_NARRAY(llati, gi%clati%dat)
    ENDIF
    
    ! for the calculation of the corners of a curvi-linear grid
    ! the interfaces must be known
    IF (l_curvilinear .AND. .NOT. l_interfaces) THEN
       status = 3021 ! for curvi-linear grids interfaces are required for 
                     ! corner calculation
       RETURN
    END IF

    IF  (gi%ulonm%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%ulonm%natts
          IF (TRIM(gi%ulonm%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%ulonm%att(i)%dat%vc)
                unit(j:j)=gi%ulonm%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_mids         = .TRUE.
       l_unstructured = .TRUE.
       CALL COPY_NARRAY(llonm, gi%ulonm%dat)
       CALL COPY_NARRAY(llatm, gi%ulatm%dat)
    END IF
    IF (gi%uloni%xtype /= NULL_XTYPE) THEN
       DO i = 1, gi%uloni%natts
          IF (TRIM(gi%uloni%att(i)%name) == 'units') THEN
             DO j=1,SIZE(gi%uloni%att(i)%dat%vc)
                unit(j:j)=gi%uloni%att(i)%dat%vc(j)
             END DO
             EXIT
          ENDIF
       END DO
       l_interfaces   = .TRUE.
       l_unstructured = .TRUE.
       CALL COPY_NARRAY(lloni, gi%uloni%dat)
       CALL COPY_NARRAY(llati, gi%ulati%dat)
    ENDIF
    ! for unstructured grids the corners must be known
    IF (l_unstructured .AND. .NOT. l_interfaces) THEN
       status = 3021 ! for unstructured grids corners (here defined as 
                     ! interfaces) are required
       RETURN
    END IF
    ENDIF

   IF (ASSOCIATED(llonm%vr)) THEN
       CALL DOUBLE_NARRAY(llonm)
    END IF
    IF(ASSOCIATED(llatm%vr)) THEN
       CALL DOUBLE_NARRAY(llatm)
    END IF
    IF (ASSOCIATED(lloni%vr)) THEN
       CALL DOUBLE_NARRAY(lloni)
    END IF
    IF (ASSOCIATED(llati%vr)) THEN
       CALL DOUBLE_NARRAY(llati)
    END IF

    IF (QDEF_NCVAR(ggrid%imask).AND. l_unstructured) THEN
       CALL COPY_NARRAY(vdmask,ggrid%imask%dat)
       CALL DOUBLE_NARRAY(vdmask)
       len = SIZE(vdmask%vd)
       IF (PRODUCT(llonm%dim) /= len .OR. &
            PRODUCT(llatm%dim) /= len )  THEN
            status = 3025 ! 'SCRIP: mask and lon/lat arrays differ in longitude
            RETURN
       END IF
       DO i = 1, SIZE(llonm%vd)
          llonm%vd(i) = llonm%vd(i) * vdmask%vd(i)
          llatm%vd(i) = llatm%vd(i) * vdmask%vd(i)
       END DO
       DO n = 1, ggrid%corners
          DO i = 1, SIZE(llonm%vd)
             lloni%vd(i + (n-1)*len) =  lloni%vd(i + (n-1)*len) * vdmask%vd(i)
             llati%vd(i + (n-1)*len) =  llati%vd(i + (n-1)*len) * vdmask%vd(i)
          END DO         
       END DO
    END IF

    IF (llsetr) THEN
       CALL SET_RANGES(status)
       IF (status /= 0) RETURN
    ENDIF
    CALL debug('grid_trafo_scrp','LON_RANGES Lower bound '//llab &
         , lon_ranges(1), substr)
    ! determine ...
    !  ... grid size
    !  .... a) determine dimension length of longitude (locdims(1)) and 
    !          latitude (locdims(2)
    if_unstruc1: IF (.NOT. l_unstructured) THEN
       if_curvi1: IF (.NOT. l_curvilinear) THEN
          if_mids1a: IF (l_mids) THEN 
             IF (gi%lonm%ndims /= 1 ) THEN
                status = 3004 ! ARRAY DIMENSIONS FOR REGULAR GRID SHOULD BE 1
                RETURN
             ENDIF
             IF (gi%lonm%ndims /= gi%latm%ndims) THEN
                status = 3001 ! ARRAY DIMENSIONS FOR LON and LAT are not conform
                RETURN
             ENDIF
             locdims(1) = SIZE(llonm%vd)
             locdims(2) = SIZE(llatm%vd)
          ELSE IF (l_interfaces) THEN
             IF (gi%loni%ndims /= 1 ) THEN
                status = 3004 ! ARRAY DIMENSIONS FOR REGULAR GRID SHOULD BE 1
                RETURN
             ENDIF
             IF (gi%loni%ndims /= gi%lati%ndims) THEN
                status = 3001 ! ARRAY DIMENSIONS FOR LON and LAT are not conform
                RETURN
             ENDIF
             locdims(1) = SIZE(lloni%vd)-1
             locdims(2) = SIZE(llati%vd)-1
          ELSE
             status = 3005 ! NEITHER MIDPOINTS NOR INTERFACES ARE DEFINED
             RETURN
          ENDIF if_mids1a
       ELSE ! if_curvi1
          if_mids1b: IF (l_mids) THEN 
             IF (gi%clonm%ndims /= 2 ) THEN
                ! ARRAY DIMENSIONS FOR CURVILINEAR GRID SHOULD BE 2
                status = 3006
                RETURN
             ENDIF
             IF (gi%clonm%ndims /= gi%clatm%ndims) THEN
                status = 3001 ! ARRAY DIMENSIONS FOR LON and LAT are not conform
                RETURN
             ENDIF
             locdims(1) = gi%clonm%dim(1)%len 
             locdims(2) = gi%clonm%dim(2)%len 
          ELSE IF (l_interfaces) THEN
             IF (gi%cloni%ndims /= 2 ) THEN
                ! ARRAY DIMENSIONS FOR CURVILINEAR GRID SHOULD BE 2
                status = 3004
                RETURN
             ENDIF
             IF (gi%cloni%ndims /= gi%clati%ndims) THEN
                status = 3001 ! ARRAY DIMENSIONS FOR LON and LAT are not conform
                RETURN
             ENDIF
             locdims(1) = gi%cloni%dim(1)%len -1 
             locdims(2) = gi%cloni%dim(2)%len -1 
          ELSE
             status = 3005 ! NEITHER MIDPOINTS NOR INTERFACES ARE DEFINED
             RETURN
          ENDIF if_mids1b
       END IF if_curvi1
       gsize = locdims(1)*locdims(2)
    ELSE ! if_unstruc1 
       locdims(1) = gi%ulonm%dim(1)%len
       locdims(2) = gi%ulonm%dim(2)%len
       gsize = locdims(1)*locdims(2)
    END IF if_unstruc1 

    rank = 2 
    !  ... lmask
    ALLOCATE(lmask(gsize))
    lmask(:) = .TRUE. ! all points take part in remapping 
    i = gsize 
    
    IF (QDEF_NCVAR(gi%imask)) THEN
       DO n = 1, SIZE(gi%imask%dat%vi)
          IF (gi%imask%dat%vi(n) == 0) THEN
             lmask(n) = .FALSE.
             i = i - 1 
          END IF
       END DO
    END IF

    if_unstruc2: IF (.NOT. l_unstructured) THEN
       ! centered /corners of longitude /latitude
       ALLOCATE(clon(gsize))
       ALLOCATE(clat(gsize))
       ALLOCATE(corlon(4,gsize))
       ALLOCATE(corlat(4,gsize))
       ALLOCATE(tcorlon(4,gsize))
       ALLOCATE(tcorlat(4,gsize))
       clon    = 0._dp  
       clat    = 0._dp
       corlon  = 0._dp
       corlat  = 0._dp
       tcorlon = 0._dp
       tcorlat = 0._dp
       
       ! determine conversion factor degree => radians
       IF (TRIM(unit) == 'radians') THEN
          fac = 1._dp
       ELSE
          ! assume degrees
          fac =DTR
       END IF

       if_curvi2: IF (.NOT. l_curvilinear) THEN
          if_mids2a: IF (l_mids) THEN
             DO i = 1, locdims(1)
                DO j = 1, locdims(2)
                   n = (j-1) * locdims(1) + i
                   clon(n) = llonm%vd(i) * fac 
                   clat(n) = llatm%vd(j) * fac
                END DO
             END DO
          ELSE ! if_mids2a
             DO i = 1, locdims(1)
                DO j = 1, locdims(2)
                   n = (j-1) * locdims(1) + i
                   clon(n) = (lloni%vd(i)+lloni%vd(i+1))/2._dp * fac
                   clat(n) = (llati%vd(j)+llati%vd(j+1))/2._dp * fac
                END DO
             END DO
          END IF if_mids2a
       ELSE ! if_curvi2
          if_mids2b: IF (l_mids) THEN
             !vdim(1) = gi%clonm%dim(1)%len
             !vdim(2) = gi%clonm%dim(2)%len
             DO n = 1, gsize
                clon(n) = llonm%vd(n) * fac 
                clat(n) = llatm%vd(n) * fac
             END DO
          ELSE IF (l_interfaces) THEN
             status = 3023 ! for curvi-linear grids mids are required for 
                           ! center calculation
             RETURN
          END IF if_mids2b
       END IF if_curvi2

       ! Determine corner: 
       if_curvi3: IF (.NOT. l_curvilinear) THEN
          if_intf3a: IF (l_interfaces) THEN

             ! ... a) interfaces are provided:
             IF (llati%vd(1) <  llati%vd(2)) THEN
                jll = 0
                jur = 1
             ELSE
                jll = 1
                jur = 0
             END IF
             DO j = 1, locdims(2)
                DO i = 1, locdims(1)
                   n = (j-1) * locdims(1) + i
                   ! counter clockwise start lower left
                   ! lower left corner
                   tcorlon(1,n) = lloni%vd(i)  
                   tcorlat(1,n) = llati%vd(j+jll) 
                   ! lower right corner
                   tcorlon(2,n) = lloni%vd(i+1) 
                   tcorlat(2,n) = tcorlat(1,n)
                   ! upper right corner
                   tcorlon(3,n) = tcorlon(2,n) 
                   tcorlat(3,n) = llati%vd(j+jur) 
                   ! upper left corner
                   tcorlon(4,n) = tcorlon(1,n)
                   tcorlat(4,n) = tcorlat(3,n)
                END DO
             END DO

          ELSE ! if_intf3a
             
             ! ... b) calculate corners from centers 
             !           (assuming equidistance for dlon)
             
             dlon = ABS(llonm%vd(2) - llonm%vd(1))
             DO j = 1, locdims(2)
                IF (j == 1) THEN
                   dlat1 = ABS(llatm%vd(2) - llatm%vd(1))
                ELSE
                   dlat1 = ABS(llatm%vd(j) - llatm%vd(j-1))
                ENDIF
                IF (j == locdims(2)) THEN
                   dlat2 = ABS(llatm%vd(j) - llatm%vd(j-1))
                ELSE
                   dlat2 = ABS(llatm%vd(j+1) - llatm%vd(j))
                ENDIF
                DO i = 1, locdims(1)
                   n = (j-1) * locdims(1) + i
                   ! lower left corner
                   tcorlon(1,n) = (llonm%vd(i) - 0.5_dp * dlon)  
                   tcorlat(1,n) = (llatm%vd(j) - 0.5_dp * dlat1) 
                   ! lower right corner
                   tcorlon(2,n) = (llonm%vd(i) + 0.5_dp * dlon)   
                   tcorlat(2,n) = tcorlat(1,n)
                   ! upper right corner
                   tcorlon(3,n) = tcorlon(2,n)   
                   tcorlat(3,n) = llatm%vd(j) + 0.5_dp * dlat2 
                   ! upper left corner
                   tcorlon(4,n) = tcorlon(1,n)
                   tcorlat(4,n) = tcorlat(3,n)
                   jur = 10
                   if_rgempty: IF (gi%ranges(2,1) /= RGEMPTY .AND. &
                        gi%ranges(2,2) /= RGEMPTY) THEN
                      IF (j == 1) THEN
                         IF (llatm%vd(1) < llatm%vd(2)) THEN
                            tcorlat(1,n) = MINVAL(gi%ranges(2,:))
                            tcorlat(2,n) = MINVAL(gi%ranges(2,:))
                            jur = 5
                         ELSE
                            tcorlat(3,n) = MAXVAL(gi%ranges(2,:))
                            tcorlat(4,n) = MAXVAL(gi%ranges(2,:))
                            jur = 6
                         END IF
                      ELSE IF (j == locdims(2)) THEN
                         IF (llatm%vd(1) < llatm%vd(2)) THEN
                            tcorlat(3,n) = MAXVAL(gi%ranges(2,:))
                            tcorlat(4,n) = MAXVAL(gi%ranges(2,:))
                            jur = 7
                         ELSE
                            tcorlat(1,n) = MINVAL(gi%ranges(2,:))
                            tcorlat(2,n) = MINVAL(gi%ranges(2,:))
                            jur = 8
                         END IF
                      END IF
                   ENDIF if_rgempty
               END DO
             END DO
             
          ENDIF if_intf3a

          DO n= 1, SIZE(tcorlon,2)

             DO ix = 1,4
                IF (tcorlon(ix,n) >= lon_ranges(2)) &
                     tcorlon(ix,n) = tcorlon(ix,n) - 360._dp
                IF (tcorlon(ix,n) <lon_ranges(1)) &
                     tcorlon(ix,n) = tcorlon(ix,n) + 360._dp
             END DO

             ! degree => radiant
             corlon(:,n) = tcorlon(:,n) * fac
             corlat(:,n) = tcorlat(:,n) * fac
             
             DO ix = 1,4
                IF (corlat(ix,n) > pi / 2._dp)   corlat(ix,n) = pi / 2._dp
                IF (corlat(ix,n) < - pi / 2._dp) corlat(ix,n) = - pi / 2._dp
             END DO
          END DO
          
       ELSE ! ifcurvi3

          ! ************************************************************
          ! corner calculation for curvi-linear grids
          ! ************************************************************
          
          if_intf3b: IF (l_interfaces) THEN

             ! ... a) interfaces are provided:
             vdim(1) = gi%cloni%dim(1)%len - 1 
             vdim(2) = gi%cloni%dim(2)%len - 1
             
             DO n = 1, SIZE(tcorlon,2)
                ! get element vector for mids
                CALL ELEMENT(vdim,n,ivec) 
                ! calculate position in interface array
                iul = n  + vdim(1) + 1 + ivec(2) - 1
                iur = n  + vdim(1) + 1 + ivec(2)    
                ilr = n                + ivec(2)
                ill = n                + ivec(2) - 1

                IF (llati%vd(ill) > llati%vd(iul)) THEN
                   ! switch upper and lower
                   ! enforce coutner-clockwise
                   ilr = iul
                   iul = n  + ivec(2)
                ENDIF
                ! lower left corner
                tcorlon(1,n) = lloni%vd(ill) 
                tcorlat(1,n) = llati%vd(ill) 
                ! lower right corner
                tcorlon(2,n) = lloni%vd(ilr) 
                tcorlat(2,n) = llati%vd(ilr) 
                ! upper right corner
                tcorlon(3,n) = lloni%vd(iur) 
                tcorlat(3,n) = llati%vd(iur) 
                ! upper left corner
                tcorlon(4,n) = lloni%vd(iul) 
                tcorlat(4,n) = llati%vd(iul) 
                
                DEALLOCATE(ivec, STAT=status)
                NULLIFY(ivec)
             END DO
             DO n= 1, SIZE(tcorlon,2)
                
                DO ix = 1,4
                   IF (tcorlon(ix,n) >= lon_ranges(2)) &
                        tcorlon(ix,n) = tcorlon(ix,n) - 360._dp
                   IF (tcorlon(ix,n) < lon_ranges(1)) &
                        tcorlon(ix,n) = tcorlon(ix,n) + 360._dp
                END DO
                
                ! degree => radiant
                corlon(:,n) = tcorlon(:,n) * fac
                corlat(:,n) = tcorlat(:,n) * fac
                
                DO ix = 1,4
                   IF (corlat(ix,n) > pi / 2._dp)   corlat(ix,n) = pi / 2._dp
                   IF (corlat(ix,n) < - pi / 2._dp) corlat(ix,n) = - pi / 2._dp
                END DO
             END DO
          ELSE !if_intf3b

             ! ... b) calculate corners from centers and dlon 
             ! I think, that is not possible for curvilinear grids
             ! with reasonable effort whoever likes to implement it ....
             ! until then:
             ! This does not work properly  => ERROR

             status = 3021
             RETURN
          END IF if_intf3b

       END IF if_curvi3

    ELSE ! if_unstruc2
       ALLOCATE(clon(gsize))
       ALLOCATE(clat(gsize))
       ALLOCATE(corlon(corners,gsize))
       ALLOCATE(corlat(corners,gsize))
       ALLOCATE(tcorlon(corners,gsize))
       ALLOCATE(tcorlat(corners,gsize))
       clon    = 0._dp  
       clat    = 0._dp
       corlon  = 0._dp
       corlat  = 0._dp
       tcorlon = 0._dp
       tcorlat = 0._dp
       
       ! determine conversion factor degree => radians
       IF (TRIM(unit) == 'radians') THEN
          fac = 1._dp
       ELSE
          IF (corners == 3) THEN
             fac = DTR_icon
          ELSE
             ! assume degrees
             fac = DTR
          END IF
       END IF
       
       vdim(1)  = gsize !gi%uloni%dim(1)%len * gi%uloni%dim(2)%len
       vdim(2)  = gi%uloni%dim(3)%len

       DO n = 1, gsize
          clon(n) = llonm%vd(n) * fac 
          clat(n) = llatm%vd(n) * fac
          DO lc =  1, corners
             tcorlon(lc,n) = lloni%vd(n + vdim(1)*(lc-1))
             tcorlat(lc,n) = llati%vd(n + vdim(1)*(lc-1))

             IF (tcorlon(lc,n) > lon_ranges(2)) &
                  tcorlon(lc,n) = tcorlon(lc,n) - 360._dp
             IF (tcorlon(lc,n) < lon_ranges(1) .AND. tcorlon(lc,n) > -360._dp) &
                  tcorlon(lc,n) = tcorlon(lc,n) + 360._dp
          END DO
          
          ! degree => radiant
          corlon(:,n) = tcorlon(:,n) * fac
          corlat(:,n) = tcorlat(:,n) * fac
       END DO
       
    END IF if_unstruc2

  CALL DEFINE_SCRIPGRID(status, rank, gsize, corners     &
       , locdims(1:2), lmask, clon, clat, corlon, corlat &
       , id=sid, grid=gs ,pgrid=pgs)
  IF (status /= 0 .AND. status /= 1) RETURN
  
  IF (PRESENT(sgrid)) THEN
     CALL COPY_SCRIPGRID(sgrid,gs)
  END IF
  IF (PRESENT(ID)) THEN
     ID = SID
  ENDIF
  IF (PRESENT(psgrid)) THEN
     psgrid => pgs
  ENDIF

   ! DEALLOCATE TERMPORARY ARRAYS
  CALL INIT_GEOHYBGRID(gi)
  CALL INIT_SCRIPGRID(gs)

  CALL INIT_NARRAY(llonm)
  CALL INIT_NARRAY(lloni)
  CALL INIT_NARRAY(llatm)
  CALL INIT_NARRAY(llati)

  NULLIFY(pgs)

  IF (ALLOCATED(lmask))    DEALLOCATE(lmask)
  IF (ALLOCATED(clon))     DEALLOCATE(clon)
  IF (ALLOCATED(clat))     DEALLOCATE(clat)
  IF (ALLOCATED(corlon))   DEALLOCATE(corlon)
  IF (ALLOCATED(corlat))   DEALLOCATE(corlat)
  IF (ALLOCATED(tcorlon))  DEALLOCATE(tcorlon)
  IF (ALLOCATED(tcorlat))  DEALLOCATE(tcorlat)

  STATUS = 0

CONTAINS
  
  SUBROUTINE SET_RANGES(status)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    ! LOCAL
    REAL(dp) :: minv
    REAL(dp) :: maxv

    minv = -9999._dp
    maxv = -9999._dp

    ! exclude masked point for range calculation
    IF (ASSOCIATED(lloni%vd)) THEN
       minv = MINVAL(lloni%vd)
       maxv = MAXVAL(lloni%vd)
    END IF
    IF (.NOT. l_curvilinear .AND. .NOT. l_unstructured) THEN
       IF (ASSOCIATED(llonm%vd)) THEN
          minv = MINVAL(llonm%vd)
          maxv = MAXVAL(llonm%vd)
       END IF
    ENDIF
    IF (minv < -9360._dp .OR. maxv < -9360._dp) THEN
       status = 3024 ! no lonm loni associated
       RETURN
    END IF
    IF (minv < -180._dp) THEN
       write (0,*) 'ERROR 3018', minv, maxv
       status = 3018 ! longitude range to small < -180.
       RETURN
    END IF
    IF (maxv > 360._dp) THEN
       write (0,*) 'ERROR 3019', minv, maxv
       status = 3019 ! longitude range to large > 360.
       RETURN
    END IF
    IF (maxv - minv > 360.0_dp) THEN
       write (0,*) 'ERROR 3020', minv, maxv
       status = 3020 ! longitude range to wide > 360.
       RETURN
    END IF

    IF (l_unstructured) THEN
       IF (maxv < 0._dp) THEN
          lon_ranges(1) = -180.0_dp
          lon_ranges(2) = 180.0_dp
       ELSE 
          lon_ranges(1) = 0.0_dp
          lon_ranges(2) = 360.0_dp
       ENDIF
    ELSE
       IF (minv < 0._dp) THEN
          lon_ranges(1) = -180.0_dp
          lon_ranges(2) = 180.0_dp
       ELSE 
          lon_ranges(1) = 0.0_dp
          lon_ranges(2) = 360.0_dp
       ENDIF
    END IF
    status = 0

  END SUBROUTINE SET_RANGES

  END SUBROUTINE GEOHYB2SCRIPGRID
  ! ---------------------------------------------------------------------------

  ! ===========================================================================
  SUBROUTINE DEFINE_SCRIPGRID (status, rank, size, corners &
       , dims, lmask, clon, clat, corlon, corlat, id, grid, pgrid)

    USE messy_main_grid_netcdf, ONLY: INIT_NCDIM

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)           :: status
    INTEGER, INTENT(IN)            :: rank
    INTEGER, INTENT(IN)            :: size
    INTEGER, INTENT(IN)            :: corners
    ! FIELD in
    INTEGER,  DIMENSION(:),   INTENT(IN) :: dims
    LOGICAL,  DIMENSION(:),   INTENT(IN) :: lmask
    REAL(dp), DIMENSION(:),   INTENT(IN) :: clon
    REAL(dp), DIMENSION(:),   INTENT(IN) :: clat
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: corlon
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: corlat

    INTEGER,            INTENT(  OUT), OPTIONAL :: ID
    TYPE(t_scrip_grid), INTENT(INOUT), OPTIONAL :: grid
    TYPE(t_scrip_grid), POINTER,       OPTIONAL :: pgrid

    ! LOCAL
    INTEGER                          :: ii  ! loop index
    TYPE(t_scrip_grid_list), POINTER :: gi     => NULL()
    TYPE(t_scrip_grid_list), POINTER :: ge     => NULL()
    TYPE(t_scrip_grid),      POINTER :: lgrid  => NULL()
    INTEGER                          :: cid

    status = 99

    ! GOTO END OF LIST
    gi => SCRIPGRIDLIST
    DO
       IF (.NOT. ASSOCIATED(gi)) EXIT
       lgrid => gi%this
       CALL COMPARE_TO_SCRIPGRID(cid, lgrid, rank, size, corners &
            , dims, lmask, clon, clat, corlon, corlat)
       IF (cid > 0) THEN
          IF (PRESENT(ID))    ID   =  cid
          IF (PRESENT(grid)) THEN
             CALL COPY_SCRIPGRID(grid,gi%this)
          ENDIF
          IF (PRESENT(pgrid)) THEN
             pgrid => gi%this
          END IF

          ! CLEAN 
          NULLIFY(lgrid)
          NULLIFY(gi)
          NULLIFY(ge)

          status = 01
          RETURN
       ENDIF

       ge => gi
       gi => gi%next
    END DO

    ALLOCATE(gi)
    NULLIFY(gi%next)
    IF (.NOT. ASSOCIATED(SCRIPGRIDLIST)) THEN
       SCRIPGRIDLIST => gi          ! SET POINTER TO FIRST GRID
    ELSE
       ge%next => gi                ! SET NEXT POINTER OF LAST GRID
       !                            ! TO NEW GRID
    END IF

    ! SET VALUES
    ! COUNT AND SET ID
    NSGRID     = NSGRID + 1
    gi%this%id = NSGRID

    ! DIMENSION GRID ... 
    gi%this%size    = size
    gi%this%rank    = rank
    gi%this%corners = corners
    ! ... and ALLOCATE and fill entries:
    ALLOCATE(gi%this%dim(rank))
    DO ii=1, rank
       CALL INIT_NCDIM(gi%this%dim(ii))
       gi%this%dim(ii)%len = dims(ii)
    END DO
    ! TODO: avoid hardcoding !!
    gi%this%dim(1)%name=TRIM('lon')
    gi%this%dim(2)%name=TRIM('lat')
    ALLOCATE(gi%this%lmask(size))
    gi%this%lmask = lmask
    ALLOCATE(gi%this%center_lon(size))
    gi%this%center_lon = clon
    ALLOCATE(gi%this%center_lat(size))
    gi%this%center_lat = clat
    ALLOCATE(gi%this%corner_lon(corners,size))
    gi%this%corner_lon = corlon
    ALLOCATE(gi%this%corner_lat(corners,size))
    gi%this%corner_lat = corlat

    IF (PRESENT(ID))    ID   =  gi%this%id
    IF (PRESENT(grid)) THEN
       CALL COPY_SCRIPGRID(grid, gi%this)
    ENDIF
    IF (PRESENT(pgrid)) THEN
       pgrid => gi%this
    END IF

    ! CLEAN 
    NULLIFY(lgrid)
    NULLIFY(gi)
    NULLIFY(ge)

    status = 0

  END SUBROUTINE DEFINE_SCRIPGRID
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE COMPARE_TO_SCRIPGRID(id, grid, rank, size, corners &
       , dims, lmask, clon, clat, corlon, corlat)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)             :: id
    TYPE(t_scrip_grid), POINTER      :: grid
    INTEGER, INTENT(IN)              :: rank
    INTEGER, INTENT(IN)              :: size
    INTEGER, INTENT(IN)              :: corners
    ! FIELD in
    INTEGER,  DIMENSION(rank),         INTENT(IN) :: dims
    LOGICAL,  DIMENSION(size),         INTENT(IN) :: lmask
    REAL(dp), DIMENSION(size),         INTENT(IN) :: clon
    REAL(dp), DIMENSION(size),         INTENT(IN) :: clat
    REAL(dp), DIMENSION(corners,size), INTENT(IN) :: corlon
    REAL(dp), DIMENSION(corners,size), INTENT(IN) :: corlat

    ! LOCAL
    INTEGER :: i,j,ir 

    ! INITIALIZE
    id = -99

    IF (grid%size    /= size)    RETURN
    IF (grid%corners /= corners) RETURN
    IF (grid%rank    /= rank)    RETURN

    DO ir = 1, rank 
       IF (grid%dim(ir)%len /= dims(ir)) RETURN
    ENDDO

    DO i = 1, size
       IF (grid%lmask(i) .NEQV. lmask(i))     RETURN 
       IF (grid%center_lon(i) /= clon(i))     RETURN
       IF (grid%center_lat(i) /= clat(i))     RETURN

       DO j = 1, corners
          IF (grid%corner_lon(j,i) /= corlon(j,i)) RETURN
          IF (grid%corner_lat(j,i) /= corlat(j,i)) RETURN
       END DO
    END DO

    ID = grid%id

  END SUBROUTINE COMPARE_TO_SCRIPGRID
  ! ===========================================================================
    
  ! ===========================================================================
  SUBROUTINE INIT_SCRIPGRID(grid)

    USE messy_main_grid_netcdf, ONLY: INIT_NCDIM

    IMPLICIT NONE

    ! I/O
    TYPE(t_scrip_grid), INTENT(INOUT):: grid

    ! LOCAL
    INTEGER   :: i

    grid%id      = -99
    grid%size    = 0
    grid%corners = 4
    grid%rank    = 0

    IF (ASSOCIATED(grid%dim)) THEN
       DO i=1,SIZE(grid%dim)
          CALL INIT_NCDIM(grid%dim(i))
       END DO
       DEALLOCATE(grid%dim)
       NULLIFY(grid%dim)
    END IF

    IF (ASSOCIATED(grid%lmask)) THEN
       DEALLOCATE(grid%lmask)
       NULLIFY(grid%lmask)
    END IF

    IF (ASSOCIATED(grid%center_lon)) THEN
       DEALLOCATE(grid%center_lon)
       NULLIFY(grid%center_lon)
    END IF
    IF (ASSOCIATED(grid%center_lat)) THEN
       DEALLOCATE(grid%center_lat)
       NULLIFY(grid%center_lat)
    END IF
    IF (ASSOCIATED(grid%corner_lon)) THEN
       DEALLOCATE(grid%corner_lon)
       NULLIFY(grid%corner_lon)
    END IF
    IF (ASSOCIATED(grid%corner_lat)) THEN
       DEALLOCATE(grid%corner_lat)
       NULLIFY(grid%corner_lat)
    END IF
  END SUBROUTINE INIT_SCRIPGRID
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE COPY_SCRIPGRID(dgrid, sgrid)

    USE messy_main_grid_netcdf, ONLY: COPY_NCDIM

    IMPLICIT NONE

    ! I/O
    TYPE(t_scrip_grid), INTENT(IN)   :: sgrid
    TYPE(t_scrip_grid), INTENT(INOUT):: dgrid

    ! LOCAL
    INTEGER   :: i

    CALL INIT_SCRIPGRID(dgrid)
    dgrid%id      = sgrid%id
    dgrid%size    = sgrid%size
    dgrid%corners = sgrid%corners
    dgrid%rank    = sgrid%rank

    ALLOCATE(dgrid%dim(SIZE(sgrid%dim)))
    DO i=1,SIZE(sgrid%dim)
       CALL COPY_NCDIM(dgrid%dim(i),sgrid%dim(i))
    END DO

    ALLOCATE(dgrid%lmask(SIZE(sgrid%lmask)))
    dgrid%lmask = sgrid%lmask

    ALLOCATE(dgrid%center_lon (SIZE(sgrid%center_lon)))
    dgrid%center_lon = sgrid%center_lon
    ALLOCATE(dgrid%center_lat (SIZE(sgrid%center_lat)))
    dgrid%center_lat = sgrid%center_lat
    ALLOCATE(dgrid%corner_lon(&
         SIZE(sgrid%corner_lon,1),SIZE(sgrid%corner_lon,2) ))
    dgrid%corner_lon = sgrid%corner_lon
    ALLOCATE(dgrid%corner_lat(&
         SIZE(sgrid%corner_lat,1),SIZE(sgrid%corner_lat,2) ))
    dgrid%corner_lat = sgrid%corner_lat

  END SUBROUTINE COPY_SCRIPGRID
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE INIT_SCRIPDATA(sdata)

    IMPLICIT NONE

    TYPE(t_scrip_data), INTENT(INOUT) :: sdata

    sdata%id   = 0
    sdata%map_type = map_type_conserv
    sdata%norm_opt = norm_opt_none
    sdata%luse_sgrd_area = .FALSE.
    sdata%luse_dgrd_area = .FALSE.

    CALL INIT_SCRIP_WEIGHTS(sdata%wghts)
    CALL INIT_SCRIPGRID(sdata%sgrd)
    CALL INIT_SCRIPGRID(sdata%dgrd)

  END SUBROUTINE INIT_SCRIPDATA
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE DEFINE_SCRIPDATA(status, SDAT_ID, maptype, normopt &
       , luse_area, dgrd, sgrd, garea2, PSD)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(OUT) :: SDAT_ID
    INTEGER,          INTENT(IN)  :: maptype
    INTEGER,          INTENT(IN)  :: normopt
    LOGICAL,          INTENT(IN)  :: luse_area

    TYPE(t_scrip_grid), POINTER :: dgrd
    TYPE(t_scrip_grid), POINTER :: sgrd

    REAL(dp),      DIMENSION(:), POINTER  :: garea2
    TYPE(t_scrip_data), POINTER, OPTIONAL :: PSD

    ! LOCAL
    TYPE(t_scrip_data_list), POINTER :: gi    => NULL()
    TYPE(t_scrip_data_list), POINTER :: ge    => NULL()
    TYPE(t_scrip_data),      POINTER :: ldata => NULL()
    INTEGER                          :: cid ! file id returned by compare

    status = 3999

    IF (PRESENT(PSD)) NULLIFY(PSD)

    ! GOTO END OF LIST
    gi => SCRIPDATALIST
    DO
       IF (.NOT. ASSOCIATED(gi)) EXIT
       ldata => gi%this
       CALL COMPARE_SCRIPDATA(cid, ldata, maptype, normopt, luse_area &
            , dgrd%id, sgrd%id)
       IF (cid > 0) THEN
          SDAT_ID =  cid
          IF (PRESENT(PSD)) PSD => gi%this

          ! CLEAN 
          NULLIFY(ldata)
          NULLIFY(gi)
          NULLIFY(ge)

          status = 01
          RETURN
       ENDIF

       ge => gi
       gi => gi%next
    END DO

    ! ADD NEW
    ALLOCATE(gi)
    NULLIFY(gi%next)
    IF (.NOT. ASSOCIATED(SCRIPDATALIST)) THEN
       SCRIPDATALIST => gi          ! SET POINTER TO FIRST DATA
    ELSE
       ge%next => gi                ! SET NEXT POINTER OF LAST DATA
       !                            ! TO NEW DATA
    END IF

    ! SET VALUES
    ! COUNT AND SET ID
    NSDATA     = NSDATA + 1
    gi%this%id = NSDATA

    ! SET REMAP OPTIONS
    gi%this%map_type       = maptype
    gi%this%norm_opt       = normopt
    gi%this%luse_dgrd_area = luse_area
    gi%this%luse_sgrd_area = .FALSE. ! TODO: calc source grid area set to true

    ! ... define grids
    gi%this%sgrd => sgrd
    gi%this%dgrd => dgrd

    IF (ASSOCIATED(garea2)) THEN
       ALLOCATE(gi%this%wghts%dstarea(gi%this%dgrd%size))
       gi%this%wghts%dstarea = garea2
       DEALLOCATE(garea2)
    ENDIF

    SDAT_ID = gi%this%id

    IF (PRESENT(PSD)) PSD => gi%this

    ! CLEAN 
    NULLIFY(ldata)
    NULLIFY(gi)
    NULLIFY(ge)

    status = 0

  END SUBROUTINE DEFINE_SCRIPDATA
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE COMPARE_SCRIPDATA(id, data, maptype, normopt, luse_area, did, sid)
    
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)             :: id
    TYPE(t_scrip_data), POINTER      :: data
    INTEGER, INTENT(IN)              :: maptype
    INTEGER, INTENT(IN)              :: normopt
    LOGICAL, INTENT(IN)              :: luse_area
    INTEGER, INTENT(IN)              :: did
    INTEGER, INTENT(IN)              :: sid

    ! INITIALIZE
    ID = -99

    IF (data%map_type   /= maptype)   RETURN
    IF (data%norm_opt   /= normopt)   RETURN
    IF (data%dgrd%id    /= did)       RETURN
    IF (data%sgrd%id    /= sid)       RETURN
    IF (data%luse_dgrd_area  .neqv. luse_area) RETURN

    ID = data%id

  END SUBROUTINE COMPARE_SCRIPDATA
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE LOCATE_SCRIPDATA(status, ID, pdata)
    
    IMPLICIT NONE
    
    ! I/O   
    INTEGER, INTENT(OUT)                  :: status
    INTEGER, INTENT(IN)                   :: ID
    TYPE(t_scrip_data),  POINTER          :: pdata      
    
    ! LOCAL
    TYPE(t_scrip_data_list), POINTER :: gi => NULL()
    TYPE(t_scrip_data_list), POINTER :: ge => NULL()
    LOGICAL                          :: lexists = .FALSE.
    
    lexists = .FALSE.
    status  = 1
    IF (ASSOCIATED(pdata)) THEN
       CALL INIT_SCRIPDATA(pdata)
    ENDIF
    NULLIFY(pdata)

    ! CHECKS
    IF (id <= 0) THEN
       status = 3015  ! INVALID DATA ID
       RETURN
    END IF
    !
    IF (.NOT. ASSOCIATED(SCRIPDATALIST)) THEN
       status = 3009  ! DATA (ID) DOES NOT EXIST
       RETURN
    END IF
    
    gi => SCRIPDATALIST
    DO
       IF (.NOT. ASSOCIATED(gi)) EXIT
       IF (id == gi%this%id) THEN
          lexists = .TRUE.
          EXIT
       END IF
       ge => gi
       gi => ge%next
    END DO
    
    IF (lexists) THEN
       pdata => gi%this
    ELSE
       status = 3013  ! DATA DOES NOT EXIST
       RETURN
    END IF
    
    status = 0
      
  END SUBROUTINE LOCATE_SCRIPDATA
  ! ===========================================================================

  ! ========================================================================
  SUBROUTINE CALC_SCRIPDATA(status, igrid, ogrid, RGT, SCRIP_ID, oarea,PSD &
       , norm_opt_in, map_type_in, label)

    USE messy_main_constants_mem,    ONLY: r_e => radius_earth
    USE messy_main_tools,            ONLY: INT2STR
   
    USE messy_main_grid,             ONLY: T_GEOHYBGRID
    USE messy_main_grid_trafo,       ONLY: RG_INT, RG_EXT, RG_IDX, RG_IXF
    USE messy_main_grid_tools,       ONLY: ALINE_ARRAY

    IMPLICIT NONE

    ! I/O
    INTEGER,            INTENT(OUT) :: status
    TYPE(t_geohybgrid), INTENT(IN)  :: igrid
    TYPE(t_geohybgrid), INTENT(IN)  :: ogrid
    INTEGER,            INTENT(IN)  :: RGT(:)  ! regridding type

    INTEGER,            INTENT(OUT) :: SCRIP_ID
    ! OPTIONAL
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: oarea
    TYPE(t_scrip_data),       POINTER, OPTIONAL :: PSD
    INTEGER,                           OPTIONAL :: norm_opt_in
    INTEGER,                           OPTIONAL :: map_type_in
    CHARACTER(LEN=6),                  OPTIONAL :: label

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'calc_scripdata'
    TYPE(t_scrip_grid),     POINTER :: psigrid  => NULL()  ! scrip input grid
    INTEGER                         :: sigID     ! ID scrip input grid
    TYPE(t_scrip_grid),     POINTER :: psogrid => NULL()   ! scrip output grid
    TYPE(t_scrip_grid)              :: sogrid    ! scrip output grid
    INTEGER                         :: sogID     ! ID scrip output grid
    INTEGER                         :: map_type
    INTEGER                         :: norm_opt
    REAL(dp), DIMENSION(:), POINTER :: gridarea => NULL()
    LOGICAL                         :: luse_area 
    INTEGER                         :: i        ! loop index
    INTEGER, SAVE                   :: CALLS = 0
    CHARACTER(LEN=4)                :: callstr   = ''
    CHARACTER(LEN=6)                :: labstr   = ''

    IF (PRESENT(label)) THEN
       labstr = label
    ELSE
       labstr = 'label_'
    END IF

    status = 3999

    luse_area = .FALSE.
    NULLIFY(psigrid)
    NULLIFY(psogrid)
    IF (PRESENT(PSD)) NULLIFY(PSD)

    CALLS =CALLS + 1

    ! 1. Determine SCRIP output grid
    CALL GEOHYB2SCRIPGRID(status, ogrid, psgrid=psogrid, sgrid=sogrid,ID=sogID &
         , l_set_ranges=.TRUE.,llab=labstr//'G01')
    IF (status /= 0) RETURN

    ! 2. Determine SCRIP input  grid
    CALL GEOHYB2SCRIPGRID(status, igrid, psgrid=psigrid,ID=sigID &
         ,llab=labstr//'G02')

    IF (status /= 0) RETURN

    ! 3. Determine SCRIP switches norm_opt and map_type
    !    TODO : FIND OUT, WHICH norm_opt / map_type combination
    !           corresponds to INT / EXT / IDX / IFX
    !           make map_type switchabel ?
    DO i=2, SIZE(RGT)
       IF (RGT(i) /= RGT(i-1)) THEN
          status = 3007 ! SCRIP can only handle equal regrid types per data set
          RETURN
       END IF
    END DO
    SELECT CASE (RGT(1))
    CASE(RG_INT, RG_IDX, RG_IXF)
       norm_opt = norm_opt_frcarea 
       map_type = map_type_conserv 
       ! norm_opt = norm_opt_dstarea
       ! norm_opt = norm_opt_none  
       !
       ! map_type = map_type_bilinear ! yes 
       ! map_type = map_type_bicubic  ! yes 
       ! map_type = map_type_distwgt  
    CASE(RG_EXT)
       norm_opt = norm_opt_none
       map_type = map_type_conserv 
    END SELECT

    ! OVERWRITE SETTINGs IF EXTERNALLY FORCED
    IF (PRESENT(norm_opt_in)) norm_opt = norm_opt_in
    IF (PRESENT(map_type_in)) map_type = map_type_in

    IF (PRESENT(oarea)) THEN ! op_mm_20131001
       IF (ASSOCIATED(oarea)) THEN
          CALL ALINE_ARRAY(status, oarea, gridarea)
          IF (status /= 0) RETURN

          ! convert unit of destination grid area
          gridarea  = gridarea / (r_e**2)
          luse_area = .TRUE.
       ENDIF
    ELSE
       NULLIFY(gridarea)
    ENDIF

    ! 4. Define SCRIP_DATA (return ID)
    CALL INT2STR(callstr, CALLS, '0') 

    CALL DEFINE_SCRIPDATA( status, SCRIP_ID, map_type, norm_opt &
                       , luse_area, psogrid, psigrid, gridarea, PSD)
    
    IF (status /= 0 .AND. status /= 01) RETURN

    CALL INIT_SCRIPGRID(sogrid)

    NULLIFY(psogrid)
    NULLIFY(psigrid)

    IF (ASSOCIATED(gridarea)) DEALLOCATE(gridarea)

  END SUBROUTINE CALC_SCRIPDATA
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE INIT_SCRIP_WEIGHTS(wght)

    IMPLICIT NONE

    TYPE(t_scrip_weights), INTENT(INOUT) :: wght

    wght%num_links = 0
    wght%ltimedependent =.FALSE.
    IF (ASSOCIATED(wght%weights)) THEN
       DEALLOCATE(wght%weights)
       NULLIFY(wght%weights)
    ENDIF
    
    IF (ASSOCIATED(wght%dstfrac)) THEN
       DEALLOCATE(wght%dstfrac)
       NULLIFY(wght%dstfrac)
    ENDIF
    IF (ASSOCIATED(wght%dstarea)) THEN
       DEALLOCATE(wght%dstarea)
       NULLIFY(wght%dstarea)
    ENDIF
    IF (ASSOCIATED(wght%srcadd)) THEN
       DEALLOCATE(wght%srcadd)
       NULLIFY(wght%srcadd)
    ENDIF
    IF (ASSOCIATED(wght%dstadd)) THEN
       DEALLOCATE(wght%dstadd)
       NULLIFY(wght%dstadd)
    ENDIF

  END SUBROUTINE INIT_SCRIP_WEIGHTS
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE CALC_SCRIP_WEIGHTS(status, psd)

    USE messy_main_grid_trafo_scrp_base,ONLY: remap_vars,   remap_conserv   &
                                           , remap_bilin,   remap_bicub     &
                                           , remap_distwgt, remap_init      &
                                           , remap_dealloc, bounds_calc     &
                                           , grid2_frac, grid1_mask         &
                                           , grid2_mask, grid2_area_in      &
                                           , grid2_area                     &
                                           , grid1_center_lon &
                                           , grid1_center_lat &
                                           , grid2_center_lon &
                                           , grid2_center_lat &
                                           , grid1_corner_lon &
                                           , grid1_corner_lat &
                                           , grid2_corner_lon &
                                           , grid2_corner_lat &
                                           , num_links_map1                  &
                                           , num_wts, wts_map1, map_type     &
                                           , norm_opt, luse_grid_centers     &
                                           , grid1_add_map1, grid2_add_map1  &
                                           , north_thresh, south_thresh      &
                                           , babystep

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)         :: status
    TYPE(t_scrip_data), POINTER  :: psd

    status = 3999

    ! ********************************************************
    ! ********************************************************
    ! NOW THE CALCULATION OF THE INTERPOLATION WEIGHTS starts
    ! ********************************************************
    ! ********************************************************

    babystep = 0.001_dp

    ! set mask
    grid1_mask => PSD%sgrd%lmask
    grid2_mask => PSD%dgrd%lmask

    ! set center_lxx
    grid1_center_lon => PSD%sgrd%center_lon
    grid1_center_lat => PSD%sgrd%center_lat
    grid2_center_lon => PSD%dgrd%center_lon
    grid2_center_lat => PSD%dgrd%center_lat

    ! .. corners
    grid1_corner_lon => PSD%sgrd%corner_lon
    grid1_corner_lat => PSD%sgrd%corner_lat
    grid2_corner_lon => PSD%dgrd%corner_lon
    grid2_corner_lat => PSD%dgrd%corner_lat

    ! ... area
    IF (ASSOCIATED(PSD%wghts%dstarea)) grid2_area_in =>  PSD%wghts%dstarea 
    ! grid1_area_in not needed

    ! ... fraction
    !  ... not yet defined here,  copy after calculation of weights

    norm_opt = PSD%norm_opt
    
    north_thresh = MAX(MAXVAL(grid1_corner_lat),MAXVAL(grid2_corner_lat))
    south_thresh = MIN(MINVAL(grid1_corner_lat),MINVAL(grid2_corner_lat))

    CALL remap_init(PSD%sgrd%corners,PSD%sgrd%dim(1)%len,PSD%sgrd%dim(2)%len &
         ,PSD%dgrd%corners,PSD%dgrd%dim(1)%len,PSD%dgrd%dim(2)%len)

    IF (PSD%map_type == map_type_conserv) THEN
       luse_grid_centers = .false.
    ELSE
       luse_grid_centers = .true.
    END IF
    CALL bounds_calc

    map_type = PSD%map_type
    CALL remap_vars(1)

    SELECT CASE(PSD%map_type)
    CASE(map_type_conserv)
!       luse_grid_centers = .false.
       CALL remap_conserv(status)
       IF (status /= 0) RETURN
    CASE(map_type_bilinear)
!       luse_grid_centers = .true.
       CALL remap_bilin 
    CASE(map_type_bicubic)
!       luse_grid_centers = .true.
       CALL remap_bicub
    CASE(map_type_distwgt)
!       luse_grid_centers = .true.
       CALL remap_distwgt 
    CASE DEFAULT
       status = 3012 ! map_type not known
       RETURN
    END SELECT

    CALL remap_vars(2)

    ! SECURE WEIGHTS
    
    PSD%wghts%num_links = num_links_map1   ! number of links
    PSD%wghts%num_wgts  = num_wts          ! number of weights
    ! TODO IF  copy anyway make weights fourdimensional for channel object
    !      same for frac
    ALLOCATE(PSD%wghts%weights(num_links_map1))
    PSD%wghts%weights = wts_map1(1,:)

    ALLOCATE(PSD%wghts%srcadd(num_links_map1))
    PSD%wghts%srcadd = grid1_add_map1
    ALLOCATE(PSD%wghts%dstadd(num_links_map1))
    PSD%wghts%dstadd = grid2_add_map1

    ! ... fraction 
    ALLOCATE(PSD%wghts%dstfrac(&
         PSD%dgrd%dim(1)%len * PSD%dgrd%dim(2)%len))
    PSD%wghts%dstfrac = grid2_frac

    IF (.NOT. ASSOCIATED( PSD%wghts%dstarea)) THEN
       ALLOCATE( PSD%wghts%dstarea(SIZE(grid2_area)))
       PSD%wghts%dstarea = grid2_area
    ENDIF

    CALL remap_dealloc

    ! CLEAN
    NULLIFY(grid1_mask)
    NULLIFY(grid2_mask)
    NULLIFY(grid1_center_lon)
    NULLIFY(grid1_center_lat)
    NULLIFY(grid2_center_lon)
    NULLIFY(grid2_center_lat)
    NULLIFY(grid1_corner_lon)
    NULLIFY(grid1_corner_lat)
    NULLIFY(grid2_corner_lon)
    NULLIFY(grid2_corner_lat)
    NULLIFY(grid2_area_in)

    status = 0

  END SUBROUTINE CALC_SCRIP_WEIGHTS
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE APPLY_SCRIP_WEIGHTS(status, pvari, pvaro, PSDID)

    USE messy_main_grid_netcdf,          ONLY: t_ncvar, VTYPE_DOUBLE &
                                             , VTYPE_REAL            &
                                             , VTYPE_INT, VTYPE_BYTE
    USE messy_main_grid_trafo_scrp_base, ONLY: zero

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT)                 :: status
    ! INTENT(IN)
    TYPE(t_ncvar), DIMENSION(:), POINTER :: pvari  ! variable to be interpolated
    INTEGER, INTENT(IN)                  :: PSDID
    ! INTENT(INOUT)
    TYPE(t_ncvar), DIMENSION(:), POINTER :: pvaro  ! interpolated fields

    ! LOCAL
    INTEGER                     :: nvars ! number of interpolated variables
    INTEGER                     :: nfree ! number of free dimensions
    INTEGER                     :: ifree ! loop index for loop over free dim.s
    INTEGER                     :: i     ! loop index for loop over variables
    INTEGER                     :: j     ! loop index for loop over num_links
    INTEGER                     :: dshft,sshft ! index shift for free dimensions
    TYPE(t_scrip_data), POINTER :: PSD   => NULL() ! SCRIP DATA 
    
    status = 3999
    nvars = SIZE(pvari)

    ! GET DATA
    CALL LOCATE_SCRIPDATA(status, PSDID, PSD)
    IF (status /= 0) RETURN

    DO i=1,nvars
       ! 1. Determine size of interpolated variable:
       ! - SCRIP only interpolates in the horizontal,
       !    thus all other dimensions are invariant
       ! => deduce number (nfree) of invariant dimensions from input variable,
       ! i.e., overall size devided by horizontal size (= SCRIP SIZE)
       ! 2. allocate destination variable accordingly    
       nfree = 1
       IF (pvari(i)%dat%type == VTYPE_REAL) THEN
          nfree = SIZE(pvari(i)%dat%vr) / PSD%sgrd%size
          pvaro(i)%dat%type = VTYPE_REAL
          ALLOCATE(pvaro(i)%dat%vr(nfree * PSD%dgrd%size))
          pvaro(i)%dat%vr(:) = 0._sp
       ELSE IF (pvari(i)%dat%type == VTYPE_DOUBLE) THEN
          nfree = SIZE(pvari(i)%dat%vd) / PSD%sgrd%size
          pvaro(i)%dat%type = VTYPE_DOUBLE
          ALLOCATE(pvaro(i)%dat%vd(nfree * PSD%dgrd%size))
          pvaro(i)%dat%vd(:) = 0._dp
       ELSE IF (pvari(i)%dat%type == VTYPE_INT) THEN
          nfree = SIZE(pvari(i)%dat%vi) / PSD%sgrd%size
          pvaro(i)%dat%type = VTYPE_DOUBLE
          ALLOCATE(pvaro(i)%dat%vd(nfree * PSD%dgrd%size))
          pvaro(i)%dat%vd(:) = 0._dp
       ELSE IF (pvari(i)%dat%type == VTYPE_BYTE) THEN
          nfree = SIZE(pvari(i)%dat%vb) / PSD%sgrd%size
          pvaro(i)%dat%type = VTYPE_DOUBLE
          ALLOCATE(pvaro(i)%dat%vd(nfree * PSD%dgrd%size))
          pvaro(i)%dat%vd(:) = 0._dp
       ELSE
          !CASE(VTYPE_CHAR)
          !CASE(VTYPE_UNDEF)
          !          CALL RGMSG('SCRIP_apply weights', RGMLE, &
          !               'INTERPOLATION IS ONLY POSSIBLE FOR DATA', .false.)
          !          CALL RGMSG('SCRIP_apply weights', RGMLE, &
          !               'OF TYPE REAL OR DOUBLE PRECISION !', .false.)
          status = 3014 ! VTYPE NOT IMPLEMENTED
          RETURN
       ENDIF

       ! 3. INTERPOLATION LOOP (loop nfree times)
       do_free: DO ifree = 0, nfree-1
          ! shift for each free dimension
          dshft = ifree * PSD%dgrd%size
          sshft = ifree * PSD%sgrd%size
          SELECT CASE(PSD%norm_opt)
          CASE(norm_opt_frcarea)
             IF (pvari(i)%dat%type == VTYPE_REAL) THEN
                DO j = 1, PSD%wghts%num_links
                  pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j)) =      &
                        pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j))  &  
                        + REAL(PSD%wghts%weights(j),sp) *           &
                        pvari(i)%dat%vr(sshft+PSD%wghts%srcadd(j))
                 END DO
             ELSE IF (pvari(i)%dat%type == VTYPE_DOUBLE) THEN
                DO j = 1, PSD%wghts%num_links
                   pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j)) =     &
                        pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j))  &  
                        + PSD%wghts%weights(j)  *                   &
                        pvari(i)%dat%vd(sshft+PSD%wghts%srcadd(j))
               END DO
             ENDIF
          CASE(norm_opt_dstarea)
             IF (pvari(i)%dat%type == VTYPE_REAL) THEN
                DO j = 1, PSD%wghts%num_links
                   IF (PSD%wghts%dstfrac(PSD%wghts%dstadd(j)) /= zero) THEN
                      pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j)) =     &
                           pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j))  &  
                           + ( REAL(PSD%wghts%weights(j),sp) *         &
                           pvari(i)%dat%vr(sshft+PSD%wghts%srcadd(j))) &
                           / REAL(PSD%wghts%dstfrac(PSD%wghts%dstadd(j)),sp)
                   ELSE
                      pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j)) = 0.
                   ENDIF
                END DO
             ELSE IF (pvari(i)%dat%type == VTYPE_DOUBLE) THEN
                DO j = 1, PSD%wghts%num_links
                   IF (PSD%wghts%dstfrac(PSD%wghts%dstadd(j)) /= zero) THEN
                      pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j)) =     &
                           pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j))  &  
                           + ( PSD%wghts%weights(j) *                  &
                           pvari(i)%dat%vd(sshft+PSD%wghts%srcadd(j))) &
                           / PSD%wghts%dstfrac(PSD%wghts%dstadd(j))
                  ELSE
                      pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j)) = 0.
                   ENDIF
                END DO
             ENDIF
          CASE(norm_opt_none)
             IF (pvari(i)%dat%type == VTYPE_REAL) THEN
                DO j = 1, PSD%wghts%num_links
                   IF (PSD%wghts%dstfrac(PSD%wghts%dstadd(j)) /= zero) THEN
                      pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j)) =       &
                           pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j))    &
                           + ( REAL(PSD%wghts%weights(j), sp) *          &
                           pvari(i)%dat%vr(sshft+PSD%wghts%srcadd(j)))   &
                           / REAL(PSD%wghts%dstfrac(PSD%wghts%dstadd(j)) &
                           *PSD%wghts%dstarea(PSD%wghts%dstadd(j)),sp) 
                   ELSE
                      pvaro(i)%dat%vr(dshft+PSD%wghts%dstadd(j)) = 0.
                   ENDIF
                END DO
             ELSE IF (pvari(i)%dat%type == VTYPE_DOUBLE) THEN             
                DO j = 1, PSD%wghts%num_links
                   IF (PSD%wghts%dstfrac(PSD%wghts%dstadd(j)) /= zero) THEN
                      pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j)) =     &
                           pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j))  &  
                           + (PSD%wghts%weights(j) *                   &
                           pvari(i)%dat%vd(sshft+PSD%wghts%srcadd(j))) &
                           / (PSD%wghts%dstfrac(PSD%wghts%dstadd(j))   &
                           *PSD%wghts%dstarea(PSD%wghts%dstadd(j))) 
   
                   ELSE
                      pvaro(i)%dat%vd(dshft+PSD%wghts%dstadd(j)) = 0._dp
                   ENDIF
             END DO
             ENDIF
          CASE DEFAULT
             ! normalize option not implemented
             status = 3030
             RETURN
          END SELECT
       END DO do_free
    END DO

    status = 0
    NULLIFY(PSD)

  END SUBROUTINE APPLY_SCRIP_WEIGHTS
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE SCRIP_CONTROL (status, SCRIP_ID, igrid, ogrid, RG_TYPE, lint &
                          , invar, var, grid, llrgz, lfirsto        )
    
    USE messy_main_grid,        ONLY: t_geohybgrid, COPY_GEOHYBGRID          &
                                    , INIT_GEOHYBGRID, EXPORT_GEOHYBGRID
    USE messy_main_grid_netcdf, ONLY: ERRMSG                                 &
                                    , RGMSG, RGMLVL, RGMLIC, RGMLVM, RGMLVMC &
                                    , RGMLW, RGMLWC, GRD_MAXSTRLEN           &
                                    , MAXFRAC2IDX_NCVAR, IDX2FRAC_NCVAR      &
                                    , t_ncvar, INIT_NCVAR, COPY_NCVAR        &
                                    , EXPORT_NCVAR, INIT_NCATT, CAT_NARRAY   &
                                    , INIT_NARRAY, t_ncatt, EXPORT_NCATT     &
                                    , IMPORT_NCATT, NF90_GLOBAL, NF90_CHAR   &
                                    , VTYPE_CHAR, VTYPE_UNDEF, QDEF_NCVAR

    USE messy_main_grid_trafo,           ONLY: RG_IDX, RG_IXF
    USE messy_main_grid_trafo_scrp_base, ONLY: SCRIPVERSION
    USE messy_main_grid_trafo,           ONLY: CHECK_NCVAR_ON_GEOHYBGRID  &
                                             , BALANCE_GEOHYBGRID_NCVAR   &
                                             , BALANCE_GEOHYBGRID_TIME    &
                                             , PACK_GEOHYBGRID_NCVAR      &
                                             , SWITCH_GEOHYBGRID          &
                                             , COMPLETE_GEOHYBGRID        &
                                             , CHECK_GEOHYBGRID

    IMPLICIT NONE

    INTEGER, INTENT(OUT)                      :: status
    INTEGER, INTENT(IN)                       :: SCRIP_ID
    TYPE(t_geohybgrid), INTENT(IN)            :: igrid
    TYPE(t_geohybgrid), INTENT(IN)            :: ogrid
    
    INTEGER, INTENT(IN)                       :: RG_TYPE(:)  ! regrid type
    LOGICAL, INTENT(IN)                       :: lint        ! input time ?
    TYPE (t_ncvar), DIMENSION(:), POINTER     :: invar ! list of input  variables
    TYPE (t_ncvar), DIMENSION(:), POINTER     ::   var ! list of output variables
    TYPE(t_geohybgrid), INTENT(INOUT), OPTIONAL :: grid
    LOGICAL, INTENT(IN)              , OPTIONAL :: llrgz
    LOGICAL, INTENT(IN)              , OPTIONAL :: lfirsto  ! first output step

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER            :: substr = 'SCRIP_CONTROL'
    ! number of interpolated variables
    INTEGER                                :: nvars  
    INTEGER                                :: i     ! loop indices
    TYPE (t_ncvar)                         :: qvari ! temporal variable
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: pvari => NULL() ! temporal variable (packed)
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: svaro => NULL() ! temporal variable (sorted)
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: pvaro => NULL() ! temporal variable (packed)
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: xovar => NULL() ! list of output variables
    ! GEOHYBRID-GRIDS
    TYPE (t_ncatt), SAVE           :: rggatt   ! RG global attribute
    TYPE(t_geohybgrid)             :: gi, gg   ! input and output geohybgrid
    TYPE(t_geohybgrid)             :: ggout
    TYPE(t_geohybgrid)             :: gis      ! sorted and index input-grid
    LOGICAL                        :: lok      ! result OK ?
    TYPE(t_geohybgrid)             :: ggs      ! sorted and index grid
    INTEGER, DIMENSION(:), POINTER :: dims => NULL() ! order of grid dimensions
    INTEGER                        :: axes(3)  ! grid axes of variable
    LOGICAL                        :: lrgz

    INTEGER                        :: ix
    CHARACTER(LEN=GRD_MAXSTRLEN)   :: timename
    INTEGER                        :: ixf(2)

    timename = ' '

    status = 3999
    IF (PRESENT(llrgz)) THEN
       lrgz = llrgz
    ELSE
       lrgz = .FALSE.
    ENDIF

    ! INITIALIZE OUTPUT PARAMETER
    IF (ASSOCIATED(var)) THEN
       DO i=1, SIZE(var)
          CALL INIT_NCVAR(var(i))
       END DO
       DEALLOCATE(var, STAT=status)
       CALL ERRMSG(substr,status,1)
    END IF
    NULLIFY(var)
    
    ! INIT
    NULLIFY(dims)
    NULLIFY(xovar)
    NULLIFY(svaro)
    NULLIFY(pvaro)
    NULLIFY(pvari)

    CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')
    CALL RGMSG(substr, RGMLVL, &
         'START SCRIP INTERPOLATION PROCEDURE (SCRIP VERSION '//&
       &TRIM(SCRIPVERSION)//')' )
    CALL RGMSG(substr, RGMLVL, &
         ' <Implementation: Astrid Kerkweg, Uni Mainz, November 2012>' )
    CALL RGMSG(substr, RGMLVL, &
         '-------------------------------------------------------')

    CALL COPY_GEOHYBGRID(gi, igrid)
    CALL COPY_GEOHYBGRID(gg, ogrid)
    CALL COPY_GEOHYBGRID(ggout, ogrid) ! define output grid
    
    ! CHECK INPUT GRID
    CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
    CALL CHECK_GEOHYBGRID(gi)
    CALL RGMSG(substr, RGMLIC, ' ... O.K.')

    ! CHECK OUTPUT GRID
    CALL RGMSG(substr, RGMLIC, ' ... SWITCHING ...')
    CALL SWITCH_GEOHYBGRID(gg, .TRUE., .TRUE., .FALSE.)
    CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
    CALL CHECK_GEOHYBGRID(gg)
    CALL RGMSG(substr, RGMLIC, ' ... O.K.')

    ! TIME BALANCING (TODO ???)
    CALL RGMSG(substr, RGMLVM, '>>> BALANCING TIME AXIS ...')
    CALL BALANCE_GEOHYBGRID_TIME(gi, gg, lint)
    CALL RGMSG(substr, RGMLVM, '<<< ... END BALANCING TIME AXIS !')

    CALL COPY_GEOHYBGRID(gis,gi)
    CALL COPY_GEOHYBGRID(ggs,gg)

    CALL RGMSG(substr, RGMLVM, '>>> COMPLETING INPUT GRID ...')
    CALL COMPLETE_GEOHYBGRID(gis)
    CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING INPUT GRID !')

    CALL RGMSG(substr, RGMLVM, '>>> COMPLETING OUTPUT GRID ...')
    CALL COMPLETE_GEOHYBGRID(ggs)
    CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING OUTPUT GRID !')

    ! BALANCE GRID
    !CALL RGMSG(substr, RGMLVM, '>>> BALANCING INPUT/OUTPUT GRID ...')
    !CALL BALANCE_GEOHYBGRID(gis, ggs)
    !CALL BALANCE_GEOHYBGRID(gix, ggx)
    !CALL RGMSG(substr, RGMLVM, '<<< END BALANCING INPUT/OUTPUT GRID !')
  
    nvars = SIZE(invar)

    CALL RGMSG(substr, RGMLVL, 'SCRIP INTERPOLATION MODE ...')

    ! ALLOCATE SPACE
    ALLOCATE(xovar(nvars))
    ALLOCATE(svaro(nvars))
    ALLOCATE(pvaro(nvars))
    ALLOCATE(pvari(nvars))

    DO i=1, nvars  ! LOOP OVER VARIABLE NAMES

       CALL RGMSG(substr, RGMLVL, &
            'VARIABLE '''//TRIM(invar(i)%name)//''' ...')
       IF ((RG_TYPE(i) == RG_IDX).OR.(RG_TYPE(i) == RG_IXF)) THEN
          ! get the number of index classes from the variable attributes,
          ! if given
          ixf(:) = 0
          DO ix=1,invar(i)%natts
             IF (TRIM(invar(i)%att(ix)%name) == 'mmig_ixf_range') THEN
                ixf(1:2) = INT(invar(i)%att(ix)%dat%vi(1:2))
             END IF
          END DO
          CALL RGMSG(substr, RGMLVMC, ' ... INDEX FRACTIONS ...')
          IF (ANY(ixf /= 0)) THEN
             CALL IDX2FRAC_NCVAR(invar(i), qvari, ixf=ixf)
          ELSE
             CALL IDX2FRAC_NCVAR(invar(i), qvari)
          END IF
          CALL RGMSG(substr, RGMLVMC, ' ... ->'''//TRIM(qvari%name)//'''')
       ELSE
          CALL COPY_NCVAR(qvari, invar(i))
       END IF

       ! GET ORDER INFORMATION
       CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
       CALL CHECK_NCVAR_ON_GEOHYBGRID(qvari, gis, dims, axes, lok)

       ! BALANCE OUTPUT VARIABLE
       CALL RGMSG(substr, RGMLIC, ' ... BALANCING ...')

       CALL BALANCE_GEOHYBGRID_NCVAR(qvari, axes, ggs, xovar(i), .FALSE.)

       ! SORT VARIABLE ACCORDING TO GRID
       !CALL RGMSG(substr, RGMLIC, ' ... SORTING ...')
       !CALL SORT_GEOHYBGRID_NCVAR('SCRIP 1', qvari, gix, axes, svari)

       CALL BALANCE_GEOHYBGRID_NCVAR(qvari, axes, ggs, svaro(i), .FALSE.)
 
       ! PACK VARIABLE
       CALL RGMSG(substr, RGMLIC, ' ... PACKING ...')
       CALL PACK_GEOHYBGRID_NCVAR(qvari, dims, axes, pvari(i), .FALSE.)
       CALL BALANCE_GEOHYBGRID_NCVAR(pvari(i), axes, ggs, pvaro(i), .FALSE.)

       CALL INIT_NARRAY(pvaro(i)%dat)
       ! CLEAN
       CALL INIT_NCVAR(qvari)

       ! DIMS and AXES NOT NEEDED ANYMPORE
       DEALLOCATE(dims, STAT=status)
       CALL ERRMSG(substr,status,8)
       NULLIFY(dims)
       axes(:) = 0

       CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
    END DO
    !***********************************************

    ! INTERPOLATION PROCEDURE
    CALL APPLY_SCRIP_WEIGHTS(status, pvari, pvaro, SCRIP_ID)
    IF (status /= 0) RETURN

    DO i = 1, nvars
       CALL RGMSG(substr, RGMLVL, &
            'VARIABLE '''//TRIM(xovar(i)%name)//''' ...')
       ! CHECK N-ARRAY
       CALL CHECK_NCVAR_ON_GEOHYBGRID(svaro(i), ggs, dims, axes, lok)
       ! UNPACK VARIABLE
       CALL RGMSG(substr, RGMLIC, ' ... UN-PACKING ...')
       CALL PACK_GEOHYBGRID_NCVAR(pvaro(i), dims, axes, svaro(i), .true.)

       !  UN-SORT
       !CALL RGMSG(substr, RGMLIC, ' ... UN-SORTING ...')
       !CALL SORT_GEOHYBGRID_NCVAR('SCRIP 1',svaro(i), ggx, axes, qvaro, .true.)

       !  UN-IDX
       CALL INIT_NCVAR(xovar(i))
       IF (RG_TYPE(i) == RG_IDX) THEN
          CALL RGMSG(substr, RGMLVMC, ' ... MAXIMUM INDEX FRACTION ...')
          CALL MAXFRAC2IDX_NCVAR(svaro(i), xovar(i))
          xovar(i)%name = TRIM(invar(i)%name)
          CALL RGMSG(substr, RGMLVMC, &
               ' ... ->'''//TRIM(xovar(i)%name)//'''')
       ELSE
          CALL COPY_NCVAR(xovar(i), svaro(i))
       ENDIF
 
       ! CLEAN
       CALL INIT_NCVAR(pvaro(i))
       CALL INIT_NCVAR(pvari(i))
       CALL INIT_NCVAR(svaro(i))
       DEALLOCATE(dims, STAT=status)
       CALL ERRMSG(substr,status,12)
       NULLIFY(dims)
       axes(:) = 0

       CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
 
    END DO

    !***********************************************

    IF (nvars > 0) THEN
       DEALLOCATE(svaro, pvaro, pvari, STAT=status)
       CALL ERRMSG(substr,status,15)
       NULLIFY(svaro)
       NULLIFY(pvaro)
       NULLIFY(pvari)
    END IF
       
    ! OUTPUT GRID AND VARIABLES
    IF (TRIM(gg%file) /= '' .AND. .NOT. lrgz) THEN
       IF (QDEF_NCVAR(gg%timem)) CALL COPY_NCVAR(ggout%timem, gg%timem)
       IF (QDEF_NCVAR(gg%timei)) CALL COPY_NCVAR(ggout%timei, gg%timei)
       ggout%t = gg%t 

       CALL EXPORT_GEOHYBGRID(ggout)
       ! EXPORT GLOBL ATTRIBUTE, ONLY AT FIRST TIME STEP
       IF (PRESENT(lfirsto)) THEN
          IF (lfirsto) THEN
             
             ! GET OLD ATTRIBUTE
             CALL IMPORT_NCATT(rggatt, varid=NF90_GLOBAL   &
                  ,attname = 'RG_HISTORY'     &
                  ,file = TRIM(gg%file), lnostop=.true.)
             ! CHECK IF EMPTY OR CHAR-ATT
             IF ((rggatt%dat%type == VTYPE_UNDEF).OR.  &
                  (rggatt%dat%type == VTYPE_CHAR)) THEN
                ! APPEND NEW DATA
                CALL CAT_NARRAY(rggatt%dat, gg%att%dat)
                rggatt%xtype = NF90_CHAR
                rggatt%len   = rggatt%dat%dim(1)
             ELSE  ! EXISTS WITH NON-CHAR TYPE
                CALL RGMSG(substr, RGMLW, &
                     'ATTRIBUTE '''//TRIM(rggatt%name)//'''')
                CALL RGMSG(substr, RGMLWC, &
                     'EXISTS ALREADY, BUT IS NOT OF TYPE CHAR !')
             END IF
             CALL EXPORT_NCATT(rggatt, file=TRIM(gg%file)  &
                  ,clobber=.true.)
             ! EXPORT VARIABLES
          END IF
       END IF
       DO i=1, nvars
          DO ix = 1, xovar(1)%ndims
             IF (xovar(i)%dim(ix)%fuid) THEN
                IF (TRIM(timename) /= '') THEN
                   xovar(i)%dim(ix)%name = TRIM(timename)
                ELSE
                   timename = TRIM(xovar(i)%dim(ix)%name)
                END IF
                EXIT
             END IF
          END DO
          xovar(i)%ustep = ggout%t
          CALL EXPORT_NCVAR(xovar(i), file=TRIM(gg%file))
       END DO
       CALL INIT_NCATT(rggatt)
    END IF
    CALL INIT_GEOHYBGRID(ggout)
         
    ! RETURN VALUES TO SUBROUTINE CALL
    CALL RGMSG(substr, RGMLVM, '... RETURNING INTERPOLATED DATA')
    
    ALLOCATE(var(nvars), STAT=status)
    CALL ERRMSG(substr,status,16)
    DO i=1, nvars
       ! assign time
       xovar(i)%ustep = gg%t
       CALL COPY_NCVAR(var(i), xovar(i))
    END DO
    IF (PRESENT(grid)) THEN
       CALL RGMSG(substr, RGMLVM, '... RETURNING OUTPUT GRID')
       CALL COPY_GEOHYBGRID(grid, ggs)
       ! ADD vertical axis of input grid
       CALL COPY_NCVAR(grid%hyam, gi%hyam)
       CALL COPY_NCVAR(grid%hyai, gi%hyai)
       CALL COPY_NCVAR(grid%hybm, gi%hybm)
       CALL COPY_NCVAR(grid%hybi, gi%hybi)
    END IF

    ! CLEAN
    DO i=1,nvars
       CALL INIT_NCVAR(xovar(i))
    END DO
    DEALLOCATE(xovar, STAT=status)
    CALL ERRMSG(substr,status,17)
    NULLIFY(xovar)
     !
    CALL RGMSG(substr, RGMLVL, '... DONE (INTERPOLATION MODE) !')
    CALL RGMSG(substr, RGMLVL, '-------------------------------')

    ! CLEAN
    ! TODO: DEALLOCATION OF VARIABLES IN ALL LEVELS of IMPORT /INTERPOLATION
    !
    CALL INIT_GEOHYBGRID(gis)
    CALL INIT_GEOHYBGRID(gg)
    CALL INIT_GEOHYBGRID(ggs)
    CALL INIT_GEOHYBGRID(gi)

    status = 0

  END SUBROUTINE SCRIP_CONTROL
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE INTERPOL_GEOHYBGRID_PS(status, gi, go, PSDID)

    USE messy_main_grid,           ONLY: t_geohybgrid, INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,    ONLY: t_ncvar, INIT_NCVAR, COPY_NCVAR &
                                       , ERRMSG &
                                       , RGMLE, RGMSG, QDEF_NCVAR
    USE messy_main_grid_trafo,     ONLY: SORT_GEOHYBGRID, COMPLETE_GEOHYBGRID&
                                       , BALANCE_GEOHYBGRID        &
                                       , CHECK_NCVAR_ON_GEOHYBGRID &
                                       , PACK_GEOHYBGRID_NCVAR     &
                                       , BALANCE_GEOHYBGRID_NCVAR

    IMPLICIT NONE

    ! I/O
    INTEGER,                  INTENT(OUT)   :: status
    TYPE (t_geohybgrid),      INTENT(IN)    :: gi
    TYPE (t_geohybgrid),      INTENT(INOUT) :: go
    INTEGER,                  INTENT(IN)    :: PSDID

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'INTERPOL_GEOHYBGRID_PS'
    TYPE (t_geohybgrid)                   :: gih, goh, gihs, gohs
    TYPE (t_geohybgrid)                   :: gix, gox
    INTEGER,        DIMENSION(:), POINTER :: dims => NULL()! order of grid dimensions
    INTEGER                               :: axes(3)   ! grid axes of variable
    TYPE (t_ncvar)                        :: psi
    TYPE (t_ncvar)                        :: psis, psos
    TYPE (t_ncvar), DIMENSION(:), POINTER :: psip => NULL(), psop => NULL()
    LOGICAL                               :: ok

    ! INIT
    NULLIFY(dims)

    ! CREATE 2-D GRIDS
    CALL INIT_GEOHYBGRID(gih)
    CALL INIT_GEOHYBGRID(goh)
    !
    gih%file = TRIM(gi%file)
    gih%t    = gi%t
    !
    goh%file = TRIM(go%file)
    goh%t    = go%t

    !
    IF (QDEF_NCVAR(gi%lonm)) THEN
       CALL COPY_NCVAR(gih%lonm, gi%lonm)
       CALL COPY_NCVAR(gih%loni, gi%loni)
       !
       CALL COPY_NCVAR(gih%latm, gi%latm)
       CALL COPY_NCVAR(gih%lati, gi%lati)
       !
    ELSEIF (QDEF_NCVAR(gi%clonm)) THEN
       CALL COPY_NCVAR(gih%lonm, gi%clonm)
       CALL COPY_NCVAR(gih%loni, gi%cloni)
       !
       CALL COPY_NCVAR(gih%latm, gi%clatm)
       CALL COPY_NCVAR(gih%lati, gi%clati)
       !
    ENDIF

    IF (QDEF_NCVAR(gi%ulonm)) THEN
       CALL COPY_NCVAR(gih%ulonm, gi%ulonm)
       CALL COPY_NCVAR(gih%uloni, gi%uloni)
       !
       CALL COPY_NCVAR(gih%ulatm, gi%ulatm)
       CALL COPY_NCVAR(gih%ulati, gi%ulati)
       !
    ENDIF

    IF (QDEF_NCVAR(go%lonm)) THEN
       CALL COPY_NCVAR(goh%latm, go%latm)
       CALL COPY_NCVAR(goh%lati, go%lati)
       !
       CALL COPY_NCVAR(goh%lonm, go%lonm)
       CALL COPY_NCVAR(goh%loni, go%loni)
    ELSEIF (QDEF_NCVAR(go%clonm)) THEN
       CALL COPY_NCVAR(goh%latm, go%clatm)
       CALL COPY_NCVAR(goh%lati, go%clati)
       !
       CALL COPY_NCVAR(goh%lonm, go%clonm)
       CALL COPY_NCVAR(goh%loni, go%cloni)
    ENDIF

    IF (QDEF_NCVAR(go%ulonm)) THEN
       CALL COPY_NCVAR(goh%latm, go%ulatm)
       CALL COPY_NCVAR(goh%lati, go%ulati)
       !
       CALL COPY_NCVAR(goh%lonm, go%ulonm)
       CALL COPY_NCVAR(goh%loni, go%uloni)
    ENDIF

    CALL COPY_NCVAR(gih%timem, gi%timem)
    CALL COPY_NCVAR(gih%timei, gi%timei)
    !
    CALL COPY_NCVAR(goh%timem, go%timem)
    CALL COPY_NCVAR(goh%timei, go%timei)

    CALL SORT_GEOHYBGRID(gih, gihs, gix)
    CALL SORT_GEOHYBGRID(goh, gohs, gox)

    CALL COMPLETE_GEOHYBGRID(gihs, gix)
    CALL COMPLETE_GEOHYBGRID(gohs, gox)

    CALL BALANCE_GEOHYBGRID(gihs, gohs)
    CALL BALANCE_GEOHYBGRID(gix,  gox)

    CALL CHECK_NCVAR_ON_GEOHYBGRID(gi%ps, gihs, dims, axes, ok)
    IF (.NOT.ok) THEN
       CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM (I) !')
    END IF

    CALL COPY_NCVAR(psis, gi%ps)
    CALL BALANCE_GEOHYBGRID_NCVAR(psis, axes, gohs, psos)

    ALLOCATE(psip(1),psop(1))

    CALL PACK_GEOHYBGRID_NCVAR(psis, dims, axes, psip(1))
    CALL BALANCE_GEOHYBGRID_NCVAR(psip(1), axes, gohs, psop(1))

    DEALLOCATE(dims, STAT=status)
    CALL ERRMSG(substr,status,1)
    NULLIFY(dims)

    CALL apply_scrip_weights(status, psip, psop, PSDID)
    IF (status /=0) RETURN

    CALL CHECK_NCVAR_ON_GEOHYBGRID(psos, gohs, dims, axes, ok)

    IF (.NOT.ok) THEN
       CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM (O) !')
    END IF

    CALL PACK_GEOHYBGRID_NCVAR(psop(1), dims, axes, psos, .true.)

    CALL COPY_NCVAR(go%ps, psos)

    ! CLEAN UP
    DEALLOCATE(dims, STAT=status)
    CALL ERRMSG(substr,status,2)
    NULLIFY(dims)

    CALL INIT_NCVAR(psi)
    CALL INIT_NCVAR(psis)
    CALL INIT_NCVAR(psos)
    CALL INIT_NCVAR(psip(1))
    CALL INIT_NCVAR(psop(1))
    DEALLOCATE(psop,psip)
    CALL INIT_GEOHYBGRID(gih)
    CALL INIT_GEOHYBGRID(goh)
    CALL INIT_GEOHYBGRID(gihs)
    CALL INIT_GEOHYBGRID(gohs)
    CALL INIT_GEOHYBGRID(gix)
    CALL INIT_GEOHYBGRID(gox)

  END SUBROUTINE INTERPOL_GEOHYBGRID_PS
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE BALANCE_CURVILINEAR_PS(status, gi, go, PSDID)

    USE messy_main_grid,        ONLY: t_geohybgrid
    USE messy_main_grid_netcdf, ONLY: RGMSG, RGMLE, RGMLEC, RGMLI &
                                    , QDEF_NCVAR, COPY_NCVAR

    IMPLICIT NONE

    ! I/O
    INTEGER,                    INTENT(OUT)   :: status
    TYPE (t_geohybgrid),        INTENT(INOUT) :: gi, go
    INTEGER, INTENT(IN)                       :: PSDID 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_CURVILINEAR_PS'
    LOGICAL :: err
    LOGICAL :: linv   ! PRE-REGRID AVAILABLE go%ps

    status = 3999

    err = (.NOT.QDEF_NCVAR(gi%ps)).AND.(.NOT.QDEF_NCVAR(go%ps))
    err = err .AND. ((QDEF_NCVAR(gi%hybi).OR.QDEF_NCVAR(gi%hybm)) .OR. &
         (QDEF_NCVAR(go%hybi).OR.QDEF_NCVAR(go%hybm)))

    IF (err) THEN
       CALL RGMSG(substr, RGMLE, &
            'HYBRID-B-COEFFICIENTS NEED I_PS AND/OR G_PS IN NAMELIST !')
    END IF

    err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(go%p0))
    err = err .AND. ((QDEF_NCVAR(gi%hyai).OR.QDEF_NCVAR(gi%hyam)) .OR. &
         (QDEF_NCVAR(go%hyai).OR.QDEF_NCVAR(go%hyam)))

    IF (err) THEN
       CALL RGMSG(substr, RGMLE, &
            'HYBRID-A-COEFFICIENTS NEED I_P0 AND/OR G_P0 IN NAMELIST !')
    END IF

    ! RETURN, IF 2D
    err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(gi%ps)).AND.     &
         (.NOT.QDEF_NCVAR(gi%hyai)).AND.(.NOT.QDEF_NCVAR(gi%hybi)).AND. &
         (.NOT.QDEF_NCVAR(gi%hyam)).AND.(.NOT.QDEF_NCVAR(gi%hybm))
    IF (err) THEN
       CALL RGMSG(substr, RGMLI, &
            'INPUT GRID IS 2-D! NO SURFACE PRESSURE REGRIDDING REQUIRED!')
       status = 0
       RETURN
    END IF

    ! BALANCE PO IN ANY CASE
    IF (QDEF_NCVAR(gi%p0).AND.(.NOT.QDEF_NCVAR(go%p0))) THEN
       CALL COPY_NCVAR(go%p0, gi%p0)
    ELSE
       IF (QDEF_NCVAR(go%p0).AND.(.NOT.QDEF_NCVAR(gi%p0))) THEN
          CALL COPY_NCVAR(gi%p0, go%p0)
       END IF
    END IF

    ! NOW PS ...

    ! CASE A: gi%ps AND go%ps BOTH AVAILABLE
    ! ADJUST TIME AXIS
    ! NOTE: THE TIME BALANCING (INCLUDING THAT FOR THE SURFACE PRESSURE
    !       TIME DIMENSION) IS PERFORMED IN
    !       SUBROUTINE BALANCE_GEOHYBGRID_TIME !

    ! CASE B: gi%ps AND go%ps BOTH UNAVAILABLE
    ! NOTHING TO DO !

    IF (.NOT.QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
       RETURN
    END IF

    ! SURFACE PRESSURE NEEDS TO BE PRE-REGRIDDED, IF NOT AVAILABLE
    ! CASE C: go%ps AVAILABLE, BUT ON WRONG HORIZONTAL GRID
    !
    linv = QDEF_NCVAR(go%ps).AND.( &
         ((.NOT.QDEF_NCVAR(go%clonm)).AND.(.NOT.QDEF_NCVAR(go%cloni)) .AND. &
         (.NOT.QDEF_NCVAR(go%ulonm)).AND.(.NOT.QDEF_NCVAR(go%uloni)) .AND. &
         (.NOT.QDEF_NCVAR(go%lonm )).AND.(.NOT.QDEF_NCVAR(go%loni ))) .OR. &
         ((.NOT.QDEF_NCVAR(go%clatm)).AND.(.NOT.QDEF_NCVAR(go%clati)) .AND. &
         (.NOT.QDEF_NCVAR(go%ulatm)).AND.(.NOT.QDEF_NCVAR(go%ulati)) .AND. &
         (.NOT.QDEF_NCVAR(go%latm )).AND.(.NOT.QDEF_NCVAR(go%lati )))      &
         )
    !
    IF (QDEF_NCVAR(go%ps).AND.linv) THEN
       CALL RGMSG(substr, RGMLE,  &
            'REGRIDDING 3-D DISTRIBUTIONS', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'ONTO A DESTINATION SURFACE PRESSURE COORDINATE', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'USING (AN) INVARIANT HORIZONTAL DIMENSION(S)', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'IS NOT POSSIBLE DUE TO A LACK OF INFORMATION', .false.)
       CALL RGMSG(substr, RGMLEC, &
            '(HORIZONTAL DESTINATION GRID)', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'FOR PRE-REGRIDDING THE DESTINATION SURFACE PRESSURE !', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'PLEASE PERFORM 2-D PRE-REGRIDDING OF SURFACE PRESSURE', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'IN SEPARATE STEP! (SCRP)')
    END IF

    ! CASE D: gi%ps XOR go%ps NOT AVAILABLE
    IF (QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
       CALL INTERPOL_GEOHYBGRID_PS(status, gi, go, PSDID)
       IF (status /= 0) RETURN
    ELSE  IF (QDEF_NCVAR(go%ps).AND.(.NOT.QDEF_NCVAR(gi%ps))) THEN
       CALL INTERPOL_GEOHYBGRID_PS(status, go, gi, PSDID)
       IF (status /= 0) RETURN
    END IF

    status = 0

  END SUBROUTINE BALANCE_CURVILINEAR_PS
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE CONSTRUCT_INPUT_SURF_PRESSURE(status, gi, gips, PSDID &
                                                 , go, RGT, lint)

    USE messy_main_grid,         ONLY: t_geohybgrid
    USE messy_main_grid_netcdf,  ONLY: QDEF_NCVAR, INIT_NARRAY, ADD_NCATT &
                                     , ERRMSG, NULL_VARID, NF90_FLOAT     &
                                     , VTYPE_REAL, POSITION, INIT_NCVAR   &
                                     , t_ncvar, COPY_NCVAR, RGMSG, RGMLVM 

    IMPLICIT NONE

    INTRINSIC :: TRIM

    TYPE(t_geohybgrid), INTENT(IN)    :: gi
    TYPE(t_geohybgrid), INTENT(INOUT) :: gips
    INTEGER           , INTENT(OUT)   :: status
    INTEGER,            INTENT(IN)    :: PSDID 
    TYPE(t_geohybgrid), INTENT(IN)    :: go
    INTEGER, DIMENSION(:),  POINTER   :: RGT  ! regridding type
    LOGICAL,            INTENT(IN)    :: lint 
    ! LOCAL
    TYPE (t_ncvar), DIMENSION(:), POINTER :: rvar => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar => NULL() ! list of variables
    CHARACTER(LEN=*), PARAMETER :: substr = 'CONSTRUCT_INPUT_SURF_PRESSURE'

    ! LOCAL
    INTEGER :: i, j

    status = 3999
    CALL RGMSG(substr, RGMLVM, '>>> IN CONSTRUCT_SURF_PRESSURE  ...')
    
    IF ( .NOT. QDEF_NCVAR(gi%ps)) THEN !RETURN
       CALL RGMSG(substr, RGMLVM, '>>> IN CONSTRUCT_SURF_PRESSURE NOT DEF ...')
       
       CALL INIT_NCVAR(gips%ps)
       gips%ps%name  = TRIM('ps')
       gips%ps%id    = NULL_VARID
       gips%ps%xtype = NF90_FLOAT
       ! ... dimensions
       gips%ps%ndims = 2
       ALLOCATE(gips%ps%dim(gips%ps%ndims), STAT=status)
       CALL ERRMSG('INTERFACE_GEOHYBGRID',status,5)
       IF (QDEF_NCVAR(gips%lonm)) THEN
          gips%ps%dim(1) = gips%lonm%dim(1)
          gips%ps%dim(2) = gips%latm%dim(1)
       ELSEIF (QDEF_NCVAR(gips%clonm)) THEN
          gips%ps%dim(1) = gips%clonm%dim(1)
          gips%ps%dim(2) = gips%clonm%dim(2)
       ELSE IF (QDEF_NCVAR(gips%ulonm)) THEN
          gips%ps%dim(1) = gips%ulonm%dim(1)
          IF (SIZE(gips%ulonm%dim)>1) THEN
             gips%ps%dim(2) = gips%ulonm%dim(2)
          ELSE
             gips%ps%dim(2)%name = 'singleton'
             gips%ps%dim(2)%len  = 1 
          END IF
       END IF
    
       ! ... data
       CALL INIT_NARRAY(gips%ps%dat, gips%ps%ndims      &
            ,(/gips%ps%dim(1)%len, gips%ps%dim(2)%len/) &
            ,VTYPE_REAL)
       DO i=1, gips%ps%dim(1)%len
          DO j=1, gips%ps%dim(2)%len
             gips%ps%dat%vr(POSITION((/gips%ps%dim(1)%len, gips%ps%dim(2)%len/)&
                  ,(/i,j/)))  = 101325.0_sp
          END DO
       END DO
    
       ! ... attributes
       CALL ADD_NCATT(gips%ps, 'long_name', vs='constants surface pressure')
       CALL ADD_NCATT(gips%ps, 'units', vs='Pa')
    ELSE

       CALL RGMSG(substr, RGMLVM, '>>> REGRID SURFACE PRESSURE  ...')
       ALLOCATE(rvar(1))
       CALL COPY_NCVAR(rvar(1),gi%ps)
       CALL SCRIP_CONTROL(status, PSDID, gi, go, RGT, lint &
                  , rvar, ovar, llrgz=.TRUE., lfirsto=.FALSE.) 
       CALL COPY_NCVAR(gips%ps, ovar(1))
       CALL INIT_NCVAR(ovar(1))
       CALL INIT_NCVAR(rvar(1))
       DEALLOCATE(ovar) ; NULLIFY(ovar)
       DEALLOCATE(rvar) ; NULLIFY(rvar)
       CALL RGMSG(substr, RGMLVM, '... REGRID SURFACE PRESSURE  <<<')

    END IF

    CALL COPY_NCVAR(gips%p0, gi%p0)
    CALL RGMSG(substr, RGMLVM, ' END CONSTRUCT_SURF_PRESSURE  >>>')
 
    status = 0

  END SUBROUTINE CONSTRUCT_INPUT_SURF_PRESSURE
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE INTERPOL_3D_PRESSURE(status, gi, gipres, PSDID, go, RGT, lint)

    USE messy_main_grid,         ONLY: t_geohybgrid, COPY_GEOHYBGRID      &
                                     , INIT_GEOHYBGRID
    USE messy_main_grid_netcdf,  ONLY: QDEF_NCVAR, INIT_NCVAR   &
                                     , t_ncvar, COPY_NCVAR, RGMSG, RGMLVM &
                                     , double_narray

    IMPLICIT NONE

    TYPE(t_geohybgrid), INTENT(IN)    :: gi
    TYPE(t_geohybgrid), INTENT(INOUT) :: gipres
    INTEGER           , INTENT(OUT)   :: status
    INTEGER,            INTENT(IN)    :: PSDID 
    TYPE(t_geohybgrid), INTENT(IN)    :: go
    INTEGER, DIMENSION(:),  POINTER   :: RGT  ! regridding type
    LOGICAL,            INTENT(IN)    :: lint 
    ! LOCAL
    TYPE (t_ncvar), DIMENSION(:), POINTER :: rvar => NULL() ! list of variables
    TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar => NULL() ! list of variables
    CHARACTER(LEN=*), PARAMETER :: substr = 'INTERPOL_3D_PRESSURE'

    ! LOCAL
    TYPE(t_geohybgrid) :: gint

    status = 3999

    CALL RGMSG(substr, RGMLVM, '>>> HOR. REMAP 3D PRESSURE  ...')
    ALLOCATE(rvar(1))
    IF (QDEF_NCVAR(gi%pressi)) THEN
       CALL COPY_NCVAR(rvar(1),gi%pressi)
    ELSE
       CALL COPY_NCVAR(rvar(1),gi%pressm)
    END IF
    CALL DOUBLE_NARRAY(rvar(1)%dat)

    CALL SCRIP_CONTROL(status, PSDID, gi, go, RGT, lint &
         , rvar, ovar, gint, llrgz=.FALSE., lfirsto=.FALSE.)
    CALL DOUBLE_NARRAY(ovar(1)%dat)

    IF (QDEF_NCVAR(gi%pressi)) THEN
       CALL COPY_NCVAR(gint%pressi, ovar(1))
    ELSE
       CALL COPY_NCVAR(gint%pressm, ovar(1))
    END IF
    gint%ranges = gipres%ranges
    ! done in COPY_GEOHYBGRID CALL INIT_GEOHYBGRID(gipres)
    CALL COPY_GEOHYBGRID(gipres, gint)

    CALL INIT_NCVAR(ovar(1))
    CALL INIT_NCVAR(rvar(1))
    CALL INIT_GEOHYBGRID(gint)
    DEALLOCATE(ovar) ; NULLIFY(ovar)
    DEALLOCATE(rvar) ; NULLIFY(rvar)
    CALL RGMSG(substr, RGMLVM, '... HOR. REMAP 3D PRESSURE  <<<')

    status = 0

  END SUBROUTINE INTERPOL_3D_PRESSURE
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE CLEAN_SCRIPGRID_LIST(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,  INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_scrip_grid_list),      POINTER :: psgi => NULL(), psge => NULL()
    TYPE(t_scrip_grid),           POINTER :: grid => NULL()

    status = 0

    psgi => SCRIPGRIDLIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(psgi)) EXIT

       grid => psgi%this

       psge => psgi
       psgi => psgi%next

       ! -------------------------------------------
       ! MEMORY: DATA
       CALL INIT_SCRIPGRID(grid)
       DEALLOCATE(psge)
       NULLIFY(psge)
       !
       ! COUNT
       NSGRID = NSGRID - 1
       ! -------------------------------------------

    END DO channel_loop

    NULLIFY(SCRIPGRIDLIST)

    IF (NSGRID /= 0) THEN
       status = 3045 ! INTERNAL SCRIP GRID COUNT ERROR
       RETURN
    END IF

    status = 0

  END SUBROUTINE CLEAN_SCRIPGRID_LIST
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE CLEAN_SCRIPDATA_LIST(status)

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! I/O
    INTEGER,  INTENT(OUT) :: status

    ! LOCAL
    TYPE(t_scrip_data_list), POINTER :: psdi => NULL(), psde => NULL()
    TYPE(t_scrip_data),      POINTER :: data => NULL()

    status = 0

    psdi => SCRIPDATALIST
    channel_loop: DO
       IF (.NOT. ASSOCIATED(psdi)) EXIT

       data => psdi%this

       psde => psdi
       psdi => psdi%next

       ! -------------------------------------------
       ! MEMORY: DATA
       CALL INIT_SCRIPDATA(data)
       DEALLOCATE(psde)
       NULLIFY(psde)
       !
       ! COUNT
       NSDATA = NSDATA - 1
       ! -------------------------------------------

    END DO channel_loop

    NULLIFY(SCRIPDATALIST)

    IF (NSDATA /= 0) THEN
       status = 3040 ! INTERNAL SCRIP DATA COUNT ERROR
       RETURN
    END IF

    status = 0

  END SUBROUTINE CLEAN_SCRIPDATA_LIST
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax, gi, lp &
       , PSDID, go, RGT, lint)


    USE messy_main_grid_trafo,   ONLY: RNGADJ_ARRAY, RNGADJ_NARRAY       &
                                     , H2PSIG_3D, IMMI_NARRAY, IMMI_NCVAR
    USE messy_main_grid_tools,   ONLY: RGTOOL_CONVERT, DEALINE_ARRAY
    USE messy_main_grid_netcdf,  ONLY: RGMSG, DOUBLE_NARRAY, QDEF_NCVAR &
                                     , RGMLI, RGMLE,NULL_XTYPE
    USE messy_main_grid,         ONLY: COPY_GEOHYBGRID, t_geohybgrid

    IMPLICIT NONE
    INTRINSIC :: PRESENT

    INTEGER           , INTENT(OUT)             :: status
    INTEGER           , INTENT(IN)              :: xsize,ysize
    REAL(dp), DIMENSION(:,:,:), POINTER         :: vax
    TYPE(t_geohybgrid), INTENT(IN)              :: gi
    LOGICAL           , INTENT(IN)              :: lp ! pressure axis
    INTEGER,            INTENT(IN),    OPTIONAL :: PSDID 
    TYPE(t_geohybgrid), INTENT(IN),    OPTIONAL :: go
    INTEGER, DIMENSION(:),  POINTER,   OPTIONAL :: RGT  ! regridding type
    LOGICAL,            INTENT(IN),    OPTIONAL :: lint 

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = ' CONSTRUCT_VERTICAL_AXIS'
    REAL(dp), DIMENSION(:,:,:,:), POINTER :: axdat => NULL()
    REAL(dp), DIMENSION(:,:),   POINTER   :: zps   => NULL()
    INTEGER                               :: nx, ny, nz, zdim
    LOGICAL                               :: linterpol
    TYPE(t_geohybgrid)                    :: gi2

    ! INITIALIZE
    status = 5555

    CALL RGMSG(substr, RGMLI, 'CONSTRUCT_VERTICAL_AXIS ===>')
    ! use presence of all optional parameters as indicator, if interpolation
    ! of pressure field is required
    linterpol =  PRESENT(PSDID) .AND. PRESENT(go) &
         .AND. PRESENT(RGT) .AND. PRESENT(lint)

    CALL COPY_GEOHYBGRID(gi2,gi)

    IF (QDEF_NCVAR(gi2%hyam) .AND. .NOT. QDEF_NCVAR(gi2%hyai))  THEN
       ! CALCULATE HYAI / HYBI from hyam /hybm
       CALL RGMSG(substr, RGMLI,' CALCULATE  hyai FROM  hyam ...')
       CALL IMMI_NARRAY(gi2%hyai%dat, gi2%hyam%dat)
       CALL IMMI_NCVAR(gi2%hyai, gi2%hyam)
    END IF

    IF  (QDEF_NCVAR(gi2%hybm) .AND. .NOT. QDEF_NCVAR(gi2%hybi)) THEN
       ! CALCULATE HYAI / HYBI from hyam /hybm
       CALL RGMSG(substr, RGMLI,' CALCULATE  hybi FROM  hybm ...')
       CALL IMMI_NARRAY(gi2%hybi%dat, gi2%hybm%dat)
       CALL IMMI_NCVAR(gi2%hybi, gi2%hybm)
    END IF

    ! a1) 3d pressure field on interfaces exists
    ifpressi: IF (QDEF_NCVAR(gi2%pressi)) THEN
       CALL RGMSG(substr, RGMLI, &
            'Grid provides interface 3D pressure field')

       IF (linterpol) THEN
          ! Horizontally interpolate pressure grid
          CALL INTERPOL_3D_PRESSURE(status, gi, gi2, PSDID, go, RGT, lint)
          IF (status /= 0) RETURN
       END IF

       CALL RGTOOL_CONVERT(gi2%pressi, axdat, gi2, order='xyzn') 
       zdim = SIZE(axdat,3) 
       ALLOCATE(VAX(xsize,ysize,zdim))
       vax(:,:,:) = axdat(:,:,:,1)
       DEALLOCATE(axdat)
       NULLIFY(axdat)

       DO nx = 1,xsize
          DO ny = 1, ysize
             ! ADJUST RANGES
             CALL RNGADJ_ARRAY(vax(nx,ny,:), gi2%ranges(3,:))
          END DO
       END DO

       IF (.NOT. lp) THEN
          DO nx = 1,xsize
             DO ny = 1, ysize
                vax(nx,ny,:) = vax(nx,ny,:) / MAXVAL(vax(nx,ny,:))
             END DO
          END DO
       END IF
       ! a2) 3d pressure field on mids exists
    ELSE IF (QDEF_NCVAR(gi2%pressm)) THEN
       CALL RGMSG(substr, RGMLI, &
            'Grid provides mids 3D pressure field')
       ! Horizontally interpolate pressure grid
       IF (linterpol) THEN
          CALL INTERPOL_3D_PRESSURE(status, gi,gi2, PSDID, go, RGT, lint)
          IF (status /= 0) RETURN
       END IF
       CALL DOUBLE_NARRAY(gi2%pressm%dat)
       CALL RGTOOL_CONVERT(gi2%pressm, axdat,gi2,order='xyzn') 
       ! FOR THE vertical INTERPOLATION INTERFACES ARE REQUIRED
       zdim = gi2%pressm%dim(3)%len+1
       ALLOCATE(VAX(xsize,ysize,zdim))
       DO nz = 2, zdim-1
          vax(:,:,nz) = &
               (axdat(:,:,nz-1,1)+axdat(:,:,nz,1))/2._dp
       END DO
       vax(:,:,1) = 2._dp *axdat(:,:,1,1) - vax(:,:,2)
       vax(:,:,zdim) =  &
            2._dp * axdat(:,:,zdim-1,1)  - vax(:,:,zdim-1)
       DEALLOCATE(axdat)
       NULLIFY(axdat)          
       ! ADJUST RANGES 
       DO nx = 1,xsize
          DO ny = 1, ysize
             CALL RNGADJ_ARRAY(vax(nx,ny,:), gi%ranges(3,:))
          END DO
       END DO
       IF (.NOT. lp) THEN
          DO nx = 1,xsize
             DO ny = 1, ysize
                vax(nx,ny,:) = vax(nx,ny,:) / MAXVAL(vax(nx,ny,:))
             END DO
          END DO
       END IF
       ! a3) construct pressure axis via p0,ps, hyai, hybi
    ELSE IF (QDEF_NCVAR(gi2%hyai).OR. QDEF_NCVAR(gi2%hybi)) THEN
       CALL RGMSG(substr, RGMLI, &
            'Grid static pressure field defined by hybrid coefficients')
       intype: IF (gi%ps%xtype /= NULL_XTYPE) THEN
          CALL DOUBLE_NARRAY(gi2%ps%dat)
          psfield: IF (SIZE(gi2%ps%dat%vd) > 1 ) THEN
             IF (linterpol) THEN
                ! PS still needs to be horizontally interpolated
                CALL CONSTRUCT_INPUT_SURF_PRESSURE(status, gi &
                     , gi2, PSDID, go, RGT, lint)
                IF (status /= 0) RETURN
                CALL DOUBLE_NARRAY(gi2%ps%dat)
             END IF
             CALL DEALINE_ARRAY(status, xsize,ysize &
                  , gi2%ps%dat%vd, zps)
             IF (status /= 0) THEN
                write (0,*) 'DIMENSION MISMATCH1', xsize,ysize &
                     , gi2%ps%dim(1)%len, gi2%ps%dim(2)%len
                CALL RGMSG(substr, RGMLE, 'DIMENSION MISMATCH1')
             END IF
          ELSE ! psfield
             ALLOCATE(zps(xsize,ysize))
             zps = gi2%ps%dat%vd(1)
          END IF psfield
       ELSE !intype
          ALLOCATE(zps(xsize,ysize))
          zps = -999.99_dp
       END IF intype
       CALL RGMSG(substr, RGMLI, 'CALCULATE BOUNDS2!')

       ! APPLY RANGE ADJUSTMENT
       IF (QDEF_NCVAR(gi2%hyai)) CALL RNGADJ_NARRAY( &
            gi2%hyai%dat,gi2%ranges(3,:),.FALSE.)
       IF (QDEF_NCVAR(gi2%hybi)) CALL RNGADJ_NARRAY( &
            gi2%hybi%dat,gi2%ranges(4,:))
       CALL H2PSIG_3D(vax,xsize, ysize,gi2%hyai%dat,gi2%hybi%dat &
            ,zps,gi2%p0%dat,lp)
       CALL RGMSG(substr, RGMLI, 'SIZE VAX: ',SIZE(VAX,3), ' ')
       DEALLOCATE(zps)
       NULLIFY(zps)
    ELSE ! ifpress
       CALL RGMSG(substr, RGMLE, &
            'Vertical regridding required, but hyai, hybi '//&
            &'for in-grid  missing!')
    END IF ifpressi
    CALL RGMSG(substr, RGMLI, '<==CONSTRUCT_VERTICAL_AXIS')

    status = 0

  END SUBROUTINE CONSTRUCT_VERTICAL_AXIS
  ! ====================================================================

! ****************************************************************************
END MODULE MESSY_MAIN_GRID_TRAFO_SCRP
! ****************************************************************************
