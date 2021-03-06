! ***********************************************************************
MODULE messy_lnox_dahl2000
! ***********************************************************************

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(LEN=*), PARAMETER :: submodstr = 'dahl2000'
  CHARACTER(LEN=*), PARAMETER :: submodver = '1.0'

  REAL(dp), PUBLIC ::   mg_range(100) = 0._dp
  REAL(dp), PUBLIC ::    v_range(100) = 0._dp
  REAL(dp), PUBLIC ::       v_gr(100) = 0._dp
  REAL(dp), PUBLIC ::    gr_diam(100) = 0._dp
  REAL(dp), PUBLIC :: rho_charge(100) = 0._dp
  REAL(dp), PUBLIC ::    delta_q(100) = 0._dp

  TYPE cluster_info
     INTEGER  ::   label
     INTEGER  ::   pixels
     INTEGER  ::   cent_pos_x
     INTEGER  ::   cent_pos_y
     INTEGER  ::   cent_pos_z
     INTEGER  ::   vert_extent
     INTEGER  ::   top_x 
     INTEGER  ::   top_y                              !   -
     INTEGER  ::   top_z                              !   - 
     REAL(dp) ::   bottom_height                      ! [m  ]
     REAL(dp) ::   top_height                         ! [m  ]
     REAL(dp) ::   average                            ! [depends]
  END TYPE cluster_info
  PUBLIC :: cluster_info

  TYPE capacitor
     INTEGER  ::   cap_label 
     ! upper plate:
     INTEGER  ::   top_label 
     INTEGER  ::   top_pixels         
     REAL(dp) ::   top_area             ! [m^2]
     REAL(dp) ::   top_depth            ! [m  ]
     INTEGER  ::   top_cent_z
     ! capacitor properties
     REAL(dp) ::   total_height         ! [m  ]
     REAL(dp) ::   separation           ! [m  ] SFC distance
     REAL(dp) ::   plate_distance       ! [m  ] center distance 
     ! (or 2D plate dist) 
     REAL(dp) ::   breakdown_alt        ! [m] between plates 
     INTEGER  ::   x_pos                ! horizontal coordinates
     INTEGER  ::   y_pos                                     

     ! lower plate:
     INTEGER  ::   bot_label
     INTEGER  ::   bot_pixels
     REAL(dp) ::   bot_area                ! [m^2]
     REAL(dp) ::   bot_depth               ! [m  ]
     INTEGER  ::   bot_cent_z

     REAL(dp) ::   graupel                 ! [kg/kg]  
     REAL(dp) ::   riming                  ! [kg/kg/s] 
  END TYPE capacitor
  PUBLIC :: capacitor

  INTEGER, PUBLIC, PARAMETER :: max_number_cap = 5000     ! max number of CBs

  TYPE(capacitor), PUBLIC :: info_cap(max_number_cap) 

  INTEGER, PUBLIC, ALLOCATABLE :: global_j(:)    ! output fron i/j_global

  ! Array containing vertical centroid 
  INTEGER, PUBLIC  ::  zrs(max_number_cap)       
  REAL (DP), PUBLIC, PARAMETER ::         &
       graupel_mass_correction = 1.2_dp,  &
       storm_width_correction  = 0.25_dp

  REAL(dp), PUBLIC, PARAMETER ::   &
       tcr         = 263.0_dp,      & ! Charge-reversal temperature  
       eps         = 8.854E-12_dp,  & ! Permittivity of air
       gamma       = 0.9_dp,        & ! how much lightning contributes to disch
       facgrdiam   = 1.5873E-4_dp,  &
       facgr       = 422.0_dp,      &  
       facrhoc     = 1.4E-10_dp,    & ! slope of grar-rho_ch relationship
       facqg       = 0.049_dp,      & 
       vfac_grar   = 6.0E9_dp,      &
       vfac_qg     = 12.0E9_dp,     &
       efac1       = -201.736_dp,   &
       efac2       = 1.E3 * 8.4_dp    ! ALT in km

  PUBLIC :: dahl_2010_initialize
  PUBLIC :: bin_field
  PUBLIC :: proplab
  PUBLIC :: two_elements
  PUBLIC :: capacitor_details

CONTAINS

  !============================================================================
  SUBROUTINE dahl_2010_initialize

    IMPLICIT NONE

    INTRINSIC :: EXP, LOG, REAL
    ! LOCAL
    INTEGER :: n_array(100)
    INTEGER :: i          


    !--------------------------------------------------------------------------
    ! Section 1: THE PARAMETER - VARIABLES
    !            
    !            Define how terminal fall velocity of graupel, space-charge
    !            density in the charging current, and the charge transferred in
    !            a flash (i.e., the parameterized variables) vary with 
    !            the parameters.
    !
    ! NOTE:
    !  In this implementation, a lookup-table style is used; a more efficient
    !  description is:
    !
    !  delta_q(V) = 25 * (1 - EXP(0.067 - 0.027 * V)), where V is the charge 
    !                                                  volume in km^3
    !  
    !  rho_ch(mg) = 4.467E-10 + 3.067E-9 * mg, where mg is the graupel mass
    !                                          in g/m^3     
    !  diam(mg)   = 1.833E-3 + 3.333E-3 * mg
    !
    !  If mg GT 3 g/m^3, then rho_ch(mg=3) and diam(mg=3) are to be used.
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! Section 1.1.: THE PARAMETERS
    ! Set up linearly-increasing arrays that contain possible ranges of 
    ! parameters
    ! used to determine charging current, channel-propagation depths, 
    ! and the like.
    !--------------------------------------------------------------------------

    ! Create array containing "independent variables"
    ! -----------------------------------------------

    DO i = 1, 100
       n_array(i)  = i
    ENDDO

    ! Range of graupel content
    ! ------------------------

    mg_range =  0.5_dp * (0.1_dp + 0.06_dp * REAL(n_array,dp))

    ! Range of volume
    ! ---------------

    ! volume ranges from [2.5, 297.025] km^3

    v_range = 1.0E9_dp * (2.5_dp + 2.975_dp * REAL(n_array - 1,dp)) 

    !--------------------------------------------------------------------------
    ! Section 1.2: THE PARAMETERIZED VARIABLES
    !
    !   Define how variables used by the flash-rate equation depend on
    !   parameters defined in Sec. 1.1.
    !--------------------------------------------------------------------------

    gr_diam = 0.002_dp + 0.0001_dp * REAL(n_array,dp) 

    v_gr    = facgr * EXP(0.89_dp * LOG(gr_diam))

    ! charge density in the charging current in C/m3
    ! ----------------------------------------------

    rho_charge = &
         2.0_dp * 0.3E-9 + 2.0_dp * (5.0E-9 - 0.4E-9) / 100. * REAL(n_array,dp)

    ! charge per flash in C
    ! ---------------------

    delta_q = 25.0_dp * (1.0_dp - EXP(-0.08_dp * REAL(n_array,dp)))

  END SUBROUTINE dahl_2010_initialize
  !============================================================================

  !============================================================================
  ELEMENTAL FUNCTION bin_field (field, binmar) RESULT(bin_field_result) 

    !--------------------------------------------------------------------------
    ! Description:
    !
    !  This internal procedure transforms the input field into a binary one,
    !  which
    !  may subsequently be labeled.  The field contains 0's and -1's, based on
    !  bnmar which determines, beyond which values the field is filtered. 
    !
    ! INPUT ARGUMENTS
    ! ----------------
    !
    !   field :   The field which is "binarized" and labeled
    !
    !   binmar:   Binarization margin: A filter, beyond which the values are
    !             set to zero, or -1, respectively.  This is required because
    !             the field which is to be labeled may only contain -1 and 
    !             zeros.
    !
    ! OUTPUT
    ! ------
    !
    !   A field that contains consecutive numbers as unique labels where the 
    !   input field contained isolated clusters/blobs of -1's.
    !
    !   TO-DO:
    !   ------
    !
    ! * Build in switch whether filter should single out values LT or GT binmar
    !   use "OPTIONAL" attribute (Check with PRESENT). 
    !   Also LT or LE (LT, LE) switches are needed
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    ! Function input arguments
    !--------------------------

    REAL(DP), INTENT(IN) :: binmar     ! binarization margin
    REAL(DP), INTENT(IN) :: field      ! input field

    ! Function result
    ! ---------------

    INTEGER  :: bin_field_result  ! output field

    ! Make binary field (0, -1) of the original field based on the margin
    ! "binmar".

    bin_field_result = -1
    IF (field <= binmar) bin_field_result = 0

    !--------------------------------------------------------------------------
    ! End function BIN_FIELD
    !--------------------------------------------------------------------------
  END FUNCTION bin_field
  !============================================================================

  !============================================================================
  ! SUBROUTINE PROPLAB(n, rr, prop, status)
  !----------------------------------------------------------------------------
  ! Purpose and method:
  !
  ! Internal function to iteratively find the end of a pointer path directed
  ! by the CSIZE array.
  !
  ! For a given position, scan through csize and follow the pointer path,
  ! starting from the input index rr.  The path is followed until a positive
  ! csize element is found.  The index pointing to the first positive element
  ! that is encountered is output as result of the function.  The trivial
  ! case of the input index already pointing at the end of a pointer path
  ! is included.
  !
  ! INPUT
  !
  !   n : csize array
  !   rr: index to rr'th element of csize - start of the pointer path
  !
  ! OUTPUT
  !
  !   prop = PROPLAB(csize, t):
  !   status = error status
  !
  !   Cluster label that the input index' cluster is associated with
  !----------------------------------------------------------------------------

  SUBROUTINE PROPLAB(n, rr, prop, status)

    IMPLICIT NONE

    !--------------------------------------------------------------------------
    ! Local variables/dummy arguments
    !--------------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: n(:)          ! csize; assumed shape from input
    INTEGER, INTENT(IN)  :: rr            ! input variables:
    !     rr: csize element
    !     (start of pointer path)
    INTEGER, INTENT(OUT) ::  prop         ! the end of the pointer path
    INTEGER, INTENT(OUT) ::  status       ! error status
    !  LOCAL
    INTEGER ::  pp                       ! loop variables
    INTEGER ::  t1, abs_t1               ! temporary/utility variables

    INTRINSIC :: ABS
    !--------------------------------------------------------------------------
    ! Do a simple iteration
    !--------------------------------------------------------------------------

    status = 0

    pp = 0
    t1 = n(rr)

    IF (n(rr) > 0 ) THEN
       prop = rr                      ! Already at end of the path

    ELSEIF (n(rr) < 0) THEN
       DO                             ! open loop (to EXIT label)
          pp = 1 + pp
          abs_t1 = ABS(t1)
          IF (n(abs_t1) > 0) THEN
             prop = abs_t1
             EXIT
          ENDIF

          t1 = n(abs_t1)

          ! Terminate iteration after 1000 steps

          IF (pp .EQ. 1000) THEN

             status = 700
             ! Some debugging information ...

             print *, &   
                  'PROPLAB: END OF POINTER PATH NOT FOUND AFTER 1000 ITERATIONS ... ABORTING'
             !   print *, my_cart_id, 'Time step: ', ntstep
             print *, ' CSIZE index, element: ', abs_t1, n(abs_t1)

             ! Terminate program
             RETURN

             ! CALL MPI_ABORT (MPI_COMM_WORLD, 100, ierr)
          ENDIF
       ENDDO
    ELSEIF (n(rr) == 0) THEN  
       prop = 0
       ! When csize element is zero (hinting at upstream error)
       print *, 'PROPLAB ERROR: REFERENCED CSIZE ELEMENT IS ZERO ... ABORTING'
       print *, 'Erroneous cluster label:', rr
       !  print *, my_cart_id, 'Time step: ', ntstep
       !  CALL MPI_ABORT (MPI_COMM_WORLD, 100, ierr)
       status = 710
       RETURN
    ENDIF

    !--------------------------------------------------------------------------
    ! End subroutine PROPLAB
    !--------------------------------------------------------------------------

  END SUBROUTINE PROPLAB
  ! ***********************************************************************

  !============================================================================

  !----------------------------------------------------------------------------
  !
  ! SUBROUTINE TWO_ELEMENTS:
  !
  ! PURPOSE
  !
  ! Internal procedure to determine, whether two of three elements
  ! contained in a three-element, 1D array, are identical.
  ! If so, a logical flag will be set to .TRUE.; in addition,
  ! the position of the value not equaling the two other elements will be
  ! determined, as well as the value of this element and the value of the
  ! two equal elements.  This should limit the number of IF-constructs in the
  ! calling subroutine.
  !
  ! INPUT
  !
  ! * elements
  !  A three-element, 1D array
  !
  ! OUTPUT
  !
  ! * eflag:
  !   A logical flag revealing whether or not precisely two elements are equal
  !
  !   --> Two of the three elements are equal:
  !       eflag = .TRUE.
  !   --> All elements being equal, or all elements being different:
  !       eflag = .FALSE.
  !
  ! * eposition:
  !   The position of the element that is not equal than the others.  Set to
  !   -999 if there is no uniquely different element (eflag == false).
  !
  ! * evalue:
  !   Value of the element being not equal to the two others. Set to
  !   -999 if there is no uniquely different element (eflag == false).
  !
  ! * ovalue:
  !   Value of the other two elements.  Set to -999 if there is no uniquely
  !   different element (eflag == false).
  !
  !----------------------------------------------------------------------------

  SUBROUTINE TWO_ELEMENTS(elements, eflag, eposition, evalue, ovalue)

    IMPLICIT NONE

    !--------------------------------------------------------------------------
    ! Dummy arguments
    !--------------------------------------------------------------------------

    ! Input arguments
    ! ---------------

    INTEGER, INTENT(IN) :: elements(3)

    ! Output arguments
    ! -----------------

    INTEGER, INTENT(OUT) :: &
         eposition,                             &
         evalue,                                &
         ovalue

    LOGICAL ::  eflag

    !--------------------------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------------------------

    INTEGER       :: e1, e2, e3

    !==========================================================================
    ! Start subroutine two_elements
    !==========================================================================

    e1 = elements(1)
    e2 = elements(2)
    e3 = elements(3)

    ! Abort if all equal or all different

    IF ( (e1 == e2 .AND. e2 == e3) .OR. (e1 /= e2 .AND. &
         e2 /= e3 .AND. e1 /=e3) ) THEN

       eflag     = .FALSE.
       ovalue    = -999
       evalue    = -999
       eposition = -999

       ! Continue if only two of the three elements are equal
       ! Time-consuming IF-constructs are used in lack of a better idea ...

    ELSEIF (e1 == e2 .AND. e2 /= e3 .OR.          &
         e1 == e3 .AND. e2 /= e3 .OR.          &
         e2 == e3 .AND. e1 /= e3) THEN

       eflag = .TRUE.

       ! Find which is the uniquely-valued element, as well as its position
       ! Not too elegant, and not debugged yet

       IF (e1 == e2) THEN
          eposition = 3
          evalue    = elements(3)
          ovalue    = elements(2)
          IF ( ovalue .NE. elements(1) ) print *, 'TWO_ELEMENTS ERROR - POINT 1'
       ELSEIF (e1 == e3) THEN
          eposition = 2
          evalue    = elements(2)
          ovalue    = elements(1)
          IF ( ovalue .NE. elements(3) ) print *, 'TWO_ELEMENTS ERROR - POINT 2'
       ELSEIF (e2 == e3) THEN
          eposition = 1
          evalue    = elements(1)
          ovalue    = elements(2)
          IF ( ovalue .NE. elements(3) ) print *, 'TWO_ELEMENTS ERROR - POINT 3'
       ENDIF

    ENDIF

    !==========================================================================
    ! End subroutine two_elements
    !==========================================================================

  END SUBROUTINE TWO_ELEMENTS
  !============================================================================

  !============================================================================
  ! Subroutine capacitor_details
  !============================================================================

  !----------------------------------------------------------------------------
  ! Description:
  !   This internal module procedure calculates the area as well as the depth
  !   of the overlapping clusters.  Also, the separatio distancen is 
  !   calculated. 
  !
  ! Method:
  !  The horizontal area: Calculated based on level of the cluster centroid.
  !  A horizontal slice through the centroid is considered and the pixels in
  !  this slice are counted.  Since the area size is known, the area of the
  !  entire cluster follows trivially.  The rotated grid results in squares
  !  quite closely, but still, dlon is scaled with R*COS(phi).  phi is taken
  !  at the centroid of the respective cell whose area is calculated.
  !  
  !  The vertical depth: Calculated by considering a vertical 1D slab 
  !  through the centroid.  The number of occupied levels is counted, and
  !  since the vertical grid spacing is known, the depth of the cluster 
  !  follows.  This is achieved with the help of the HHL field. 
  !
  ! Note:
  !  This is serial code, run on processor number zero.
  !----------------------------------------------------------------------------

  SUBROUTINE capacitor_details (zinfo_top, zinfo_bottom,   &
       global_top, global_bottom, &
       znumber_overlaps, zhfl, ke     &
       , ie_tot, je_tot, ke_tot, dlon, dlat &
       , startlat_tot, lprint)


    USE messy_main_constants_mem, ONLY: r_earth => radius_earth, pi
    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    IMPLICIT NONE

    ! 0D  variables
    ! -------------

    INTEGER,  INTENT(IN) :: znumber_overlaps   ! number of overlapping clusters
    INTEGER,  INTENT(IN) :: ie_tot, je_tot, ke_tot, ke
    REAL(DP), INTENT(IN) :: dlon, dlat, startlat_tot

    ! 1D structures
    ! -------------

    TYPE(cluster_info), INTENT(IN)  ::  &
         zinfo_top   (:),      &      ! cluster statistics of top cluster
         zinfo_bottom(:)              ! cluster statistics of bottom cluster 

    !TYPE(capacitor), INTENT(INOUT)  ::  &
    !  zinfo_capa(:)               ! capacitor-info structure

    ! 3D fields
    ! ---------

    INTEGER, INTENT(IN)  ::  &
         global_top   (ie_tot, je_tot, ke_tot),  &  ! global upper field
         global_bottom(ie_tot, je_tot, ke_tot)      ! global bottom field

    ! height of full levels
    REAL (DP),  INTENT(IN) :: zhfl(ie_tot, je_tot, ke_tot)

    ! 2D fields
    ! ---------

    ! 1D fields
    ! ---------

    !-------------------------------------------------------------------------
    ! Local variables
    !-------------------------------------------------------------------------

    LOGICAL,  INTENT(IN) :: lprint
    ! Local parameters
    ! ----------------

    ! None

    ! Local scalars
    ! -------------

    REAL (DP)  :: &
         zlat,           &          ! lat/lon at cell centroid 
         dx, dy     !         &          ! grid spacing

    INTEGER  ::  &
         i, j, k, ii,        &  ! loop variables
         zbot_label,                 &  ! label of lower cluster
         zcount,                     &  ! counter variable
         lzrs,                       &  ! scalar version of zrs
         xpos, ypos,                 &  ! position of cell centroid
         ztop_label,                 &  ! label of top cluster
         zvdist,                     &  ! distance of plates
         bot_cent                       ! INT of bottom centroid height

    ! Local dynamic 1D arrays
    ! -----------------------

    INTEGER, ALLOCATABLE  ::  &
         area_bot_pixels(:),  & ! pixels in horizontal slice
         summit_bottom(:),    & ! top height of bottom cluster 
         base_bottom(:),      & ! base height of bottom cluster
         summit_top(:),       & ! top height of upper cluster 
         base_top(:)            ! base height of upper cluster

    ! Local 2D arrays
    ! ---------------

    INTEGER  ::   utility_field_bot(ie_tot, je_tot) 

    !REAL (KIND=ireals) ::                  &
    !  qg_slice_glob_temp(ie_tot, je_tot),  &
    !  qg_slice_glob     (ie_tot, je_tot),  & 
    !  qg_slice_local    (ie,     je    )

    ! Local 3D arrays
    ! ---------------

    ! None

    INTRINSIC :: COS, INT, REAL
    !-------------------------------------------------------------------------
    ! Allocate dynamic arrays
    !-------------------------------------------------------------------------

    ALLOCATE(area_bot_pixels(znumber_overlaps))
    ALLOCATE(summit_bottom  (znumber_overlaps))
    ALLOCATE(base_bottom    (znumber_overlaps))
    ALLOCATE(summit_top     (znumber_overlaps))
    ALLOCATE(base_top       (znumber_overlaps))

    !=========================================================================
    ! Start routine
    !=========================================================================

    ! Initialize arrays

    summit_bottom = 0
    base_bottom   = 0
    summit_top    = 0
    base_top      = 0

    !-------------------------------------------------------------------------
    ! Section 1: Find area and depth of the "lower plate"
    !-------------------------------------------------------------------------

    ! Count pixels of slab through centroid

    zrs = 0

    DO ii = 1, znumber_overlaps

       ! Calculate area: Count pixels
       ! ----------------------------

       ztop_label = info_cap(ii) % top_label
       zbot_label = info_cap(ii) % bot_label
       zrs(ii)    = zinfo_bottom(zbot_label) % cent_pos_z  ! Module variable 
       lzrs       = zrs(ii)
       zcount     = 0
       utility_field_bot(:,:) = global_bottom(:,:,lzrs)

       DO j = 1, je_tot
          DO i = 1, ie_tot
             IF (utility_field_bot(i,j) == zbot_label) THEN
                zcount = zcount + 1
             ENDIF
          ENDDO
       ENDDO

       area_bot_pixels(ii) = zcount

       ! calculate the area
       ! ------------------

       xpos = zinfo_bottom(zbot_label) % cent_pos_x
       ypos = zinfo_bottom(zbot_label) % cent_pos_y

       zlat = startlat_tot + dlat * REAL( (ypos - 1) ,dp)

       dx = r_earth * COS(zlat * pi/180.0_dp) * dlon * pi/180.0_dp
       dy = r_earth * dlat * pi/180.0_dp    

       ! Witdh scaling factor  

       info_cap(ii) % bot_area = &
            storm_width_correction * REAL(area_bot_pixels(ii),dp) * dx * dy

       ! Calculate depth (lower plate)
       ! -----------------------------

       downward_loop: DO k = 1, ke
          IF (global_bottom(xpos, ypos, k) == zbot_label) THEN
             summit_bottom(ii) = INT(zhfl(xpos, ypos, k) )
             EXIT downward_loop
          ENDIF
       ENDDO downward_loop

       upward_loop: DO k = ke, 1, -1   
          IF (global_bottom(xpos, ypos, k) == zbot_label) THEN
             base_bottom(ii) = INT(zhfl(xpos, ypos, k) )
             EXIT upward_loop
          ENDIF
       ENDDO upward_loop

       info_cap(ii) % bot_depth = REAL(summit_bottom(ii) - base_bottom(ii),dp)

       !-----------------------------------------------------------------------
       ! Section 2: Calculate area and depth of upper plate 
       !
       !  A projection of the lower plate onto the upper plate is used.  In
       !  this case, the projection is one-to-one, i.e., the area is equal
       !  to that of the lower plate.  The upper-plate depth is calculated based
       !  on the centroid position of the lower plate.
       !-----------------------------------------------------------------------

       info_cap(ii) % top_area   = info_cap(ii) % bot_area 
       info_cap(ii) % top_pixels = area_bot_pixels(ii)

       ! Calculate depth (upper plate)
       ! -----------------------------

       loop_down: DO k = 1, ke
          IF (global_top(xpos, ypos, k) == ztop_label) THEN
             summit_top(ii) = INT(zhfl(xpos, ypos, k))
             EXIT loop_down
          ENDIF
       ENDDO loop_down

       info_cap(ii) % total_height = REAL(summit_top(ii),dp)

       loop_up: DO k = ke, 1, -1
          IF (global_top(xpos, ypos, k) == ztop_label) THEN
             base_top(ii) = INT(zhfl(xpos, ypos, k))
             EXIT loop_up 
          ENDIF
       ENDDO loop_up

       !! Added later
       info_cap(ii) % top_depth = REAL(summit_top(ii) - base_top(ii),dp)

       ! Calculate vertical separation distance
       ! --------------------------------------

       zvdist = base_top(ii) - summit_bottom(ii)

       ! Vertical separation distance is equal to the centroid distance.

       info_cap(ii) % separation = REAL(zvdist,dp)

       IF (zvdist <= 0) THEN   ! Overlap

          ! Round off 2nd term but keep it real

          bot_cent = INT(zhfl(xpos, ypos, zinfo_bottom(ii) % cent_pos_z))

          info_cap(ii) % plate_distance = &
               REAL(summit_bottom(ii) - bot_cent,dp)

          info_cap(ii) % separation = 0.0_dp

       ELSEIF (zvdist > 0) THEN

          info_cap(ii) % plate_distance =                &
               zhfl(xpos, ypos, zinfo_top(ii) % cent_pos_z) - &
               zhfl(xpos, ypos, zinfo_bottom(ii) % cent_pos_z)

       ENDIF

       ! IF the distance still turns out to be (smaller than) zero, 
       ! set it to 500 m     

       IF (info_cap(ii) % plate_distance <= 0.0) THEN
          info_cap(ii) % plate_distance = 500.0_dp
          IF (lprint) THEN
             !WRITE (*,*) 'ATTENTION: Separation distance error: ' &
             !, ntstep, ii, &
             WRITE (*,*) 'ATTENTION: Separation distance error: ', ii, &
                  info_cap(ii) % plate_distance 
          ENDIF
       ENDIF

       ! One-pixel cells have depth zero and cause an exception farther above.
       ! The following avoids that.

       IF (info_cap(ii) % top_depth < 1.E-5_dp) THEN
          info_cap(ii) % top_depth = 1.E-2_dp
       ENDIF
       IF (info_cap(ii) % bot_depth < 1.E-5_dp) THEN
          info_cap(ii) % bot_depth = 1.E-2_dp
       ENDIF

    ENDDO ! ii-loop through clusters

    !-------------------------------------------------------------------------
    ! Deallocate dynamic arrays
    !-------------------------------------------------------------------------

    DEALLOCATE(area_bot_pixels)
    DEALLOCATE(summit_bottom  )
    DEALLOCATE(base_bottom    )
    DEALLOCATE(summit_top     )
    DEALLOCATE(base_top       )

    !=========================================================================
    ! End subroutine capacitor_details
    !=========================================================================

  END SUBROUTINE capacitor_details
  !============================================================================

! ***********************************************************************
END MODULE messy_lnox_dahl2000
! ***********************************************************************
