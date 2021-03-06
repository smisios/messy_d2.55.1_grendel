MODULE messy_main_grid_trafo_scrp_base

  !***********************************************************************
  !
  ! This submodel contains the SCRIP tools that are used in MESSy 
  ! Author:  Pozzer Andrea, October 2007
  !
  !***********************************************************************
  !     Copyright (c) 1997, 1998 the Regents of the University of 
  !       California.
  !
  !     This software and ancillary information (herein called software) 
  !     called SCRIP is made available under the terms described here.  
  !     The software has been approved for release with associated 
  !     LA-CC Number 98-45.
  !
  !     Unless otherwise indicated, this software has been authored
  !     by an employee or employees of the University of California,
  !     operator of the Los Alamos National Laboratory under Contract
  !     No. W-7405-ENG-36 with the U.S. Department of Energy.  The U.S.
  !     Government has rights to use, reproduce, and distribute this
  !     software.  The public may copy and use this software without
  !     charge, provided that this Notice and any statement of authorship
  !     are reproduced on all copies.  Neither the Government nor the
  !     University makes any warranty, express or implied, or assumes
  !     any liability or responsibility for the use of this software.
  !
  !     If software is modified to produce derivative works, such modified
  !     software should be clearly marked, so as not to confuse it with 
  !     the version available from Los Alamos National Laboratory.
  !
  !***********************************************************************
  !
  ! USAGE (see messy_a2o_e5.f90 as example) :
  !
  !   * give some variables values (generally from namelist):
  !    
  !       Check SCRIP documentation for the meaning!
  ! 
  !       VARIABLE      | DEFAULT VALUE
  !       num_neighbors =  4 
  !       north_thresh  =  1.45 
  !       south_thresh  = -2.00
  !       max_subseg    =  100000
  !       max_iter      =  100
  !       converge      =  1E-10
  !
  !
  !   * CALCULATIONS :
  !
  !     grid1_mask =>       1D source grid mask
  !     grid1_center_lat => 1D source grid latitude fo grid centers
  !     grid1_center_lon => 1D source grid longitude fo grid centers
  !     grid1_area_in =>    1D source grid of area (in Steradian!) 
  !     grid1_corner_lat => 2D source grid corners latitude (:,numberofcorners) counterclockwise!
  !     grid1_corner_lon => 2D source grid corners longitude (:,numberofcorners) counterclockwise!
  !     grid2_mask =>       1D destination grid mask
  !     grid2_center_lat => 1D destination grid latitude fo grid centers
  !     grid2_center_lon => 1D destination grid longitude fo grid centers
  !     grid2_area_in =>    1D destination grid of area (in Steradian!) 
  !     grid2_corner_lat => 2D destination grid corners latitude (:,numberofcorners) counterclockwise!
  !     grid2_corner_lon => 2D destination grid corners longitude (:,numberofcorners) counterclockwise!
  !
  !     luse_grid_centers = .false.(conservative) /.true. (other transformations)
  !     map_type = map_type_conserv / map_type_bilinear / map_type_bicubic / map_type_distwgt
  !     norm_opt = norm_opt_none / norm_opt_frcarea / norm_opt_dstarea 
  !
  !     call remap_init(  number of corners in source grid,        &
  !                       grid boxes x direction source grid,      &
  !                       grid boxes y direction source grid,      & 
  !                       number of corners in destination grid,   &
  !                       grid boxes x direction destination grid, &
  !                       grid boxes y direction destination grid) 
  !     call bounds_calc
  !     call remap_vars(1)
  !     call remap_conserv
  !     call remap_vars(2)
  !     HERE CHECK THE GRID MASK FOR MISSING VALUES !!!
  !     HERE SAVE THE RESULTS ON LOCAL ARRAYS!!!
  !     call remap_dealloc
  !
  !
  ! N.B. lat-lon must be in radiants!!!
  !
  !***********************************************************************

  USE messy_main_constants_mem,  ONLY: DP, RTD, pi

  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

! um_ak_20121031+
  CHARACTER(LEN=*), PARAMETER :: scripversion = '3.1??'
! um_ak_20121031-

!mz_ap_20080714
      integer , parameter :: &
           atm    = 1,       &
           oces   = 2,       &
           oceu   = 3,       &
           ocev   = 4
           ! see store_link_cnsrv
      integer , dimension(:,:), allocatable, save :: &
              link_add1,  & ! min,max link add to restrict search
              link_add2     ! min,max link add to restrict search
      
!mz_ap_20080714+
      logical , save :: first_call = .true.
!mz_ap_20080714-

      integer , parameter :: &
           norm_opt_none    = 1,            &
           norm_opt_dstarea = 2,            &
           norm_opt_frcarea = 3

      integer , parameter :: &
           map_type_conserv  = 1,           &
           map_type_bilinear = 2,           &
           map_type_bicubic  = 3,           &
           map_type_distwgt  = 4

      integer , save :: &
           max_links_map1,  &! current size of link arrays
           num_links_map1,  &! actual number of links for remapping
           max_links_map2,  &! current size of link arrays
           num_links_map2,  &! actual number of links for remapping
           num_maps,        &! num of remappings for this grid pair
           num_wts,         &! num of weights used in remapping
           map_type,        &! identifier for remapping method
           norm_opt,        &! option for normalization (conserv only)
           resize_increment  ! default amount to increase array size

      integer , dimension(:), allocatable, save :: &
            grid1_add_map1, & ! grid1 address for each link in mapping 1
            grid2_add_map1, & ! grid2 address for each link in mapping 1
            grid1_add_map2, & ! grid1 address for each link in mapping 2
            grid2_add_map2  ! grid2 address for each link in mapping 2

      real (DP), dimension(:,:), allocatable, save :: &
           wts_map1, & ! map weights for each link (num_wts,max_links)
           wts_map2   ! map weights for each link (num_wts,max_links)

!-----------------------------------------------------------------------
! constant!
!-----------------------------------------------------------------------

      integer, parameter :: char_len = 80

      real (DP), parameter :: zero   = 0.0_dp,           &     
                              one    = 1.0_dp,           & 
                              two    = 2.0_dp,           & 
                              three  = 3.0_dp,           &      
                              four   = 4.0_dp,           &       
                              five   = 5.0_dp,           &       
                              half   = 0.5_dp,           &         
                              quart  = 0.25_dp,          &       
                              bignum = 1.e+20_dp,        &     
                              tiny   = 1.e-14_dp,        &   
!                              pi     = 3.14159265359_dp, &        
                              pi2    = two*pi,           &    
                              pih    = half*pi            

!-----------------------------------------------------------------------
!
!     variables that describe each grid
!
!-----------------------------------------------------------------------

      integer , save ::                    &
                   grid1_size, grid2_size, & ! total points on each grid
                   grid1_rank, grid2_rank, & ! rank of each grid
                   grid1_corners, grid2_corners ! number of corners
                                                ! for each grid cell

      integer , dimension(:), allocatable, save ::  &
                   grid1_dims, grid2_dims  ! size of each grid dimension

      character(char_len), save ::  &
                   grid1_name, grid2_name  ! name for each grid

      character(char_len), save :: & 
                  grid1_units,      & ! units for grid coords (degs/radians)
                  grid2_units         ! units for grid coords

      real (DP), parameter :: deg2rad = pi/180._dp ! conversion for deg to rads

!-----------------------------------------------------------------------
!
!     grid coordinates and masks
!
!-----------------------------------------------------------------------

      logical , dimension(:), pointer, save :: grid1_mask => NULL()  ! flag which cells participate
      logical , dimension(:), pointer, save :: grid2_mask => NULL()  ! flag which cells participate

      real (DP), dimension(:), pointer, save :: grid1_center_lat  => NULL()! lat/lon coordinates for
      real (DP), dimension(:), pointer, save :: grid1_center_lon  => NULL()! each grid center in radians
      real (DP), dimension(:), pointer, save :: grid2_center_lat  => NULL() 
      real (DP), dimension(:), pointer, save :: grid2_center_lon  => NULL()
      real (DP), dimension(:), pointer, save :: grid1_area        => NULL()! tot area of each grid1 cell
      real (DP), dimension(:), pointer, save :: grid2_area        => NULL()! tot area of each grid2 cell
      real (DP), dimension(:), pointer, save :: grid1_area_in     => NULL()! area of grid1 cell from file
      real (DP), dimension(:), pointer, save :: grid2_area_in     => NULL()! area of grid2 cell from file
      real (DP), dimension(:), pointer, save :: grid1_frac        => NULL()! fractional area of grid cells
      real (DP), dimension(:), pointer, save :: grid2_frac        => NULL()! participating in remapping

      real (DP), dimension(:,:), pointer, save :: grid1_corner_lat => NULL() ! lat/lon coordinates for
      real (DP), dimension(:,:), pointer, save :: grid1_corner_lon => NULL() ! each grid corner in radians
      real (DP), dimension(:,:), pointer, save :: grid2_corner_lat => NULL()
      real (DP), dimension(:,:), pointer, save :: grid2_corner_lon => NULL()

      logical , save ::  &
                   luse_grid_centers, & ! use centers for bounding boxes
                   luse_grid1_area,   & ! use area from grid file
                   luse_grid2_area   ! use area from grid file

      real (DP), dimension(:,:), allocatable, save :: &
                  grid1_bound_box,  & ! lat/lon bounding box for use
                  grid2_bound_box     ! in restricting grid searches

! op_pj_20130703+
!-----------------------------------------------------------------------
      real (DP), save :: babystep = 0.0001_dp
!-----------------------------------------------------------------------
! op_pj_20130703-

!-----------------------------------------------------------------------
!
!     module variables (for distance_weight)
!
!-----------------------------------------------------------------------

      integer , save :: & 
           num_neighbors = 4 ! num nearest neighbors to interpolate from

      real (DP), dimension(:), allocatable, save :: &
           coslat, sinlat, & ! cosine, sine of grid lats (for distance)
           coslon, sinlon, & ! cosine, sine of grid lons (for distance)
           wgtstmp           ! an array to hold the link weight

!-----------------------------------------------------------------------
!
!     module variables (for remap_conserv)
!
!-----------------------------------------------------------------------

      integer, DIMENSION(:), allocatable :: grid2_overlap ! overlapping points

      integer ,  save :: num_srch_cells ! num cells in restricted search arrays

      integer , dimension(:), allocatable, save ::  &
              srch_add  ! global address of cells in srch arrays

      real (DP), save ::      &
           north_thresh = 1.45_dp, & ! threshold for coord transf.
           south_thresh =-2.00_dp    ! threshold for coord transf.

      integer , save :: & 
! mz_bk_20090728+
!!$             max_subseg = 10000 ! max number of subsegments per segment
!!$                                ! to prevent infinite loop
             max_subseg = 1000000 ! max number of subsegments per segment
                                ! to prevent infinite loop
! mz_bk_20090728-

      real (DP), dimension(:,:), allocatable, save ::  &
           srch_corner_lat,  & ! lat of each corner of srch cells
           srch_corner_lon     ! lon of each corner of srch cells

      logical :: l_ver = .FALSE. !verbose
!      logical :: l_ver = .TRUE. !verbose

!-----------------------------------------------------------------------
!
!     module variables (for bicubic/bilinear)
!
!-----------------------------------------------------------------------

      integer , save :: &
          max_iter = 100   ! max iteration count for i,j iteration

      real (DP), save :: &
           converge = 1.e-10_dp  ! convergence criterion

!-----------------------------------------------------------------------


   PUBLIC :: remap_init         ! initialization remap variables 
   PUBLIC :: remap_vars         ! allocation pointers for calculations
   PUBLIC :: remap_conserv      ! routines for conservative remap
   PUBLIC :: remap_distwgt      ! routines for dist-weight remap
   PUBLIC :: remap_bilin        ! routines for bilinear interp
   PUBLIC :: remap_bicub        ! routines for bicubic  interp
   PUBLIC :: remap_dealloc      ! deallocation pointers


CONTAINS

!***********************************************************************

      subroutine remap_init(corners_grid1,dim1_grid1,dim2_grid1 &
           ,corners_grid2,dim1_grid2,dim2_grid2)


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: corners_grid1
      INTEGER, INTENT(IN) :: corners_grid2
      INTEGER, INTENT(IN) :: dim1_grid1
      INTEGER, INTENT(IN) :: dim2_grid1
      INTEGER, INTENT(IN) :: dim1_grid2
      INTEGER, INTENT(IN) :: dim2_grid2

    
       num_maps=1
      
       ALLOCATE(grid1_dims(2),grid2_dims(2))
       grid1_corners = corners_grid1
       grid2_corners = corners_grid2
       ! we are always using the calculated area 
       ! however script calculate it!
       luse_grid1_area =.TRUE.
       ! um_ak_20130301+
       IF (ASSOCIATED(grid2_area_in)) THEN
          luse_grid2_area =.TRUE.
       ELSE
          luse_grid2_area = .FALSE.
       ENDIF
       ! um_ak_20130301-
       grid1_rank = 2
       grid2_rank = 2
       grid1_size = dim1_grid1*dim2_grid1
       grid2_size = dim1_grid2*dim2_grid2
       grid1_dims(1)=dim1_grid1
       grid1_dims(2)=dim2_grid1
       grid2_dims(1)=dim1_grid2
       grid2_dims(2)=dim2_grid2

       allocate( grid1_frac      (grid1_size),     &
                 grid2_frac      (grid2_size),     &
                 grid1_area      (grid1_size),     &
                 grid2_area      (grid2_size),     &
                 grid1_bound_box (4, grid1_size),  & 
                 grid2_bound_box (4, grid2_size))

       ! Nullify array needed for calculations 

       grid1_frac(:) = 0.0_dp
       grid2_frac(:) = 0.0_dp
       grid1_area(:) = 0.0_dp
       grid2_area(:) = 0.0_dp
    
    
      end subroutine remap_init

!***********************************************************************

      subroutine remap_vars(phase)
    
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: phase
    
    !-----------------------------------------------------------------------
    !
    !     this routine initializes some variables and provides an initial
    !     allocation of arrays (fairly large so frequent resizing 
    !     unnecessary).
    !
    !-----------------------------------------------------------------------
    
    
          select case(phase)
          case(1)
          !-----------------------------------------------------------------------
          !
          !     determine the number of weights
          !
          !-----------------------------------------------------------------------
    
          select case (map_type)
          case(map_type_conserv)
            num_wts = 3
          case(map_type_bilinear)
            num_wts = 1
          case(map_type_bicubic)
            num_wts = 4
          case(map_type_distwgt)
            num_wts = 1
          end select

          !-----------------------------------------------------------------------
          !
          !     initialize num_links and set max_links to four times the largest 
          !     of the destination grid sizes initially (can be changed later).
          !     set a default resize increment to increase the size of link
          !     arrays if the number of links exceeds the initial size
          !   
          !-----------------------------------------------------------------------
    
          num_links_map1 = 0
          max_links_map1 = 4*grid2_size
!          if (num_maps > 1) then
!            num_links_map2 = 0
!            max_links_map1 = max(4*grid1_size,4*grid2_size)
!            max_links_map2 = max_links_map1
!          endif
    
          resize_increment = INT(0.1*max(grid1_size,grid2_size))
          
          !-----------------------------------------------------------------------
          !
          !     allocate address and weight arrays for mapping 1
          !   
          !-----------------------------------------------------------------------
          
          allocate (grid1_add_map1(max_links_map1), &
                   grid2_add_map1(max_links_map1),  &
                   wts_map1(num_wts, max_links_map1))
                   grid1_add_map1=1
                   grid2_add_map1=1
                   wts_map1=0._dp
    

          case(2)
          !-----------------------------------------------------------------------
          !
          !     reduce size of remapping arrays and then write remapping info
          !     to a file.
          !
          !-----------------------------------------------------------------------
          if (num_links_map1 /= max_links_map1) then
            call resize_remap_vars(1, num_links_map1-max_links_map1)
          endif
!          if ((num_maps > 1) .and. (num_links_map2 /= max_links_map2)) then
!            call resize_remap_vars(2, num_links_map2-max_links_map2)
!          endif
    
    
          end select
    
    
      end subroutine remap_vars

!***********************************************************************

      subroutine bounds_calc
    
        IMPLICIT NONE
    
          integer  :: n                  &! loop counter
!!$             , nele                      &! element loop counter
             , i,j                       &! logical 2d addresses
             , ip1,jp1                   &
             , n_add, e_add, ne_add      &
             , nx, ny
    
         real (dp), dimension(4) :: &
           tmp_lats, tmp_lons  ! temps for computing bounding boxes
    
    !-----------------------------------------------------------------------
    !
    !     compute bounding boxes for restricting future grid searches
    !
    !-----------------------------------------------------------------------
    
          if (.not. luse_grid_centers) then

            grid1_bound_box(1,:) = minval(grid1_corner_lat, DIM=1)
            grid1_bound_box(2,:) = maxval(grid1_corner_lat, DIM=1)
            grid1_bound_box(3,:) = minval(grid1_corner_lon, DIM=1)
            grid1_bound_box(4,:) = maxval(grid1_corner_lon, DIM=1)
    

            grid2_bound_box(1,:) = minval(grid2_corner_lat, DIM=1)
            grid2_bound_box(2,:) = maxval(grid2_corner_lat, DIM=1)
            grid2_bound_box(3,:) = minval(grid2_corner_lon, DIM=1)
            grid2_bound_box(4,:) = maxval(grid2_corner_lon, DIM=1)
            
          else
    
            nx = grid1_dims(1)
            ny = grid1_dims(2)
    
            do n=1,grid1_size
    
              !*** find N,S and NE points to this grid point
    
              j = (n - 1)/nx +1
              i = n - (j-1)*nx
    
              if (i < nx) then
                ip1 = i + 1
              else
                !*** assume cyclic
                ip1 = 1
                !*** but if it is not, correct
                e_add = (j - 1)*nx + ip1
                if (abs(grid1_center_lat(e_add) - & 
                       grid1_center_lat(n   )) > pih) then
                  ip1 = i
                endif
              endif
    
              if (j < ny) then
                jp1 = j+1
              else
                !*** assume cyclic
                jp1 = 1
                !*** but if it is not, correct
                n_add = (jp1 - 1)*nx + i
                if (abs(grid1_center_lat(n_add) - & 
                       grid1_center_lat(n   )) > pih) then
                  jp1 = j
                endif
              endif
    
              n_add = (jp1 - 1)*nx + i
              e_add = (j - 1)*nx + ip1
              ne_add = (jp1 - 1)*nx + ip1
    
              !*** find N,S and NE lat/lon coords and check bounding box
    
              tmp_lats(1) = grid1_center_lat(n)
              tmp_lats(2) = grid1_center_lat(e_add)
              tmp_lats(3) = grid1_center_lat(ne_add)
              tmp_lats(4) = grid1_center_lat(n_add)
    
              tmp_lons(1) = grid1_center_lon(n)
              tmp_lons(2) = grid1_center_lon(e_add)
              tmp_lons(3) = grid1_center_lon(ne_add)
              tmp_lons(4) = grid1_center_lon(n_add)
    
              grid1_bound_box(1,n) = minval(tmp_lats)
              grid1_bound_box(2,n) = maxval(tmp_lats)
              grid1_bound_box(3,n) = minval(tmp_lons)
              grid1_bound_box(4,n) = maxval(tmp_lons)
            end do
    
            nx = grid2_dims(1)
            ny = grid2_dims(2)
    
            do n=1,grid2_size
    
              !*** find N,S and NE points to this grid point
    
              j = (n - 1)/nx +1
              i = n - (j-1)*nx
    
              if (i < nx) then
                ip1 = i + 1
              else
                !*** assume cyclic
                ip1 = 1
                !*** but if it is not, correct
                e_add = (j - 1)*nx + ip1
                if (abs(grid2_center_lat(e_add) - & 
                       grid2_center_lat(n   )) > pih) then
                  ip1 = i
                endif
              endif
    
              if (j < ny) then
                jp1 = j+1
              else
                !*** assume cyclic
                jp1 = 1
                !*** but if it is not, correct
                n_add = (jp1 - 1)*nx + i
                if (abs(grid2_center_lat(n_add) - & 
                       grid2_center_lat(n   )) > pih) then
                  jp1 = j
                endif
              endif
    
              n_add = (jp1 - 1)*nx + i
              e_add = (j - 1)*nx + ip1
              ne_add = (jp1 - 1)*nx + ip1
    
              !*** find N,S and NE lat/lon coords and check bounding box
    
              tmp_lats(1) = grid2_center_lat(n)
              tmp_lats(2) = grid2_center_lat(e_add)
              tmp_lats(3) = grid2_center_lat(ne_add)
              tmp_lats(4) = grid2_center_lat(n_add)
    
              tmp_lons(1) = grid2_center_lon(n)
              tmp_lons(2) = grid2_center_lon(e_add)
              tmp_lons(3) = grid2_center_lon(ne_add)
              tmp_lons(4) = grid2_center_lon(n_add)
    
              grid2_bound_box(1,n) = minval(tmp_lats)
              grid2_bound_box(2,n) = maxval(tmp_lats)
              grid2_bound_box(3,n) = minval(tmp_lons)
              grid2_bound_box(4,n) = maxval(tmp_lons)
            end do
    
          endif
    
          where (abs(grid1_bound_box(4,:) - grid1_bound_box(3,:)) > pi)
            grid1_bound_box(3,:) = zero
            grid1_bound_box(4,:) = pi2
          end where
    
          where (abs(grid2_bound_box(4,:) - grid2_bound_box(3,:)) > pi)
            grid2_bound_box(3,:) = zero
            grid2_bound_box(4,:) = pi2
          end where
    
          !***
          !*** try to check for cells that overlap poles
          !***
    
          where (grid1_center_lat > grid1_bound_box(2,:)) &
           grid1_bound_box(2,:) = pih
    
          where (grid1_center_lat < grid1_bound_box(1,:)) &
           grid1_bound_box(1,:) = -pih
    
          where (grid2_center_lat > grid2_bound_box(2,:)) &
           grid2_bound_box(2,:) = pih
    
          where (grid2_center_lat < grid2_bound_box(1,:)) &
           grid2_bound_box(1,:) = -pih
    
    
      end subroutine bounds_calc

!***********************************************************************

      subroutine resize_remap_vars(nmap, increment)

!-----------------------------------------------------------------------
!
!     this routine resizes remapping arrays by increasing(decreasing)
!     the max_links by increment
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
           nmap,      & ! identifies which mapping array to resize
           increment    ! the number of links to add(subtract) to arrays

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: &
         ierr,    &  ! error flag
         mxlinks     ! size of link arrays

      integer , dimension(:), allocatable :: &
        add1_tmp, & ! temp array for resizing address arrays
        add2_tmp    ! temp array for resizing address arrays

      real (DP), dimension(:,:), allocatable :: &
          wts_tmp   ! temp array for resizing weight arrays

!-----------------------------------------------------------------------
!
!     resize map 1 arrays if required.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(grid1_add_map1)
        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), & 
                  wts_tmp(num_wts,mxlinks))

        add1_tmp = 1
        add2_tmp = 1
        wts_tmp  = 0._dp
        
        add1_tmp = grid1_add_map1
        add2_tmp = grid2_add_map1
        wts_tmp  = wts_map1
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (grid1_add_map1, grid2_add_map1, wts_map1)
        max_links_map1 = mxlinks + increment
        allocate (grid1_add_map1(max_links_map1), &
                  grid2_add_map1(max_links_map1), &
                  wts_map1(num_wts,max_links_map1))

        grid1_add_map1=1
        grid2_add_map1=1
        wts_map1=0._dp

        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***

        mxlinks = min(mxlinks, max_links_map1)
        grid1_add_map1(1:mxlinks) = add1_tmp (1:mxlinks)
        grid2_add_map1(1:mxlinks) = add2_tmp (1:mxlinks)
        wts_map1    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
        deallocate(add1_tmp, add2_tmp, wts_tmp)

!-----------------------------------------------------------------------
!
!     resize map 2 arrays if required.
!
!-----------------------------------------------------------------------

      case(2)

        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(grid1_add_map2)
        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), & 
                 wts_tmp(num_wts,mxlinks),stat=ierr)
        if (ierr .ne. 0) then
          print *,'error allocating temps in resize: ',ierr
          stop
        endif

        add1_tmp = grid1_add_map2
        add2_tmp = grid2_add_map2
        wts_tmp  = wts_map2
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (grid1_add_map2, grid2_add_map2, wts_map2)
        max_links_map2 = mxlinks + increment
        allocate (grid1_add_map2(max_links_map2), &
                  grid2_add_map2(max_links_map2), &
                  wts_map2(num_wts,max_links_map2),stat=ierr)
        if (ierr .ne. 0) then
          print *,'error allocating new arrays in resize: ',ierr
          stop
        endif


        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***

        mxlinks = min(mxlinks, max_links_map2)
        grid1_add_map2(1:mxlinks) = add1_tmp (1:mxlinks)
        grid2_add_map2(1:mxlinks) = add2_tmp (1:mxlinks)
        wts_map2    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
        deallocate(add1_tmp, add2_tmp, wts_tmp)

      end select

!-----------------------------------------------------------------------

      end subroutine resize_remap_vars

!***********************************************************************

      subroutine remap_conserv(status)

!-----------------------------------------------------------------------
!
!     this routine traces the perimeters of every grid cell on each
!     grid checking for intersections with the other grid and computing
!     line integrals for each subsegment.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

        INTEGER, INTENT(OUT), OPTIONAL :: status


      integer  ::        & 
             grid1_add,  &  ! current linear address for grid1 cell
             grid2_add,  &  ! current linear address for grid2 cell
             min_add,    &  ! addresses for restricting search of
             max_add,    &  !   destination grid
             n, nwgt,    &  ! generic counters
             corner,     &  ! corner of cell that segment starts from
             next_corn,  &  ! corner of cell that segment ends on
             num_subseg     ! number of subsegments 

      logical  ::     & 
             lcoinc,  & ! flag for coincident segments
             lrevers, & ! flag for reversing direction of segment
             lbegin   ! flag for first integration of a segment

      logical , dimension(:), allocatable :: &
              srch_mask   ! mask for restricting searches

      real (DP) ::  &
          intrsct_lat, intrsct_lon,       & ! lat/lon of next intersect
          beglat, endlat, beglon, endlon, & ! endpoints of current seg.
          norm_factor,                    & ! factor for normalizing wts
          delta                             ! precision 

      real (DP), dimension(:), allocatable :: &
            grid2_centroid_lat, grid2_centroid_lon, & ! centroid coords
            grid1_centroid_lat, grid1_centroid_lon    ! on each grid

      real (DP), dimension(2) :: begseg ! begin lat/lon for
                                                   ! full segment

      real (DP), dimension(6) :: weights ! local wgt array

!-----------------------------------------------------------------------
!
!     initialize centroid arrays
!
!-----------------------------------------------------------------------
      IF (PRESENT(status)) THEN
         status = 0
      END IF

      allocate( grid1_centroid_lat(grid1_size), &
                grid1_centroid_lon(grid1_size), &
                grid2_centroid_lat(grid2_size), &
                grid2_centroid_lon(grid2_size))

      grid1_centroid_lat = zero
      grid1_centroid_lon = zero
      grid2_centroid_lat = zero
      grid2_centroid_lon = zero

!-----------------------------------------------------------------------
!
!     integrate around each cell on grid1
!
!-----------------------------------------------------------------------

      allocate(srch_mask(grid2_size))

! Check overlapping point of the source grid

      delta = epsilon(1.)
      do grid1_add = 1,grid1_size
        IF (grid1_mask(grid1_add)) THEN
            DO n = grid1_add, grid1_size
              IF (n == grid1_add) cycle
              IF ((ABS(grid1_center_lon(grid1_add)-  &
                  grid1_center_lon(n))<delta).and.   &
                  (ABS(grid1_center_lat(grid1_add)-  &
                  grid1_center_lat(n))<delta)) THEN
                  grid1_mask(n) = .false.
                  exit
              END IF              
            END DO
        END IF
      END DO

! Check overlapping point of the target grid
      allocate(grid2_overlap(grid2_size))
      grid2_overlap = -1
      delta = epsilon(1.)
      do grid2_add = 1,grid2_size
        IF ((grid2_overlap(grid2_add)==-1).and.(grid2_mask(grid2_add))) &
            THEN
            DO n = grid2_add, grid2_size
              IF (n == grid2_add) cycle
              IF ((ABS(grid2_center_lon(grid2_add)-   &
                  grid2_center_lon(n))<delta).and.    & 
                  (ABS(grid2_center_lat(grid2_add)-   &
                  grid2_center_lat(n))<delta)) THEN
                  grid2_overlap(n) = grid2_add
                  grid2_overlap(grid2_add) = n            
                  grid2_mask(n) = .false.
                  exit
              END IF              
            END DO
        ELSE
            grid2_overlap(grid2_add) = -1
        END IF

      END DO

      IF (l_ver) print *,'grid1 sweep '
      do grid1_add = 1,grid1_size

!        call timer_start(1)
        min_add = 1
        max_add = grid2_size

        !***
        !*** further restrict searches using bounding boxes
        !***

        num_srch_cells = 0
        do grid2_add = min_add,max_add
          srch_mask(grid2_add) = (grid2_bound_box(1,grid2_add) <=      &  
                                  grid1_bound_box(2,grid1_add)) .and.  & 
                                 (grid2_bound_box(2,grid2_add) >=      & 
                                  grid1_bound_box(1,grid1_add)) .and.  & 
                                 (grid2_bound_box(3,grid2_add) <=      & 
                                  grid1_bound_box(4,grid1_add)) .and.  & 
                                 (grid2_bound_box(4,grid2_add) >=      & 
                                  grid1_bound_box(3,grid1_add))

          if (srch_mask(grid2_add)) num_srch_cells = num_srch_cells+1
        end do

        !***
        !*** create search arrays
        !***

        allocate(srch_add(num_srch_cells),                     &
                srch_corner_lat(grid2_corners,num_srch_cells), &
                srch_corner_lon(grid2_corners,num_srch_cells))

        n = 0
        gather1: do grid2_add = min_add,max_add
          if (srch_mask(grid2_add)) then
            n = n+1
            srch_add(n) = grid2_add
            srch_corner_lat(:,n) = grid2_corner_lat(:,grid2_add)
            srch_corner_lon(:,n) = grid2_corner_lon(:,grid2_add)
          endif
        end do gather1
!        call timer_stop(1)

        !***
        !*** integrate around this cell
        !***

        do corner = 1,grid1_corners
          next_corn = mod(corner,grid1_corners) + 1

          !***
          !*** define endpoints of the current segment
          !***

          beglat = grid1_corner_lat(corner,grid1_add)
          beglon = grid1_corner_lon(corner,grid1_add)
          endlat = grid1_corner_lat(next_corn,grid1_add)
          endlon = grid1_corner_lon(next_corn,grid1_add)
          lrevers = .false.

          !***
          !*** to ensure exact path taken during both
          !*** sweeps, always integrate segments in the same 
          !*** direction (SW to NE).
          !***

          if ((endlat < beglat) .or. &
              (endlat == beglat .and. endlon < beglon)) then 
            beglat = grid1_corner_lat(next_corn,grid1_add)
            beglon = grid1_corner_lon(next_corn,grid1_add)
            endlat = grid1_corner_lat(corner,grid1_add)
            endlon = grid1_corner_lon(corner,grid1_add)
            lrevers = .true.
          endif

          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.
          num_subseg = 0

          !***
          !*** if this is a constant-longitude segment, skip the rest 
          !*** since the line integral contribution will be zero.
          !***

          if (endlon /= beglon) then

          !***
          !*** integrate along this segment, detecting intersections 
          !*** and computing the line integral for each sub-segment
          !***
          
          !um_ch_20150629+
          !do while (beglat /= endlat .or. beglon /= endlon)         
          do while ((ABS(beglat - endlat) > 0.01*tiny) .or. &
               (ABS(beglon - endlon) > 0.01*tiny)) 
          !um_ch_20150629-
                
            !***
            !*** prevent infinite loops if integration gets stuck
            !*** near cell or threshold boundary
            !***

            num_subseg = num_subseg + 1
            if (num_subseg > max_subseg) then
               write (0,*) 'SCRIP: exceeded1 ', grid1_add,beglat*RTD , beglon*RTD ,  endlat*RTD,  endlon  *RTD 
               write (0,*) 'SCRIP: exceeded2 ', grid1_add,grid1_corner_lat(corner,grid1_add)*RTD, grid1_corner_lon(corner,grid1_add)*RTD, grid1_corner_lat(next_corn,grid1_add)*RTD, grid1_corner_lon(next_corn,grid1_add)*RTD

              IF (PRESENT(status)) THEN
                 status = 3998
                 RETURN
              ELSE
                 stop 'integration stalled: num_subseg exceeded limit'
              END IF
            endif

            !***
            !*** find next intersection of this segment with a grid
            !*** line on grid 2.
            !***

 !            call timer_start(2)
            call intersection(grid2_add,intrsct_lat,intrsct_lon,lcoinc, &
                             beglat, beglon, endlat, endlon, begseg,   &
                             lbegin, lrevers)
!            call timer_stop(2)
            lbegin = .false.

            !***
            !*** compute line integral for this subsegment.
            !***

!            call timer_start(3)
            if (grid2_add /= 0) then
              call line_integral(weights, num_wts,                       &
                               beglon, intrsct_lon, beglat, intrsct_lat, &
                               grid1_center_lat(grid1_add),              & 
                               grid1_center_lon(grid1_add),              &
                               grid2_center_lat(grid2_add),              & 
                               grid2_center_lon(grid2_add))
            else
              call line_integral(weights, num_wts,                       &
                               beglon, intrsct_lon, beglat, intrsct_lat, &
                               grid1_center_lat(grid1_add),              & 
                               grid1_center_lon(grid1_add),              &
                               grid1_center_lat(grid1_add),              & 
                               grid1_center_lon(grid1_add))
            endif
!            call timer_stop(3)

            !***
            !*** if integrating in reverse order, change
            !*** sign of weights
            !***

            if (lrevers) then
              weights = -weights
            endif

            !***
            !*** store the appropriate addresses and weights. 
            !*** also add contributions to cell areas and centroids.
            !***

            !if (grid1_add == 119247) then
            !  print *,grid1_add,grid2_add,corner,weights(1)
            !  print *,grid1_corner_lat(:,grid1_add)
            !  print *,grid1_corner_lon(:,grid1_add)
            !  print *,grid2_corner_lat(:,grid2_add)
            !  print *,grid2_corner_lon(:,grid2_add)
            !  print *,beglat,beglon,intrsct_lat,intrsct_lon
            !endif

            if (grid2_add /= 0) then
              if (grid1_mask(grid1_add)) then
!                call timer_start(4)
                call store_link_cnsrv(grid1_add, grid2_add, weights)
                IF (grid2_overlap(grid2_add)/=-1) then
                    call store_link_cnsrv(grid1_add,            & 
                        grid2_overlap(grid2_add), weights)
                ENDIF
!                call timer_stop(4)
                grid1_frac(grid1_add) = grid1_frac(grid1_add) + & 
                                       weights(1)
                grid2_frac(grid2_add) = grid2_frac(grid2_add) + & 
                                       weights(num_wts+1)
                IF (grid2_overlap(grid2_add)/=-1)               &
                    grid2_frac(grid2_overlap(grid2_add)) =       &
                    grid2_frac(grid2_overlap(grid2_add)) +       &
                    weights(num_wts+1)
              endif

            endif


            grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
            grid1_centroid_lat(grid1_add) =            & 
            grid1_centroid_lat(grid1_add) + weights(2) 
            grid1_centroid_lon(grid1_add) =            & 
            grid1_centroid_lon(grid1_add) + weights(3)  

            !***
            !*** reset beglat and beglon for next subsegment.
            !***

            beglat = intrsct_lat
            beglon = intrsct_lon
          end do ! do while

          endif

          !***
          !*** end of segment
          !***

        end do

        !***
        !*** finished with this cell: deallocate search array and
        !*** start on next cell

        deallocate(srch_add, srch_corner_lat, srch_corner_lon)

      end do

      deallocate(srch_mask)


!-----------------------------------------------------------------------
!
!     integrate around each cell on grid2
!
!-----------------------------------------------------------------------

      allocate(srch_mask(grid1_size))

      IF (l_ver) print *,'grid2 sweep '
      do grid2_add = 1,grid2_size


!        call timer_start(5)
        min_add = 1
        max_add = grid1_size

        !***
        !*** further restrict searches using bounding boxes
        !***

        num_srch_cells = 0
        do grid1_add = min_add, max_add
          srch_mask(grid1_add) = (grid1_bound_box(1,grid1_add) <=     & 
                                  grid2_bound_box(2,grid2_add)) .and. &
                                 (grid1_bound_box(2,grid1_add) >=     & 
                                  grid2_bound_box(1,grid2_add)) .and. &
                                 (grid1_bound_box(3,grid1_add) <=     & 
                                  grid2_bound_box(4,grid2_add)) .and. & 
                                 (grid1_bound_box(4,grid1_add) >=     & 
                                  grid2_bound_box(3,grid2_add))

          if (srch_mask(grid1_add)) num_srch_cells = num_srch_cells+1
        end do

        allocate(srch_add(num_srch_cells),  &
                srch_corner_lat(grid1_corners,num_srch_cells), &
                srch_corner_lon(grid1_corners,num_srch_cells))

        n = 0
        gather2: do grid1_add = min_add,max_add
          if (srch_mask(grid1_add)) then
            n = n+1
            srch_add(n) = grid1_add
            srch_corner_lat(:,n) = grid1_corner_lat(:,grid1_add)
            srch_corner_lon(:,n) = grid1_corner_lon(:,grid1_add)
          endif
        end do gather2
!        call timer_stop(5)

        !***
        !*** integrate around this cell
        !***

        do corner = 1,grid2_corners
          next_corn = mod(corner,grid2_corners) + 1

          beglat = grid2_corner_lat(corner,grid2_add)
          beglon = grid2_corner_lon(corner,grid2_add)
          endlat = grid2_corner_lat(next_corn,grid2_add)
          endlon = grid2_corner_lon(next_corn,grid2_add)
          lrevers = .false.

          !***
          !*** to ensure exact path taken during both
          !*** sweeps, always integrate in the same direction
          !***

          if ((endlat < beglat) .or. &
              (endlat == beglat .and. endlon < beglon)) then 
            beglat = grid2_corner_lat(next_corn,grid2_add)
            beglon = grid2_corner_lon(next_corn,grid2_add)
            endlat = grid2_corner_lat(corner,grid2_add)
            endlon = grid2_corner_lon(corner,grid2_add)
            lrevers = .true.
          endif

          begseg(1) = beglat
          begseg(2) = beglon
          lbegin = .true.

          !***
          !*** if this is a constant-longitude segment, skip the rest 
          !*** since the line integral contribution will be zero.
          !***

          if (endlon /= beglon) then
          num_subseg = 0

          !***
          !*** integrate along this segment, detecting intersections 
          !*** and computing the line integral for each sub-segment
          !***
          
          !um_ch_20150629+
          !do while (beglat /= endlat .or. beglon /= endlon)         
          do while ((ABS(beglat - endlat) > 0.01*tiny) .or. &
               (ABS(beglon - endlon) > 0.01*tiny)) 
          !um_ch_20150629-
             
            !***
            !*** prevent infinite loops if integration gets stuck
            !*** near cell or threshold boundary
            !***

            num_subseg = num_subseg + 1
            if (num_subseg > max_subseg) then
               write (0,*) 'SCRIP1: exceeded1 ', grid2_add,beglat*RTD , beglon*RTD ,  endlat*RTD,  endlon*RTD   
               write (0,*) 'SCRIP1: exceeded2 ', grid2_add,grid2_corner_lat(corner,grid2_add)*RTD, grid2_corner_lon(corner,grid2_add)*RTD, grid2_corner_lat(next_corn,grid2_add)*RTD, grid2_corner_lon(next_corn,grid2_add)*RTD
              IF (PRESENT(status)) THEN
                 status = 3998
                 RETURN
              ELSE
                 stop 'integration stalled: num_subseg exceeded limit'
              ENDIF
            endif

            !***
            !*** find next intersection of this segment with a line 
            !*** on grid 2.
            !***

!            call timer_start(6)
            call intersection(grid1_add,intrsct_lat,intrsct_lon,lcoinc, &
                             beglat, beglon, endlat, endlon, begseg,    &
                             lbegin, lrevers)
!            call timer_stop(6)
            lbegin = .false.

            !***
            !*** compute line integral for this subsegment.
            !***

!            call timer_start(7)
            if (grid1_add /= 0) then
              call line_integral(weights, num_wts,                       &  
                               beglon, intrsct_lon, beglat, intrsct_lat, & 
                               grid1_center_lat(grid1_add),              & 
                               grid1_center_lon(grid1_add),              &  
                               grid2_center_lat(grid2_add),              & 
                               grid2_center_lon(grid2_add))
            else
              call line_integral(weights, num_wts,                       &
                               beglon, intrsct_lon, beglat, intrsct_lat, &
                               grid2_center_lat(grid2_add),              & 
                               grid2_center_lon(grid2_add),              &
                               grid2_center_lat(grid2_add),              &
                               grid2_center_lon(grid2_add))
            endif
!            call timer_stop(7)

            if (lrevers) then
              weights = -weights
            endif

            !***
            !*** store the appropriate addresses and weights. 
            !*** also add contributions to cell areas and centroids.
            !*** if there is a coincidence, do not store weights
            !*** because they have been captured in the previous loop.
            !*** the grid1 mask is the master mask
            !***

            !if (grid1_add == 119247) then
            !  print *,grid1_add,grid2_add,corner,weights(1)
            !  print *,grid1_corner_lat(:,grid1_add)
            !  print *,grid1_corner_lon(:,grid1_add)
            !  print *,grid2_corner_lat(:,grid2_add)
            !  print *,grid2_corner_lon(:,grid2_add)
            !  print *,beglat,beglon,intrsct_lat,intrsct_lon
            !endif

            if (.not. lcoinc .and. grid1_add /= 0) then
              if (grid1_mask(grid1_add)) then
!                call timer_start(8)
                call store_link_cnsrv(grid1_add, grid2_add, weights)
!                call timer_stop(8)
                grid1_frac(grid1_add) = grid1_frac(grid1_add) + & 
                                        weights(1)
                grid2_frac(grid2_add) = grid2_frac(grid2_add) + &
                                        weights(num_wts+1)
             endif

            endif          

            grid2_area(grid2_add) = grid2_area(grid2_add) +    & 
                                            weights(num_wts+1)
            grid2_centroid_lat(grid2_add) =                    & 
            grid2_centroid_lat(grid2_add) + weights(num_wts+2)
            grid2_centroid_lon(grid2_add) =                    &
            grid2_centroid_lon(grid2_add) + weights(num_wts+3)

            !***
            !*** reset beglat and beglon for next subsegment.
            !***

            beglat = intrsct_lat
            beglon = intrsct_lon
          end do ! while

          endif

          !***
          !*** end of segment
          !***

        end do ! corner

        !***
        !*** finished with this cell: deallocate search array and
        !*** start on next cell

        deallocate(srch_add, srch_corner_lat, srch_corner_lon)

      end do ! grid2_add


      deallocate(grid2_overlap)
      deallocate(srch_mask)

!-----------------------------------------------------------------------
!
!     correct for situations where N/S pole not explicitly included in
!     grid (i.e. as a grid corner point). if pole is missing from only
!     one grid, need to correct only the area and centroid of that 
!     grid.  if missing from both, do complete weight calculation.
!
!-----------------------------------------------------------------------

      !*** North Pole
      weights(1) =  pi2
      weights(2) =  pi*pi
      weights(3) =  zero
      weights(4) =  pi2
      weights(5) =  pi*pi
      weights(6) =  zero

      grid1_add = 0
      pole_loop1: do n=1,grid1_size
        if (grid1_area(n) < -three*pih .and. &
            grid1_center_lat(n) > zero) then
          grid1_add = n
          exit pole_loop1
        endif
      end do pole_loop1

      grid2_add = 0
      pole_loop2: do n=1,grid2_size
        if (grid2_area(n) < -three*pih .and. &
            grid2_center_lat(n) > zero) then
          grid2_add = n
          exit pole_loop2
        endif
      end do pole_loop2

      if (grid1_add /=0) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) =            & 
        grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) =            &
        grid1_centroid_lon(grid1_add) + weights(3)
      endif

      if (grid2_add /=0) then
        grid2_area(grid2_add) = grid2_area(grid2_add) +    & 
                                        weights(num_wts+1)
        grid2_centroid_lat(grid2_add) =                    &
        grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) =                    &
        grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      endif

      if (grid1_add /= 0 .and. grid2_add /=0) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)
        grid1_frac(grid1_add) = grid1_frac(grid1_add) +   &
                                weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) +   &
                                weights(num_wts+1)
     endif

      !*** South Pole
      weights(1) =  pi2
      weights(2) = -pi*pi
      weights(3) =  zero
      weights(4) =  pi2
      weights(5) = -pi*pi
      weights(6) =  zero

      grid1_add = 0
      pole_loop3: do n=1,grid1_size
        if (grid1_area(n) < -three*pih .and.  &
            grid1_center_lat(n) < zero) then
          grid1_add = n
          exit pole_loop3
        endif
      end do pole_loop3

      grid2_add = 0
      pole_loop4: do n=1,grid2_size
        if (grid2_area(n) < -three*pih .and. &
            grid2_center_lat(n) < zero) then
          grid2_add = n
          exit pole_loop4
        endif
      end do pole_loop4

      if (grid1_add /=0) then
        grid1_area(grid1_add) = grid1_area(grid1_add) + weights(1)
        grid1_centroid_lat(grid1_add) =            &
        grid1_centroid_lat(grid1_add) + weights(2)
        grid1_centroid_lon(grid1_add) =            &
        grid1_centroid_lon(grid1_add) + weights(3)
      endif

      if (grid2_add /=0) then
        grid2_area(grid2_add) = grid2_area(grid2_add) +    & 
                                        weights(num_wts+1)
        grid2_centroid_lat(grid2_add) =                    &
        grid2_centroid_lat(grid2_add) + weights(num_wts+2)
        grid2_centroid_lon(grid2_add) =                    &
        grid2_centroid_lon(grid2_add) + weights(num_wts+3)
      endif

      if (grid1_add /= 0 .and. grid2_add /=0) then
        call store_link_cnsrv(grid1_add, grid2_add, weights)

        grid1_frac(grid1_add) = grid1_frac(grid1_add) + & 
                                weights(1)
        grid2_frac(grid2_add) = grid2_frac(grid2_add) + &
                                weights(num_wts+1)
     endif

!-----------------------------------------------------------------------
!
!     finish centroid computation
!
!-----------------------------------------------------------------------

      where (grid1_area /= zero)
        grid1_centroid_lat = grid1_centroid_lat/grid1_area
        grid1_centroid_lon = grid1_centroid_lon/grid1_area
      end where

      where (grid2_area /= zero)
        grid2_centroid_lat = grid2_centroid_lat/grid2_area
        grid2_centroid_lon = grid2_centroid_lon/grid2_area
      end where

!-----------------------------------------------------------------------
!
!     include centroids in weights and normalize using destination
!     area if requested
!
!-----------------------------------------------------------------------

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        do nwgt=1,num_wts
          weights(        nwgt) = wts_map1(nwgt,n)
!          if (num_maps > 1) then
!            weights(num_wts+nwgt) = wts_map2(nwgt,n)
!          endif
        end do

        select case(norm_opt)
        case (norm_opt_dstarea)
          if (grid2_area(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = one/grid2_area_in(grid2_add)
            else
              norm_factor = one/grid2_area(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_frcarea)
           if (grid2_frac(grid2_add) /= zero) then
            if (luse_grid2_area) then
              norm_factor = grid2_area(grid2_add)/    &
                           (grid2_frac(grid2_add)*    &
                            grid2_area_in(grid2_add))
            else
              norm_factor = one/grid2_frac(grid2_add)
            endif
          else
            norm_factor = zero
          endif
        case (norm_opt_none)
          norm_factor = one
        end select

        wts_map1(1,n) =  weights(1)*norm_factor
        wts_map1(2,n) = (weights(2) - weights(1)*                    &
                                    grid1_centroid_lat(grid1_add))*  &
                                    norm_factor
        wts_map1(3,n) = (weights(3) - weights(1)*                    &
                                    grid1_centroid_lon(grid1_add))*  &
                                    norm_factor
!        if (num_maps > 1) then
!          select case(norm_opt)
!          case (norm_opt_dstarea)
!            if (grid1_area(grid1_add) /= zero) then
!              if (luse_grid1_area) then
!                norm_factor = one/grid1_area_in(grid1_add)
!              else
!                norm_factor = one/grid1_area(grid1_add)
!              endif
!            else
!              norm_factor = zero
!            endif
!          case (norm_opt_frcarea)
!            if (grid1_frac(grid1_add) /= zero) then
!              if (luse_grid1_area) then
!                norm_factor = grid1_area(grid1_add)/    &
!                             (grid1_frac(grid1_add)*    &
!                              grid1_area_in(grid1_add))
!              else
!                norm_factor = one/grid1_frac(grid1_add)
!              endif
!            else
!              norm_factor = zero
!            endif
!          case (norm_opt_none)
!            norm_factor = one
!          end select
!
!          wts_map2(1,n) =  weights(num_wts+1)*norm_factor
!          wts_map2(2,n) = (weights(num_wts+2) - weights(num_wts+1)*   &
!                                      grid2_centroid_lat(grid2_add))* &
!                                      norm_factor
!          wts_map2(3,n) = (weights(num_wts+3) - weights(num_wts+1)*   &
!                                      grid2_centroid_lon(grid2_add))* &
!                                      norm_factor
!        endif

      end do

      if (l_ver) WRITE(*,*) ' Total number of links : ',num_links_map1

      where (grid1_area /= zero) grid1_frac = grid1_frac/grid1_area
      where (grid2_area /= zero) grid2_frac = grid2_frac/grid2_area

!-----------------------------------------------------------------------
!
!     perform some error checking on final weights
!
!-----------------------------------------------------------------------

      grid2_centroid_lat = zero
      grid2_centroid_lon = zero

      do n=1,grid1_size
        if (grid1_area(n) < -.01) then
                IF (l_ver) print *,'Grid 1 area error: ',n,grid1_area(n)
        endif
        if (grid1_centroid_lat(n) < -pih-.01 .or. &
            grid1_centroid_lat(n) >  pih+.01) then
            IF (l_ver) print *,'Grid 1 centroid lat error: ',n,grid1_centroid_lat(n)
        endif
        grid1_centroid_lat(n) = zero
        grid1_centroid_lon(n) = zero
      end do

      do n=1,grid2_size
        if (grid2_area(n) < -.01) then
                IF (l_ver) print *,'Grid 2 area error: ',n,grid2_area(n)
        endif
        if (grid2_centroid_lat(n) < -pih-.01 .or. &
            grid2_centroid_lat(n) >  pih+.01) then
            IF (l_ver) print *,'Grid 2 centroid lat error: ',n,grid2_centroid_lat(n)
        endif
        grid2_centroid_lat(n) = zero
        grid2_centroid_lon(n) = zero
      end do

      do n=1,num_links_map1
        grid1_add = grid1_add_map1(n)
        grid2_add = grid2_add_map1(n)
        
        ! um_ak_20160128+
        if (wts_map1(1,n) < -.01) then
!                IF (l_ver) print *,'Map 1 weight < 0 ',grid1_add,grid2_add,wts_map1(1,n), 'LON1', grid1_corner_lon(:,grid1_add),'LAT1', grid1_corner_lat(:,grid1_add),'LON2', grid2_corner_lon(:,grid2_add),'LAT2', grid2_corner_lat(:,grid2_add)
                IF (l_ver) print *,'Map 1 weight < 0 Deg error ',grid1_add &
                     ,grid2_add,wts_map1(1,n) &
                     , 'LON1', grid1_corner_lon(:,grid1_add)*RTD &
                     , 'LAT1', grid1_corner_lat(:,grid1_add)*RTD &
                     , 'LON2', grid2_corner_lon(:,grid2_add)*RTD &
                     , 'LAT2', grid2_corner_lat(:,grid2_add)*RTD
        endif
        if (norm_opt /= norm_opt_none .and. wts_map1(1,n) > 1.01) then
!                IF (l_ver) print *,'Map 1 weight > 1 ',grid1_add,grid2_add,wts_map1(1,n), 'LON1', grid1_corner_lon(:,grid1_add),'LAT1', grid1_corner_lat(:,grid1_add),'LON2', grid2_corner_lon(:,grid2_add),'LAT2', grid2_corner_lat(:,grid2_add)
                IF (l_ver) print *,'Map 1 weight > 1 Deg error ', n, grid1_add &
                     ,grid2_add,wts_map1(1,n) &
                     , 'LON1', grid1_corner_lon(:,grid1_add)*RTD &
                     , 'LAT1', grid1_corner_lat(:,grid1_add)*RTD &
                     , 'LON2', grid2_corner_lon(:,grid2_add)*RTD &
                     , 'LAT2', grid2_corner_lat(:,grid2_add)*RTD
        endif
!        if (wts_map1(1,n) < -.01) then
!                IF (l_ver) print *,'Map 1 weight < 0 ',grid1_add,grid2_add,wts_map1(1,n)
!        endif
!        if (norm_opt /= norm_opt_none .and. wts_map1(1,n) > 1.01) then
!                IF (l_ver) print *,'Map 1 weight > 1 ',grid1_add,grid2_add,wts_map1(1,n)
!        endif
        ! um_ak_20160128-
        grid2_centroid_lat(grid2_add) = &
        grid2_centroid_lat(grid2_add) + wts_map1(1,n)

!        if (num_maps > 1) then
!          if (wts_map2(1,n) < -.01) then
!                  IF (l_ver) print *,'Map 2 weight < 0 ',grid1_add,grid2_add, &
!                                        wts_map2(1,n)
!          endif
!          if (norm_opt /= norm_opt_none .and. wts_map2(1,n) > 1.01) then
!                  IF (l_ver) print *,'Map 2 weight < 0 ',grid1_add,grid2_add, &
!                                        wts_map2(1,n)
!          endif
!          grid1_centroid_lat(grid1_add) =  &
!          grid1_centroid_lat(grid1_add) + wts_map2(1,n)
!        endif
      end do

      do n=1,grid2_size
        grid2_add = n
        select case(norm_opt)
        case (norm_opt_dstarea)
          norm_factor = grid2_frac(grid2_add)
        case (norm_opt_frcarea)
          norm_factor = one
        case (norm_opt_none)
          if (luse_grid2_area) then
            norm_factor = grid2_area_in(grid2_add)
          else
            norm_factor = grid2_area(grid2_add)
          endif
        end select
        if (abs(grid2_centroid_lat(grid2_add)-norm_factor) > .01) then
! um_ak_20130318+
           IF (l_ver) print *,'Error: sum of wts for map1 ',grid2_add, &
!!$           IF (l_ver .AND.&
!!$                ( grid2_centroid_lat(grid2_add) > 0._dp .OR. &
!!$                grid2_centroid_lat(grid2_add) < 0._dp)) &
!!$                print *,'Error: sum of wts for map1 ',grid2_add, &
! um_ak_20130318-
                  grid2_centroid_lat(grid2_add),norm_factor
        endif
      end do



!      if (num_maps > 1) then
!        do n=1,grid1_size
!          grid1_add = n
!          select case(norm_opt)
!          case (norm_opt_dstarea)
!            norm_factor = grid1_frac(grid1_add)
!          case (norm_opt_frcarea)
!            norm_factor = one
!          case (norm_opt_none)
!            if (luse_grid1_area) then
!              norm_factor = grid1_area_in(grid1_add)
!            else
!              norm_factor = grid1_area(grid1_add)
!            endif
!          end select
!          if (abs(grid1_centroid_lat(grid1_add)-norm_factor) > .01) then
!            print *,'Error: sum of wts for map2 ',grid1_add, &
!                   grid1_centroid_lat(grid1_add),norm_factor
!          endif
!        end do
!      endif
!-----------------------------------------------------------------------

!mz_ap_20080722+
      deallocate( grid1_centroid_lat, &
                  grid1_centroid_lon, &
                  grid2_centroid_lat, &
                  grid2_centroid_lon)
      IF (ALLOCATED(link_add1)) deallocate( link_add1,link_add2 )   
       first_call = .true.
!mz_ap_20080722-

      end subroutine remap_conserv

!***********************************************************************

      subroutine intersection(location,intrsct_lat,intrsct_lon,lcoinc, &
                              beglat, beglon, endlat, endlon, begseg,  &
                              lbegin, lrevers)

!-----------------------------------------------------------------------
!
!     this routine finds the next intersection of a destination grid 
!     line with the line segment given by beglon, endlon, etc.
!     a coincidence flag is returned if the segment is entirely 
!     coincident with an ocean grid line.  the cells in which to search
!     for an intersection must have already been restricted in the
!     calling routine.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      logical , intent(in) :: &
          lbegin, & ! flag for first integration along this segment
          lrevers  ! flag whether segment integrated in reverse

      real (DP), intent(in) :: &
           beglat, beglon,     & ! beginning lat/lon endpoints for segment
           endlat, endlon        ! ending    lat/lon endpoints for segment

      real (DP), dimension(2), intent(inout) :: & 
           begseg ! begin lat/lon of full segment

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer , intent(out) :: &
              location  ! address in destination array containing this
                        ! segment

      logical , intent(out) :: &
              lcoinc    ! flag segments which are entirely coincident
                        ! with a grid line

      real (DP), intent(out) :: &
           intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n, next_n, cell, srch_corners !!$, pole_loc

      integer , save :: & 
           last_loc  ! save location when crossing threshold

      logical  :: & 
           loutside  ! flags points outside grid

      logical , save ::  &
           lthresh = .false.  ! flags segments crossing threshold bndy

      real (DP) ::             &
           lon1, lon2,         & ! local longitude variables for segment
           lat1, lat2,         & ! local latitude  variables for segment
           grdlon1, grdlon2,   & ! local longitude variables for grid cell
           grdlat1, grdlat2,   & ! local latitude  variables for grid cell
           vec1_lat, vec1_lon, & ! vectors and cross products used
           vec2_lat, vec2_lon, & ! during grid search
           cross_product,      & 
           eps, offset,        & ! small offset away from intersect
           s1, s2, determ,     & ! variables used for linear solve to
           mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection

      real (DP), save ::  &
           intrsct_lat_off, intrsct_lon_off ! lat/lon coords offset 
                                            ! for next search

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon

      if (num_srch_cells == 0) return

      if (beglat > north_thresh .or. beglat < south_thresh) then

        if (lthresh) location = last_loc
        call pole_intersection(location,                              &
                     intrsct_lat,intrsct_lon,lcoinc,lthresh,          &
                     beglat, beglon, endlat, endlon, begseg, lrevers)
        if (lthresh) then
          last_loc = location
          intrsct_lat_off = intrsct_lat
          intrsct_lon_off = intrsct_lon
        endif
       return

      endif

      loutside = .false.
      if (lbegin) then
        lat1 = beglat
        lon1 = beglon
      else
        lat1 = intrsct_lat_off
        lon1 = intrsct_lon_off
      endif
      lat2 = endlat
      lon2 = endlon
      if ((lon2-lon1) > three*pih) then
        lon2 = lon2 - pi2
      else if ((lon2-lon1) < -three*pih) then
        lon2 = lon2 + pi2
      endif
      s1 = zero

!-----------------------------------------------------------------------
!
!     search for location of this segment in ocean grid using cross
!     product method to determine whether a point is enclosed by a cell
!
!-----------------------------------------------------------------------

!      call timer_start(12)
      srch_corners = size(srch_corner_lat,DIM=1)

      srch_loop: do

        !***
        !*** if last segment crossed threshold, use that location
        !***

        if (lthresh) then
          do cell=1,num_srch_cells
            if (srch_add(cell) == last_loc) then
              location = last_loc
              eps = tiny
              exit srch_loop
            endif
          end do
        endif

        !***
        !*** otherwise normal search algorithm
        !***

        cell_loop: do cell=1,num_srch_cells
          corner_loop: do n=1,srch_corners
            next_n = MOD(n,srch_corners) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each cell side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** positive, the point is contained in the cell.
            !***

            vec1_lat = srch_corner_lat(next_n,cell) - & 
                       srch_corner_lat(n     ,cell)
            vec1_lon = srch_corner_lon(next_n,cell) - &
                       srch_corner_lon(n     ,cell)
            vec2_lat = lat1 - srch_corner_lat(n,cell)
            vec2_lon = lon1 - srch_corner_lon(n,cell)

            !***
            !*** if endpoint coincident with vertex, offset
            !*** the endpoint
            !***

            if (vec2_lat == 0 .and. vec2_lon == 0) then
              lat1 = lat1 + 1.d-10*(lat2-lat1)
              lon1 = lon1 + 1.d-10*(lon2-lon1)
              vec2_lat = lat1 - srch_corner_lat(n,cell)
              vec2_lon = lon1 - srch_corner_lon(n,cell)
            endif

            !***
            !*** check for 0,2pi crossings
            !***

            if (vec1_lon >  pi) then
              vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -pi) then
              vec1_lon = vec1_lon + pi2
            endif
            if (vec2_lon >  pi) then
              vec2_lon = vec2_lon - pi2
            else if (vec2_lon < -pi) then
              vec2_lon = vec2_lon + pi2
            endif

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat

            !***
            !*** if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the side is degenerate
            !***   (zero length).  if degenerate, set the cross 
            !***   product to a positive number.  otherwise perform 
            !***   another cross product between the side and the 
            !***   segment itself. 
            !*** if this cross product is also zero, the line is 
            !***   coincident with the cell boundary - perform the 
            !***   dot product and only choose the cell if the dot 
            !***   product is positive (parallel vs anti-parallel).
            !***

            if (cross_product == zero) then
              if (vec1_lat /= zero .or. vec1_lon /= zero) then
                vec2_lat = lat2 - lat1
                vec2_lon = lon2 - lon1

                if (vec2_lon >  pi) then
                  vec2_lon = vec2_lon - pi2
                else if (vec2_lon < -pi) then
                  vec2_lon = vec2_lon + pi2
                endif

                cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat
              else
                cross_product = one
              endif

              if (cross_product == zero) then
                lcoinc = .true.
                cross_product = vec1_lon*vec2_lon + vec1_lat*vec2_lat
                if (lrevers) cross_product = -cross_product
              endif
            endif

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (cross_product < zero) exit corner_loop

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > srch_corners) then
            location = srch_add(cell)

            !***
            !*** if the beginning of this segment was outside the
            !*** grid, invert the segment so the intersection found
            !*** will be the first intersection with the grid
            !***

            if (loutside) then
              lat2 = beglat
              lon2 = beglon
              location = 0
              eps  = -tiny
            else
              eps  = tiny
            endif

            exit srch_loop
          endif

          !***
          !*** otherwise move on to next cell
          !***

        end do cell_loop

        !***
        !*** if still no cell found, the point lies outside the grid.
        !***   take some baby steps along the segment to see if any
        !***   part of the segment lies inside the grid.  
        !***

        loutside = .true.
        !mz_ap_20090317+
        ! NB: we decreased the babystep values:
        ! for some resolution we need better resolution
        ! to see if the grid intersect....
! op_pj_20130703+
!!$        !s1 = s1 + 0.001_dp
!!$        s1 = s1 + 0.0001_dp
        s1 = s1 + babystep
! op_pj_20130703+
        lat1 = beglat + s1*(endlat - beglat)
        lon1 = beglon + s1*(lon2   - beglon)

        !***
        !*** reached the end of the segment and still outside the grid
        !*** return no intersection
        !***

        if (s1 >= one) return

      end do srch_loop
!      call timer_stop(12)

!-----------------------------------------------------------------------
!
!     now that a cell is found, search for the next intersection.
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

!      call timer_start(13)
      intrsct_loop: do n=1,srch_corners
        next_n = mod(n,srch_corners) + 1

        grdlon1 = srch_corner_lon(n     ,cell)
        grdlon2 = srch_corner_lon(next_n,cell)
        grdlat1 = srch_corner_lat(n     ,cell)
        grdlat2 = srch_corner_lat(next_n,cell)

        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = lat2 - lat1
        mat2 = grdlat1 - grdlat2
        mat3 = lon2 - lon1
        mat4 = grdlon1 - grdlon2
        rhs1 = grdlat1 - lat1
        rhs2 = grdlon1 - lon1

        if (mat3 >  pi) then
          mat3 = mat3 - pi2
        else if (mat3 < -pi) then
          mat3 = mat3 + pi2
        endif
        if (mat4 >  pi) then
          mat4 = mat4 - pi2
        else if (mat4 < -pi) then
          mat4 = mat4 + pi2
        endif
        if (rhs2 >  pi) then
          rhs2 = rhs2 - pi2
        else if (rhs2 < -pi) then
          rhs2 = rhs2 + pi2
        endif

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident.  coincidences were detected 
        !***   above so do nothing.
        !*** if the determinant is non-zero, solve for the linear 
        !***   parameters s for the intersection point on each line 
        !***   segment.
        !*** if 0<s1,s2<1 then the segment intersects with this side.
        !***   return the point of intersection (adding a small
        !***   number so the intersection is off the grid line).
        !***

        if (abs(determ) > 1.e-30) then

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zero .and. s2 <= one .and.  &
              s1 >  zero .and. s1 <= one) then

            !***
            !*** recompute intersection based on full segment
            !*** so intersections are consistent for both sweeps
            !***

            if (.not. loutside) then
              mat1 = lat2 - begseg(1)
              mat3 = lon2 - begseg(2)
              rhs1 = grdlat1 - begseg(1)
              rhs2 = grdlon1 - begseg(2)
            else
              mat1 = begseg(1) - endlat
              mat3 = begseg(2) - endlon
              rhs1 = grdlat1 - endlat
              rhs2 = grdlon1 - endlon
            endif

            if (mat3 >  pi) then
              mat3 = mat3 - pi2
            else if (mat3 < -pi) then
              mat3 = mat3 + pi2
            endif
            if (rhs2 >  pi) then
              rhs2 = rhs2 - pi2
            else if (rhs2 < -pi) then
              rhs2 = rhs2 + pi2
            endif

            determ = mat1*mat4 - mat2*mat3

            !***
            !*** sometimes due to roundoff, the previous 
            !*** determinant is non-zero, but the lines
            !*** are actually coincident.  if this is the
            !*** case, skip the rest.
            !***

            if (determ /= zero) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ

              offset = s1 + eps/determ
              if (offset > one) offset = one

              if (.not. loutside) then
                intrsct_lat = begseg(1) + mat1*s1
                intrsct_lon = begseg(2) + mat3*s1
                intrsct_lat_off = begseg(1) + mat1*offset
                intrsct_lon_off = begseg(2) + mat3*offset
              else
                intrsct_lat = endlat + mat1*s1
                intrsct_lon = endlon + mat3*s1
                intrsct_lat_off = endlat + mat1*offset
                intrsct_lon_off = endlon + mat3*offset
              endif
              exit intrsct_loop
            endif

          endif
        endif

        !***
        !*** no intersection this side, move on to next side
        !***

      end do intrsct_loop
!      call timer_stop(13)

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude.  only check if this was not a
!     threshold segment since sometimes coordinate transform can end
!     up on other side of threshold again.
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat < north_thresh .or. intrsct_lat > south_thresh) &
            lthresh = .false.
      else if (lat1 > zero .and. intrsct_lat > north_thresh) then
        intrsct_lat = north_thresh + tiny
        intrsct_lat_off = north_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      else if (lat1 < zero .and. intrsct_lat < south_thresh) then
        intrsct_lat = south_thresh - tiny
        intrsct_lat_off = south_thresh + eps*mat1
        s1 = (intrsct_lat - begseg(1))/mat1
        intrsct_lon     = begseg(2) + s1*mat3
        intrsct_lon_off = begseg(2) + (s1+eps)*mat3
        last_loc = location
        lthresh = .true.
      endif

!-----------------------------------------------------------------------

      end subroutine intersection

!***********************************************************************

      subroutine pole_intersection(location,                            &
                       intrsct_lat,intrsct_lon,lcoinc,lthresh,          &
                       beglat, beglon, endlat, endlon, begseg, lrevers)

!-----------------------------------------------------------------------
!
!     this routine is identical to the intersection routine except
!     that a coordinate transformation (using a Lambert azimuthal
!     equivalent projection) is performed to treat polar cells more
!     accurately.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real (DP), intent(in) :: & 
           beglat, beglon,     & ! beginning lat/lon endpoints for segment
           endlat, endlon        ! ending    lat/lon endpoints for segment

      real (DP), dimension(2), intent(inout) :: & 
           begseg ! begin lat/lon of full segment

      logical , intent(in) :: &
              lrevers   ! flag true if segment integrated in reverse

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      integer , intent(inout) :: &
             location  ! address in destination array containing this
                       ! segment -- also may contain last location on
                       ! entry

      logical , intent(out) :: &
              lcoinc    ! flag segment coincident with grid line

      logical , intent(inout) :: &
              lthresh   ! flag segment crossing threshold boundary

      real (DP), intent(out) :: &
           intrsct_lat, intrsct_lon ! lat/lon coords of next intersect.

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n, next_n, cell, srch_corners !!$, pole_loc

      logical  :: loutside ! flags points outside grid

      real (DP) :: pi4, rns,  & ! north/south conversion
          x1, x2,            & ! local x variables for segment
          y1, y2,            & ! local y variables for segment
          begx, begy,        & ! beginning x,y variables for segment
          endx, endy,        & ! beginning x,y variables for segment
          begsegx, begsegy,  & ! beginning x,y variables for segment
          grdx1, grdx2,      & ! local x variables for grid cell
          grdy1, grdy2,      & ! local y variables for grid cell
          vec1_y, vec1_x,    & ! vectors and cross products used
          vec2_y, vec2_x,    & ! during grid search
          cross_product, eps,& ! eps=small offset away from intersect
          s1, s2, determ,    & ! variables used for linear solve to
          mat1, mat2, mat3, mat4, rhs1, rhs2  ! find intersection

      real (DP), dimension(:,:), allocatable :: &
          srch_corner_x, & ! x of each corner of srch cells
          srch_corner_y   ! y of each corner of srch cells

      !***
      !*** save last intersection to avoid roundoff during coord
      !*** transformation
      !***

      logical , save :: luse_last = .false.

      real (DP), save :: & 
           intrsct_x, intrsct_y  ! x,y for intersection

      !***
      !*** variables necessary if segment manages to hit pole
      !***

      integer , save :: & 
          avoid_pole_count = 0  ! count attempts to avoid pole

      real (DP), save :: &
          avoid_pole_offset = tiny  ! endpoint offset to avoid pole

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      if (.not. lthresh) location = 0
      lcoinc = .false.
      intrsct_lat = endlat
      intrsct_lon = endlon

      loutside = .false.
      s1 = zero

!-----------------------------------------------------------------------
!
!     convert coordinates
!
!-----------------------------------------------------------------------

      allocate(srch_corner_x(size(srch_corner_lat,DIM=1), &
                             size(srch_corner_lat,DIM=2)),&
               srch_corner_y(size(srch_corner_lat,DIM=1), &
                             size(srch_corner_lat,DIM=2)))

      if (beglat > zero) then
        pi4 = quart*pi
        rns = one
      else
        pi4 = -quart*pi
        rns = -one
      endif

      if (luse_last) then
        x1 = intrsct_x
        y1 = intrsct_y
      else
        x1 = rns*two*sin(pi4 - half*beglat)*cos(beglon)
        y1 =     two*sin(pi4 - half*beglat)*sin(beglon)
        luse_last = .true.
      endif
      x2 = rns*two*sin(pi4 - half*endlat)*cos(endlon)
      y2 =     two*sin(pi4 - half*endlat)*sin(endlon)
      srch_corner_x = rns*two*sin(pi4 - half*srch_corner_lat)* &
                              cos(srch_corner_lon)
      srch_corner_y =     two*sin(pi4 - half*srch_corner_lat)* &
                              sin(srch_corner_lon)

      begx = x1
      begy = y1
      endx = x2
      endy = y2
      begsegx = rns*two*sin(pi4 - half*begseg(1))*cos(begseg(2))
      begsegy =     two*sin(pi4 - half*begseg(1))*sin(begseg(2))
      intrsct_x = endx
      intrsct_y = endy

!-----------------------------------------------------------------------
!
!     search for location of this segment in ocean grid using cross
!     product method to determine whether a point is enclosed by a cell
!
!-----------------------------------------------------------------------

!      call timer_start(12)
      srch_corners = size(srch_corner_lat,DIM=1)
      srch_loop: do

        !***
        !*** if last segment crossed threshold, use that location
        !***

        if (lthresh) then
          do cell=1,num_srch_cells
            if (srch_add(cell) == location) then
              eps = tiny
              exit srch_loop
            endif
          end do
        endif

        !***
        !*** otherwise normal search algorithm
        !***

        cell_loop: do cell=1,num_srch_cells
          corner_loop: do n=1,srch_corners
            next_n = MOD(n,srch_corners) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each cell side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** positive, the point is contained in the cell.
            !***

            vec1_x = srch_corner_x(next_n,cell) - & 
                     srch_corner_x(n     ,cell)
            vec1_y = srch_corner_y(next_n,cell) - &
                     srch_corner_y(n     ,cell)
            vec2_x = x1 - srch_corner_x(n,cell)
            vec2_y = y1 - srch_corner_y(n,cell)

            !***
            !*** if endpoint coincident with vertex, offset
            !*** the endpoint
            !***

            if (vec2_x == 0 .and. vec2_y == 0) then
              x1 = x1 + 1.d-10*(x2-x1)
              y1 = y1 + 1.d-10*(y2-y1)
              vec2_x = x1 - srch_corner_x(n,cell)
              vec2_y = y1 - srch_corner_y(n,cell)
            endif

            cross_product = vec1_x*vec2_y - vec2_x*vec1_y

            !***
            !*** if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the length of a side
            !***   is zero.  if the length is zero set det > 0.
            !***   otherwise, perform another cross 
            !***   product between the side and the segment itself. 
            !*** if this cross product is also zero, the line is 
            !***   coincident with the cell boundary - perform the 
            !***   dot product and only choose the cell if the dot 
            !***   product is positive (parallel vs anti-parallel).
            !***

            if (cross_product == zero) then
              if (vec1_x /= zero .or. vec1_y /= 0) then
                vec2_x = x2 - x1
                vec2_y = y2 - y1
                cross_product = vec1_x*vec2_y - vec2_x*vec1_y
              else
                cross_product = one
              endif

              if (cross_product == zero) then
                lcoinc = .true.
                cross_product = vec1_x*vec2_x + vec1_y*vec2_y
                if (lrevers) cross_product = -cross_product
              endif
            endif

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (cross_product < zero) exit corner_loop

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > srch_corners) then
            location = srch_add(cell)

            !***
            !*** if the beginning of this segment was outside the
            !*** grid, invert the segment so the intersection found
            !*** will be the first intersection with the grid
            !***

            if (loutside) then
              x2 = begx
              y2 = begy
              location = 0
              eps  = -tiny
            else
              eps  = tiny
            endif

            exit srch_loop
          endif

          !***
          !*** otherwise move on to next cell
          !***

        end do cell_loop

        !***
        !*** if no cell found, the point lies outside the grid.
        !***   take some baby steps along the segment to see if any
        !***   part of the segment lies inside the grid.  
        !***

        loutside = .true.
        s1 = s1 + 0.001_dp
        x1 = begx + s1*(x2 - begx)
        y1 = begy + s1*(y2 - begy)

        !***
        !*** reached the end of the segment and still outside the grid
        !*** return no intersection
        !***

        if (s1 >= one) then
          deallocate(srch_corner_x, srch_corner_y)
          luse_last = .false.
          return
        endif

      end do srch_loop
!      call timer_stop(12)

!-----------------------------------------------------------------------
!
!     now that a cell is found, search for the next intersection.
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

!      call timer_start(13)
      intrsct_loop: do n=1,srch_corners
        next_n = mod(n,srch_corners) + 1

        grdy1 = srch_corner_y(n     ,cell)
        grdy2 = srch_corner_y(next_n,cell)
        grdx1 = srch_corner_x(n     ,cell)
        grdx2 = srch_corner_x(next_n,cell)

        !***
        !*** set up linear system to solve for intersection
        !***

        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident or one segment has zero length.  
        !***   coincidences were detected above so do nothing.
        !*** if the determinant is non-zero, solve for the linear 
        !***   parameters s for the intersection point on each line 
        !***   segment.
        !*** if 0<s1,s2<1 then the segment intersects with this side.
        !***   return the point of intersection (adding a small
        !***   number so the intersection is off the grid line).
        !***

        if (abs(determ) > 1.e-30) then

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zero .and. s2 <= one .and. &
              s1 >  zero .and. s1 <= one) then

            !***
            !*** recompute intersection using entire segment
            !*** for consistency between sweeps
            !***

            if (.not. loutside) then
              mat1 = x2 - begsegx
              mat3 = y2 - begsegy
              rhs1 = grdx1 - begsegx
              rhs2 = grdy1 - begsegy
            else 
              mat1 = x2 - endx
              mat3 = y2 - endy
              rhs1 = grdx1 - endx
              rhs2 = grdy1 - endy
            endif

            determ = mat1*mat4 - mat2*mat3

            !***
            !*** sometimes due to roundoff, the previous 
            !*** determinant is non-zero, but the lines
            !*** are actually coincident.  if this is the
            !*** case, skip the rest.
            !***

            if (determ /= zero) then
              s1 = (rhs1*mat4 - mat2*rhs2)/determ
              s2 = (mat1*rhs2 - rhs1*mat3)/determ

              if (.not. loutside) then
                intrsct_x = begsegx + s1*mat1
                intrsct_y = begsegy + s1*mat3
              else 
                intrsct_x = endx + s1*mat1
                intrsct_y = endy + s1*mat3
              endif

              !***
              !*** convert back to lat/lon coordinates
              !***

              intrsct_lon = rns*atan2(intrsct_y,intrsct_x)
              if (intrsct_lon < zero) & 
                intrsct_lon = intrsct_lon + pi2

              if (abs(intrsct_x) > 1.d-10) then
                intrsct_lat = (pi4 - & 
                  asin(rns*half*intrsct_x/cos(intrsct_lon)))*two
              else if (abs(intrsct_y) > 1.d-10) then
                intrsct_lat = (pi4 - & 
                  asin(half*intrsct_y/sin(intrsct_lon)))*two
              else
                intrsct_lat = two*pi4
              endif

              !***
              !*** add offset in transformed space for next pass.
              !***

              if (s1 - eps/determ < one) then
                intrsct_x = intrsct_x - mat1*(eps/determ)
                intrsct_y = intrsct_y - mat3*(eps/determ)
              else
                if (.not. loutside) then
                  intrsct_x = endx
                  intrsct_y = endy
                  intrsct_lat = endlat
                  intrsct_lon = endlon
                else 
                  intrsct_x = begsegx
                  intrsct_y = begsegy
                  intrsct_lat = begseg(1)
                  intrsct_lon = begseg(2)
                endif
              endif

              exit intrsct_loop
            endif
          endif
        endif

        !***
        !*** no intersection this side, move on to next side
        !***

      end do intrsct_loop
!      call timer_stop(13)

      deallocate(srch_corner_x, srch_corner_y)

!-----------------------------------------------------------------------
!
!     if segment manages to cross over pole, shift the beginning 
!     endpoint in order to avoid hitting pole directly
!     (it is ok for endpoint to be pole point)
!
!-----------------------------------------------------------------------

      if (abs(intrsct_x) < 1.e-10 .and. abs(intrsct_y) < 1.e-10 .and. &
          (endx /= zero .and. endy /=0)) then
        if (avoid_pole_count > 2) then
           avoid_pole_count = 0
           avoid_pole_offset = 10.*avoid_pole_offset
        endif

        cross_product = begsegx*(endy-begsegy) - begsegy*(endx-begsegx)
        intrsct_lat = begseg(1)
        if (cross_product*intrsct_lat > zero) then
          intrsct_lon = beglon    + avoid_pole_offset
          begseg(2)   = begseg(2) + avoid_pole_offset
        else
          intrsct_lon = beglon    - avoid_pole_offset
          begseg(2)   = begseg(2) - avoid_pole_offset
        endif

        avoid_pole_count = avoid_pole_count + 1
        luse_last = .false.
      else
        avoid_pole_count = 0
        avoid_pole_offset = tiny
      endif

!-----------------------------------------------------------------------
!
!     if the segment crosses a pole threshold, reset the intersection
!     to be the threshold latitude and do not reuse x,y intersect
!     on next entry.  only check if did not cross threshold last
!     time - sometimes the coordinate transformation can place a
!     segment on the other side of the threshold again
!
!-----------------------------------------------------------------------

      if (lthresh) then
        if (intrsct_lat > north_thresh .or. intrsct_lat < south_thresh) &
          lthresh = .false.
      else if (beglat > zero .and. intrsct_lat < north_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
        intrsct_lat = north_thresh - tiny
        s1 = (north_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      else if (beglat < zero .and. intrsct_lat > south_thresh) then
        mat4 = endlat - begseg(1)
        mat3 = endlon - begseg(2)
        if (mat3 >  pi) mat3 = mat3 - pi2
        if (mat3 < -pi) mat3 = mat3 + pi2
        intrsct_lat = south_thresh + tiny
        s1 = (south_thresh - begseg(1))/mat4
        intrsct_lon = begseg(2) + s1*mat3
        luse_last = .false.
        lthresh = .true.
      endif

      !***
      !*** if reached end of segment, do not use x,y intersect 
      !*** on next entry
      !***

      if (intrsct_lat == endlat .and. intrsct_lon == endlon) then
        luse_last = .false.
      endif

!-----------------------------------------------------------------------

      end subroutine pole_intersection

!***********************************************************************

      subroutine line_integral(weights, num_wts,               & 
                             in_phi1, in_phi2, theta1, theta2, &
                             grid1_lat, grid1_lon, grid2_lat, grid2_lon)

!-----------------------------------------------------------------------
!
!     this routine computes the line integral of the flux function 
!     that results in the interpolation weights.  the line is defined
!     by the input lat/lon of the endpoints.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
              num_wts  ! number of weights to compute

      real (DP), intent(in) ::   & 
           in_phi1, in_phi2,     & ! longitude endpoints for the segment
           theta1, theta2,       & ! latitude  endpoints for the segment
           grid1_lat, grid1_lon, & ! reference coordinates for each
           grid2_lat, grid2_lon    ! grid (to ensure correct 0,2pi interv.

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real (DP), dimension(2*num_wts), intent(out) :: &
           weights   ! line integral contribution to weights

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real (DP) :: dphi, sinth1, sinth2, costh1, costh2, fac, &
                              phi1, phi2 !!$, phidiff1, phidiff2, sinint
      real (DP) :: f1, f2, fint

!-----------------------------------------------------------------------
!
!     weights for the general case based on a trapezoidal approx to
!     the integrals.
!
!-----------------------------------------------------------------------

      sinth1 = SIN(theta1)
      sinth2 = SIN(theta2)
      costh1 = COS(theta1)
      costh2 = COS(theta2)

      dphi = in_phi1 - in_phi2
      if (dphi >  pi) then
        dphi = dphi - pi2
      else if (dphi < -pi) then
        dphi = dphi + pi2
      endif
      dphi = half*dphi

!-----------------------------------------------------------------------
!
!     the first weight is the area overlap integral. the second and
!     fourth are second-order latitude gradient weights.
!
!-----------------------------------------------------------------------

      weights(        1) = dphi*(sinth1 + sinth2)
      weights(num_wts+1) = dphi*(sinth1 + sinth2)
      weights(        2) = dphi*(costh1 + costh2 + (theta1*sinth1 + &
                                                    theta2*sinth2))
      weights(num_wts+2) = dphi*(costh1 + costh2 + (theta1*sinth1 + &
                                                    theta2*sinth2))

!-----------------------------------------------------------------------
!
!     the third and fifth weights are for the second-order phi gradient
!     component.  must be careful of longitude range.
!
!-----------------------------------------------------------------------

      f1 = half*(costh1*sinth1 + theta1)
      f2 = half*(costh2*sinth2 + theta2)

      phi1 = in_phi1 - grid1_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid1_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(3) = half*phi1*(phi1-fac)*f1 - &
                     half*phi2*(phi2+fac)*f2 + &
                     half*fac*(phi1+phi2)*fint
      endif

      phi1 = in_phi1 - grid2_lon
      if (phi1 >  pi) then
        phi1 = phi1 - pi2
      else if (phi1 < -pi) then
        phi1 = phi1 + pi2
      endif

      phi2 = in_phi2 - grid2_lon
      if (phi2 >  pi) then
        phi2 = phi2 - pi2
      else if (phi2 < -pi) then
        phi2 = phi2 + pi2
      endif

      if ((phi2-phi1) <  pi .and. (phi2-phi1) > -pi) then
        weights(num_wts+3) = dphi*(phi1*f1 + phi2*f2)
      else
        if (phi1 > zero) then
          fac = pi
        else
          fac = -pi
        endif
        fint = f1 + (f2-f1)*(fac-phi1)/abs(dphi)
        weights(num_wts+3) = half*phi1*(phi1-fac)*f1 - &
                             half*phi2*(phi2+fac)*f2 + &
                             half*fac*(phi1+phi2)*fint
      endif

!-----------------------------------------------------------------------

      end subroutine line_integral

!***********************************************************************

      subroutine store_link_cnsrv(add1, add2, weights)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for this link in
!     the appropriate address and weight arrays and resizes those
!     arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
              add1,  & ! address on grid1
              add2     ! address on grid2

      real (DP), dimension(:), intent(in) :: &
              weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: nlink, min_link, max_link ! link index

!mz_ap_20080722+
      ! move as global variable so that we can deallocate them
      !integer , dimension(:,:), allocatable, save :: &
      !        link_add1,  & ! min,max link add to restrict search
      !        link_add2     ! min,max link add to restrict search
      !
      !logical , save :: first_call = .true.
!mz_ap_20080722-

!-----------------------------------------------------------------------
!
!     if all weights are zero, do not bother storing the link
!
!-----------------------------------------------------------------------

      if (all(weights == zero)) return

!-----------------------------------------------------------------------
!
!     restrict the range of links to search for existing links
!
!-----------------------------------------------------------------------

      if (first_call) then
        allocate(link_add1(2,grid1_size), link_add2(2,grid2_size))
        link_add1 = 0
        link_add2 = 0
        first_call = .false.
        min_link = 1
        max_link = 0
      else
        min_link = min(link_add1(1,add1),link_add2(1,add2))
        max_link = max(link_add1(2,add1),link_add2(2,add2))
        if (min_link == 0) then
          min_link = 1
          max_link = 0
        endif
      endif

!-----------------------------------------------------------------------
!
!     if the link already exists, add the weight to the current weight
!     arrays
!
!-----------------------------------------------------------------------

      do nlink=min_link,max_link
        if (add1 == grid1_add_map1(nlink)) then
        if (add2 == grid2_add_map1(nlink)) then

          wts_map1(:,nlink) = wts_map1(:,nlink) + weights(1:num_wts)
!          if (num_maps == 2) then
!            wts_map2(:,nlink) = wts_map2(:,nlink) + & 
!                                        weights(num_wts+1:2*num_wts)
!          endif
          return

        endif
        endif
      end do

!-----------------------------------------------------------------------
!
!     if the link does not yet exist, increment number of links and 
!     check to see if remap arrays need to be increased to accomodate 
!     the new link.  then store the link.
!
!-----------------------------------------------------------------------

      num_links_map1  = num_links_map1 + 1
      if (num_links_map1 > max_links_map1) & 
         call resize_remap_vars(1,resize_increment)

      grid1_add_map1(num_links_map1) = add1
      grid2_add_map1(num_links_map1) = add2
      wts_map1    (:,num_links_map1) = weights(1:num_wts)

!      if (num_maps > 1) then
!        num_links_map2  = num_links_map2 + 1
!        if (num_links_map2 > max_links_map2) & 
!          call resize_remap_vars(2,resize_increment)
!
!        grid1_add_map2(num_links_map2) = add1
!        grid2_add_map2(num_links_map2) = add2
!        wts_map2    (:,num_links_map2) = weights(num_wts+1:2*num_wts)
!      endif

      if (link_add1(1,add1) == 0) link_add1(1,add1) = num_links_map1
      if (link_add2(1,add2) == 0) link_add2(1,add2) = num_links_map1
      link_add1(2,add1) = num_links_map1
      link_add2(2,add2) = num_links_map1


!-----------------------------------------------------------------------

      end subroutine store_link_cnsrv

!***********************************************************************

      subroutine remap_distwgt

!-----------------------------------------------------------------------
!
!     this routine computes the inverse-distance weights for a
!     nearest-neighbor interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      logical , dimension(num_neighbors) :: &
           nbr_mask        ! mask at nearest neighbors

      integer  :: n, &
           dst_add,        & ! destination address
           nmap              ! index of current map being computed

      integer , dimension(num_neighbors) :: &
           nbr_add         ! source address at nearest neighbors

      real (DP), dimension(num_neighbors) :: &
           nbr_dist        ! angular distance four nearest neighbors

      real (DP) :: &
           coslat_dst,        & ! cos(lat) of destination grid point
           coslon_dst,        & ! cos(lon) of destination grid point
           sinlat_dst,        & ! sin(lat) of destination grid point
           sinlon_dst,        & ! sin(lon) of destination grid point
           dist_tot             ! sum of neighbor distances (for normalizing)

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1

      !***
      !*** allocate wgtstmp to be consistent with store_link interface
      !***

      allocate (wgtstmp(num_wts))

      !***
      !*** compute cos, sin of lat/lon on source grid for distance
      !*** calculations
      !***

      allocate (coslat(grid1_size), coslon(grid1_size), &
                sinlat(grid1_size), sinlon(grid1_size))

      coslat = cos(grid1_center_lat)
      coslon = cos(grid1_center_lon)
      sinlat = sin(grid1_center_lat)
      sinlon = sin(grid1_center_lon)

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        if (.not. grid2_mask(dst_add)) cycle grid_loop1

        coslat_dst = cos(grid2_center_lat(dst_add))
        coslon_dst = cos(grid2_center_lon(dst_add))
        sinlat_dst = sin(grid2_center_lat(dst_add))
        sinlon_dst = sin(grid2_center_lon(dst_add))

        !***
        !*** find nearest grid points on source grid and
        !*** distances to each point
        !***

        call grid_search_nbr(nbr_add, nbr_dist,         & 
                             grid2_center_lat(dst_add), &
                             grid2_center_lon(dst_add), &
                             coslat_dst, coslon_dst,    &
                             sinlat_dst, sinlon_dst)

        !***
        !*** compute weights based on inverse distance
        !*** if mask is false, eliminate those points
        !***

        dist_tot = zero
        do n=1,num_neighbors
          if (grid1_mask(nbr_add(n))) then
            nbr_dist(n) = one/nbr_dist(n)
            dist_tot = dist_tot + nbr_dist(n)
            nbr_mask(n) = .true.
          else
            nbr_mask(n) = .false.
          endif
        end do

        !***
        !*** normalize weights and store the link
        !***

        do n=1,num_neighbors
          if (nbr_mask(n)) then
            wgtstmp(1) = nbr_dist(n)/dist_tot
            call store_link_nbr(nbr_add(n), dst_add, wgtstmp, nmap)
            grid2_frac(dst_add) = one
          endif
        end do

      end do grid_loop1

      deallocate (coslat, coslon, sinlat, sinlon)

!!-----------------------------------------------------------------------
!!
!!     compute mappings from grid2 to grid1 if necessary
!!
!!-----------------------------------------------------------------------

!      if (num_maps > 1) then
!
!      nmap = 2
!
!      !***
!      !*** compute cos, sin of lat/lon on source grid for distance
!      !*** calculations
!      !***
!
!      allocate (coslat(grid2_size), coslon(grid2_size), &
!                sinlat(grid2_size), sinlon(grid2_size))
!
!      coslat = cos(grid2_center_lat)
!      coslon = cos(grid2_center_lon)
!      sinlat = sin(grid2_center_lat)
!      sinlon = sin(grid2_center_lon)
!
!      !***
!      !*** loop over destination grid 
!      !***
!
!      grid_loop2: do dst_add = 1, grid1_size
!
!        if (.not. grid1_mask(dst_add)) cycle grid_loop2
!
!        coslat_dst = cos(grid1_center_lat(dst_add))
!        coslon_dst = cos(grid1_center_lon(dst_add))
!        sinlat_dst = sin(grid1_center_lat(dst_add))
!        sinlon_dst = sin(grid1_center_lon(dst_add))
!
!        !***
!        !*** find four nearest grid points on source grid and
!        !*** distances to each point
!        !***
!
!        call grid_search_nbr(nbr_add, nbr_dist,         &
!                             grid1_center_lat(dst_add), &
!                             grid1_center_lon(dst_add), &
!                             coslat_dst, coslon_dst,    & 
!                             sinlat_dst, sinlon_dst)
!
!        !***
!        !*** compute weights based on inverse distance
!        !*** if mask is false, eliminate those points
!        !***
!
!        dist_tot = zero
!        do n=1,num_neighbors
!          if (grid2_mask(nbr_add(n))) then
!            nbr_dist(n) = one/nbr_dist(n)
!            dist_tot = dist_tot + nbr_dist(n)
!            nbr_mask(n) = .true.
!          else
!            nbr_mask(n) = .false.
!          endif
!        end do
!
!        !***
!        !*** normalize weights and store the link
!        !***
!
!        do n=1,num_neighbors
!          if (nbr_mask(n)) then
!            wgtstmp(1) = nbr_dist(n)/dist_tot
!            call store_link_nbr(dst_add, nbr_add(n), wgtstmp, nmap)
!            grid1_frac(dst_add) = one
!          endif
!        end do
!
!      end do grid_loop2
!
!      deallocate (coslat, coslon, sinlat, sinlon)
!
!      endif

      deallocate(wgtstmp)

!-----------------------------------------------------------------------

      end subroutine remap_distwgt

!***********************************************************************

      subroutine grid_search_nbr(nbr_add, nbr_dist, plat, plon,      &  
                     coslat_dst, coslon_dst, sinlat_dst, sinlon_dst)

!-----------------------------------------------------------------------
!
!     this routine finds the closest num_neighbor points to a search 
!     point and computes a distance to each of the neighbors.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer , dimension(num_neighbors), intent(out) :: &
             nbr_add  ! address of each of the closest points

      real (DP), dimension(num_neighbors), intent(out) ::  &
             nbr_dist ! distance to each of the closest points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      real (DP), intent(in) :: &
              plat,         & ! latitude  of the search point
              plon,         & ! longitude of the search point
              coslat_dst,   & ! cos(lat)  of the search point
              coslon_dst,   & ! cos(lon)  of the search point
              sinlat_dst,   & ! sin(lat)  of the search point
              sinlon_dst      ! sin(lon)  of the search point

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

!!$      integer  :: n, nmax, nadd, nchk, & ! dummy indices
!!$             min_add, max_add, nm1, np1, i, j, ip1, im1, jp1, jm1
      integer  :: n, nadd, nchk, & ! dummy indices
             min_add, max_add

      real (DP) :: &
              distance      ! angular distance

!-----------------------------------------------------------------------
!
!     loop over source grid and find nearest neighbors
!
!-----------------------------------------------------------------------

      !***
      !*** initialize distance and address arrays
      !***

      min_add = 1
      max_add = grid1_size
      nbr_add = 0
      nbr_dist = bignum

      do nadd=min_add,max_add

        !***
        !*** find distance to this point
        !***

!mz_ap_20081006
! if the source mask is masked... then skip it
        if (.not.grid1_mask(nadd)) CYCLE
! this will ensure that we compute all the points with
! grid2_mask = TRUE
        

        distance = acos(sinlat_dst*sinlat(nadd) + &
                        coslat_dst*coslat(nadd)*  &
                       (coslon_dst*coslon(nadd) + &
                        sinlon_dst*sinlon(nadd)) )

        !***
        !*** store the address and distance if this is one of the
        !*** smallest four so far
        !***

        check_loop: do nchk=1,num_neighbors
          if (distance .lt. nbr_dist(nchk)) then
            do n=num_neighbors,nchk+1,-1
              nbr_add(n) = nbr_add(n-1)
              nbr_dist(n) = nbr_dist(n-1)
            end do
            nbr_add(nchk) = nadd
            nbr_dist(nchk) = distance
            exit check_loop
          endif
        end do check_loop

      end do

!-----------------------------------------------------------------------

      end subroutine grid_search_nbr 

!***********************************************************************

      subroutine store_link_nbr(add1, add2, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for this link in
!     the appropriate address and weight arrays and resizes those
!     arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
              add1,   &! address on grid1
              add2,   &! address on grid2
              nmap     ! identifies which direction for mapping

      real (DP), dimension(:), intent(in) :: &
              weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_map1  = num_links_map1 + 1

        if (num_links_map1 > max_links_map1) & 
           call resize_remap_vars(1,resize_increment)

        grid1_add_map1(num_links_map1) = add1
        grid2_add_map1(num_links_map1) = add2
        wts_map1    (:,num_links_map1) = weights

      case(2)

        num_links_map2  = num_links_map2 + 1

        if (num_links_map2 > max_links_map2) &
           call resize_remap_vars(2,resize_increment)

        grid1_add_map2(num_links_map2) = add1
        grid2_add_map2(num_links_map2) = add2
        wts_map2    (:,num_links_map2) = weights

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_nbr

!***********************************************************************

!***********************************************************************

      subroutine remap_bicub

!-----------------------------------------------------------------------
!
!     this routine computes the weights for a bicubic interpolation.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local parameters
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n,icount, &
           dst_add,     &   ! destination address
           iter,        &   ! iteration counter
           nmap             ! index of current map being computed

      integer , dimension(4) :: & 
           src_add         ! address for the four source points

      real (DP), dimension(4)  :: &
           src_lats,       & ! latitudes  of four bilinear corners
           src_lons          ! longitudes of four bilinear corners

      real (DP), dimension(4,4)  :: &
           wgts            ! bicubic weights for four corners

      real (DP) ::      &
          plat, plon,              & ! lat/lon coords of destination point
          iguess, jguess,          & ! current guess for bilinear coordinate
!!$          thguess, phguess,        & ! current guess for lat/lon coordinate
          deli, delj,              & ! corrections to i,j
          dth1, dth2, dth3,        & ! some latitude  differences
          dph1, dph2, dph3,        & ! some longitude differences
          dthp, dphp,              & ! difference between point and sw corner
          mat1, mat2, mat3, mat4,  & ! matrix elements
          determinant,             & ! matrix determinant
          sum_wgts !!$,                & ! sum of weights for normalization
!!$          w1,w2,w3,w4,w5,w6,w7,w8, & ! 16 bicubic weight functions
!!$          w9,w10,w11,w12,w13,w14,w15,w16

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1
      if (grid1_rank /= 2) then
        stop 'Can not do bicubic interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        if (.not. grid2_mask(dst_add)) cycle grid_loop1

        plat = grid2_center_lat(dst_add)
        plon = grid2_center_lon(dst_add)

!-----------------------------------------------------------------------
!
!       find nearest square of grid points on source grid
!
!-----------------------------------------------------------------------

        call grid_search_bicub(src_add, src_lats, src_lons,         & 
                              plat, plon, grid1_dims,               &
                              grid1_center_lat, grid1_center_lon,   &
                              grid1_bound_box)

        !***
        !*** check to see if points are land points
        !***

        if (src_add(1) > 0) then
          do n=1,4
            if (.not. grid1_mask(src_add(n))) src_add(1) = 0
          end do
        endif

!-----------------------------------------------------------------------
!
!       if point found, find local i,j coordinates for weights
!
!-----------------------------------------------------------------------

        if (src_add(1) > 0) then

          grid2_frac(dst_add) = one

          !***
          !*** iterate to find i,j for bicubic approximation
          !***

          dth1 = src_lats(2) - src_lats(1)
          dth2 = src_lats(4) - src_lats(1)
          dth3 = src_lats(3) - src_lats(2) - dth2

          dph1 = src_lons(2) - src_lons(1)
          dph2 = src_lons(4) - src_lons(1)
          dph3 = src_lons(3) - src_lons(2)

          if (dph1 >  three*pih) dph1 = dph1 - pi2
          if (dph2 >  three*pih) dph2 = dph2 - pi2
          if (dph3 >  three*pih) dph3 = dph3 - pi2
          if (dph1 < -three*pih) dph1 = dph1 + pi2
          if (dph2 < -three*pih) dph2 = dph2 + pi2
          if (dph3 < -three*pih) dph3 = dph3 + pi2

          dph3 = dph3 - dph2

          iguess = half
          jguess = half

          iter_loop1: do iter=1,max_iter

            dthp = plat - src_lats(1) - dth1*iguess -       &
                          dth2*jguess - dth3*iguess*jguess
            dphp = plon - src_lons(1)

            if (dphp >  three*pih) dphp = dphp - pi2
            if (dphp < -three*pih) dphp = dphp + pi2

            dphp = dphp - dph1*iguess - dph2*jguess -       &
                          dph3*iguess*jguess

            mat1 = dth1 + dth3*jguess
            mat2 = dth2 + dth3*iguess
            mat3 = dph1 + dph3*jguess
            mat4 = dph2 + dph3*iguess

            determinant = mat1*mat4 - mat2*mat3

            deli = (dthp*mat4 - mat2*dphp)/determinant
            delj = (mat1*dphp - dthp*mat3)/determinant

            if (abs(deli) < converge .and.                 & 
                abs(delj) < converge) exit iter_loop1

            iguess = iguess + deli
            jguess = jguess + delj

          end do iter_loop1

          if (iter <= max_iter) then

!-----------------------------------------------------------------------
!
!           successfully found i,j - compute weights
!
!-----------------------------------------------------------------------

            wgts(1,1) = (one - jguess**2*(three-two*jguess))* &
                        (one - iguess**2*(three-two*iguess))
            wgts(1,2) = (one - jguess**2*(three-two*jguess))* &
                               iguess**2*(three-two*iguess)
            wgts(1,3) =        jguess**2*(three-two*jguess)*  &
                               iguess**2*(three-two*iguess)
            wgts(1,4) =        jguess**2*(three-two*jguess)*  &
                        (one - iguess**2*(three-two*iguess))
            wgts(2,1) = (one - jguess**2*(three-two*jguess))* &
                               iguess*(iguess-one)**2
            wgts(2,2) = (one - jguess**2*(three-two*jguess))* &
                               iguess**2*(iguess-one)
            wgts(2,3) =        jguess**2*(three-two*jguess)*  &
                               iguess**2*(iguess-one) 
            wgts(2,4) =        jguess**2*(three-two*jguess)*  &
                               iguess*(iguess-one)**2
            wgts(3,1) =        jguess*(jguess-one)**2*        &
                        (one - iguess**2*(three-two*iguess))
            wgts(3,2) =        jguess*(jguess-one)**2*        &
                               iguess**2*(three-two*iguess)
            wgts(3,3) =        jguess**2*(jguess-one)*        &
                               iguess**2*(three-two*iguess)
            wgts(3,4) =        jguess**2*(jguess-one)*        &
                        (one - iguess**2*(three-two*iguess))
            wgts(4,1) =        iguess*(iguess-one)**2*        &
                               jguess*(jguess-one)**2
            wgts(4,2) =        iguess**2*(iguess-one)*        &
                               jguess*(jguess-one)**2
            wgts(4,3) =        iguess**2*(iguess-one)*        &
                               jguess**2*(jguess-one)
            wgts(4,4) =        iguess*(iguess-one)**2*        &
                               jguess**2*(jguess-one)

            call store_link_bicub(dst_add, src_add, wgts, nmap)

          else
            stop 'Iteration for i,j exceed max iteration count'
          endif

!-----------------------------------------------------------------------
!
!       search for bilinear failed - use a distance-weighted
!       average instead (this is typically near the pole)
!
!-----------------------------------------------------------------------

        else if (src_add(1) < 0) then

          src_add = abs(src_add)

          icount = 0
          do n=1,4
            if (grid1_mask(src_add(n))) then
              icount = icount + 1
            else
              src_lats(n) = zero
            endif
          end do

          if (icount > 0) then
            !*** renormalize weights

            sum_wgts = sum(src_lats)
            wgts(1,1) = src_lats(1)/sum_wgts
            wgts(1,2) = src_lats(2)/sum_wgts
            wgts(1,3) = src_lats(3)/sum_wgts
            wgts(1,4) = src_lats(4)/sum_wgts
            wgts(2:4,:) = zero

            grid2_frac(dst_add) = one
            call store_link_bicub(dst_add, src_add, wgts, nmap)
          endif

        endif
      end do grid_loop1

!!-----------------------------------------------------------------------
!!
!!     compute mappings from grid2 to grid1 if necessary
!!
!!-----------------------------------------------------------------------

!      if (num_maps > 1) then
!
!      nmap = 2
!      if (grid2_rank /= 2) then
!        stop 'Can not do bicubic interpolation when grid_rank /= 2'
!      endif
!
!      !***
!      !*** loop over destination grid 
!      !***
!
!      grid_loop2: do dst_add = 1, grid1_size
!
!        if (.not. grid1_mask(dst_add)) cycle grid_loop2
!
!        plat = grid1_center_lat(dst_add)
!        plon = grid1_center_lon(dst_add)
!
!        !***
!        !*** find nearest square of grid points on source grid
!        !***
!
!        call grid_search_bicub(src_add, src_lats, src_lons,         & 
!                               plat, plon, grid2_dims,              &
!                               grid2_center_lat, grid2_center_lon,  &
!                               grid2_bound_box)
!
!        !***
!        !*** check to see if points are land points
!        !***
!
!        if (src_add(1) > 0) then
!          do n=1,4
!            if (.not. grid2_mask(src_add(n))) src_add(1) = 0
!          end do
!        endif
!
!        !***
!        !*** if point found, find i,j coordinates for weights
!        !***
!
!        if (src_add(1) > 0) then
!
!          grid1_frac(dst_add) = one
!
!          !***
!          !*** iterate to find i,j for bilinear approximation
!          !***
!
!          dth1 = src_lats(2) - src_lats(1)
!          dth2 = src_lats(4) - src_lats(1)
!          dth3 = src_lats(3) - src_lats(2) - dth2
!
!          dph1 = src_lons(2) - src_lons(1)
!          dph2 = src_lons(4) - src_lons(1)
!          dph3 = src_lons(3) - src_lons(2)
!
!          if (dph1 >  pi) dph1 = dph1 - pi2
!          if (dph2 >  pi) dph2 = dph2 - pi2
!          if (dph3 >  pi) dph3 = dph3 - pi2
!          if (dph1 < -pi) dph1 = dph1 + pi2
!          if (dph2 < -pi) dph2 = dph2 + pi2
!          if (dph3 < -pi) dph3 = dph3 + pi2
!
!          dph3 = dph3 - dph2
!
!          iguess = zero
!          jguess = zero
!
!          iter_loop2: do iter=1,max_iter
!
!            dthp = plat - src_lats(1) - dth1*iguess -      &
!                          dth2*jguess - dth3*iguess*jguess
!            dphp = plon - src_lons(1)
!
!            if (dphp >  pi) dphp = dphp - pi2
!            if (dphp < -pi) dphp = dphp + pi2
!
!            dphp = dphp - dph1*iguess - dph2*jguess -      & 
!                          dph3*iguess*jguess
!
!            mat1 = dth1 + dth3*jguess
!            mat2 = dth2 + dth3*iguess
!            mat3 = dph1 + dph3*jguess
!            mat4 = dph2 + dph3*iguess
!
!            determinant = mat1*mat4 - mat2*mat3
!
!            deli = (dthp*mat4 - mat2*dphp)/determinant
!            delj = (mat1*dphp - dthp*mat3)/determinant
!
!            if (abs(deli) < converge .and.                 &
!                abs(delj) < converge) exit iter_loop2
!
!            iguess = iguess + deli
!            jguess = jguess + delj
!
!          end do iter_loop2
!
!          if (iter <= max_iter) then
!
!            !***
!            !*** successfully found i,j - compute weights
!            !***
!
!            wgts(1,1) = (one - jguess**2*(three-two*jguess))* &
!                        (one - iguess**2*(three-two*iguess))   
!            wgts(1,2) = (one - jguess**2*(three-two*jguess))* &
!                               iguess**2*(three-two*iguess)    
!            wgts(1,3) =        jguess**2*(three-two*jguess)*  &
!                               iguess**2*(three-two*iguess)    
!            wgts(1,4) =        jguess**2*(three-two*jguess)*  &
!                        (one - iguess**2*(three-two*iguess))   
!            wgts(2,1) = (one - jguess**2*(three-two*jguess))* &
!                               iguess*(iguess-one)**2          
!            wgts(2,2) = (one - jguess**2*(three-two*jguess))* &
!                               iguess**2*(iguess-one)          
!            wgts(2,3) =        jguess**2*(three-two*jguess)*  &
!                               iguess**2*(iguess-one)          
!            wgts(2,4) =        jguess**2*(three-two*jguess)*  &
!                               iguess*(iguess-one)**2          
!            wgts(3,1) =        jguess*(jguess-one)**2*        &
!                        (one - iguess**2*(three-two*iguess))   
!            wgts(3,2) =        jguess*(jguess-one)**2*        &
!                               iguess**2*(three-two*iguess)    
!            wgts(3,3) =        jguess**2*(jguess-one)*        &
!                               iguess**2*(three-two*iguess)    
!            wgts(3,4) =        jguess**2*(jguess-one)*        &
!                        (one - iguess**2*(three-two*iguess))   
!            wgts(4,1) =        iguess*(iguess-one)**2*        &
!                               jguess*(jguess-one)**2          
!            wgts(4,2) =        iguess**2*(iguess-one)*        &
!                               jguess*(jguess-one)**2          
!            wgts(4,3) =        iguess**2*(iguess-one)*        &
!                               jguess**2*(jguess-one)          
!            wgts(4,4) =        iguess*(iguess-one)**2*        &
!                               jguess**2*(jguess-one)
!
!            call store_link_bicub(dst_add, src_add, wgts, nmap)
!
!          else
!            stop 'Iteration for i,j exceed max iteration count'
!          endif
!
!        !***
!        !*** search for bilinear failed - us a distance-weighted
!        !*** average instead
!        !***
!
!        else if (src_add(1) < 0) then
!
!          src_add = abs(src_add)
!
!          icount = 0
!          do n=1,4
!            if (grid2_mask(src_add(n))) then
!              icount = icount + 1
!            else
!              src_lats(n) = zero
!            endif
!          end do
!
!          if (icount > 0) then
!            !*** renormalize weights
!
!            sum_wgts = sum(src_lats)
!            wgts(1,1) = src_lats(1)/sum_wgts
!            wgts(1,2) = src_lats(2)/sum_wgts
!            wgts(1,3) = src_lats(3)/sum_wgts
!            wgts(1,4) = src_lats(4)/sum_wgts
!            wgts(2:4,:) = zero
!
!            grid1_frac(dst_add) = one
!            call store_link_bicub(dst_add, src_add, wgts, nmap)
!          endif
!
!        endif
!      end do grid_loop2
!
!      endif ! nmap=2

!-----------------------------------------------------------------------

      end subroutine remap_bicub

!***********************************************************************

      subroutine grid_search_bicub(src_add, src_lats, src_lons,   &
                                  plat, plon, src_grid_dims,      &
                                  src_center_lat, src_center_lon, &
                                  src_bound_box)

!-----------------------------------------------------------------------
!
!     this routine finds the location of the search point plat, plon
!     in the source grid and returns the corners needed for a bicubic
!     interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer , dimension(4), intent(out) :: &
              src_add  ! address of each corner point enclosing P

      real (DP), dimension(4), intent(out) :: &
              src_lats, & ! latitudes  of the four corner points
              src_lons    ! longitudes of the four corner points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      real (DP), intent(in) :: &
            plat,   & ! latitude  of the search point
            plon      ! longitude of the search point

      integer , dimension(2), intent(in) :: &
             src_grid_dims  ! size of each src grid dimension

      real (DP), dimension(:), intent(in) :: &
              src_center_lat, & ! latitude  of each src grid center 
              src_center_lon    ! longitude of each src grid center

      real (DP), dimension(:,:), intent(in) :: &
              src_bound_box   ! bounding box for src grid search

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n, next_n, srch_add, &  ! dummy indices
          nx, ny,                      &  ! dimensions of src grid
          min_add, max_add,            &  ! addresses for restricting search
          i, j, jp1, ip1, n_add, e_add, ne_add  ! addresses

      real (DP) ::   &! vectors for cross-product check
            vec1_lat, vec1_lon,                                    &
            vec2_lat, vec2_lon, cross_product, cross_product_last, &
            coslat_dst, sinlat_dst, coslon_dst, sinlon_dst,        &
            dist_min, distance ! for computing dist-weighted avg


      src_add = 0

      max_add = size(src_center_lat)
      min_add = 1
 
!-----------------------------------------------------------------------
!
!     now perform a more detailed search 
!
!-----------------------------------------------------------------------

      nx = src_grid_dims(1)
      ny = src_grid_dims(2)

      srch_loop: do srch_add = min_add,max_add

        if (plat <= src_bound_box(2,srch_add) .and.  &
            plat >= src_bound_box(1,srch_add) .and.  &
            plon <= src_bound_box(4,srch_add) .and.  &
            plon >= src_bound_box(3,srch_add)) then

          !***
          !*** we are within bounding box so get really serious
          !***

          !*** find N,S and NE points to this grid point

          j = (srch_add - 1)/nx +1
          i = srch_add - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            ip1 = 1
          endif

          if (j < ny) then
            jp1 = j+1
          else
            jp1 = 1
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          !***
          !*** find N,S and NE lat/lon coords and check bounding box
          !***

          src_lats(1) = src_center_lat(srch_add)
          src_lats(2) = src_center_lat(e_add)
          src_lats(3) = src_center_lat(ne_add)
          src_lats(4) = src_center_lat(n_add)

          src_lons(1) = src_center_lon(srch_add)
          src_lons(2) = src_center_lon(e_add)
          src_lons(3) = src_center_lon(ne_add)
          src_lons(4) = src_center_lon(n_add)

          !***
          !*** for consistency, we must make sure all lons are in
          !*** same 2pi interval
          !***

          vec1_lon = src_lons(1) - plon
          if (vec1_lon > pi) then
            src_lons(1) = src_lons(1) - pi2
          else if (vec1_lon < -pi) then
            src_lons(1) = src_lons(1) + pi2
          endif
          do n=2,4
            vec1_lon = src_lons(n) - src_lons(1)
            if (vec1_lon > pi) then
              src_lons(n) = src_lons(n) - pi2
            else if (vec1_lon < -pi) then
              src_lons(n) = src_lons(n) + pi2
            endif
          end do

          corner_loop: do n=1,4
            next_n = MOD(n,4) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each box side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** same sign, the point is contained in the box.
            !***

            vec1_lat = src_lats(next_n) - src_lats(n)
            vec1_lon = src_lons(next_n) - src_lons(n)
            vec2_lat = plat - src_lats(n)
            vec2_lon = plon - src_lons(n)

            !***
            !*** check for 0,2pi crossings
            !***

            if (vec1_lon >  three*pih) then
              vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -three*pih) then
              vec1_lon = vec1_lon + pi2
            endif
            if (vec2_lon >  three*pih) then
              vec2_lon = vec2_lon - pi2
            else if (vec2_lon < -three*pih) then
              vec2_lon = vec2_lon + pi2
            endif

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (n==1) cross_product_last = cross_product
            if (cross_product*cross_product_last < zero) then
              exit corner_loop
            else
              cross_product_last = cross_product
            endif

          end do corner_loop

          !***
          !*** if cross products all positive, we found the location
          !***

          if (n > 4) then
            src_add(1) = srch_add
            src_add(2) = e_add
            src_add(3) = ne_add
            src_add(4) = n_add

            return
          endif

          !***
          !*** otherwise move on to next cell
          !***

        endif !bounding box check
      end do srch_loop

      !***
      !*** if no cell found, point is likely either in a box that
      !*** straddles either pole or is outside the grid.  fall back
      !*** to a distance-weighted average of the four closest
      !*** points.  go ahead and compute weights here, but store
      !*** in src_lats and return -add to prevent the parent
      !*** routine from computing bilinear weights
      !***

      coslat_dst = cos(plat)
      sinlat_dst = sin(plat)
      coslon_dst = cos(plon)
      sinlon_dst = sin(plon)

      dist_min = bignum
      src_lats = bignum
      do srch_add = min_add,max_add
        distance = acos(coslat_dst*cos(src_center_lat(srch_add))*  &
                       (coslon_dst*cos(src_center_lon(srch_add)) + &
                        sinlon_dst*sin(src_center_lon(srch_add)))+ &
                        sinlat_dst*sin(src_center_lat(srch_add)))

        if (distance < dist_min) then
          sort_loop: do n=1,4
            if (distance < src_lats(n)) then
              do i=4,n+1,-1
                src_add (i) = src_add (i-1)
                src_lats(i) = src_lats(i-1)
              end do
              src_add (n) = -srch_add
              src_lats(n) = distance
              dist_min = src_lats(4)
              exit sort_loop
            endif
          end do sort_loop
        endif
      end do

      src_lons = one/(src_lats + tiny)
      distance = sum(src_lons)
      src_lats = src_lons/distance

!-----------------------------------------------------------------------

      end subroutine grid_search_bicub 

!***********************************************************************

      subroutine store_link_bicub(dst_add, src_add, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for four links 
!     associated with one destination point in the appropriate address 
!     and weight arrays and resizes those arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
             dst_add,         & ! address on destination grid
             nmap               ! identifies which direction for mapping

      integer , dimension(4), intent(in) :: &
             src_add   ! addresses on source grid

      real (DP), dimension(4,4), intent(in) :: &
             weights ! array of remapping weights for these links

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: n, & ! dummy index
             num_links_old            ! placeholder for old link number

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_old  = num_links_map1
        num_links_map1 = num_links_old + 4

        if (num_links_map1 > max_links_map1)               & 
           call resize_remap_vars(1,resize_increment)

        do n=1,4
          grid1_add_map1(num_links_old+n) = src_add(n)
          grid2_add_map1(num_links_old+n) = dst_add
          wts_map1    (:,num_links_old+n) = weights(:,n)
        end do

      case(2)

        num_links_old  = num_links_map2
        num_links_map2 = num_links_old + 4

        if (num_links_map2 > max_links_map2)               & 
           call resize_remap_vars(2,resize_increment)

        do n=1,4
          grid1_add_map2(num_links_old+n) = dst_add
          grid2_add_map2(num_links_old+n) = src_add(n)
          wts_map2    (:,num_links_old+n) = weights(:,n)
        end do

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_bicub

!***********************************************************************
!***********************************************************************

      subroutine remap_bilin

!-----------------------------------------------------------------------
!
!     this routine computes the weights for a bilinear interpolation.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     local parameters
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: n,icount,  & 
                dst_add,    &  ! destination address
                iter,       &  ! iteration counter
                nmap           ! index of current map being computed

      integer, dimension(4) :: src_add ! address for the four source points

      real (DP), dimension(4)  ::  & 
                 src_lats,         &  ! latitudes  of four bilinear corners
                 src_lons,         &  ! longitudes of four bilinear corners
                 wgts                 ! bilinear weights for four corners

      real (DP) ::                 & 
           plat, plon,             &  ! lat/lon coords of destination point
           iguess, jguess,         &  ! current guess for bilinear coordinate
!!$           thguess, phguess,       &  ! current guess for lat/lon coordinate
           deli, delj,             &  ! corrections to i,j
           dth1, dth2, dth3,       &  ! some latitude  differences
           dph1, dph2, dph3,       &  ! some longitude differences
           dthp, dphp,             &  ! difference between point and sw corner
           mat1, mat2, mat3, mat4, &  ! matrix elements
           determinant,            &  ! matrix determinant
           sum_wgts                   ! sum of weights for normalization

!-----------------------------------------------------------------------
!
!     compute mappings from grid1 to grid2
!
!-----------------------------------------------------------------------

      nmap = 1
      if (grid1_rank /= 2) then
        stop 'Can not do bilinear interpolation when grid_rank /= 2'
      endif

      !***
      !*** loop over destination grid 
      !***

      grid_loop1: do dst_add = 1, grid2_size

        if (.not. grid2_mask(dst_add)) cycle grid_loop1

        plat = grid2_center_lat(dst_add)
        plon = grid2_center_lon(dst_add)

        !***
        !*** find nearest square of grid points on source grid
        !***

        call grid_search_bilin(src_add, src_lats, src_lons,            & 
                               plat, plon, grid1_dims,                 &
                               grid1_center_lat, grid1_center_lon,     &
                               grid1_bound_box)

        !***
        !*** check to see if points are land points
        !***

        if (src_add(1) > 0) then
          do n=1,4
            if (.not. grid1_mask(src_add(n))) src_add(1) = 0
          end do
        endif

        !***
        !*** if point found, find local i,j coordinates for weights
        !***

        if (src_add(1) > 0) then

          grid2_frac(dst_add) = one

          !***
          !*** iterate to find i,j for bilinear approximation
          !***

          dth1 = src_lats(2) - src_lats(1)
          dth2 = src_lats(4) - src_lats(1)
          dth3 = src_lats(3) - src_lats(2) - dth2

          dph1 = src_lons(2) - src_lons(1)
          dph2 = src_lons(4) - src_lons(1)
          dph3 = src_lons(3) - src_lons(2)

          if (dph1 >  three*pih) dph1 = dph1 - pi2
          if (dph2 >  three*pih) dph2 = dph2 - pi2
          if (dph3 >  three*pih) dph3 = dph3 - pi2
          if (dph1 < -three*pih) dph1 = dph1 + pi2
          if (dph2 < -three*pih) dph2 = dph2 + pi2
          if (dph3 < -three*pih) dph3 = dph3 + pi2

          dph3 = dph3 - dph2

          iguess = half
          jguess = half

          iter_loop1: do iter=1,max_iter

            dthp = plat - src_lats(1) - dth1*iguess -      &
                          dth2*jguess - dth3*iguess*jguess
            dphp = plon - src_lons(1)

            if (dphp >  three*pih) dphp = dphp - pi2
            if (dphp < -three*pih) dphp = dphp + pi2

            dphp = dphp - dph1*iguess - dph2*jguess -      & 
                          dph3*iguess*jguess

            mat1 = dth1 + dth3*jguess
            mat2 = dth2 + dth3*iguess
            mat3 = dph1 + dph3*jguess
            mat4 = dph2 + dph3*iguess

            determinant = mat1*mat4 - mat2*mat3

            deli = (dthp*mat4 - mat2*dphp)/determinant
            delj = (mat1*dphp - dthp*mat3)/determinant

            if (abs(deli) < converge .and.                 &
                abs(delj) < converge) exit iter_loop1

            iguess = iguess + deli
            jguess = jguess + delj

          end do iter_loop1

          if (iter <= max_iter) then

            !***
            !*** successfully found i,j - compute weights
            !***

            wgts(1) = (one-iguess)*(one-jguess)
            wgts(2) = iguess*(one-jguess)
            wgts(3) = iguess*jguess
            wgts(4) = (one-iguess)*jguess

            call store_link_bilin(dst_add, src_add, wgts, nmap)

          else
            print *,'Point coords: ',plat,plon
            print *,'Dest grid lats: ',src_lats
            print *,'Dest grid lons: ',src_lons
            print *,'Dest grid addresses: ',src_add
            print *,'Current i,j : ',iguess, jguess
            stop 'Iteration for i,j exceed max iteration count'
          endif

        !***
        !*** search for bilinear failed - use a distance-weighted
        !*** average instead (this is typically near the pole)
        !***

        else if (src_add(1) < 0) then

          src_add = abs(src_add)
          icount = 0
          do n=1,4
            if (grid1_mask(src_add(n))) then
              icount = icount + 1
            else
              src_lats(n) = zero
            endif
          end do

          if (icount > 0) then
            !*** renormalize weights

            sum_wgts = sum(src_lats)
            wgts(1) = src_lats(1)/sum_wgts
            wgts(2) = src_lats(2)/sum_wgts
            wgts(3) = src_lats(3)/sum_wgts
            wgts(4) = src_lats(4)/sum_wgts

            grid2_frac(dst_add) = one
            call store_link_bilin(dst_add, src_add, wgts, nmap)
          endif

        endif
      end do grid_loop1

!!-----------------------------------------------------------------------
!!
!!     compute mappings from grid2 to grid1 if necessary
!!
!!-----------------------------------------------------------------------
!
!      if (num_maps > 1) then
!
!      nmap = 2
!      if (grid2_rank /= 2) then
!        stop 'Can not do bilinear interpolation when grid_rank /= 2'
!      endif
!
!      !***
!      !*** loop over destination grid 
!      !***
!
!      grid_loop2: do dst_add = 1, grid1_size
!
!        if (.not. grid1_mask(dst_add)) cycle grid_loop2
!
!        plat = grid1_center_lat(dst_add)
!        plon = grid1_center_lon(dst_add)
!
!        !***
!        !*** find nearest square of grid points on source grid
!        !***
!
!        call grid_search_bilin(src_add, src_lats, src_lons,         &
!                               plat, plon, grid2_dims,              &
!                               grid2_center_lat, grid2_center_lon,  &
!                               grid2_bound_box)
!
!        !***
!        !*** check to see if points are land points
!        !***
!
!        if (src_add(1) > 0) then
!          do n=1,4
!            if (.not. grid2_mask(src_add(n))) src_add(1) = 0
!          end do
!        endif
!
!        !***
!        !*** if point found, find i,j coordinates for weights
!        !***
!
!        if (src_add(1) > 0) then
!
!          grid1_frac(dst_add) = one
!
!          !***
!          !*** iterate to find i,j for bilinear approximation
!          !***
!
!          dth1 = src_lats(2) - src_lats(1)
!          dth2 = src_lats(4) - src_lats(1)
!          dth3 = src_lats(3) - src_lats(2) - dth2
!
!          dph1 = src_lons(2) - src_lons(1)
!          dph2 = src_lons(4) - src_lons(1)
!          dph3 = src_lons(3) - src_lons(2)
!
!          if (dph1 >  pi) dph1 = dph1 - pi2
!          if (dph2 >  pi) dph2 = dph2 - pi2
!          if (dph3 >  pi) dph3 = dph3 - pi2
!          if (dph1 < -pi) dph1 = dph1 + pi2
!          if (dph2 < -pi) dph2 = dph2 + pi2
!          if (dph3 < -pi) dph3 = dph3 + pi2
!
!          dph3 = dph3 - dph2
!
!          iguess = zero
!          jguess = zero
!
!          iter_loop2: do iter=1,max_iter
!
!            dthp = plat - src_lats(1) - dth1*iguess -      &
!                          dth2*jguess - dth3*iguess*jguess
!            dphp = plon - src_lons(1)
!
!            if (dphp >  pi) dphp = dphp - pi2
!            if (dphp < -pi) dphp = dphp + pi2
!
!            dphp = dphp - dph1*iguess - dph2*jguess -      & 
!                          dph3*iguess*jguess
!
!            mat1 = dth1 + dth3*jguess
!            mat2 = dth2 + dth3*iguess
!            mat3 = dph1 + dph3*jguess
!            mat4 = dph2 + dph3*iguess
!
!            determinant = mat1*mat4 - mat2*mat3
!
!            deli = (dthp*mat4 - mat2*dphp)/determinant
!            delj = (mat1*dphp - dthp*mat3)/determinant
!
!            if (abs(deli) < converge .and. & 
!               abs(delj) < converge) exit iter_loop2
!
!            iguess = iguess + deli
!            jguess = jguess + delj
!
!          end do iter_loop2
!
!          if (iter <= max_iter) then
!
!            !***
!            !*** successfully found i,j - compute weights
!            !***
!
!            wgts(1) = (one-iguess)*(one-jguess)
!            wgts(2) = iguess*(one-jguess)
!            wgts(3) = iguess*jguess
!            wgts(4) = (one-iguess)*jguess
!
!            call store_link_bilin(dst_add, src_add, wgts, nmap)
!
!          else
!            print *,'Point coords: ',plat,plon
!            print *,'Dest grid lats: ',src_lats
!            print *,'Dest grid lons: ',src_lons
!            print *,'Dest grid addresses: ',src_add
!            print *,'Current i,j : ',iguess, jguess
!            stop 'Iteration for i,j exceed max iteration count'
!          endif
!
!        !***
!        !*** search for bilinear failed - us a distance-weighted
!        !*** average instead
!        !***
!
!        else if (src_add(1) < 0) then
!
!          src_add = abs(src_add)
!          icount = 0
!          do n=1,4
!            if (grid2_mask(src_add(n))) then
!              icount = icount + 1
!            else
!              src_lats(n) = zero
!            endif
!          end do
!
!          if (icount > 0) then
!            !*** renormalize weights
!
!            sum_wgts = sum(src_lats)
!            wgts(1) = src_lats(1)/sum_wgts
!            wgts(2) = src_lats(2)/sum_wgts
!            wgts(3) = src_lats(3)/sum_wgts
!            wgts(4) = src_lats(4)/sum_wgts
!
!            grid1_frac(dst_add) = one
!            call store_link_bilin(dst_add, src_add, wgts, nmap)
!          endif
!
!        endif
!      end do grid_loop2
!
!      endif ! nmap=2
!
!-----------------------------------------------------------------------

      end subroutine remap_bilin

!***********************************************************************

      subroutine grid_search_bilin(src_add, src_lats, src_lons,     &
                                   plat, plon, src_grid_dims,       &
                                   src_center_lat, src_center_lon,  &
                                   src_grid_bound_box)

!-----------------------------------------------------------------------
!
!     this routine finds the location of the search point plat, plon
!     in the source grid and returns the corners needed for a bilinear
!     interpolation.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     output variables
!
!-----------------------------------------------------------------------

      integer , dimension(4), intent(out) :: &
              src_add  ! address of each corner point enclosing P

      real (DP), dimension(4), intent(out) :: &
              src_lats, & ! latitudes  of the four corner points
              src_lons    ! longitudes of the four corner points

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      real (DP), intent(in) :: &
              plat,   & ! latitude  of the search point
              plon      ! longitude of the search point

      integer , dimension(2), intent(in) :: &
             src_grid_dims  ! size of each src grid dimension

      real (DP), dimension(:), intent(in) :: &
              src_center_lat, & ! latitude  of each src grid center 
              src_center_lon    ! longitude of each src grid center

      real (DP), dimension(:,:), intent(in) :: &
              src_grid_bound_box ! bound box for source grid

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n, next_n, srch_add,   & ! dummy indices
          nx, ny,                        & ! dimensions of src grid
          min_add, max_add,              & ! addresses for restricting search
          i, j, jp1, ip1, n_add, e_add, ne_add  ! addresses

      real (DP) ::             & ! vectors for cross-product check
            vec1_lat, vec1_lon,                                    &
            vec2_lat, vec2_lon, cross_product, cross_product_last, &
            coslat_dst, sinlat_dst, coslon_dst, sinlon_dst,        &
            dist_min, distance ! for computing dist-weighted avg


      src_add = 0

      max_add = size(src_center_lat)
      min_add = 1
 
!-----------------------------------------------------------------------
!
!     now perform a more detailed search 
!
!-----------------------------------------------------------------------

      nx = src_grid_dims(1)
      ny = src_grid_dims(2)

      srch_loop: do srch_add = min_add,max_add

        !*** first check bounding box

        if (plat <= src_grid_bound_box(2,srch_add) .and.  & 
            plat >= src_grid_bound_box(1,srch_add) .and.  &
            plon <= src_grid_bound_box(4,srch_add) .and.  & 
            plon >= src_grid_bound_box(3,srch_add)) then

          !***
          !*** we are within bounding box so get really serious
          !***

          !*** determine neighbor addresses

          j = (srch_add - 1)/nx +1
          i = srch_add - (j-1)*nx

          if (i < nx) then
            ip1 = i + 1
          else
            ip1 = 1
          endif

          if (j < ny) then
            jp1 = j+1
          else
            jp1 = 1
          endif

          n_add = (jp1 - 1)*nx + i
          e_add = (j - 1)*nx + ip1
          ne_add = (jp1 - 1)*nx + ip1

          src_lats(1) = src_center_lat(srch_add)
          src_lats(2) = src_center_lat(e_add)
          src_lats(3) = src_center_lat(ne_add)
          src_lats(4) = src_center_lat(n_add)

          src_lons(1) = src_center_lon(srch_add)
          src_lons(2) = src_center_lon(e_add)
          src_lons(3) = src_center_lon(ne_add)
          src_lons(4) = src_center_lon(n_add)

          !***
          !*** for consistency, we must make sure all lons are in
          !*** same 2pi interval
          !***

          vec1_lon = src_lons(1) - plon
          if (vec1_lon >  pi) then
            src_lons(1) = src_lons(1) - pi2
          else if (vec1_lon < -pi) then
            src_lons(1) = src_lons(1) + pi2
          endif
          do n=2,4
            vec1_lon = src_lons(n) - src_lons(1)
            if (vec1_lon >  pi) then
              src_lons(n) = src_lons(n) - pi2
            else if (vec1_lon < -pi) then
              src_lons(n) = src_lons(n) + pi2
            endif
          end do

          corner_loop: do n=1,4
            next_n = MOD(n,4) + 1

            !***
            !*** here we take the cross product of the vector making 
            !*** up each box side with the vector formed by the vertex
            !*** and search point.  if all the cross products are 
            !*** positive, the point is contained in the box.
            !***

            vec1_lat = src_lats(next_n) - src_lats(n)
            vec1_lon = src_lons(next_n) - src_lons(n)
            vec2_lat = plat - src_lats(n)
            vec2_lon = plon - src_lons(n)

            !***
            !*** check for 0,2pi crossings
            !***

            if (vec1_lon >  three*pih) then
              vec1_lon = vec1_lon - pi2
            else if (vec1_lon < -three*pih) then
              vec1_lon = vec1_lon + pi2
            endif
            if (vec2_lon >  three*pih) then
              vec2_lon = vec2_lon - pi2
            else if (vec2_lon < -three*pih) then
              vec2_lon = vec2_lon + pi2
            endif

            cross_product = vec1_lon*vec2_lat - vec2_lon*vec1_lat

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***

            if (n == 1) cross_product_last = cross_product
            if (cross_product*cross_product_last < zero)   & 
                exit corner_loop
            cross_product_last = cross_product

          end do corner_loop

          !***
          !*** if cross products all same sign, we found the location
          !***

          if (n > 4) then
            src_add(1) = srch_add
            src_add(2) = e_add
            src_add(3) = ne_add
            src_add(4) = n_add

            return
          endif

          !***
          !*** otherwise move on to next cell
          !***

        endif !bounding box check
      end do srch_loop

      !***
      !*** if no cell found, point is likely either in a box that
      !*** straddles either pole or is outside the grid.  fall back
      !*** to a distance-weighted average of the four closest
      !*** points.  go ahead and compute weights here, but store
      !*** in src_lats and return -add to prevent the parent
      !*** routine from computing bilinear weights
      !***

      !print *,'Could not find location for ',plat,plon
      !print *,'Using nearest-neighbor average for this point'

      coslat_dst = cos(plat)
      sinlat_dst = sin(plat)
      coslon_dst = cos(plon)
      sinlon_dst = sin(plon)

      dist_min = bignum
      src_lats = bignum
      do srch_add = min_add,max_add
        distance = acos(coslat_dst*cos(src_center_lat(srch_add))*  &
                       (coslon_dst*cos(src_center_lon(srch_add)) + &
                        sinlon_dst*sin(src_center_lon(srch_add)))+ &
                        sinlat_dst*sin(src_center_lat(srch_add)))

        if (distance < dist_min) then
          sort_loop: do n=1,4
            if (distance < src_lats(n)) then
              do i=4,n+1,-1
                src_add (i) = src_add (i-1)
                src_lats(i) = src_lats(i-1)
              end do
              src_add (n) = -srch_add
              src_lats(n) = distance
              dist_min = src_lats(4)
              exit sort_loop
            endif
          end do sort_loop
        endif
      end do

      src_lons = one/(src_lats + tiny)
      distance = sum(src_lons)
      src_lats = src_lons/distance

!-----------------------------------------------------------------------

      end subroutine grid_search_bilin 

!***********************************************************************

      subroutine store_link_bilin(dst_add, src_add, weights, nmap)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for four links 
!     associated with one destination point in the appropriate address 
!     and weight arrays and resizes those arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer , intent(in) :: &
             dst_add,  & ! address on destination grid
             nmap        ! identifies which direction for mapping

      integer , dimension(4), intent(in) :: &
             src_add   ! addresses on source grid

      real (DP), dimension(4), intent(in) :: &
             weights ! array of remapping weights for these links

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer  :: n,  &     ! dummy index
             num_links_old  ! placeholder for old link number

!-----------------------------------------------------------------------
!
!     increment number of links and check to see if remap arrays need
!     to be increased to accomodate the new link.  then store the
!     link.
!
!-----------------------------------------------------------------------

      select case (nmap)
      case(1)

        num_links_old  = num_links_map1
        num_links_map1 = num_links_old + 4

        if (num_links_map1 > max_links_map1) & 
          call resize_remap_vars(1,resize_increment)

        do n=1,4
          grid1_add_map1(num_links_old+n) = src_add(n)
          grid2_add_map1(num_links_old+n) = dst_add
          wts_map1    (1,num_links_old+n) = weights(n)
        end do

      case(2)

        num_links_old  = num_links_map2
        num_links_map2 = num_links_old + 4

        if (num_links_map2 > max_links_map2) & 
         call resize_remap_vars(2,resize_increment)

        do n=1,4
          grid1_add_map2(num_links_old+n) = dst_add
          grid2_add_map2(num_links_old+n) = src_add(n)
          wts_map2    (1,num_links_old+n) = weights(n)
        end do

      end select

!-----------------------------------------------------------------------

      end subroutine store_link_bilin

!***********************************************************************

  SUBROUTINE remap_dealloc

      DEALLOCATE(grid1_add_map1, grid2_add_map1, wts_map1)

      DEALLOCATE( grid1_frac     ,  &
                 grid2_frac      ,  &
                 grid1_area      ,  &
                 grid2_area      ,  &
                 grid1_bound_box ,  & 
                 grid2_bound_box )

      DEALLOCATE(grid1_dims,grid2_dims)

  END SUBROUTINE remap_dealloc



! ***************************************************************************
END MODULE messy_main_grid_trafo_scrp_base
! ***************************************************************************
