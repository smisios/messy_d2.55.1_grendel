!******************************************************************************
!
! Modul clamssedi_triang
!
! Written by Paul Konopka, Nicole Thomas, Thomas Breuer
! IEK-7, Forschungszentrum Juelich
! Last Modified By: N.Thomas
! Last Modified On: Tue Apr 28 11:09:02 2015
!
! subroutine triangulation(status)
! subroutine get_triangles (ilev, nairparcels, max_ntriang, ntriang, status)
! subroutine get_triangles_for_particle(ipart, status)
! subroutine set_particle_kart (lon, lat, particle_kart)
! integer FUNCTION find_triang (t_act, ilev, x, ind_next, write_warnings) 
! logical FUNCTION out_of_triang (t, i, x)
! integer FUNCTION nn_index (itriang, ilev, iside) 
!
!******************************************************************************

Module messy_clamssedi_triang

CONTAINS
  
  !****************************************************************************
  ! Triangulation for a given time
  !****************************************************************************
  subroutine triangulation(status)
    
    use messy_clamssedi_global,     only: nlevs, nairparcels_lev, ntriang_lev, &
                                          triangles, max_ntriang_lev
    
    implicit none
    
    integer :: status
    integer :: ilev
        
    allocate (triangles%airparcel_indices(max_ntriang_lev, nlevs, 3))
    allocate (triangles%nparticles_nat(max_ntriang_lev, nlevs))
    allocate (triangles%nparticles_ice(max_ntriang_lev, nlevs))
!!!!! Momentan: Jeder Prozessor bearbeitet alle Level
!!!!! Ausblick: verteilte Arbeit, in dem Level auf Prozessoren verteilt werden
!!!!!           anschliessend: Daten einsammeln, zB durch Merken des Level oder airparcel Index
   
    !do ilev = rank+1, nlevs, ntasks
    do ilev = 1, nlevs
       call get_triangles (ilev, nairparcels_lev(ilev), max_ntriang_lev, ntriang_lev(ilev), status)
       if (status/=0) return
    enddo

    
  end subroutine triangulation
     
  !****************************************************************************
  ! Carry out the Delaunay triangulation (using qhull) of an ensamble
  ! of air parcels for a given level (ilev).
  ! Indicies of Delaunay triangles saved in "triangles".
  ! Number of adjecant triangles for each AP and the indicies are
  ! saved  in airparcels%ntriang(:) and airparcels%triang_ids(:,:).
  !****************************************************************************
  subroutine get_triangles (ilev, nairparcels, max_ntriang, ntriang, status)

    use messy_clamssedi_global,   only: airparcels, triangles, max_nb
    use messy_clams_global,       only: nparts, rank

    implicit none

    integer, intent(in)  :: ilev, nairparcels, max_ntriang
    integer, intent(out) :: ntriang
    integer              :: status

    integer, dimension(:), allocatable :: airparcel_ids

    real,    dimension(:), pointer :: in_points
    integer, dimension(:), pointer :: out_points 

    integer :: ipart, count, i, itr, ap
    integer :: nnb, inb
    logical :: found

    ! indices of the airparcels for the given level (ilev) in the global airparcels array
    allocate (airparcel_ids(nairparcels))

    ! allocate the input-field for qhull (in_points) containing:
    ! nairparcels - # of nodes (APs)
    ! dimension   - space-dimension for the determination of the convex hull
    !               3 = 2d triangulation (on a plane or on a unit sphere)
    !               4 = 3d triangulation in the 3d-space
    ! coor        - cartesian coordinates of the APs
    allocate (in_points(3 * nairparcels + 2))
    in_points(1) = real(nairparcels)
    in_points(2) = 3.

    count = 1
    do ipart = 1, nparts
       if (airparcels%ilev(ipart) == ilev) then
          airparcel_ids(count) = ipart     
          in_points(3*count:3*count+2) = airparcels%coor(:,ipart)
          count = count + 1
       endif
    end do

    ! allocate output field for qhull (out_points)
    allocate (out_points(3*3*nairparcels+1))
    out_points = 0

    ! run qhull
    call qhull_dll (in_points, out_points)

    ! check if the allocated memory was large enough
    ! to read completely the qhull output
    ntriang = out_points(1)
    if (ntriang >= max_ntriang) then
       write(*,*) 'mem_alloc_size for qhull-call has been to small (change it)'
       status = 1
       return
    endif

    ! save output in a field containing only the indicies
    ! of the Voronoi triangles (tetrahedras), i.e. triangle_s
    do i = 1, ntriang
       ! Indices in global airparcels array 
       triangles%airparcel_indices(i,ilev,:) = airparcel_ids(out_points(3*i-1:3*i+1) + 1) ! C to Fortran
       ! set initial values (previously in function get_triangles_for_tboxes)
       triangles%nparticles_nat(i,ilev) = 0
       triangles%nparticles_ice(i,ilev) = 0
    end do

    deallocate (in_points)
    deallocate (out_points)
    deallocate (airparcel_ids)

    ! determine neighbouring triangles for each APs
    do itr = 1, ntriang
       do ap = 1, 3
          found = .false.
          inb = 1
          do while (inb <= airparcels%ntriang(triangles%airparcel_indices(itr,ilev,ap)) .and. &
               .not. found) 
             ! check if current triangle index (itr) is already existing in 
             ! the list of neighbours of the current airparcel (ap)
             if (airparcels%triang_ids(inb,triangles%airparcel_indices(itr,ilev,ap)) == itr) & 
                  found = .true.
             inb = inb + 1
          enddo
          if (.not. found) then
             nnb = airparcels%ntriang(triangles%airparcel_indices(itr,ilev,ap))
             if (nnb < max_nb) then 
                ! add current triangle index to the list of neighbours of the current airparcel
                airparcels%ntriang(triangles%airparcel_indices(itr,ilev,ap)) = nnb + 1
                airparcels%triang_ids(nnb+1,triangles%airparcel_indices(itr,ilev,ap)) = itr
             else
                write (*,*) 'number of neighbours >',max_nb,'!!!'
                status = 1
                return
             endif
          endif
       enddo
    enddo

  end subroutine get_triangles

  !****************************************************************************
  ! - For each particle: get indices of triangles and levels containing 
  !                   this particle on lower and upper level
  ! - For each triangle: get number of particles within triangle
  !****************************************************************************
  subroutine get_triangles_for_particle(ipart, ipart_local, status)
    
    use messy_clamssedi_global,   only: particles, triangles, nlevs, &
                                        lev_grid, ntriang_lev
    use messy_clams_global,       only: prec, nparts, rank


    implicit none

    integer, intent(in) :: ipart, ipart_local
    real(prec), dimension(3) :: particle_kart
    integer :: ilev, lev_ind_down, lev_ind_up
    integer :: t_act
    integer :: status
    
    do ilev = 1, nlevs+1
       if (lev_grid(ilev) <= particles%lev(ipart) .and.  &
            particles%lev(ipart) < lev_grid(ilev+1)) then
          
          ! determine the indices of the levels below and above the particle
          if (ilev == 1) then
             lev_ind_down = 1
             lev_ind_up = 1
          elseif (ilev == nlevs+1) then
             lev_ind_down = nlevs
             lev_ind_up = nlevs
          else
             lev_ind_down = ilev - 1
             lev_ind_up = ilev
          end if
          
          particles%lev_down(ipart_local) = lev_ind_down
          particles%lev_up(ipart_local) = lev_ind_up
          
          call set_particle_kart(particles%lon(ipart), particles%lat(ipart), particle_kart)
          !---------------------------------------------------------------------------------
          ! determine the indices of the triangles on the lower resp. upper 
          ! level surrounding the current particle
          ! -> particles%tr_ind_down(ipart_local), particles%tr_ind_up(ipart_local)
          !---------------------------------------------------------------------------------
          t_act = int((ntriang_lev(lev_ind_down)-1)/2)     
          if (t_act <= 0) then
             write(*,*)'warning get_triangles_for_particle: ',t_act, ilev
             t_act=1
          endif
          particles%tr_ind_down(ipart_local) = find_triang(t_act, lev_ind_down, &
               particle_kart, write_warnings=.false.)
          if (particles%tr_ind_down(ipart_local) == -1) then
             write (*,*) 'in sub. get_triangles_for_particles:'
             write (*,*) 'Triangle on the level below the current particle not found:'
             write (*,*) 'rank, ipart, ilev=',rank,ipart,ilev
             write (*,*) 'particle: lon,lat=',particles%lon(ipart),particles%lat(ipart)
             status = 1 
             return
          endif
             
          t_act = int((ntriang_lev(lev_ind_up)-1)/2)
          particles%tr_ind_up(ipart_local) = find_triang(t_act, lev_ind_up, &
               particle_kart, write_warnings=.false.)
          if (particles%tr_ind_up(ipart_local) == -1) then
             write (*,*) 'in sub. get_triangles_for_particles:'
             write (*,*) 'Triangle on the level above the current particle not found:'
             write (*,*) 'ipart, ilev=',ipart,ilev
             write (*,*) 'particle: lon,lat=',particles%lon(ipart),particles%lat(ipart)
             status = 1
             return
          endif
             
          EXIT ! continue with next particle
       endif
       
    enddo
    
  end subroutine get_triangles_for_particle

  !****************************************************************************
  ! convert lon and lat to Cartesian coordinates on a unit sphere
  !****************************************************************************
  subroutine set_particle_kart (lon, lat, particle_kart)

    use messy_clams_global, only: prec, pi
    
    implicit none

    real(prec), intent(in)                :: lon, lat
    real(prec), dimension(3), intent(out) :: particle_kart

    real  :: lon_rad, lat_rad
    
    !!!! lon and lat are already radiant
    !lon_rad = lon * pi / 180.          ! transform from deg to rad
    !lat_rad = (90. - lat) * pi / 180.  ! use spherical coordinates
                                                    !   0. < lons < 2*PI
                                                    !   0. < lats < PI
    lon_rad = lon
    lat_rad = (pi/2. - lat)
    particle_kart(1) = cos(lon_rad)*sin(lat_rad) ! transf from spheric
    particle_kart(2) = sin(lon_rad)*sin(lat_rad) ! coord. to  kart.
    particle_kart(3) = cos(lat_rad)              ! coord on a unit sphere (r=1)

  end subroutine set_particle_kart

  !****************************************************************************
  ! PURPOSE: The walking triangle algorithm (Sambridge et al. 1995)
  !   Input: t_act -  first triangle for the search
  !          vert_s -  vertices of the triangles
  !          triangle_s - Delaunay triangles
  !          x - point where the triangle t with x in t has to be found
  !  Output: find_triang - triangle containing x
  !****************************************************************************
  integer FUNCTION find_triang (t_act, ilev, x, write_warnings)

    use messy_clams_global, only: radius_earth, prec, eps, rank

    use messy_clamssedi_global, only: max_ntriang_lev, triangles, airparcels
    
    implicit none

    integer                               :: t_act
    integer, intent(in)                   :: ilev
    real(prec), dimension(3), intent(in)  :: x
    logical, intent(in), optional         :: write_warnings

    logical              :: warnings_out, found
    integer              :: nsteps, t_last, t_last2
    integer              :: i, i_in, iside, kside
    integer              :: ind(1)
    real                 :: dist_max
    real, parameter      :: t_r = 1000.  ! criterium for the triangle containing a given AP, r_max < t_r
    real, dimension(3,3) :: t
    real, dimension(3)   :: dist
    
    if (present(write_warnings)) then
       warnings_out = write_warnings
    else
       warnings_out = .true.
    endif
    
    found = .false.
    nsteps = 0
    t_last = -1
    t_last2 = -1
    
    ! write(*,*) 'input t_act=',t_act, ' ilev=',ilev,' x=', x
    ! write(*,*) ''

    do while (.not. found .and. nsteps < max_ntriang_lev)
       
       nsteps = nsteps + 1
       
       t(1,:) = airparcels%coor(:,triangles%airparcel_indices(t_act,ilev,1))
       t(2,:) = airparcels%coor(:,triangles%airparcel_indices(t_act,ilev,2))
       t(3,:) = airparcels%coor(:,triangles%airparcel_indices(t_act,ilev,3))
       
       ! check, if x is vertex of the current triangle
       if ((ABS(t(1,1)-x(1))<eps .AND. ABS(t(1,2)-x(2))<eps .AND. ABS(t(1,3)-x(3))<eps) .OR. &
            (ABS(t(2,1)-x(1))<eps .AND. ABS(t(2,2)-x(2))<eps .AND. ABS(t(2,3)-x(3))<eps) .OR. &
            (ABS(t(3,1)-x(1))<eps .AND. ABS(t(3,2)-x(2))<eps .AND. ABS(t(3,3)-x(3))<eps)) then
          found = .true.
          write(*,*) 'particle is vertex of the current triangle'
       else
          i_in = 0
          do iside = 1, 3
             if (out_of_triang(t, iside, x)) then
                t_last2 = t_last
                t_last = t_act
                t_act = nn_index(t_act, ilev, iside)
                if(t_act == -1) then
                   find_triang = -1
                   return
                endif

                ! write(*,*) 'nsteps,t_last2,t_last,t_act', nsteps,t_last2,t_last,t_act
                ! write(*,*) 't',t
                ! write(*,*) ''
                
                ! Attention: return to previously tested triangle !!!
                ! verify the other sides of the last triangle
                if (t_act == t_last2) then
                   if (warnings_out) then
                      write(*,*) 'Warnung triang...t=',t
                      write(*,*) 'Warnung triang...x=',x
                      write(*,*) 'Warnung triang...iside=',iside
                      write(*,*) 't_last2, t_last, t_act=',t_last2,t_last,t_act
                   endif
                   
                   ! iside is the recently tested side in t_last (not in t_act) !!!
                   ! in t_last i_in=iside side is already tested !
                   i_in = iside
                   t(1,:) = airparcels%coor(:,triangles%airparcel_indices(t_last,ilev,1))
                   t(2,:) = airparcels%coor(:,triangles%airparcel_indices(t_last,ilev,2))
                   t(3,:) = airparcels%coor(:,triangles%airparcel_indices(t_last,ilev,3))
                   do kside = iside+1, 3
                      if (out_of_triang(t, kside, x)) then
                         t_act = nn_index(t_last, ilev, kside)
                         if(t_act == -1) then
                            find_triang = -1
                            return
                         endif
                         if (warnings_out) write(*,*) 'neues t_act=',t_act
                      else
                         i_in = i_in + 1
                      endif
                   enddo
                endif
                exit
             
             else
                i_in = i_in + 1
             endif
             
          end do
          
          if (i_in == 3) found = .true.
          
       endif
       
    enddo
    
    do i=1, 3
       dist(i) = sqrt((t(i,1)-x(1))**2 + (t(i,2)-x(2))**2 + (t(i,3)-x(3))**2)
    end do
    dist_max = radius_earth/1000. * minval(dist)
    
    if (.not. found) then
       print *, 'Triangle not found !!!'
       print *, 'Distance=', radius_earth/1000.*dist
    endif

    if (dist_max > t_r) then
       print *, 'Wrong triangle, dist_max larger than', t_r, 'km'
       print *, 'Max dist: ', dist_max
       print *, 't_act', t_act
       ! find_triang = -1
    endif
    
    find_triang = t_act
   
    
  end function find_triang

  !****************************************************************************
  ! PURPOSE decide if you are outside of the triangle t (out_of_triang = 1)
  !         out_of_triang is set to 0 if x and the triangle point facing
  !         side i are on the same side
  !   Input: t - a given triangle
  !          i - considered triangle side
  !          x - point where decision, left or right of i should be made
  !  Output: out_of_triang =
  !          0 - x and i-point are on different sides
  !          1 - x and i-point are on the same sides
  !****************************************************************************
  logical FUNCTION out_of_triang (t, i, x)
    
    use messy_clams_global, only: prec
    
    implicit none
    
    real, dimension(3,3)      :: t
    real(prec), dimension(3)  :: x
    integer                   :: i
    
    real               :: hh
    real, dimension(3) :: a, b, c, n
    real, parameter    :: t_tol = 1E-20  ! tolerance for determination of the triangles
    
    a = t(MOD(i,3)+1,:)
    b = t(MOD(i+1,3)+1,:)
    c = t(i,:)
    n = (/ a(2)*b(3)-a(3)*b(2),  a(3)*b(1)-a(1)*b(3),  a(1)*b(2)-a(2)*b(1) /)
    
    hh = ((a(1)-c(1))*n(1) + (a(2)-c(2))*n(2) + (a(3)-c(3))*n(3)) *  &
         ((a(1)-x(1))*n(1) + (a(2)-x(2))*n(2) + (a(3)-x(3))*n(3))
    
    out_of_triang = .false.
    if (hh < -t_tol)  out_of_triang = .true.
    
  END FUNCTION out_of_triang

  !****************************************************************************
  ! PURPOSE determine the index
  ! of the adjecant triangle to the side i of the given triangle itriang
  !   Input: itriang - considered triangle
  !          iside   - considered side
  !          vert_s - list of adjecant triangles
  !          triangle_s - Delaunay triangles
  !  Output: nn_index - index of the adjecant triangle
  !****************************************************************************
  integer FUNCTION nn_index (itriang, ilev, iside)

    use messy_clamssedi_global,   only: max_ntriang_lev, triangles, airparcels

    implicit none

    integer, intent(in) :: itriang, ilev, iside

    integer :: status
    integer :: va, vb, count
    integer :: cc = -1

    if ((itriang > max_ntriang_lev) .OR. (iside > 3)) then
       write (*,*) 'Wrong input in nn_index'
       nn_index = -1
       return
    endif

    va = triangles%airparcel_indices(itriang,ilev,MOD(iside,3)+1)
    vb = triangles%airparcel_indices(itriang,ilev,MOD(iside+1,3)+1)

    do count = 1, airparcels%ntriang(va)
       cc = airparcels%triang_ids(count,va)
       if ((va == triangles%airparcel_indices(cc,ilev,1) .OR. &
            va == triangles%airparcel_indices(cc,ilev,2) .OR. &
            va == triangles%airparcel_indices(cc,ilev,3)) .AND. &
           (vb == triangles%airparcel_indices(cc,ilev,1) .OR. &
            vb == triangles%airparcel_indices(cc,ilev,2) .OR. &
            vb == triangles%airparcel_indices(cc,ilev,3)) .AND. cc /= itriang)  exit
    end do

    nn_index = cc

  end FUNCTION nn_index

End Module messy_clamssedi_triang
