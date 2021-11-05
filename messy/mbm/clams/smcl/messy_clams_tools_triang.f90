!******************************************************************************
! File utils/src/lib_triangulate.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! P. Konopka, N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Mon Nov 25 13:39:23 2019
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version. This program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more
! details. You should have received a copy of the GNU General Public
! License along with this program; if not, see <https://www.gnu.org/licenses/>.
!
!******************************************************************************
!
! Module lib_triangulate
!
!   Subroutine  triang_qhull  (vert_coor, vert_nb, vert_ind, triangles, ntriang)
!   logical function out_of_triangle (t, i, x)
!   integer function get_triangle_index (itriang, i, vert_nb, vert_ind, triangles)
!   integer function find_triangle (t_act,vert_coor,vert_nb,vert_ind,triangles, &
!                                   x, max_dist, write_warnings)
!
!******************************************************************************

Module messy_clams_tools_triang

CONTAINS

  !****************************************************************************
  ! PURPOSE: Carry out the Delaunay triangulation (using qhull) of an ensamble
  !          of air parcels (vert_s) and determine only those triangulation
  !          properties which are important for the interpolation procedure
  !          clams2grid:
  !          Indicies of Delaunay triangles saved in triangles
  !          Number of adjecant triangles for each AP and the indices are 
  !          saved  in vert_nb (number) and vert_ind (indices)
  !****************************************************************************
  
  Subroutine  triang_qhull  (status, vert_coor, vert_nb, vert_ind, triangles, ntriang)

    use messy_clams_global, only: prec

    implicit none

    real(prec), dimension(:,:), pointer :: vert_coor
    integer,    dimension(:),   pointer :: vert_nb
    integer,    dimension(:,:), pointer :: vert_ind
    integer,    dimension(:,:), pointer :: triangles
    integer                             :: ntriang
    integer                             :: status

!!!!! points muss REAL bleiben wg. QHULL !?!
    real,    dimension(:), pointer :: in_points
    integer, dimension(:), pointer :: out_points 

    integer :: nparts, nnb, i, itr, inb, k
    logical :: found

    status = 0 ! no error

    nparts = size(vert_nb)

    ! allocate the input-field for qhull (in_points) containing:
    ! nparts     - # of nodes (APs)
    ! dimension  - space-dimension for the determination of the convex hull 
    !              3 = 2d triangulation (on a plane or on a unit sphere)
    !              4 = 3d triangulation in the 3d-space
    ! coor       - cartesian coordinates of the APs
    allocate (in_points(3*nparts+2))
    in_points(1) = real(nparts)
    in_points(2) = 3.
    
    do i = 1, nparts
       in_points(3*i:3*i+2) = vert_coor(1:3,i)
    end do

    ! allocate output field for qhull (out_points)
    allocate (out_points(3*3*nparts+1))
    out_points = 0

    ! run qhull
    call qhull_dll (in_points, out_points)

    ! check if the allocated memory was large enough 
    ! to read completely the qhull output
    ntriang=out_points(1)
    if (ntriang >= 3*nparts) then 
       write(*,*) 'mem_alloc_size for qhull-call has been to small (change it) !!!'
       status = -1
       return
    endif
    
    ! save output in a field containing only the indicies 
    ! of the Voronoi triangles (tetrahedras), i.e. triangles
    allocate (triangles(3,ntriang))
    do i = 1, ntriang
       triangles(1:3,i) = out_points(3*i-1:3*i+1) + 1 ! C to Fortran
    end do

    deallocate (in_points)
    deallocate (out_points)

    vert_nb(:) = 0

    ! determine neighbouring triangles for each APs
    do itr = 1, ntriang
       do k = 1, 3 
          found = .false.
          inb = 1
          do while (inb <= vert_nb(triangles(k,itr)) .and. .not. found) 
             if (vert_ind(inb,triangles(k,itr)) == itr) found = .true.
             inb = inb+1
          end do
          if (.not. found) then
             nnb = vert_nb(triangles(k,itr))
             vert_nb(triangles(k,itr)) = nnb+1
             vert_ind(nnb+1,triangles(k,itr)) = itr
          end if
       end do
    end do

    !write (*,*) 'Triangulation is ready !'
    !write (*,*) 'Number of triangles: ', ntriang

  end Subroutine triang_qhull

  !****************************************************************************
  ! PURPOSE decide if you are outside of the triangle t (out_of_triangle = 1) 
  !         out_of_triangle is set to 0 if x and the triangle point facing
  !         side i are on the same side
  !   Input: t - a given triangle
  !          i - considered triangle side
  !          x - point where decision, left or right of i should be made
  !  Output: out_of_triangle =
  !          0 - x and i-point are on different sides 
  !          1 - x and i-point are on the same sides
  !****************************************************************************

  logical function out_of_triangle (t, i, x)

    USE messy_clams_global, ONLY: prec

    implicit none

    real(prec), dimension(3,3) :: t
    real(prec), dimension(3)   :: x
    integer                    :: i
    
    real(PREC)               :: hh
    real(PREC), dimension(3) :: a, b, c, n
    real(PREC), parameter    :: t_tol=1E-20  ! tolerance for determination of the triangles

    a = t(MOD(i,3)+1,:)
    b = t(MOD(i+1,3)+1,:)
    c = t(i,:)      
    n = (/ a(2)*b(3)-a(3)*b(2),  a(3)*b(1)-a(1)*b(3),  a(1)*b(2)-a(2)*b(1) /)   

    hh = ((a(1)-c(1))*n(1) + (a(2)-c(2))*n(2) + (a(3)-c(3))*n(3)) *  &
         ((a(1)-x(1))*n(1) + (a(2)-x(2))*n(2) + (a(3)-x(3))*n(3))

    out_of_triangle = .false.
    if (hh < -t_tol)  out_of_triangle = .true.

  end function out_of_triangle

  !****************************************************************************
  ! PURPOSE determine the index 
  ! of the adjecant triangle to the side i of the given triangle itriang
  !   Input: itriang - considered triangle
  !          i - considered side 
  !          vert_s - list of adjecant triangles
  !          triangles - Delaunay triangles
  !  Output: get_triangle_index - index of the adjecant triangle
  !****************************************************************************

  integer function get_triangle_index (itriang, i, vert_nb, vert_ind, triangles, ntriang)

    implicit none

    integer                          :: itriang, i
    integer, dimension(:),   pointer :: vert_nb
    integer, dimension(:,:), pointer :: vert_ind
    integer, dimension(:,:), pointer :: triangles
    integer                          :: ntriang

    integer :: va, vb, cc, count

    if ((itriang > ntriang) .OR. (i > 3)) then  
       write (*,*) 'Wrong input in get_triangle_index'
       stop
    endif

    va = triangles(MOD(i,3)+1,itriang)
    vb = triangles(MOD(i+1,3)+1,itriang)
    
    do count = 1, vert_nb(va) 
       cc = vert_ind(count,va)
       if ((va == triangles(1,cc) .OR. &
            va == triangles(2,cc) .OR. &
            va == triangles(3,cc)) .AND. &
           (vb == triangles(1,cc) .OR. &
            vb == triangles(2,cc) .OR. &
            vb == triangles(3,cc)) .AND. cc /= itriang)  exit
    end do

    get_triangle_index = cc

  end function get_triangle_index

  !****************************************************************************
  ! PURPOSE: The walking triangle algorithm (Sambridge et al. 1995)
  !   Input: t_act -  first triangle for the search
  !          vert_s -  vertices of the triangles
  !          triangles - Delaunay triangles
  !          x - point where the triangle t with x in t has to be found
  !  Output: find_triangle - triangle containing x
  !****************************************************************************

  integer function find_triangle (t_act, vert_coor, vert_nb, vert_ind, &
                                triangles, ntriang, x, max_dist, write_warnings)

    USE messy_clams_global, ONLY: prec

    implicit none

    integer                             :: t_act
    real(prec), dimension(:,:), pointer :: vert_coor
    integer,    dimension(:),   pointer :: vert_nb
    integer,    dimension(:,:), pointer :: vert_ind
    integer,    dimension(:,:), pointer :: triangles
    integer                             :: ntriang
    real(prec), dimension(3)            :: x
    real(prec), intent(in), optional    :: max_dist
    logical,    intent(in), optional    :: write_warnings
  
    real(prec)                 :: r_earth, dist_max
    real(prec), dimension(3)   :: dist
    real(prec), dimension(3,3) :: t
    integer                    :: i, i_in, iside, t_last, t_last2, kside, nsteps
    logical                    :: found, ctrl_out

    real(PREC), parameter :: eps = 1E-6

    real(PREC) :: t_r    ! criterium for the triangle containing a given AP, r_max < t_r

    if (present(write_warnings)) then
       ctrl_out = write_warnings
    else
       ctrl_out = .true.
    endif

    if (present(max_dist)) then
       t_r = max_dist
    else
       t_r = 1000.
    endif

    r_earth = 6378.1690  ! in km

    found = .false.
    nsteps = 0
    t_last=-1
    t_last2=-1

    do while (.not. found .and. nsteps<ntriang)

       nsteps = nsteps+1

       t(1,:) = vert_coor(:,triangles(1,t_act))
       t(2,:) = vert_coor(:,triangles(2,t_act))
       t(3,:) = vert_coor(:,triangles(3,t_act))

       ! Ueberpruefe, ob x einem der Eckpunkte der akt. Triangel entspricht
       if ((ABS(t(1,1)-x(1))<eps .AND. ABS(t(1,2)-x(2))<eps .AND. ABS(t(1,3)-x(3))<eps) .OR. &
            (ABS(t(2,1)-x(1))<eps .AND. ABS(t(2,2)-x(2))<eps .AND. ABS(t(2,3)-x(3))<eps) .OR. &
            (ABS(t(3,1)-x(1))<eps .AND. ABS(t(3,2)-x(2))<eps .AND. ABS(t(3,3)-x(3))<eps)) then
          !write (*,*) 'Eckpunkt !!!'
          found = .true.
       else

          i_in = 0
          do iside=1, 3 
             if (out_of_triangle(t, iside, x)) then  
                t_last2 = t_last
                t_last=t_act
                t_act = get_triangle_index(t_act,iside, vert_nb, vert_ind, triangles, ntriang)
                
                ! Achtung: Ruecksprung in zuvor ueberpruefte Triangel !!!
                ! Ueberpruefe die anderen Seiten der letzten Triangel 
                if (t_act == t_last2) then
           
                   if (ctrl_out) then
                      write(*,*) 'Warnung triang...t=',t
                      write(*,*) 'Warnung triang...x=',x
                      write(*,*) 'Warnung triang...iside=',iside
                      write(*,*) 't_last2, t_last, t_act=',t_last2,t_last,t_act
                   endif

                   ! iside ist die zuletzt ueberpruefte Seite in t_last (nicht in t_act) !!!
                   ! in t_last wurden bereits i_in=iside Seiten ueberprueft !
                   i_in = iside
                   t(1,:) = vert_coor(:,triangles(1,t_last))
                   t(2,:) = vert_coor(:,triangles(2,t_last))
                   t(3,:) = vert_coor(:,triangles(3,t_last))
                   do kside=iside+1,3
                      if (out_of_triangle(t,kside,x)) then
                         t_act = get_triangle_index(t_last,kside,vert_nb, vert_ind,triangles, ntriang)
                         if (ctrl_out) write(*,*) 'neues t_act=',t_act
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
          
    end do

    do i=1, 3  
       dist(i) = sqrt((t(i,1)-x(1))**2 + (t(i,2)-x(2))**2 + (t(i,3)-x(3))**2)
    end do
    dist_max = r_earth*minval(dist)
    ! if max_dist>=20000 specified: do not check triangle
    if (t_r >= 20000) then
       find_triangle = t_act
    else
       if (dist_max > t_r) then  
          if (ctrl_out) then
             print *, 'Wrong triangle, dist_max larger than', t_r, 'km'
             print *, 'Max dist: ', dist_max
          endif
          find_triangle = -1
       else
          find_triangle = t_act
       endif
    endif

  end function find_triangle


End Module messy_clams_tools_triang
