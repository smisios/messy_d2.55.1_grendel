!********************************************************************************!
! File clams/dynmod/source/lib_triang.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Juergen Ankenbrand, Paul Konopka
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Fri Nov 30 11:12:39 2012
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
!***********************************************************************************
!
!   Library  lib_triang  
!
!    Functions for triangulation, getting neighbours and Voronoi areas for all APs  
!                                                                                   
!    Functions:                                                                     
!                                                                                   
!       subroutine qhull(ap_s,nparts,dim,switch)                                    
!                                                                                   
!    For further informations, see inline commentars.                             
!                                                                                   
!***********************************************************************************

Module messy_clamsmix_lib_triang

contains


! Subroutine to carry out the qhull-triangulation (based on the convex hull algorithm)
! Input for qhull-call is written to file 'qhull_in.dat', output from the qhull
! call to 'qhull_out.dat' (old version)
! Now a C-interface is used
! after the triangulation, all neighbours of APs and the corresonding
! Voronoi areas (optional) are determined
! These results are written into the ap_s-structure

subroutine qhull(status,ap_s,nparts,dim,switch)

  use messy_clams_global, only: prec
  use messy_clamsmix_global, only: triangle, max_nb
  use messy_clamsmix_ap_m_access 
 
  implicit none

  integer, intent(inout)                    :: status
  type(ap_struct)                           :: ap_s
  integer, intent(in)                       :: dim
  integer,intent(in)                        :: nparts
  character                                 :: switch

  integer                                   :: i,j,k,m, dim_tr
  type(triangle), dimension(:), allocatable :: triangle_s
  real(prec),    dimension(:),  pointer     :: area    ! Voronoi area of air parcels

  integer,dimension(:),allocatable          ::zaehler
  integer,dimension(:,:),allocatable        ::hilfsfeld
  logical                                   ::found

!!!!! points muss REAL bleiben wg. QHULL !?!
  real,dimension(:),pointer                 :: points
  real(prec)                                :: voronoi
  integer,dimension(:),pointer              :: outpoints
  integer                                   :: mem_alloc_size
  external qhull_dll


  status = 0 ! no error

  ! set the dimension of the convex hull triangulation
  ap_s%dimension=dim
  if (ctrl_out) print *, ap_s%dimension , '-D qhull, Triang. for: ', n, &
                            'points with (a)ctual /(o)ld NN-relations: ', switch

  ! allocate the input-field for qhull (points) containing:
  ! nparts     - # of nodes (APs)
  ! dimension  - space-dimension for the determination of the convex hull 
  !              3 = 2d triangulation (on a plane or on a unit sphere)
  !              4 = 3d triangulation in the 3d-space
  ! coor       - cartesian coordinates of the APs
  allocate(points(3*nparts+2))
  points(1)=real(nparts)
  points(2)=real(dim)
  
  ! dependent on "switch" use actuel (a) or old (o) coordinates of the APs
  ! a - after the advection step, i.e. before the deformation
  ! o - before the advection step, i.e. after the deformation
  ! v - same as "a", used for the determination of the Voronoi areas
  select case (switch)
  case ('a') ! actual
     do i=1,nparts
        points(3*i:3*i+2)=real(coor(ap_s,i))
     end do
  case ('o') ! old
     do i=1,nparts
        points(3*i:3*i+2)=real(coor_old(ap_s,i))
     end do
  case ('v') ! actual, use for calculation of the Voronoi areas
     do i=1,nparts
        points(3*i:3*i+2)=real(coor(ap_s,i))
     end do
  case default
     print *, 'Wrong choice for switch-parameter'
     status = -1
     return
  end select
  
  ! Allocate the output field for qhull (outpoints)
  ! note: different memory size for different qhull configurations
  if (dim ==3) then 
     mem_alloc_size=3*dim*nparts
  else
     mem_alloc_size=3*dim*nparts*4
  end if
  allocate(outpoints(mem_alloc_size))
  outpoints=0

  ! run qhull
  ! for the used options see qhull_dll in unix.c in the qhull-directory
  call qhull_dll(points,outpoints)

  ! check if the allocated memory was large enough 
  ! to read completely the qhull output
  dim_tr=outpoints(1)
  if (dim_tr*dim+1 >= mem_alloc_size) then 
     write(*,*) 'mem_alloc_size for qhull-call has been to small (change it)'
     status = -1
     return
  endif

  ! save output in a field containing only the indicies 
  ! of the Voronoi triangles (tetrahedras), i.e. triangle_s
  allocate(triangle_s(dim_tr))
  do i=1, dim_tr
     triangle_s(i)%ind(1:dim)=outpoints(dim*(i-1)+2:dim*i+1) +1  !C to Fortran
  enddo
  
  deallocate(points)
  deallocate(outpoints)

  !allocate the array to store the neighbours
  allocate(hilfsfeld(nparts,max_nb))
  !allocate the array for counting the neighbours of every APs
  allocate(zaehler(nparts))
  zaehler=0
  hilfsfeld=0

  ! allocate array for calculation of the Voronoi areas 
  if (switch == 'v') then 
    allocate(area(nparts))
    area=0.0
  endif

  ! derive from the Voronoi-triangles (tetrahedras) the NNs of every APs  
  do i=1,dim_tr  !for all triangles
     ! in the NN-list (1, 2, ..., j, m,..., dimension) of ith traingle
     ! j and m are NNs. Save these NN in hilfsfeld and count the total number
     ! in zaehler
     do j=1,dim-1  
        do m=j+1,dim
           ! found=.false. mean that this NN can be put on the list 
           found=.false. 
           ! Check first if this NN is still present in the list
           ! if yes, set found=.true.
           k=1
           do while (k<=zaehler(triangle_s(i)%ind(j)) .and. .not. found)
              if (hilfsfeld(triangle_s(i)%ind(j),k)==triangle_s(i)%ind(m)) found=.true. 
              k=k+1
           end do
           if (.not. found) then 
              ! this NN can be put on the list, check first if the array list
              ! is not too small
              if ((zaehler(triangle_s(i)%ind(j))+1 >max_nb) .or.&
                   & (zaehler(triangle_s(i)%ind(m))+1 >max_nb)) then
                 ! in this case the NNs are not taken into account 
                 ! and are not put on the list of the NNs
                 ! Important for the boundary APs
                 ! print *, '# of nbs > max_nb'
                 ! status = -1
                 ! return
              else
              ! this NN is correct, put it on the list, increase the counter
              ! and calculate the Voronoi area if necessary
                 if (switch == 'v') then 
                    voronoi=area_triangle(coor(ap_s,triangle_s(i)%ind(1)), &
                                          coor(ap_s,triangle_s(i)%ind(2)), & 
                                          coor(ap_s,triangle_s(i)%ind(3)))/3.
                    area(triangle_s(i)%ind(j))=area(triangle_s(i)%ind(j))+voronoi
                    area(triangle_s(i)%ind(m))=area(triangle_s(i)%ind(m))+voronoi
                 endif
                 zaehler(triangle_s(i)%ind(j))=zaehler(triangle_s(i)%ind(j))+1
                 zaehler(triangle_s(i)%ind(m))=zaehler(triangle_s(i)%ind(m))+1
                 hilfsfeld(triangle_s(i)%ind(j),zaehler(triangle_s(i)%ind(j)))=triangle_s(i)%ind(m) ! A neighbour B
                 hilfsfeld(triangle_s(i)%ind(m),zaehler(triangle_s(i)%ind(m)))=triangle_s(i)%ind(j) ! B neighbour A
              endif
           endif
        end do
     end do
  end do

  if (switch == 'v') print *, 'Total Voronoi area/4pi= ', sum(area)/(4.*3.14159)

  ! copy the NN-informations into the ap_s-structure
  do i=1,nparts  
     call nb(ap_s,i,zaehler(i))
     call ind(ap_s,i,hilfsfeld(i,:))
  end do

  ! deallocate the memory
  deallocate(hilfsfeld,zaehler)
  deallocate(triangle_s)
!!!!! ???
  if (switch == 'v') then 
     deallocate(area)
  endif

  contains

! -----------------------------------------------------------------

  function area_triangle(a,b,c)
    ! area of a triangle described in terms of 3 3d-vectors
    ! see Bronstein, page 608, eq 4.3   

    real(prec),dimension(3), intent(in)                    :: a, b, c
    real(prec)                                             :: area_triangle

    real(prec),dimension(3)                                :: h1, h2

    h1=b-a
    h2=c-a
    area_triangle=0.5*sqrt((h1(2)*h2(3)-h1(3)*h2(2))**2+ &
                           (h1(3)*h2(1)-h1(1)*h2(3))**2+ &
                           (h1(1)*h2(2)-h1(2)*h2(1))**2)

  end function area_triangle

! -----------------------------------------------------------------

end subroutine qhull

end Module messy_clamsmix_lib_triang



