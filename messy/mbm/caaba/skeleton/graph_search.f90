!*****************************************************************************

! Dijkstra's algorithm, as described by Niemeyer and Sung, Combust.
! Flame 158, 1439–1443 (2011), doi: 10.1016/j.combustflame.2010.12.010

! Author: Kyle Niemeyer

! This file is available under the MIT License (MIT)
 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
 
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

!*****************************************************************************

module graph_search
	implicit none
	
	private
	public dijkstra, dijkstra_adj
	
	integer, parameter :: wp = kind(1.0d0)
	
contains

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
	
subroutine dijkstra (nnode, adj, targ, RAB)

	implicit none

	integer, intent(in) :: nnode, targ
	real(wp), intent(in) :: adj(:, :)

	real(wp), intent(out) :: RAB(:)
	
	integer :: i, u, v, nQ
	integer, dimension(nnode) :: Q
	real(wp) :: Rtmp
!!!!

	RAB = 0.0d0
	
	do i = 1, nnode
		Q(i) = i
	end do
	nQ = nnode
	RAB(targ) = 1.0d0

	search_loop: do while (Q(1) > 0)
		u = maxQ(nQ, nnode, Q, RAB)
		
!		if (RAB(u) < 1.0d-80) exit search_loop
		
		call removeQ (nQ, nnode, Q, u)
		nQ = nQ - 1
		
		do i = 1, nQ
			v = Q(i)
			if (adj(u,v) > 1.0d-80) then
				Rtmp = RAB(u) * adj(u,v)
				RAB(v) = DMAX1( RAB(v), Rtmp)
			end if
		end do
		
	end do search_loop

	return
end subroutine dijkstra

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

subroutine dijkstra_adj ( nnode, adj, neigh, n_neigh, targ, RAB )
! Dijkstra's algorithm using adjacency list

	implicit none

	integer, intent(in) :: nnode, targ, neigh(nnode, nnode), n_neigh(nnode)
	real(wp), intent(in) :: adj(nnode, nnode)

	real(wp), intent(out) :: RAB(nnode)
!
	integer :: i, u, v, nQ
	integer, dimension(nnode) :: Q
	real(wp) :: Rtmp

	RAB = 0.0d0
	
	do i = 1, nnode
		Q(i) = i
	end do
	nQ = nnode
	RAB(targ) = 1.0d0

	search_loop: do while (Q(1) /= 0)
		u = maxQ(nQ, nnode, Q, RAB)
!		if (RAB(u) < 1.0d-80) exit search_loop
		call removeQ (nQ, nnode, Q, u)
		nQ = nQ - 1
		
		do i = 1, n_neigh(u)
			v = neigh(u, i)
			if (adj(u,v) > 1.0d-80) then
				Rtmp = RAB(u) * adj(u,v)
				RAB(v) = DMAX1( RAB(v), Rtmp)
			end if
		end do
		
	end do search_loop

	return
end subroutine dijkstra_adj

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

subroutine removeQ (n, nnode, Q, u)
! remove element u from Q
	implicit none
	
	integer, intent(in) :: n, nnode, u
	integer, intent(inout) :: Q(nnode)
	
	integer :: i
	logical :: flag
!!!
	flag = .false.
	
	do i = 1, n - 1
		if ( (Q(i) == u) .or. (flag) ) then
			Q(i) = Q(i + 1)
			flag = .true.
		end if	
	end do
	Q(n) = 0
	
	return
end subroutine removeQ

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

integer function maxQ(n, nnode, Q, RAB)
! find element of Q with maximum RAB
	implicit none
	
	integer, intent(in) :: n, nnode, Q(nnode)
	real(wp), intent(in) :: RAB(nnode)
	
	integer :: i, tmp
	real(wp) :: maxtmp
!!!!!
	
	maxtmp = 0.0d0
	
	do i = 1, n
		if ( Q(i) /= 0 ) then
			if (RAB(Q(i)) >= maxtmp) then
				tmp = Q(i)
				maxtmp = RAB(tmp)
			end if
		end if
	end do
	
	maxQ = tmp

end function maxQ

end module graph_search
