      SUBROUTINE SUCHIJ(abrei,alaij,k,ipos,jpos,dist,wetref)
!
!     subroutine to determine i and j indices of nearest point
!     from lat and long
!
! RJ: i and j are GLOBAL indices!!!
!
!     input:
!     abrei      latitude (deg.)
!     alaij      longitude (deg.)
!     k          level (integer)
!     wetref     if land shall be considered too : 0.
!                ocean points only                 1.
!
!     output:
!     ipos       i index
!     jpos       j index
!     dist       distance from point in m
!
      use mo_kind
      USE MO_PARAM1, ONLY: ie_g, je_g
      USE MO_PARALLEL, ONLY: p_bcast, p_pe, p_io, gather
      USE MO_COMMO1, ONLY: alon_g, alat_g, weto
!      USE MO_UNITS
      USE mo_constants, ONLY: agratorad
      IMPLICIT NONE

      REAL(wp) ::  abrei,alaij,dist,wetref
      INTEGER k,ipos,jpos
      REAL(wp), ALLOCATABLE :: w_g(:,:)

      REAL(wp) :: aphi,alam,xx,yy,zz,ax,ay,az,d,alo,ala
      INTEGER i,j
!
      if (p_pe == p_io) then
      allocate(w_g(ie_g,je_g))
      else
      allocate(w_g(0,0))

      endif

      call gather(weto(:,:,k),w_g,p_io)

      if (p_pe == p_io) then

!
      aphi=agratorad*abrei
      alam=agratorad*alaij
!
      xx=cos(alam)*cos(aphi)
      yy=sin(alam)*cos(aphi)
      zz=sin(aphi)
!
      dist = 1.e20
      ipos = 0
      jpos = 0
!
!CDIR NOLOOPCHG
      do j=2,je_g-1
!CDIR NOVECTOR
       do i=2,ie_g-1
        if(w_g(i,j).gt.wetref-0.1)then
          alo=agratorad*alon_g(i,j)
          ala=agratorad*alat_g(i,j)
          ax=cos(ala)*cos(alo)
          ay=cos(ala)*sin(alo)
          az=sin(ala)
          d=(xx-ax)**2+(yy-ay)**2+(zz-az)**2
          if (d<dist) then
            dist = d
            ipos = i
            jpos = j
          endif
        endif
       enddo
      enddo

      dist=sqrt(dist)*6350000.

      write(6,*)'in suchij end: ',ipos,jpos,k,w_g(ipos,jpos)

      endif

      call p_bcast(ipos,p_io)
      call p_bcast(jpos,p_io)

      deallocate(w_g)

      return
      end
