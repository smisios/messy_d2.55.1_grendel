      SUBROUTINE SEDSHI(kpie,kpje)
!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/sedshi.f90,v $\\
!$Revision: 1.2.2.1.20.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!
!**********************************************************************
!
!**** *SEDSHI* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename ssssil(i,j,k)=sedlay(i,j,k,issssil) etc.
!     I. Kriest         *MPI-Met, HH*,   27.05.03
!     - change specific weights for opal, CaCO3, POC --> bodensed.f90
!     - include upward transport
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     .
!
!**   Interface.
!     ----------
!
!     *CALL*       *SEDSHI*
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      use mo_param1_bgc

implicit none

      INTEGER :: kpie,kpje,i,j,k,iv
      REAL(wp) :: wsed(kpie,kpje)                  ! shifting velocity for upward/downward shifts
      REAL(wp) :: fulsed(kpie,kpje)
      REAL(wp) :: sedlo,uebers
      REAL(wp) :: seddef                           ! sediment deficiency
      REAL(wp) :: spresent, buried
      REAL(wp) :: refill,frac

! DOWNWARD SHIFTING
! shift solid sediment downwards, if layer is full, i.e., if
! the volume filled by the four constituents poc, opal, caco3, clay
! is more than porsol*seddw
! the outflow of layer i is given by sedlay(i)*porsol(i)*seddw(i), it is
! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)

!e    write(0,121)orgfa,calfa,oplfa,clafa
!121   format('sediments',4e12.4)

      do k=1,ks-1

        do j=1,kpje
        do i=1,kpie
          if(bolay(i,j).gt.0._wp) then
             sedlo = orgfa*rcar*sedlay(i,j,k,issso12)                  &
     &              +calfa*sedlay(i,j,k,isssc12)                       &
     &              +oplfa*sedlay(i,j,k,issssil)                       &
     &              +clafa*sedlay(i,j,k,issster)
! "full" sediment has sedlo=1. for sedlo>1., wsed is >0.
             wsed(i, j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))          ! downward shifting velocity (?)
          endif
       enddo !end i-loop
       enddo !end j-loop

! filling downward  (accumulation)
        do iv=1,nsedtra
        do j=1,kpje
        do i=1,kpie
          if(bolay(i,j).gt.0._wp) then
            uebers=wsed(i,j)*sedlay(i,j,k,iv)                     ! 'uebersaettigung?'
            sedlay(i,j,k  ,iv)=sedlay(i,j,k  ,iv)-uebers
            sedlay(i,j,k+1,iv)=sedlay(i,j,k+1,iv)+uebers               &
     &        *(seddw(k)*porsol(k))/(seddw(k+1)*porsol(k+1))
          endif
        enddo !end i-loop
        enddo !end j-loop
        enddo !end iv-loop

      enddo !end k-loop


! store amount lost from bottom sediment layer - this is a kind of
! permanent burial in deep consolidated layer, and this stuff is
! effectively lost from the whole ocean+sediment(+atmosphere) system.
! Would have to be supplied by river runoff or simple addition e.g.
! to surface layers in the long range. Can be supplied again if a
! sediment column has a deficiency in volume.

       do j=1,kpje
       do i=1,kpie
          if(bolay(i,j).gt.0._wp) then
             sedlo = orgfa*rcar*sedlay(i,j,ks,issso12)                 &
     &              +calfa*sedlay(i,j,ks,isssc12)                      &
     &              +oplfa*sedlay(i,j,ks,issssil)                      &
     &              +clafa*sedlay(i,j,ks,issster)
             wsed(i, j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
          endif
      enddo !end i-loop
      enddo !end j-loop

      do iv=1,nsedtra
      do j=1,kpje
      do i=1,kpie
          if(bolay(i,j).gt.0._wp) then
            uebers=wsed(i,j)*sedlay(i,j,k,iv)
            sedlay(i,j,ks ,iv)=sedlay(i,j,ks ,iv)-uebers
            burial(i,j,iv)=burial(i,j,iv)+uebers*seddw(k)*porsol(k)
          endif
      enddo !end i-loop
      enddo !end j-loop
      enddo !end iv-loop

      return       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< switch off upward shift


! now the loading nowhere exceeds 1.

! digging from below in case of erosion (js and initialization with 0 or none, as currently)
! UPWARD SHIFTING
! shift solid sediment upwards, if total sediment volume is less
! than required, i.e., if the volume filled by the four constituents
! poc, opal, caco3, clay (integrated over the total sediment column)
! is less than porsol*seddw (integrated over the total sediment column)
! first, the deepest box is filled from below with total required volume;
! then, successively, the following layers are filled upwards.
! if there is not enough solid matter to fill the column, add clay. (js: so implicite initial state is all clay)

      do j=1,kpje
      do i=1,kpie
        fulsed(i,j)=0._wp
      enddo !end i-loop
      enddo !end j-loop

! determine how the total sediment column is filled
      do k=1,ks
      do j=1,kpje
      do i=1,kpie
        if(bolay(i,j).gt.0._wp) then
          sedlo=orgfa*rcar*sedlay(i,j,k,issso12)                       &
     &         +calfa*sedlay(i,j,k,isssc12)                            &
     &         +oplfa*sedlay(i,j,k,issssil)                            &
     &         +clafa*sedlay(i,j,k,issster)
          fulsed(i,j)=fulsed(i,j)+porsol(k)*seddw(k)*sedlo
        endif
      enddo !end i-loop
      enddo !end j-loop
      enddo !end k-loop

! shift the sediment deficiency from the deepest (burial)
! layer into layer ks
      do j=1,kpje
      do i=1,kpie
      if(bolay(i,j).gt.0._wp) then

! deficiency with respect to fully loaded sediment |packed in sedlay(i,j,ks) ??
! this is the volume of sediment shifted upwards from the burial layer

        seddef=solfu-fulsed(i,j)    ! 'sediment deficiency', solfu = total column inegrated solid fraction volume (bodensed)

! total volume of solid constituents in buried layer
        spresent=orgfa*rcar*burial(i,j,issso12)                        &
     &         +calfa*burial(i,j,isssc12)                              &
     &         +oplfa*burial(i,j,issssil)                              &
     &         +clafa*burial(i,j,issster)

! determine whether an additional amount of clay is needed from the burial
! layer to fill the whole sediment; I assume that there is an infinite
! supply of clay from below
        burial(i,j,issster) = burial(i,j,issster)                      &
     &                      + MAX(0._wp, seddef - spresent) / clafa

! determine new volume of buried layer
        buried=orgfa*rcar*burial(i,j,issso12)                          &
     &        +calfa*burial(i,j,isssc12)                               &
     &        +oplfa*burial(i,j,issssil)                               &
     &        +clafa*burial(i,j,issster)

! fill the deepest active sediment layer
        refill=seddef/buried
        frac = porsol(ks)*seddw(ks) !changed k to ks, ik

        sedlay(i,j,ks,issso12)=sedlay(i,j,ks,issso12)                  &
     &                        +refill*burial(i,j,issso12)/frac
        sedlay(i,j,ks,isssc12)=sedlay(i,j,ks,isssc12)                  &
     &                        +refill*burial(i,j,isssc12)/frac
        sedlay(i,j,ks,issssil)=sedlay(i,j,ks,issssil)                  &
     &                        +refill*burial(i,j,issssil)/frac
        sedlay(i,j,ks,issster)=sedlay(i,j,ks,issster)                  &
     &                        +refill*burial(i,j,issster)/frac

! account for losses in buried sediment
        burial(i,j,issso12) = burial(i,j,issso12)                      &
     &                      - refill*burial(i,j,issso12)
        burial(i,j,isssc12) = burial(i,j,isssc12)                      &
     &                      - refill*burial(i,j,isssc12)
        burial(i,j,issssil) = burial(i,j,issssil)                      &
     &                      - refill*burial(i,j,issssil)
        burial(i,j,issster) = burial(i,j,issster)                      &
     &                      - refill*burial(i,j,issster)

      endif ! bolay >0
      enddo !end i-loop
      enddo !end j-loop

!     redistribute overload of deepest layer ks to layers 2 to ks
      do  k=ks,2,-1
      do j=1,kpje
      do i=1,kpie
        if(bolay(i,j).gt.0._wp) then
          sedlo=orgfa*rcar*sedlay(i,j,k,issso12)                       &
     &         +calfa*sedlay(i,j,k,isssc12)                            &
     &         +oplfa*sedlay(i,j,k,issssil)                            &
     &         +clafa*sedlay(i,j,k,issster)
          wsed(i, j) = MAX(0._wp, (sedlo - 1._wp) / (sedlo + 1.e-10_wp))
        endif
      enddo !end i-loop
      enddo !end j-loop

      do iv=1,4
      do j=1,kpje
      do i=1,kpie
        if(bolay(i,j).gt.0._wp) then
          uebers=sedlay(i,j,k,iv)*wsed(i,j)
          frac=porsol(k)*seddw(k)/(porsol(k-1)*seddw(k-1))
          sedlay(i,j,k,iv)=sedlay(i,j,k,iv)-uebers
          sedlay(i,j,k-1,iv)=sedlay(i,j,k-1,iv)+uebers*frac               ! note k-1 here =upward shift
        endif
      enddo !end i-loop
      enddo !end j-loop
      enddo !end iv-loop

      enddo  !end k-loop


      RETURN
      END
