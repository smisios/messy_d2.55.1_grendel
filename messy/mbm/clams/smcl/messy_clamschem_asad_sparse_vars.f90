! *****************************COPYRIGHT*******************************
! 
! Copyright (c) 2008, Regents of the University of California 
! All rights reserved. 
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are 
! met: 
! 
!     * Redistributions of source code must retain the above copyright 
!       notice, this list of conditions and the following disclaimer.  
!     * Redistributions in binary form must reproduce the above 
!       copyright notice, this list of conditions and the following 
!       disclaimer in the documentation and/or other materials provided 
!       with the distribution.  
!     * Neither the name of the University of California, Irvine nor the 
!       names of its contributors may be used to endorse or promote 
!       products derived from this software without specific prior 
!       written permission. 
! 
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER 
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
! 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module holding variables and routines for sparse algebra.
!    Part of the ASAD chemical solver. Contains the following
!    routines:
!      setup_spfuljac
!      spfuljac
!      splinslv2
!      spresolv2
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
!     Method
!     ------
!     Sparse algebra works in the same way as dense algebra, namely by LU
!     decomposition (Gaussian elimination) of the linear system. In contrast
!     to dense algebra, here we keep track of non-zero matrix elements
!     and hence cut out algebraic amnipulations whose result would be zero.
!     To make this efficient, species need to be reorder, such that those
!     with the fewest non-zero matrix elements (reactions) associated with
!     them, occur first in the list, and those with most (generally OH)
!     occur last.
!
! ######################################################################
!
MODULE messy_clamschem_asad_sparse_vars

      USE messy_clamschem_asad_mod
      USE messy_clamschem_asad_dummy, ONLY: ereport
      USE messy_clams_global,         ONLY: prec

      IMPLICIT NONE


      INTEGER, PARAMETER :: maxterms = 160  ! maximum number of nonzero
                                           ! terms for individual species
      INTEGER, PARAMETER :: maxfterms = 100    ! maximum number of terms
                                         ! involving fractional  products
      INTEGER :: total
      INTEGER :: total1    ! total number of nonzero entries in Jacobian

      INTEGER :: nposterms(spfjsize_max), nnegterms(spfjsize_max),      &
                 nfracterms(spfjsize_max)

      INTEGER :: posterms(spfjsize_max, maxterms)
      INTEGER :: negterms(spfjsize_max, maxterms)
      INTEGER :: fracterms(spfjsize_max, maxfterms)
      INTEGER :: base_tracer(spfjsize_max)

      REAL(PREC) :: fraction(spfjsize_max, maxfterms)


      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE SETUP_SPFULJAC(n_points)
!
!  Routine to set up pointer arrays for sparse full jacobian.
!  Note that this routine mirrors the "dense" fuljac routine
!  but only calculates where the nonsero entries would go
!  in a full Jacobian. It then calculates the pointer arrays
!  for a sparse representation of the Jacobian.

        USE messy_clamschem_asad_dummy, ONLY : ereport
        USE messy_clamschem_asad_mod_clams, ONLY: mype, &
                                printstatus, prstatus_oper, prstatus_normal
        IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! Local variables

      INTEGER :: errcode                ! Variable passed to ereport


      INTEGER :: irj
      INTEGER :: ifamd
      INTEGER :: itrd
      INTEGER :: i
      INTEGER :: j
      INTEGER :: jc
      INTEGER :: itrcr
      INTEGER :: j3
      INTEGER :: jn
      INTEGER :: js
      INTEGER :: jl
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: is
      INTEGER :: kr
      INTEGER :: ikr
      INTEGER :: krj
      INTEGER :: ij1
      INTEGER :: itemp1
      INTEGER :: ik(jpmsp)
      INTEGER :: ij(jpmsp)
      INTEGER :: activity(jpctr)

      INTEGER, ALLOCATABLE :: p(:,:)    ! permutation matrix
      INTEGER, ALLOCATABLE :: temp(:,:) ! temporary index matrix
      INTEGER, ALLOCATABLE :: map(:,:)

      LOGICAL :: lj(jpmsp)

      CHARACTER :: ityped*2
      CHARACTER(LEN=72) :: cmessage


! index map for nonzero entries in Jacobian
!
!
      IF (.NOT. ALLOCATED(map)) ALLOCATE(map(jpctr, jpctr))
      map = 0
      DO i=1,jpctr
        map(i,i) = 1
      END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
      DO jc = 1, ntrf
        itrcr = nltrf(jc)
        DO j3 = 1,nmzjac(itrcr)
          irj = nzjac1(j3,itrcr)
!
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,jpmsp
            IF(lj(jn)) map(ij(jn),itrcr) = 1
          END DO
!
          IF(npdfr(irj,1) /= 0) THEN
            i1 = npdfr(irj,1)
            i2 = npdfr(irj,2)
            DO jn = i1, i2
              is = ntabpd(jn,1)
              map(is,itrcr) = 1
            END DO
          END IF
!
        END DO
      END DO
!
!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) and O(3P) are modelled, and that both are
!  required for the O3 loss rate.
!
      DO jc = 1, nstst
        DO j3 = 1,nmsjac(jc)
          irj = nsjac1(j3,jc)
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,jpmsp
            IF(lj(jn)) THEN
              map(ij(jn),ntro3 ) = 1
!              map(ij(jn),ntroh ) = 1
!              map(ij(jn),ntrho2) = 1
!              map(ij(jn),ntrno ) = 1
            END IF
          END DO
        END DO
      END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
      IF ( (ndepw  /=  0) .OR. (ndepd  /=  0) ) THEN
        DO js = 1, jpspec
          ifamd = moffam(js)
          itrd = madvtr(js)
          ityped = ctype(js)
          IF ( itrd /= 0 ) map(itrd,itrd) = 1
        END DO
      END IF
!
! produce pointer variables

      total = SUM(map)
      IF (mype == 0 .AND. printstatus >= prstatus_oper) WRITE(6,*)      &
        'TOTAL NUMBER OF NONZERO ENTRIES IN JACOBIAN: ',total

      is=0
      pntr2(:,:)=0
      DO i=1,jpctr
        DO j=1,jpctr
! calculate forward and backward pointers for sparse representation
          IF (map(i,j) == 1) THEN
            is = is + 1
! backward pointer
            pntr2(i,j) = is
          END IF
        END DO
      END DO

! Reorder species to minimize fill-in
      DO i=1,jpctr
        activity(i) = SUM(map(:,i)) + SUM(map(i,:))
        ro(i) = i
      END DO
      DO i=1,jpctr-1
        DO j=i+1,jpctr
          IF (activity(i) > activity(j)) THEN
! exchange i and j tracers if i is more active than j.
            itemp1 = ro(i)
            ro(i) = ro(j)
            ro(j) = itemp1
            itemp1 = activity(i)
            activity(i) = activity(j)
            activity(j) = itemp1
          END IF
        END DO
      END DO
      DO i=1,jpctr
        IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
          write(6,*) 'IN ASAD SETUP_SPFULJAC: ',&
               advt(ro(i)),' ',activity(i)
        END IF
      END DO

! reorganize pointer variable to account for varying fill-in
      IF (.NOT. ALLOCATED(p)) ALLOCATE(p(jpctr,jpctr))
      IF (.NOT. ALLOCATED(temp)) ALLOCATE(temp(jpctr,jpctr))

      p = 0
      DO i=1,jpctr
        p(i,ro(i)) = 1
      END DO

! Calculate P*A
      DO i=1,jpctr
        DO j=1,jpctr
          temp(i,j) = SUM(p(i,:) * pntr2(:,j))
        END DO
      END DO

! form P*A*P' -> A
      pntr(:,:) = 0.0
      DO i=1,jpctr
        DO j=1,jpctr
          pntr(i,j) = SUM(temp(i,:) * p(j,:))
        END DO
      END DO

      IF (ALLOCATED(p)) DEALLOCATE(p)
      IF (ALLOCATED(temp)) DEALLOCATE(temp)
      IF (ALLOCATED(map)) DEALLOCATE(map)

! calculate production and loss terms
      nposterms = 0
      nnegterms = 0
      nfracterms = 0
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
      DO jc = 1, ntrf
        itrcr = nltrf(jc)
!
        DO j3 = 1,nmzjac(itrcr)
          irj = nzjac1(j3,itrcr)
!
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)= (ij(jn) /= 0)
          END DO
          DO jn=1,2
            IF (lj(jn)) THEN
              i = pntr2(ij(jn),itrcr)
              nnegterms(i) = nnegterms(i) + 1
              IF (nnegterms(i) > maxterms) THEN
                errcode=1
                cmessage=' Increase maxterms'
                CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
              END IF
              negterms(i,nnegterms(i)) = irj
            END IF
          END DO
          IF (nfrpx(irj) == 0) THEN
            DO jn=3,jpmsp
              IF (lj(jn)) THEN
                i = pntr2(ij(jn),itrcr)
                nposterms(i) = nposterms(i) + 1
                IF (nposterms(i) > maxterms) THEN
                  errcode=2
                  cmessage=' Increase maxterms'
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                posterms(i,nposterms(i)) = irj
              END IF
            END DO
          ELSE
            DO jn=3,jpmsp
              IF (lj(jn)) THEN
                i = pntr2(ij(jn),itrcr)
                nfracterms(i) = nfracterms(i) + 1
                IF (nfracterms(i) > maxfterms) THEN
                  errcode = 3
                  cmessage=' Increase maxfterms'
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                IF (nfrpx(irj)+jn-3 > jpfrpx) THEN
                  cmessage = ' frpx array index > jpfrpx'
                  errcode = 3
                  write(6,*) cmessage
                  write(6,*) 'irj: ',irj,' nfrpx(irj): ',nfrpx(irj)
                  write(6,*) 'i: ',i,' jn: ',jn
                  CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
                END IF
                fracterms(i,nfracterms(i)) = irj
                fraction(i,nfracterms(i)) = frpx(nfrpx(irj)+jn-3)
              END IF
            END DO
          END IF
!
          IF (npdfr(irj,1) /= 0) THEN
            i1 = npdfr(irj,1)
            i2 = npdfr(irj,2)
            DO jn = i1, i2
              is = ntabpd(jn,1)
              i = pntr2(is,itrcr)
              nfracterms(i) = nfracterms(i) + 1
              IF (nfracterms(i) > maxfterms) THEN
                errcode=3
                cmessage=' Increase maxfterms'
                CALL EREPORT('SETUP_SPFULJAC',errcode,cmessage)
              END IF
              fracterms(i,nfracterms(i)) = irj
              fraction(i,nfracterms(i)) = ztabpd(jn,1)
            END DO
          END IF
!
        END DO
      END DO

! calculate base_tracer (i.e., number of base tracer that goes with each
! matrix element
      DO i=1,jpctr
        DO j=1,jpctr
          i1 = pntr2(j,i)
          IF (i1 > 0) base_tracer(i1)=i
        END DO
      END DO

! calculate number of matrix elements after factorization. Make sure it is
! smaller than the set maximum.
      total1 = total
      pntr1 = pntr
      DO kr=1,jpctr
        DO i=kr+1,jpctr
          ikr = pntr1(i,kr)
          IF (ikr > 0) THEN
            DO j=kr+1,jpctr
              krj = pntr1(kr,j)
              IF (krj > 0) THEN
                ij1 = pntr1(i,j)
! Distinguish whether matrix element is zero or not. If not, proceed
! as in dense case. If it is, create new matrix element.
                IF (ij1 <= 0) THEN
                  total1 = total1 + 1
                  IF (total1 > spfjsize_max) THEN
                    errcode=total1
                    write(6,*) 'Total1 exceeded spfjsize_max: ',total1, &
                                spfjsize_max
                    cmessage='Increase spfjsize_max'
                    CALL EREPORT                                        &
                      ('SETUP_SPFULJAC',errcode,cmessage)
                  END IF
                  pntr1(i,j) = total1
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO

      IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
        WRITE (6,*) 'Initial and final fill-in: ',total, total1
        WRITE (6,*) 'spfjsize_max is currently set to: ',spfjsize_max
        WRITE (6,*) 'Reduce spfjsize_max if much greater than total1 '
      END IF

      RETURN
      END SUBROUTINE setup_spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE spfuljac(n_points)
!
!  Routine to calculate the sparse full Jacobian
!  Filter for negatives before entering!
!
! 28/3/2006 Include derivative of SS species with respect to
! ozone (deriv) in calculation of Jacobian
!                                 Olaf Morgenstern
!
! 5/4/2006 Bug removed with indexing of dpd and dpw
!                                 Olaf Morgenstern
! 21/6/2006 Adapted from fuljac to calculate sparse
!           full Jacobian.
!                                 Olaf Morgenstern

        USE messy_clamschem_asad_dummy, ONLY : ereport

        IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! Local variables
      REAL(PREC) :: deltt
      REAL(PREC) :: fr

      CHARACTER (LEN=2) :: ityped

      LOGICAL, SAVE :: first = .TRUE.

      INTEGER :: p
      INTEGER :: irj
      INTEGER :: ifamd
      INTEGER :: itrd
      INTEGER :: i
      INTEGER :: j
      INTEGER :: jc
      INTEGER :: itrcr
      INTEGER :: j3
      INTEGER :: jn
      INTEGER :: js
      INTEGER :: jl
      INTEGER :: i1
      INTEGER :: i2
      INTEGER :: is
      INTEGER :: ij(jpmsp)
      INTEGER :: ik(jpmsp)

      LOGICAL :: lj(jpmsp)

!
      deltt=1./cdt
!
      IF (first) THEN
! Determine number and positions of nonzero elements in sparse
! full Jacobian
        CALL setup_spfuljac(n_points)
        first = .FALSE.
      END IF

      spfj = 0.

! Calculate diagonal element of Jacobian
      DO i=1,jpctr
        spfj(:,pntr2(i,i))=-deltt*f(:,i)
      END DO
!
!
!     -----------------------------------------------------------------
!           2.  Calc. full Jacobian matrix.
!               ----- ---- -------- -------
!
! Sum up positive non-fractional, negative, and fractional terms
! Count down (I think that's faster...)
!
      DO p=1,total
        DO i=nposterms(p),1,-1
          spfj(:,p) = spfj(:,p) + prk(:,posterms(p,i))
        END DO
        DO i=nnegterms(p),1,-1
          spfj(:,p) = spfj(:,p) - prk(:,negterms(p,i))
        END DO
        DO i=nfracterms(p),1,-1
          spfj(:,p) = spfj(:,p) + fraction(p,i)*prk(:,fracterms(p,i))
        END DO
      END DO

!  Go through the steady state additions to the Jacobian; currently
!  assume that only O(1D) [and O(3P)] are modelled, and that [both] are
!  required for the O3 loss rate.
!
      DO jc = 1, nstst
        DO j3 = 1,nmsjac(jc)
          irj = nsjac1(j3,jc)
          DO jn=1,jpmsp
            ij(jn)=njcoth(irj,jn)
            lj(jn)=ij(jn) /= 0
          END DO
          DO jn=1,2
            IF (lj(jn)) THEN
              p = pntr2(ij(jn),ntro3)
              spfj(:,p) = spfj(:,p) -                                   &
                prk(:,irj)*deriv(:,jc,1) ! *prkrat(jc)
            END IF
          END DO
          DO jn=3,jpmsp
            IF (lj(jn)) THEN
              p = pntr2(ij(jn),ntro3)
              spfj(:,p) = spfj(:,p) +                                   &
                prk(:,irj)*deriv(:,jc,1) ! *prkrat(jc)
            END IF
          END DO
        END DO
      END DO
!
!     -----------------------------------------------------------------
!          4.  Add deposition terms to Jacobian diagonal.
!              --- ---------- ----- -- -------- ---------
!
      IF ( (ndepw  /=  0) .OR. (ndepd  /=  0) ) THEN
        DO js = 1, jpspec
          ifamd = moffam(js)
          itrd = madvtr(js)
          ityped = ctype(js)
!
          IF ( ifamd /= 0 ) THEN
            DO jl=1,n_points
              IF ((ityped == jpfm) .OR. ((ityped == jpif) .AND.         &
              linfam(jl,itrd))) THEN
                p = pntr2(ifamd,ifamd)
                spfj(jl,p)=spfj(jl,p)                                   &
                - nodd(js)*(dpd(jl,js)+dpw(jl,js))*y(jl,js)
              END IF
            END DO
          END IF
          IF ( itrd /= 0 ) THEN
            p = pntr2(itrd,itrd)
            spfj(:,p) = spfj(:,p) - (dpd(:,js)+dpw(:,js))*y(:,js)
          END IF
        END DO
      END IF
!
!     -------------------------------------------------------------
!          5.  Jacobian elements in final form
!              -------- -------- -- ----- ----
!
      DO p=1,total
        spfj(:,p)=spfj(:,p)/f(:,base_tracer(p))   ! filter f earlier!
      END DO
!
      RETURN
      END SUBROUTINE spfuljac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE splinslv2(zb,zx,n_points,rafeps,rafbig)

        USE messy_clamschem_asad_dummy, ONLY : ereport
        USE messy_clamschem_asad_mod_clams, only: theta_field_size

        IMPLICIT NONE
!
!***** SUB -LINSLV- does simple Gaussian elimination
!*****     fj(N,N)*X(N) = B(N)  BY REDUCING THE A-MATRIX IN PLACE.
!***** THE ENTRY -RESOLV- ASSUMES THAT THE fj-MATRIX HAS BEEN PROPERLY
!*****     REDUCED AND JUST SOLVES FOR X(N).  THIS OPTION SAVES TIME
!*****     WHEN THE SYSTEM IS TO BE RESOLVED WITH A NEW B-VECTOR.
!
! 7/7/2006  Adapted to do sparse matrix algebra and species reordering
!           to improve throughput. For a typical application, this
!           reduces CPU time by at least a factor of 2.
!                                     Olaf Morgenstern
! 1/9/2006  Check diagonal elements of Jacobian whether they are close
!           to 0. Improves stability. Olaf Morgenstern
!
!
! Subroutine interface

      INTEGER, INTENT(IN) :: n_points
      REAL(PREC),    INTENT(IN) :: rafeps
      REAL(PREC),    INTENT(IN) :: rafbig

      REAL(PREC), INTENT(INOUT) :: zb(theta_field_size,jpctr)
      REAL(PREC), INTENT(OUT)   :: zx(theta_field_size,jpctr)

! Local variables
      INTEGER :: kr
      INTEGER :: jl
      INTEGER :: i
      INTEGER :: j
      INTEGER :: ikr
      INTEGER :: krj
      INTEGER :: ij

      REAL(PREC) :: zb1(theta_field_size,jpctr)
      REAL(PREC) :: zx1(theta_field_size,jpctr)
      REAL(PREC) :: pivot(theta_field_size)
      REAL(PREC) :: kfact(theta_field_size)


!
! copy pointer variables into new representations; we do not
! want to change these.

      spfj = MIN(MAX(spfj,-rafbig),rafbig)

      total1 = total
      pntr1 = pntr
      DO kr=1,jpctr
        pivot = spfj(:,pntr1(kr,kr))
        WHERE (ABS(pivot) > rafeps)
          pivot = 1./pivot
        ELSEWHERE
          pivot = rafbig
        ENDWHERE
!        PIVOT = 1./spfj(:,pntr1(kr,kr))
        DO i=kr+1,jpctr
          ikr = pntr1(i,kr)
          IF (ikr > 0) THEN
            KFACT = spfj(:,ikr)*pivot
            spfj(:,ikr) = KFACT
            DO j=kr+1,jpctr
              krj = pntr1(kr,j)
              IF (krj > 0) THEN
                ij = pntr1(i,j)
! Distinguish whether matrix element is zero or not. If not, proceed
! as in dense case. If it is, create new matrix element.
                IF (ij > 0) THEN
                  spfj(:,ij) = spfj(:,ij) - KFACT*spfj(:,krj)
                ELSE
                  total1 = total1 + 1
                  pntr1(i,j) = total1
                  spfj(:,total1) = -KFACT*spfj(:,krj)
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO

! Filter spfj

      spfj = MIN(MAX(spfj,-rafbig),rafbig)

!
!c      call resolv2(zb,zx,n_points,rafeps)!!resolv2 manually inlined below
!
! Form P*zb (permutation of RHS)
      DO i=1,jpctr
        zb1(:,i) = zb(:,ro(i))
      END DO

      DO kr=1,jpctr-1
        DO i=kr+1,jpctr
          ikr = pntr1(i,kr)
          IF (ikr > 0)                                                  &
            zb1(:,i) = zb1(:,i) - spfj(:,ikr) * zb1(:,kr)
        END DO
      END DO
!
      DO kr=jpctr,1,-1
        zx1(:,kr) = zb1(:,kr)
        DO j=kr+1,jpctr
          krj = pntr1(kr,j)
          IF (krj > 0)                                                  &
            zx1(:,kr) = zx1(:,kr) - spfj(:,krj) * zx1(:,j)
        END DO
        pivot = spfj(:,pntr1(kr,kr))
        WHERE (ABS(pivot) < rafeps) pivot = rafeps
        zx1(:,kr) = MIN(MAX(zx1(:,kr),-rafbig),rafbig)/pivot
      END DO

! Form P' zx1 = zx (result of linear system)

      DO j=1,jpctr
        zx(:,ro(j)) = zx1(:,j)
      END DO

      RETURN
      END SUBROUTINE splinslv2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE spresolv2(zb,zx,n_points,rafeps)

!            ENTRY FOR BACK SOLUTION WITH DIFFERENT B-VALUE
! 7/7/2006 Modified to perform sparse matrix processing and species
!          reordering.               Olaf Morgenstern
! 17/8/2006 Bug removed              Olaf Morgenstern

        USE messy_clamschem_asad_dummy, ONLY : ereport
        USE messy_clamschem_asad_mod_clams,  ONLY: theta_field_size
        IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
      REAL(PREC),    INTENT(IN) :: rafeps
      REAL(PREC), INTENT(INOUT) :: zb(theta_field_size,jpctr)
      REAL(PREC), INTENT(INOUT) :: zx(theta_field_size,jpctr)

! Local variables
      INTEGER :: kr
      INTEGER :: jl
      INTEGER :: i
      INTEGER :: j
      INTEGER :: krj
      INTEGER :: ikr

      REAL(PREC) :: zb1(theta_field_size, jpctr)
      REAL(PREC) :: zx1(theta_field_size, jpctr)
      REAL(PREC) :: pivot(theta_field_size)

!
! Form P*zb (permutation of RHS)
      DO i=1,jpctr
        zb1(:,i) = zb(:,ro(i))
      END DO

      DO kr=1,jpctr-1
        DO i=kr+1,jpctr
          ikr = pntr1(i,kr)
          IF (ikr > 0) zb1(:,i) = zb1(:,i) - spfj(:,ikr)*zb1(:,kr)
        END DO
      END DO
!
      DO kr=jpctr,1,-1
        zx1(:,kr) = zb1(:,kr)
        DO j=kr+1,jpctr
          krj = pntr1(kr,j)
          IF (krj > 0) zx1(:,kr) = zx1(:,kr) - spfj(:,krj)*zx1(:,j)
        END DO
        pivot = spfj(:,pntr1(kr,kr))
        WHERE(ABS(pivot) < rafeps) pivot = rafeps
        zx1(:,kr) = zx1(:,kr) / pivot
      END DO

      DO i=1,jpctr
        zx(:,ro(i)) = zx1(:,i)
      END DO

      RETURN
      END SUBROUTINE spresolv2

END MODULE messy_clamschem_asad_sparse_vars

