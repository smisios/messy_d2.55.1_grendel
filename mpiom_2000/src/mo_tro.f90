MODULE mo_tro

! module file containing the sbrs for iterative sea level solver (SOR)
  USE mo_kind, ONLY: wp
    USE mo_boundsexch, ONLY: bounds_exch

    IMPLICIT NONE

    REAL(wp) ,ALLOCATABLE :: dilcor(:,:),susalo(:,:),suwath(:,:)

  CONTAINS

    SUBROUTINE itprep

    ! rj: rewritten for mpi-parallelization nov 2003
    ! calculates some arrays needed for the iterative solver

    USE mo_param1, ONLY: ie, je, stabn
    USE mo_planetary_constants, ONLY: g
    USE mo_para2, ONLY: xx, uf, vf, uf0,vf0, ff
    USE mo_commo1, ONLY: area,dlxu,dlxv,dlyu,dlyv,deute,deuto &
                        ,conn,dt,amsuo,amsue,lwith_barotropic_stokes_drift

    INTEGER i,j

    xx(:,:) = area(:,:)
    uf(:,:) = deuto(:,:)*(dlyu(:,:)/dlxu(:,:))*g*conn*stabn*dt**2
    vf(:,:) = deute(:,:)*(dlxv(:,:)/dlyv(:,:))*g*conn*stabn*dt**2

    IF ( lwith_barotropic_stokes_drift ) THEN
      uf0(:,:) = amsuo(:,:,1)*(dlyu(:,:)/dlxu(:,:))*g*conn*stabn*dt**2
      vf0(:,:) = amsue(:,:,1)*(dlxv(:,:)/dlyv(:,:))*g*conn*stabn*dt**2
    ENDIF

    ff(:,:) = 0._wp
    DO j=2,je-1
      DO i=2,ie-1
        ff(i,j) = 1._wp / (xx(i,j)+uf(i,j)+uf(i-1,j)+vf(i,j)+vf(i,j-1))
      ENDDO
    ENDDO

    CALL bounds_exch(1,'p',ff,'mo_tro 4')

  END SUBROUTINE itprep

  FUNCTION sqrtrnd()

    REAL(wp) :: sqrtrnd
    REAL(wp), SAVE :: d=0.34679201245_wp
    REAL(wp), PARAMETER :: b=34251._wp
    REAL(wp) :: c
    c = SQRT(b*d)
    d = c - AINT(c)
    sqrtrnd = d
  END FUNCTION sqrtrnd

!!$  SUBROUTINE trotest_lk
!!$
!!$    USE mo_param1,   ONLY: ie_g, je_g
!!$    USE mo_para2,    ONLY: xx, uf, vf, ff, sorpar
!!$    USE mo_commo1,   ONLY: wetol1_g,z1o,b1o
!!$    USE mo_units,    ONLY: io_stdout
!!$    USE mo_parallel, ONLY: p_pe, p_io, nprocxy, p_recv, p_send, p_bcast, &
!!$         gather_arr,scatter
!!$    !
!!$    !=======================================================================
!!$    !     purpose :
!!$    !
!!$    !
!!$    !uwe  preparation for iterative solution of zeta field
!!$    !     optimises sorpar
!!$    !
!!$    ! rj: rewritten for mpi-parallelization nov 2003
!!$    !
!!$
!!$    IMPLICIT NONE
!!$
!!$    INTEGER, PARAMETER :: isormax = 250
!!$
!!$    REAL(wp) :: sor_omega(isormax), convergence(isormax)
!!$
!!$    REAL(wp) :: b1o_g(ie_g,je_g), z1o_g(ie_g,je_g)
!!$    REAL(wp) :: xx_g(ie_g,je_g)
!!$    REAL(wp) :: uf_g(ie_g,je_g), vf_g(ie_g,je_g), ff_g(ie_g,je_g)
!!$
!!$    INTEGER :: i, j, is, n, iss, ies, isor, kiter, isor_start, isor_end
!!$    INTEGER :: i0(1),jb
!!$    REAL(wp) :: chabas, sorpai, sumcha, zalt
!!$    REAL(wp) :: sqrtrnd, rnd
!!$
!!$    ! build global fields on io pe
!!$
!!$    CALL gather_arr(xx,xx_g,p_io)
!!$    CALL gather_arr(uf,uf_g,p_io)
!!$    CALL gather_arr(vf,vf_g,p_io)
!!$    CALL gather_arr(ff,ff_g,p_io)
!!$
!!$!    ! distribute global fields to all pe
!!$!
!!$!    CALL p_bcast(xx_g,p_io)
!!$!    CALL p_bcast(uf_g,p_io)
!!$!    CALL p_bcast(vf_g,p_io)
!!$!    CALL p_bcast(ff_g,p_io)
!!$
!!$    !uwe  initialize with random numbers
!!$    !hh random_number doesn't work with fujitsu f90!
!!$    !hh    call random_number(xhh)
!!$    !hh    b1o(i,j)=(xhh-0.5)*weto(i,j,1)
!!$    !otb  using a simple one
!!$
!!$    ! rj: we do that in the following way so that the result
!!$    !     does not depend on the decomposition
!!$    !     why is that done for every step in the original version ????
!!$    !     we do it only once here
!!$
!!$    IF (p_pe==p_io) THEN
!!$    b1o_g(:,:) = 0.0
!!$    DO j = 1, je_g
!!$       DO i = 2, ie_g-1
!!$          rnd = sqrtrnd()
!!$          IF (i >=1 .AND. i <= ie_g .AND. j >= 1 .AND. j <= je_g) THEN
!!$             b1o_g(i,j) = (rnd-0.5)*wetol1_g(i,j)*xx_g(i,j)
!!$          ENDIF
!!$       ENDDO
!!$    ENDDO
!!$    ENDIF
!!$
!!$    CALL scatter(b1o_g,b1o,p_io)
!!$
!!$!    b1o_g(ie_g,:) = b1o_g(2,:)
!!$!    b1o_g(1,:)    = b1o_g(ie_g-1,:)
!!$
!!$    jb=2
!!$    DO isor = 1, isormax
!!$       sor_omega(isor) = 1.70+isor*0.001
!!$    ENDDO
!!$
!!$    convergence(:) = 9.e99
!!$
!!$    CALL decompose (isormax, nprocxy, p_pe, isor_start, isor_end)
!!$
!!$    !$omp parallel private (isor, sorpar, sorpai, z1o_g, kiter, sumcha, i, j, is, &
!!$    !$omp                   zalt, chabas)
!!$
!!$    !$omp do
!!$    DO isor = isor_start, isor_end
!!$       !
!!$       sorpar = sor_omega(isor)
!!$       sorpai = 1.0-sorpar
!!$       !
!!$       z1o_g(:,:) = 0.
!!$       !
!!$       DO kiter = 1, 250
!!$
!!$         IF (p_pe==p_io) THEN
!!$
!!$          sumcha = 0.0
!!$
!!$          ! odd/odd and even/even z1o elements
!!$
!!$          DO j=jb,je_g-1
!!$             is = MOD(j,2)+2
!!$             DO i=is,ie_g-1,2
!!$                zalt = z1o_g(i,j)
!!$                z1o_g(i,j) = sorpar*ff_g(i,j)*(b1o_g(i,j)                      &
!!$                     +  uf_g(i,j)*z1o_g(i+1,j) + uf_g(i-1,j)*z1o_g(i-1,j)      &
!!$                     +  vf_g(i,j)*z1o_g(i,j+1) + vf_g(i,j-1)*z1o_g(i,j-1))     &
!!$                     + sorpai*z1o_g(i,j)
!!$                sumcha=sumcha+(zalt-z1o_g(i,j))**2
!!$             ENDDO
!!$          ENDDO
!!$
!!$          ENDIF
!!$
!!$          CALL scatter(z1o_g,z1o,p_io)
!!$          CALL gather_arr(z1o,z1o_g,p_io)
!!$
!!$
!!$!          z1o_g(ie_g,:) = z1o_g(2,:)
!!$!          z1o_g(1,:)    = z1o_g(ie_g-1,:)
!!$
!!$          ! odd/even and even/odd z1o elements
!!$
!!$       IF (p_pe==p_io) THEN
!!$
!!$          DO j=jb,je_g-1
!!$             is = MOD(j+1,2)+2
!!$             DO i=is,ie_g-1,2
!!$                zalt = z1o_g(i,j)
!!$                z1o_g(i,j) = sorpar*ff_g(i,j)*(b1o_g(i,j)                      &
!!$                     +  uf_g(i,j)*z1o_g(i+1,j) + uf_g(i-1,j)*z1o_g(i-1,j)      &
!!$                     +  vf_g(i,j)*z1o_g(i,j+1) + vf_g(i,j-1)*z1o_g(i,j-1))     &
!!$                     + sorpai*z1o_g(i,j)
!!$                sumcha=sumcha+(zalt-z1o_g(i,j))**2
!!$             ENDDO
!!$          ENDDO
!!$
!!$       ENDIF
!!$
!!$          CALL scatter(z1o_g,z1o,p_io)
!!$          CALL gather_arr(z1o,z1o_g,p_io)
!!$
!!$!          z1o_g(ie_g,:) = z1o_g(2,:)
!!$!          z1o_g(1,:)    = z1o_g(ie_g-1,:)
!!$
!!$          IF (kiter == 1) chabas = sumcha
!!$
!!$       ENDDO ! kiter loop
!!$
!!$       convergence(isor) = sumcha/chabas
!!$
!!$    ENDDO ! isor loop
!!$
!!$    !$omp end parallel
!!$
!!$    IF (p_pe == p_io) THEN
!!$       DO n = 0, nprocxy
!!$          IF  (p_pe == p_io) CYCLE
!!$          CALL decompose (isormax, nprocxy, n, iss, ies)
!!$          CALL p_recv(convergence(iss), n, 200, &
!!$               p_count=ies-iss+1)
!!$       ENDDO
!!$    ELSE
!!$       CALL p_send(convergence(isor_start),p_io, 200, &
!!$            p_count=isor_end-isor_start+1)
!!$    ENDIF
!!$
!!$    IF (p_pe == p_io) THEN
!!$       i0 = MINLOC(convergence)
!!$       sorpar = sor_omega(i0(1))
!!$    ENDIF
!!$
!!$    CALL p_bcast(sorpar,p_io)
!!$
!!$    WRITE(io_stdout,*) 'trotest: sor parameter = ',sorpar
!!$
!!$  CONTAINS
!!$
!!$    SUBROUTINE decompose (ns, npes, pe, istart, iend)
!!$
!!$      INTEGER, INTENT(in) :: ns       ! domain size
!!$      INTEGER, INTENT(in) :: npes     ! pes availbale for that decomposition
!!$      INTEGER, INTENT(in) :: pe       ! this pe
!!$
!!$      INTEGER, INTENT(out) :: istart  ! start index of sub domain
!!$      INTEGER, INTENT(out) :: iend    ! end index of sub domain
!!$
!!$      INTEGER :: nlocal, nremain
!!$
!!$      nlocal = ns/npes
!!$      istart = pe*nlocal+1
!!$      nremain = MOD(ns,npes)
!!$      istart = istart+MIN(pe,nremain)
!!$      IF (pe < nremain) THEN
!!$         nlocal = nlocal+1
!!$      ENDIF
!!$      iend = istart+nlocal-1
!!$      IF (iend > ns .OR. pe == npes-1) iend = ns
!!$
!!$    END SUBROUTINE decompose
!!$
!!$  END SUBROUTINE trotest_lk

  SUBROUTINE TROTEST
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS

      USE MO_PARALLEL
!
!=======================================================================
!     PURPOSE :
!
!
!UWE  PREPARATION FOR ITERATIVE SOLUTION OF ZETA FIELD
!     OPTIMISES SORPAR
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!

      INTEGER I, J, II, JJ, IS, ISOR, KITER
      REAL(wp) CHABAS, CONTRIC, CONVMIN, SORPAI, SORPMIN, SUMCHA, ZALT
      REAL(wp) RND

      sorpar = 1.75_wp

      convmin = 9.E99_wp
      SORPMIN=SORPAR
!
      DO ISOR=1,200
        sorpar = sorpar + 0.001_wp
!
!UWE  INITIALIZE WITH RANDOM NUMBERS
!hh random_number doesn't work with fujitsu f90!
!hh    CALL RANDOM_NUMBER(XHH)
!hh    B1O(I,J)=(XHH-0.5)*WETO(I,J,1)
!OtB  Using a simple one

! RJ: We do that in the following way so that the result
!     does not depend on the decomposition
!     Why is that done for every step in the original version ????
!     We do it only once here

      IF (ISOR==1) THEN
        DO JJ=1,JE_G
        DO II=2,IE_G-1
          I = II - p_ioff
          J = JJ - p_joff
          RND = SQRTRND()
          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
            b1o(i, j) = (rnd - 0.5_wp) * weto(i, j, 1) * xx(i, j)
          ENDIF
        ENDDO
        ENDDO

        CALL bounds_exch(1,'p',B1O,'mo_tro 5')
      ENDIF

      Z1O(:,:) = 0._wp

      sorpai = 1._wp - sorpar

      DO KITER=1,200

        sumcha = 0._wp

        ! odd/odd and even/even Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=2,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence ',sorpar,contric
      ENDIF

      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop

      SORPAR=SORPMIN
      WRITE(IO_STDOUT,*) 'TROTEST SORPAR=',SORPAR

    END SUBROUTINE TROTEST

    SUBROUTINE TROTEST2
      USE MO_PARAM1
      USE MO_PARA2
      USE MO_COMMO1
      USE MO_UNITS

      USE MO_PARALLEL
!
!=======================================================================
!     PURPOSE :
!
!
!UWE  PREPARATION FOR ITERATIVE SOLUTION OF ZETA FIELD
!     OPTIMISES SORPAR
!
! RJ: Rewritten for MPI-Parallelization Nov 2003
!

      INTEGER I, J, II, JJ, IS, ISOR, KITER,jb
      REAL(wp) CHABAS, CONTRIC, CONVMIN, SORPAI, SORPMIN, SUMCHA, ZALT
      REAL(wp) RND

      sorpar = 0.9_wp
!
      convmin = 9.E99_wp
      SORPMIN=SORPAR

      jb=2
      if( p_joff.eq.0 .and. lbounds_exch_tp ) jb=3


!     1) find SORPMIN between 1. and 2.0
      DO ISOR=1,10
        sorpar = sorpar + 0.1_wp

      IF (ISOR==1) THEN
        DO JJ=1,JE_G
        DO II=2,IE_G-1
          I = II - p_ioff
          J = JJ - p_joff
          RND = SQRTRND()
          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
            b1o(i, j) = (rnd - 0.5_wp) * weto(i, j, 1) * xx(i, j)
          ENDIF
        ENDDO
        ENDDO

        CALL bounds_exch(1,'p',B1O,'mo_tro 5')
      ENDIF

      z1o(:,:) = 0._wp

      sorpai = 1._wp - sorpar

      DO KITER=1,200

        sumcha = 0._wp

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence 1: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop


!     2) find SORPMIN between -0.1 and +0.1
      ! FIXME: this should probably read
      ! convmin = HUGE(convmin)
      convmin = 9.e99_wp
      sorpar = sorpmin - 0.1_wp
      SORPMIN=SORPAR

      DO ISOR=1,20
      sorpar = sorpar + 0.01_wp

!!$      IF (ISOR==1) THEN
!!$        DO JJ=1,JE_G
!!$        DO II=2,IE_G-1
!!$          I = II - p_ioff
!!$          J = JJ - p_joff
!!$          RND = SQRTRND()
!!$          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
!!$            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
!!$          ENDIF
!!$        ENDDO
!!$        ENDDO
!!$
!!$        CALL bounds_exch(1,'p',B1O,'mo_tro 5')
!!$      ENDIF

      z1o(:,:) = 0._wp

      sorpai = 1._wp - sorpar

      DO KITER=1,200

        sumcha = 0._wp

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.NE.0)THEN
         WRITE(6,*)'sorpa convergence 2: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop

!     3) find SORPMIN between -0.01 and +0.01

      ! FIXME: this should probably read
      ! convmin = HUGE(convmin)
      convmin = 9.e99_wp
      sorpar = sorpmin - 0.01_wp
      SORPMIN=SORPAR

      DO ISOR=1,20
        sorpar = sorpar + 0.001_wp

!!$      IF (ISOR==1) THEN
!!$        DO JJ=1,JE_G
!!$        DO II=2,IE_G-1
!!$          I = II - p_ioff
!!$          J = JJ - p_joff
!!$          RND = SQRTRND()
!!$          IF(I>=1 .AND. I<=IE .AND. J>=1 .AND. J<=JE) THEN
!!$            B1O(I,J)=(RND-0.5)*WETO(I,J,1)*XX(I,J)
!!$          ENDIF
!!$        ENDDO
!!$        ENDDO
!!$
!!$        CALL bounds_exch(1,'p',B1O,'mo_tro 5')
!!$      ENDIF

      z1o(:,:) = 0._wp

      sorpai = 1._wp - sorpar

      DO KITER=1,200

        sumcha = 0._wp

        ! odd/odd and even/even Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 6')

        ! odd/even and even/odd Z1O elements

        DO J=jb,JE-1
          IS = MOD(p_ioff+p_joff+j+1,2)+2
          DO I=IS,IE-1,2
            ZALT = Z1O(I,J)
            Z1O(I,J) = SORPAR*FF(I,J)*(B1O(I,J)                         &
     &                +  UF(I,J)*Z1O(I+1,J) + UF(I-1,J)*Z1O(I-1,J)      &
     &                +  VF(I,J)*Z1O(I,J+1) + VF(I,J-1)*Z1O(I,J-1))     &
     &              + SORPAI*Z1O(I,J)
            SUMCHA=SUMCHA+(ZALT-Z1O(I,J))**2
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',Z1O,'mo_tro 7')

        CALL global_sum(SUMCHA)

        IF(KITER.EQ.1) CHABAS=SUMCHA

      ENDDO ! KITER Loop

      CONTRIC=SUMCHA/CHABAS

      IF (icontro.GT.0)THEN
         WRITE(6,*)'sorpa convergence 3: ',sorpar,contric
      ENDIF


      IF(CONTRIC.LT.CONVMIN)THEN
        SORPMIN=SORPAR
        CONVMIN=CONTRIC
      ENDIF

      ENDDO ! ISOR Loop




      SORPAR=SORPMIN
      WRITE(IO_STDOUT,*) 'TROTEST SORPAR2=',SORPAR

    END SUBROUTINE TROTEST2


  SUBROUTINE troneu
    !**********************************************************************
    !
    !
    !     tttttt  rrrrr    ooo    nn   n  eeeeee  u    u
    !       tt    r    r  o   o   n n  n  e       u    u
    !       tt    rrrrr   o   o   n  n n  eeeee   u    u
    !       tt    r rr    o   o   n   nn  e       u    u
    !       tt    r   rr   ooo    n    n  eeeeee  uuuuuu
    !
    !
    !**********************************************************************
    !
    USE mo_param1, ONLY: ie,je,ie1,je1
    USE mo_para2, ONLY: ff,uf,vf,sorpar,iter_sor,rtsorpar,iter_sor_hack,rtsorpar_hack,uf0,vf0,xx
    USE mo_commo1, ONLY: weto,dt,cono,u1o,uzo,conn,dlyu,v1e,vze   &
                        ,dlxv,b1o,z1o,icontro,lbounds_exch_tp     &
                        ,lwith_barotropic_stokes_drift,deute,deuto &
                        ,deution,deutien,amsue,dduo,ddue   &
                        ,amsuo,uso,vse,surlen,surlon
    USE mo_units
    USE mo_parallel, ONLY: p_ioff,p_joff
    USE mo_mpi, ONLY: p_pe,p_io


    INTEGER i, j, is, kiter,jb
    REAL(wp) sorpai,zurr,zalt

    REAL(wp), ALLOCATABLE :: b1oo(:,:),d1o(:,:),d1u(:,:),d1v(:,:),deuton(:,:),deuten(:,:)


    !
    !=======================================================================
    !     sbr tropit
    !
    !     purpose :
    !
    !     a) iterative solution of zeta field
    !
    !
    !   iterative solution of the system
    !   z*(dx*dy+g*dt**2(hw*dyw/dxw+ho*dyo/dxo+hs*dxs/dys+hn*dxn/dyn) )
    !
    ! =            g*dt**2*(hw*zww*dyw/dxw+ho*zoo*dyo/dxo
    !                      +hs*zss*dxs/dys+hn*znn*dxn/dyn  )
    !+g*dt**3*(f*hw*(zsw-znw)+f*ho*(zno-zso)+fs*hs*(zso-zsw)+fn*hn*(znw-zno)
    !                 +b
    !
    !   where b contains the winstress, divergence of old flow and old z
    !
    !
    ! rj: rewritten for mpi-parallelization nov 2003

    IF ( lwith_barotropic_stokes_drift ) THEN

      ALLOCATE( b1oo(ie,je),d1o(ie,je),d1u(ie,je),d1v(ie,je))
      d1o = 0._wp

      DO j=2,je-1
        DO i=2,ie-1
          b1oo(i,j) = weto(i,j,1) *dt* (                                      &
               (cono * u1o(i-1,j) +  conn * uzo(i-1,j) ) * dlyu(i-1,j)        &
               - (cono * u1o(i,j)   +  conn * uzo(i,j)   ) * dlyu(i,j)        &
               + (cono * v1e(i,j)   +  conn * vze(i,j)   ) * dlxv(i,j)        &
               - (cono * v1e(i,j-1) +  conn * vze(i,j-1) ) * dlxv(i,j-1))
        ENDDO
      ENDDO

      CALL bounds_exch(1,'p',b1oo,'mo_tro 8')

    ELSE

      DO j=2,je-1
        DO i=2,ie-1
          b1o(i,j) = weto(i,j,1) *dt* (                                         &
                     (cono * u1o(i-1,j) +  conn * uzo(i-1,j) ) * dlyu(i-1,j)    &
                     - (cono * u1o(i,j)   +  conn * uzo(i,j)   ) * dlyu(i,j)    &
                     + (cono * v1e(i,j)   +  conn * vze(i,j)   ) * dlxv(i,j)    &
                     - (cono * v1e(i,j-1) +  conn * vze(i,j-1) ) * dlxv(i,j-1))
        ENDDO
      ENDDO

      CALL bounds_exch(1,'p',b1o,'mo_tro 8')

    ENDIF

!   DO j=1,je
!      DO i=1,ie
!         z1o(i,j)=0.
!      ENDDO
!   ENDDO


    jb=2
    IF( p_joff.EQ.0 .AND. lbounds_exch_tp ) jb=3

    IF (rtsorpar .GT. 0._wp) THEN
      sorpar=rtsorpar
      IF (p_pe==p_io) WRITE(0,*) 'attn: user definded sorpar: ', sorpar
    ENDIF

    sorpai = 1._wp - sorpar

    DO kiter=1,iter_sor

      IF (kiter.GT.iter_sor-iter_sor_hack) THEN
        IF (rtsorpar_hack .GT. 0._wp) THEN
          sorpar=rtsorpar_hack
          sorpai = 1._wp - sorpar
        ENDIF
      ENDIF

      zurr = 0._wp

      IF ( lwith_barotropic_stokes_drift ) THEN

        d1o = 0.9_wp * d1o + 0.1_wp * z1o

        DO j=1,je1
          DO i=1,ie1
            uf(i,j) = MAX(0._wp, uf0(i, j) * (deuto(i, j) &
                 + 0.5_wp * (d1o(i, j) + d1o(i+1, j))))
            vf(i,j) = MAX(0._wp, vf0(i,j) * (deute(i,j) &
                 + 0.5_wp *(d1o(i, j) + d1o(i, j+1))))
            d1u(i,j) = 0.5_wp * (d1o(i,j)+d1o(i+1,j))
            d1v(i,j) = 0.5_wp * (d1o(i,j)+d1o(i,j+1))
          ENDDO
        ENDDO

        DO j=2,je-1
          DO i=2,ie-1

            b1o(i,j) = b1oo(i,j)+weto(i,j,1) *dt* conn*(         &
                 (  uso(i-1,j)*d1u(i-1,j) ) * dlyu(i-1,j)        &
                 - (  uso(i,j)*d1u(i,j)   ) * dlyu(i,j)          &
                 + (  vse(i,j)*d1v(i,j)   ) * dlxv(i,j)          &
                 - (  vse(i,j-1) )*d1v(i,j-1) * dlxv(i,j-1))

          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',b1o,'mo_tro 8')
        CALL bounds_exch(1,'u+',uf,'mo_tro 8')
        CALL bounds_exch(1,'v+',vf,'mo_tro 8')

        DO j=2,je-1
          DO i=2,ie-1
            ff(i,j) = 1._wp/(xx(i, j) + uf(i, j) + uf(i-1, j) &
                 &           + vf(i, j) + vf(i, j-1))
          ENDDO
        ENDDO

        CALL bounds_exch(1,'p',ff,'mo_tro 4')

      ENDIF



      ! SOR with red/black splitting
      ! 1step: update of the "red" points

      IF (icontro.NE.0) THEN

        DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j,2)+2
          DO i=is,ie-1,2
            zalt=z1o(i,j)
            z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                     &
                 +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)       &
                 +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))      &
                 + sorpai*z1o(i,j)
            zurr=zurr+(zalt-z1o(i,j))**2
          ENDDO
        ENDDO

      ELSE


!      WRITE(0,*) p_pe, MAXVAL(ff),MAXVAL(b1o),MAXVAL(uf),MAXVAL(vf),MAXVAL(z1o)


        DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j,2)+2
          DO i=is,ie-1,2
            zalt=z1o(i,j)
            z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                     &
                 +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)       &
                 +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))      &
                 + sorpai*z1o(i,j)
          ENDDO
        ENDDO

      ENDIF

      CALL bounds_exch(1,'p',z1o,'mo_tro 9')


      IF ( lwith_barotropic_stokes_drift ) THEN

        DO j=2,je1
          DO i=2,ie1
            uf(i,j) = MAX(0._wp, &
                 uf0(i, j) * (deuto(i,j) + 0.5_wp * (d1o(i, j) + d1o(i+1, j))))
            vf(i,j) = MAX(0._wp, &
                 vf0(i, j) * (deute(i,j) + 0.5_wp * (d1o(i, j) + d1o(i, j+1))))
          ENDDO
        ENDDO

        CALL bounds_exch(1,'u+',uf,'mo_tro 8')
        CALL bounds_exch(1,'v+',vf,'mo_tro 8')
        d1o = 0.9_wp * d1o + 0.1_wp * z1o

      ENDIF


      ! 2step: update of the "black" points

      IF (icontro.NE.0) THEN

        DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
          DO i=is,ie-1,2
            zalt=z1o(i,j)
            z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                        &
                 +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)      &
                 +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))     &
                 + sorpai*z1o(i,j)
            zurr=zurr+(zalt-z1o(i,j))**2
          ENDDO
        ENDDO

      ELSE

        DO j=jb,je-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
          DO i=is,ie-1,2
            zalt=z1o(i,j)
            z1o(i,j) = sorpar*ff(i,j)*(b1o(i,j)                        &
                 +  uf(i,j)*z1o(i+1,j) + uf(i-1,j)*z1o(i-1,j)      &
                 +  vf(i,j)*z1o(i,j+1) + vf(i,j-1)*z1o(i,j-1))     &
                 + sorpai*z1o(i,j)
          ENDDO
        ENDDO

      ENDIF

      CALL bounds_exch(1,'p',z1o,'mo_tro 10')

    ENDDO

    IF (icontro.NE.0) THEN
      WRITE(0,*)'zurr',zurr,'sorpar',sorpar
    ENDIF

    IF ( lwith_barotropic_stokes_drift ) THEN

      ALLOCATE(deuton(ie,je),deuten(ie,je))

      deuton = 0._wp
      deuten = 0._wp

      !      write(io_stdout,*)'d1o-z1o',maxval(d1o-z1o),minval(d1o-z1o),zurr
      DO j=1,je1
        DO i=1,ie1
          deuton(i, j) = amsuo(i, j, 1) * &
               MAX(0._wp, deuto(i, j) + 0.5_wp * (z1o(i, j) + z1o(i+1, j)))
          deuten(i, j) = amsue(i, j, 1) * &
               MAX(0._wp, deute(i, j) + 0.5_wp * (z1o(i, j) + z1o(i, j+1)))
          surlon(i, j) = amsuo(i, j, 1) * &
               MAX(0._wp, dduo(i, j, 1) + 0.5_wp * (z1o(i, j) + z1o(i+1, j)))
          surlen(i, j) = amsue(i, j, 1) * &
               MAX(0._wp, ddue(i, j, 1) + 0.5_wp * (z1o(i, j) + z1o(i, j+1)))
        ENDDO
      ENDDO

      CALL bounds_exch(1,'u+',deuton,'mo_tro 10')
      CALL bounds_exch(1,'u+',surlon,'mo_tro 10')
      CALL bounds_exch(1,'v+',deuten,'mo_tro 10')
      CALL bounds_exch(1,'v+',surlen,'mo_tro 10')

      WHERE (deuton .GT. 0._wp) deution = 1._wp / deuton
      WHERE (deuten .GT. 0._wp) deutien = 1._wp / deuten

    ENDIF


  END SUBROUTINE troneu


  SUBROUTINE troneu2(ihalo_sor,imod)
    !**********************************************************************
    !
    !
    !     tttttt  rrrrr    ooo    nn   n  eeeeee  u    u  222222
    !       tt    r    r  o   o   n n  n  e       u    u       2
    !       tt    rrrrr   o   o   n  n n  eeeee   u    u  222222
    !       tt    r rr    o   o   n   nn  e       u    u  2
    !       tt    r   rr   ooo    n    n  eeeeee  uuuuuu  222222
    !
    !
    !**********************************************************************
    !
    USE mo_param1, ONLY: ie,je
    USE mo_para2, ONLY: ff,uf,vf,sorpar,iter_sor,rtsorpar,iter_sor_hack,rtsorpar_hack
    USE mo_commo1, ONLY: weto,dt,cono,u1o,uzo,conn,dlyu,v1e,vze   &
                        ,dlxv,b1o,z1o,icontro
    USE mo_units
    USE mo_parallel, ONLY: p_ioff, p_joff, p_pe, p_io &
#ifdef CARTESIAN_COORD
                           ,nprocx,nprocy &
#endif
                           ,sethaloN


    INTEGER i, j, is, kiter,jb,ihm,ihalo_sor,imod
    REAL(wp) sorpai,zurr,zalt(ie,je)
    REAL(wp),allocatable:: zz1o(:,:),zff(:,:),zb1o(:,:),zuf(:,:),zvf(:,:)

#ifdef CARTESIAN_COORD
    INCLUDE 'mpif.h'

    INTEGER max_dims
    PARAMETER(max_dims=2)
    INTEGER ifehler, old_rank, my_rank, p_size

    INTEGER dims(max_dims), coords(max_dims)
    LOGICAL reorder, periods(max_dims)

    INTEGER status(MPI_STATUS_SIZE)
    INTEGER request

#endif
    INTEGER new_comm
    !
    !=======================================================================
    !     sbr tropit
    !
    !     purpose :
    !
    !     a) iterative solution of zeta field
    !
    !
    !   iterative solution of the system
    !   z*(dx*dy+g*dt**2(hw*dyw/dxw+ho*dyo/dxo+hs*dxs/dys+hn*dxn/dyn) )
    !
    ! =            g*dt**2*(hw*zww*dyw/dxw+ho*zoo*dyo/dxo
    !                      +hs*zss*dxs/dys+hn*znn*dxn/dyn  )
    !+g*dt**3*(f*hw*(zsw-znw)+f*ho*(zno-zso)+fs*hs*(zso-zsw)+fn*hn*(znw-zno)
    !                 +b
    !
    !   where b contains the winstress, divergence of old flow and old z
    !
    !
    ! rj: rewritten for mpi-parallelization nov 2003


    ihm=2*(ihalo_sor-1)
    allocate(zz1o(ie+ihm,je+ihm),zff(ie+ihm,je+ihm),zb1o(ie+ihm,je+ihm),zuf(ie+ihm,je+ihm),zvf(ie+ihm,je+ihm))


    DO j=2,je-1
       DO i=2,ie-1

          b1o(i,j) = weto(i,j,1) *dt* (                                       &
                     (cono * u1o(i-1,j) +  conn * uzo(i-1,j) ) * dlyu(i-1,j)  &
                   - (cono * u1o(i,j)   +  conn * uzo(i,j)   ) * dlyu(i,j)    &
                   + (cono * v1e(i,j)   +  conn * vze(i,j)   ) * dlxv(i,j)    &
                   - (cono * v1e(i,j-1) +  conn * vze(i,j-1) ) * dlxv(i,j-1))

       ENDDO
    ENDDO

    CALL bounds_exch(1,'p',b1o,'mo_tro 8')

    DO j=1,je
       DO i=1,ie
          z1o(i, j) = 0._wp
       ENDDO
    ENDDO

    DO j=1,je+ihm
       DO i=1,ie+ihm
          zz1o(i, j) = 0._wp
       ENDDO
    ENDDO


     jb=2





#ifdef CARTESIAN_COORD

     CALL MPI_COMM_SIZE(MPI_COMM_WORLD, p_size, ifehler)
!
! Set one-dimensional cartesian topology.
!
!      write(0,*)'In CARTESIAN_COORD p_size=',p_size
!      write(0,*)'In CARTESIAN_COORD npes=',npes
      dims(1) = nprocy
      dims(2) = nprocx
      periods(1) = .TRUE.
      periods(2) = .TRUE.
      reorder = .FALSE.

      CALL MPI_CART_CREATE(MPI_COMM_WORLD, max_dims, dims, &
                          periods, reorder, new_comm, ifehler)

!      CALL MPI_COMM_RANK(MPI_COMM_WORLD, old_rank, ifehler)

      CALL MPI_COMM_RANK(new_comm, my_rank, ifehler)
!      write(0,*)'In CARTESIAN_COORD old_rank=',old_rank,'my_rank=',my_rank

#endif
!    do i=1,ie
!    do j=1,je
!       z1o(i,j)=float(i)
!       z1o(i,j)=float(j)
!   enddo
!   enddo



!    call sethalo2(z1o,zz1o,2)

!    do j=1,je+2
!    write(0,*) '2',j,zz1o(1,j),zz1o(2,j),zz1o(3,j),zz1o(4,j)
!    enddo

!   CALL bounds_exch(1,'p',z1o,'mo_tro 10')
   call sethaloN('p',z1o,zz1o,ihalo_sor,new_comm)

!    if ( p_pe == 0 ) then

!XXX
!       do j=1,je+ihm
!         write(0,66388) (zz1o(i,j),i=1,8),(zz1o(i,j),i=ie+ihm-7,ie+ihm)
! 66388    format(8(1x,F3.0),' ..',8(1x,F3.0))
!       enddo



!    do j=1,je+ihm
!    write(0,*) 'links ',p_pe,j,zz1o(1,j),zz1o(2,j),zz1o(3,j),zz1o(4,j)
!    enddo

!    do j=1,je+ihm
!    write(0,*) 'rechts',p_pe,j,zz1o(ie,j),zz1o(ie+1,j),zz1o(ie+2,j),zz1o(ie+3,j),zz1o(ie+4,j)
!    enddo
!
!    do i=1,ie+ihm
!    write(0,*) 'oben ',p_pe,i,zz1o(i,1),zz1o(i,2),zz1o(i,3),zz1o(i,4)
!    enddo
!
!    do i=1,ie+ihm
!    write(0,*) 'unten',p_pe,i,zz1o(i,je),zz1o(i,je+1),zz1o(i,je+2),zz1o(i,je+3),zz1o(i,je+4)
!    enddo
!

!   endif
!   stop

!   CALL bounds_exch(1,'p',b1o,'mo_tro 10')
   call sethaloN('p',b1o,zb1o,ihalo_sor,new_comm)

!   CALL bounds_exch(1,'p',ff,'mo_tro 10')
   call sethaloN('p',ff,zff,ihalo_sor,new_comm)

!   CALL bounds_exch(1,'u+',uf,'mo_tro 10')
   CALL sethaloN('u+',uf,zuf,ihalo_sor,new_comm)

!   CALL bounds_exch(1,'vf+',vf,'mo_tro 10')
   CALL sethaloN('v+',vf,zvf,ihalo_sor,new_comm)


!    write(0,*)'my_rank=',my_rank,'p_ioff=',p_ioff
!    write(0,*)'my_rank=',my_rank,'p_joff=',p_joff
!    write(0,*)'icontro =',icontro

   IF (rtsorpar .GT. 0._wp) THEN
     sorpar = rtsorpar
     IF (p_pe==p_io) WRITE(0,*)'attn: user definded sorpar: ',sorpar
   ENDIF

   sorpai = 1._wp - sorpar

   DO kiter=1,iter_sor

     IF (kiter.GT.iter_sor-iter_sor_hack) THEN
       IF (rtsorpar_hack .GT. 0._wp) THEN
         sorpar=rtsorpar_hack
         sorpai = 1._wp - sorpar
       ENDIF
     ENDIF


       zurr = 0._wp

       ! SOR with red/black splitting
       ! 1step: update of the "red" points

       IF ( icontro.NE.0 ) THEN
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j,2)+2
             DO i=is,ie-1,2
                zalt(i,j)=z1o(i,j)
             ENDDO
          ENDDO
        ENDIF


       DO j=jb,je+ihm-1
          is = MOD(p_ioff+p_joff+j,2)+2
!             write(0,*)'my_rank=',my_rank,'red is =',is
          DO i=is,ie+ihm-1,2
             zz1o(i,j) = sorpar*zff(i,j)*(zb1o(i,j)                      &
                  +  zuf(i,j)*zz1o(i+1,j) + zuf(i-1,j)*zz1o(i-1,j)       &
                  +  zvf(i,j)*zz1o(i,j+1) + zvf(i,j-1)*zz1o(i,j-1))      &
                  + sorpai*zz1o(i,j)
          ENDDO
       ENDDO

       if ( icontro.NE.0 ) THEN
          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo_sor-1,j+ihalo_sor-1)
             enddo
          enddo
!
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j,2)+2
             DO i=is,ie-1,2
                zurr=zurr+(zalt(i,j)-z1o(i,j))**2
             enddo
          enddo
        endif

!------------------------------ Exchange -------------------------------

       if(imod.eq.1 .and. ihalo_sor.eq.1 ) then
          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo_sor-1,j+ihalo_sor-1)
             enddo
          enddo

          CALL bounds_exch(1,'p',z1o,'mo_tro 10')
          call sethaloN('p',z1o,zz1o,ihalo_sor,new_comm)

       endif
!------------------------------ Exchange -------------------------------

       ! 2step: update of the "black" points

       if ( icontro.NE.0 ) THEN
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j+1,2)+2
             DO i=is,ie-1,2
                zalt(i,j)=z1o(i,j)
             ENDDO
          ENDDO
       endif


       DO j=jb,je+ihm-1
          is = MOD(p_ioff+p_joff+j+1,2)+2
!             write(0,*)'my_rank=',my_rank,'black is =',is
          DO i=is,ie+ihm-1,2
             zz1o(i,j) = sorpar*zff(i,j)*(zb1o(i,j)                      &
                  +  zuf(i,j)*zz1o(i+1,j) + zuf(i-1,j)*zz1o(i-1,j)       &
                  +  zvf(i,j)*zz1o(i,j+1) + zvf(i,j-1)*zz1o(i,j-1))      &
                  + sorpai*zz1o(i,j)
          ENDDO
       ENDDO

       if ( icontro.NE.0)  THEN
          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo_sor-1,j+ihalo_sor-1)
             enddo
          enddo
!
          DO j=jb,je-1
             is = MOD(p_ioff+p_joff+j+1,2)+2
             DO i=is,ie-1,2
                zurr=zurr+(zalt(i,j)-z1o(i,j))**2
             enddo
          enddo
       endif


! ATTN : these can be reduced for larger halos than 2
!       if ( (kiter.eq.1).or.( mod(kiter,2).eq.0) )  then
       if ( MOD(kiter,imod).eq.0)  then

          do i=1,ie
             do j=1,je
                z1o(i,j)=zz1o(i+ihalo_sor-1,j+ihalo_sor-1)
             enddo
          enddo
!HH       for strange reasons the next call to bounds_exch is still needed
!HH       in order to work with all combinations of ihalo_sor and imod=ihalo_sor/2

          CALL bounds_exch(1,'p',z1o,'mo_tro 10')
          call sethaloN('p',z1o,zz1o,ihalo_sor,new_comm)

       endif


    ENDDO   !kiter
    IF (icontro.NE.0) THEN
    WRITE(0,*)'zurr',zurr,'sorpar',sorpar,'p_pe',p_pe
    ENDIF

    deallocate(zz1o,zff,zb1o,zuf,zvf)

  END SUBROUTINE troneu2



  SUBROUTINE OCVTRO

      USE MO_PARAM1
      USE mo_planetary_constants, ONLY: g
      USE MO_PARALLEL

      USE MO_COMMO1
      !USE MO_UNITS, only: io_stdout

      REAL(wp) :: DTHI,DTH

      INTEGER :: I,J

      dthi = 0.5_wp * dti
      dth  = 0.5_wp * dt


      IF ( .NOT. lwith_barotropic_stokes_drift ) THEN

        DO j=1,je
          DO i=1,ie
            uzo(i,j)=uzo(i,j)*deutio(i,j)
            vze(i,j)=vze(i,j)*deutie(i,j)
          END DO
        END DO

        DO j=2,je1
          DO i=2,ie1
            uso(i,j)=uzo(i,j)                                            &
                 +g*stabn*dt*(z1o(i,j)-z1o(i+1,j))                       &
                 *amsuo(i,j,1)/dlxu(i,j)
            vse(i,j)=vze(i,j)                                            &
                 +g*stabn*dt*(z1o(i,j+1)-z1o(i,j))                       &
                 *amsue(i,j,1)/dlyv(i,j)
          END DO
        END DO

      ELSE

        DO j=2,je1
          DO i=2,ie1
            uso(i,j)=uso(i,j)                                            &
                 +g*stabn*dt*(z1o(i,j)-z1o(i+1,j))                       &
                 *amsuo(i,j,1)/dlxu(i,j)
            vse(i,j)=vse(i,j)                                            &
                 +g*stabn*dt*(z1o(i,j+1)-z1o(i,j))                       &
                 *amsue(i,j,1)/dlyv(i,j)

          END DO
        END DO

      ENDIF


!#ifdef bounds_exch_save
       CALL bounds_exch(1,'u',USO,'ocvtro 1')
       CALL bounds_exch(1,'v',VSE,'ocvtro 2')
!#endif

!#ifdef bounds_exch_save
      CALL bounds_exch(1,'p',Z1O,'ocvtro 3')
!#endif

      END SUBROUTINE OCVTRO



      SUBROUTINE update_zo

        !   USE MO_COMMO1, ONLY: ZO, Z1O
        USE MO_COMMO1, ONLY: ZO, WO, dt
#ifdef PBGC
        USE mo_param1_bgc, ONLY: nocetra, ntraad
        USE mo_carbch, ONLY: ocetra
#endif
#ifdef PBGC
        INTEGER :: l
#endif

#ifdef PBGC
        CALL dilcor_gtrf2
#endif

! emr suggest this formulation for more stable solutions

        ZO(:,:)=Zo(:,:)+WO(:,:,1)*dt

#ifdef PBGC
        DO l=ntraad+1,nocetra
          CALL dilcor_ptrf2(ocetra(1,1,1,l))
        ENDDO

#endif

      END SUBROUTINE update_zo

SUBROUTINE alloc_mem_dilcor

      USE MO_PARAM1, ONLY: ie,je

      ALLOCATE(susalo(ie,je),suwath(ie,je))
      ALLOCATE(dilcor(ie,je))


END SUBROUTINE alloc_mem_dilcor

SUBROUTINE dilcor_gtrf

      USE mo_param1,ONLY: ie,je
      USE mo_commo1,ONLY: sao,zo,ddpo,sictho        &
      ,sicsno
      USE MO_COMMOAU1,ONLY: sice
      USE mo_planetary_constants, ONLY: rhoicwa,rhosnwa
      INTEGER :: i,j

!  WRITE(0,*)'in dilcor_pre',size(susalo(:,1)),size(susalo(1,:))

      DO J=1,JE
        DO I=1,IE
           susalo(i,j)=(SAO(I,J,1)*(ZO(I,J)+ddpo(i,j,1)             &
                -SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)           &
                +SICE*SICTHO(I,J)*RHOICWA)
           suwath(i,j)=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
                -sicsno(i,j)*rhosnwa

        ENDDO
     ENDDO

END SUBROUTINE dilcor_gtrf

SUBROUTINE dilcor_gtrf2

      USE mo_param1,ONLY: ie,je
      USE mo_commo1,ONLY: zo,ddpo,sictho        &
      ,sicsno
      USE mo_planetary_constants, ONLY :rhoicwa,rhosnwa
      INTEGER :: i,j

!  WRITE(0,*)'in dilcor_pre',size(susalo(:,1)),size(susalo(1,:))

      DO J=1,JE
        DO I=1,IE
           suwath(i,j)=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
                -sicsno(i,j)*rhosnwa
        ENDDO
     ENDDO

END SUBROUTINE dilcor_gtrf2

SUBROUTINE dilcor_ptrf

  USE mo_param1,ONLY: ie,je
  USE mo_commo1,ONLY: sao,zo,ddpo,sictho        &
       ,sicsno,lweto
  USE mo_planetary_constants, ONLY :rhoicwa,rhosnwa
  USE MO_COMMOAU1,ONLY: sice

  INTEGER :: i,j
  REAL(wp) :: wathne

!  WRITE(0,*)'in dilcor_post',size(dilcor(:,1)),size(dilcor(1,:))

  dilcor(:,:) = 1._wp
  DO j=1,je
     DO i=1,ie
        IF(lweto(i,j,1))THEN

           wathne=zo(i,j)+ddpo(i,j,1)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa

           sao(i,j,1)=(susalo(i,j)-sice*sictho(i,j)*rhoicwa)/wathne

        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE dilcor_ptrf

SUBROUTINE dilcor_ptrf2(TRF)

  USE mo_param1,ONLY: ie,je
  USE mo_commo1,ONLY: zo,ddpo,sictho        &
       ,sicsno,lweto
  USE mo_planetary_constants, ONLY :rhoicwa,rhosnwa
  INTEGER :: i,j
  REAL(wp) :: wathne
  REAL(wp) :: trf(ie,je,1)

!  WRITE(0,*)'in dilcor_post',size(dilcor(:,1)),size(dilcor(1,:))

  dilcor(:,:) = 1._wp
  DO j=1,je
     DO i=1,ie
        IF(lweto(i,j,1))THEN
           wathne=zo(i,j)+ddpo(i,j,1)-sictho(i,j)*rhoicwa-sicsno(i,j)*rhosnwa
           dilcor(i,j)=suwath(i,j)/wathne
        ENDIF
     ENDDO
  ENDDO

  DO j=1,je
     DO i=1,ie
        trf(i,j,1)=trf(i,j,1)*dilcor(i,j)
     ENDDO
  ENDDO

END SUBROUTINE dilcor_ptrf2

SUBROUTINE correct_zo

      USE MO_PARAM1
      USE MO_COMMO1
      USE MO_PARALLEL

      USE MO_UNITS, ONLY: io_stdout
#ifdef PBGC
      USE MO_PARAM1_BGC,only: nocetra
      USE MO_CARBCH,only: ocetra
#endif
      INTEGER :: I,J
#ifdef PBGC
      INTEGER :: l
#endif
      REAL(wp),ALLOCATABLE  :: ZO_g(:,:), Z1O_g(:,:)

      REAL(wp) SUGG, SUZZ, SUZZ1, ZQQ, ZQQ1

      IF (p_pe==p_io)then
         ALLOCATE(ZO_g(ie_g,je_g), Z1O_g(ie_g,je_g))
      ELSE
         ALLOCATE(ZO_g(0,0), Z1O_g(0,0))
      ENDIF
#ifdef PBGC
      call dilcor_gtrf2
#endif

! SET SEA LEVEL (ZO and Z1O) BACK TO GLOBAL ZERO. Done on one processor
! to ensure reproducibility.


!      CALL gather_arr(ZO,ZO_g,p_io)
!      CALL gather_arr(Z1O,Z1O_g,p_io)
      CALL gather(ZO,ZO_g,p_io)
      CALL gather(Z1O,Z1O_g,p_io)

      IF (p_pe==p_io) THEN
         sugg = 0.0_wp
         suzz = 0.0_wp
         suzz1 = 0.0_wp
         DO  J=1,je_g
            DO  I=2,ie_g-1
               IF(lwetol1_g(i,j)) THEN
                  IF (zo_g(i, j) .LT. - 3._wp) THEN
                     WRITE(IO_STDOUT,*) 'ZO LT -3! at i=',i,'j=',j
                  ENDIF
                  SUGG=SUGG+DLXP_g(I,J)*DLYP_g(I,J)
                  SUZZ=SUZZ+DLXP_g(I,J)*DLYP_g(I,J)*ZO_g(I,J)
                  SUZZ1=SUZZ1+DLXP_g(I,J)*DLYP_g(I,J)*Z1O_g(I,J)
               ENDIF
            ENDDO
         ENDDO
         ZQQ=SUZZ/SUGG
         ZQQ1=SUZZ1/SUGG
         WRITE(IO_STDOUT,*) ' MEAN ZETA ', ZQQ,ZQQ1
         DO J=1,je_g
            DO I=1,ie_g
               ZO_g(I,J)=(ZO_g(I,J)-ZQQ)*WETOl1_g(I,J)
               Z1O_g(I,J)=(Z1O_g(I,J)-ZQQ1)*WETOL1_g(I,J)
            ENDDO
         ENDDO
      END IF
      CALL scatter(ZO_g,ZO,p_io)
      CALL scatter(Z1O_g,Z1O,p_io)
      CALL p_bcast(ZQQ,p_io)
      CALL bounds_exch(1,'p',ZO,'mpiom 32')
      CALL bounds_exch(1,'p',Z1O,'mpiom 33')

#ifdef PBGC
      do l=1,nocetra
         call dilcor_ptrf2(ocetra(1,1,1,l))
      enddo
#endif/*PBGC*/

      deallocate(ZO_g,Z1O_g)

END SUBROUTINE correct_zo

END MODULE mo_tro






