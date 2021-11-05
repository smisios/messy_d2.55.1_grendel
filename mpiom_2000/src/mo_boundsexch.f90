MODULE MO_BOUNDSEXCH

  USE mo_kind,   ONLY: dp, wp
  USE mo_param1, ONLY: kep, ke, ie, je, ie_g, je_g, icycli
  USE mo_mpi
  USE mo_parallel
  USE mo_commo1, ONLY : lbounds_exch_tp,lmpitype,lnonblock
  use mpi, only: mpi_proc_null

  IMPLICIT NONE

  REAL(dp), PRIVATE :: t2d=0.0_dp, t3d=0.0_dp, ts, te
  INTEGER,  PRIVATE :: n2d=0, n3d=0

  logical :: lold

!   send/recv buffers for the 2d case
  REAL(wp), ALLOCATABLE :: xyr1(:),xyr2(:)         ! recv buffers
  REAL(wp), ALLOCATABLE :: xys1(:),xys2(:)         ! send buffers

  REAL(wp), ALLOCATABLE :: buffer(:,:,:)        ! send/recv buffer for northern boundary


  INTERFACE bounds_exch

    MODULE PROCEDURE bounds_exch_halo_2d
    MODULE PROCEDURE bounds_exch_halo_3d

  END INTERFACE

CONTAINS

  !>
  !! Allocate auxiliary variables for bounds exchange.
  !!
  SUBROUTINE alloc_mem_bounds_exch_halo

    INTEGER :: i_up,i_lo,i_size,j_size

    IF( .NOT. lmpitype  .OR. nprocx == 1 .OR. p_pe == nprocxy .OR. lbounds_exch_tp ) THEN   !use mpitypes/sendrecv buffers

      j_size=je*(ke+1)*ihalo_max

      i_up=ie+ihalo_max-1
      i_lo=2-ihalo_max
      i_size=MAX(((i_up-i_lo+1)*(ke+1)*ihalo_max+1), j_size)


      ALLOCATE ( xyr1(i_size), xyr2(i_size), xys1(i_size), xys2(i_size) )

      xyr1(:) = 0._wp
      xyr2(:) = 0._wp
      xys1(:) = 0._wp
      xys2(:) = 0._wp


   ENDIF

    IF ( lbounds_exch_tp ) THEN  !   send/recv buffer for northern boundary

      i_up=ie+ihalo_max-1
      i_lo=2-ihalo_max

      ALLOCATE ( buffer(i_lo:i_up,1:ihalo_max,1:ke+1) )

      buffer(:,:,:) = 0._wp

    ENDIF

  END SUBROUTINE alloc_mem_bounds_exch_halo

  !>
  !! Exchange boundary rows for 2D variables, according to grid mode.
  !!
  SUBROUTINE bounds_exch_halo_2d(ihalo,ttt,a0,text)


#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

    ! Exchanges boundaries of 2D arrays
    ! adopted to support tripolar type grid (hh, 01/2006)
    ! adopted to larger halos (hh, 09/2009)

#ifdef bounds_exch_check
    INTEGER nn,jj
    REAL(wp), ALLOCATABLE :: aa0(:,:),aa1(:,:)
#endif

    INTEGER :: ihalo ! halo width
    INTEGER :: i_up  ! upper i-index
    INTEGER :: i_lo  ! lower i-index
    INTEGER :: j_up  ! upper j-index
    INTEGER :: j_lo  ! lower j-index

    INTEGER :: itype
    INTEGER :: xup,xlo,yup,ylo,i_size,j_size

    REAL(wp), INTENT(INOUT) :: a0(:,:)
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER nm, np, ii


    INTEGER  :: is, ir, il, ne, ih                                   ! for northern boundary
    REAL(wp) :: multip

    ts = p_time()

#ifdef _PROFILE
    CALL trace_start ('bounds_exch_2d total', 2)
    CALL trace_start ('bounds_exch_2d ew/ns', 3)
#endif

!   calculate field bounds incl. halos
    i_up=ie+ihalo-1
    i_lo=2-ihalo
    j_up=je+ihalo-1
    j_lo=2-ihalo
    i_size=i_up-i_lo+1
    j_size=j_up-j_lo+1

#ifdef bounds_exch_check
    ALLOCATE(aa0(i_lo:i_up,j_lo:j_up),aa1(i_lo:i_up,j_lo:j_up))
    aa0(:,:) = a0(:,:)
#endif

!    lold=lmpitype
!    lmpitype=.false.

    ! x-direction

    IF( ( .NOT. lmpitype ) .OR. nprocx == 1 .OR. p_pe == nprocxy ) THEN       ! prepare sendrecv buffers if needed

      DO ii=1,ihalo

        xlo=((ii-1)*je)+1
        xup=ii*je

        xys1(xlo:xup) = a0(1+ii,1:je)          ! prepare left send buffer
        xys2(xlo:xup) = a0(ie-1-ihalo+ii,1:je) ! prepare right send buffer

      ENDDO

    ENDIF

    ! only do mpi communication if needed
    IF ( nprocx > 1 .AND. p_pe < nprocxy ) THEN                        ! mpi communication needed

      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

#ifdef NOLAND

      IF (iproc(p_pe,1).NE.0) THEN ! do not communicate if neighbour p_pe is pure land patch

        IF (iproc(nm,1).EQ.0) nm=MPI_PROC_NULL
        IF (iproc(np,1).EQ.0) np=MPI_PROC_NULL

      ELSE                          !local p_pe is land only

        nm=MPI_PROC_NULL
        np=MPI_PROC_NULL

      ENDIF

#endif

      IF ( .NOT. lmpitype ) THEN            ! use sendrecv buffers

        xup=ihalo*je

        IF ( lnonblock ) THEN               ! use nonblocking communication

          CALL p_isend(xys1(1:xup),nm,1) ! send to left
          CALL p_irecv(xyr2(1:xup),nm,2) ! recv from left
          CALL p_isend(xys2(1:xup),np,2) ! send to right
          CALL p_irecv(xyr1(1:xup),np,1) ! recv from right

          CALL p_wait

        ELSE                                ! use blocking communication

          CALL p_sendrecv(xys1(1:xup),nm,xyr1(1:xup),np,1) ! send to left, recv from right
          CALL p_sendrecv(xys2(1:xup),np,xyr2(1:xup),nm,2) ! send to right, recv from left

        ENDIF

      ELSE                                             ! use mpitypes

        itype=ew_boundary(ihalo)                       ! help variable

        IF ( lnonblock ) THEN                          ! use nonblocking communication

          CALL p_isend(itype,a0(2,1),nm,1)             ! send to left
          CALL p_irecv(itype,a0(2-ihalo,1),nm,2)       ! recv from right
          CALL p_isend(itype,a0(ie-ihalo,1),np,2)      ! send to right
          CALL p_irecv(itype,a0(ie,1),np,1)            ! recv from left

          CALL p_wait

        ELSE                                           ! use blocking communication

          CALL p_sendrecv(itype,a0(2,1),nm,a0(ie,1),np,1)             ! send to left, recv from right
          CALL p_sendrecv(itype,a0(ie-ihalo,1),np,a0(2-ihalo,1),nm,2) ! send to right, recv from left

        ENDIF


      ENDIF   !use mpitypes/sendrecv buffers

    ELSE                      ! nprocx == 1 => no mpi communication needed

      xup=ihalo*je

      xyr1(1:xup)=xys1(1:xup) ! copy left send buffer to right recv buffer
      xyr2(1:xup)=xys2(1:xup) ! copy right send buffer to left recv buffer

    ENDIF

    ! update halos in x direction ; needed for (.not. lmpitype)
    ! also needed if  nprocx == 1 or in serial test mode p_pe == nprocxy

    IF( ( .NOT. lmpitype ) .OR. nprocx == 1 .OR. p_pe == nprocxy ) THEN   !use sendrecv buffers

      DO ii=1,ihalo

        xlo=((ii-1)*je)+1
        xup=ii*je

        IF( icycli /= 0 .OR. ( .NOT. have_g_is ) ) THEN  ! update left halo

          a0(1-ihalo+ii,1:je)=xyr2(xlo:xup)

        ENDIF

        IF( icycli /= 0 .OR. ( .NOT. have_g_ie ) ) THEN  ! update right halo

          a0(ie-1+ii,1:je)=xyr1(xlo:xup)

        ENDIF


      ENDDO


    ENDIF !use sendrecv buffers


!    lmpitype=lold

    ! y-direction

!    lold=lmpitype
!    lmpitype=.false.

    ! Note that there is no action required if nprocy==1

    IF ( nprocy > 1 .AND. p_pe < nprocxy ) THEN

      IF ( .NOT. lmpitype ) THEN        ! use sendrecv buffers

        DO ii=1,ihalo

          ylo=(ii-1)*i_size+1
          yup=ii*i_size

          xys1(ylo:yup) = a0(i_lo:i_up,ii+1)    ! prepare upper send buffer
          xys2(ylo:yup) = a0(i_lo:i_up,je-ii)   ! prepare lower send bufffer

        ENDDO

      ENDIF

      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions

      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      IF ( have_g_js ) nm=MPI_PROC_NULL
      IF ( have_g_je ) np=MPI_PROC_NULL


#ifdef NOLAND

      IF (iproc(p_pe,1).NE.0) THEN ! do not communicate if neighbour p_pe is pure land patch

        IF (iproc(nm,1).EQ.0) nm=MPI_PROC_NULL
        IF (iproc(np,1).EQ.0) np=MPI_PROC_NULL

      ELSE                          !local p_pe is land only

        nm=MPI_PROC_NULL
        np=MPI_PROC_NULL

      ENDIF

#endif

      IF ( .NOT. lmpitype ) THEN                !use sendrecv buffers

        yup=ihalo*i_size

        IF ( lnonblock ) THEN                ! use nonblocking communication

          CALL p_isend(xys1(1:yup),nm,1) ! send to up
          CALL p_irecv(xyr2(1:yup),nm,2) ! recv from up
          CALL p_isend(xys2(1:yup),np,2) ! send to down
          CALL p_irecv(xyr1(1:yup),np,1) ! recv from down
          CALL p_wait

        ELSE                               ! use blocking communication

          CALL p_sendrecv(xys1(1:yup),nm,xyr1(1:yup),np,1) ! send to up , recv from down
          CALL p_sendrecv(xys2(1:yup),np,xyr2(1:yup),nm,2) ! send to down , recv from up

        ENDIF


        DO ii=1,ihalo

          ylo=(ii-1)*i_size+1
          yup=ii*i_size

          IF ( .NOT. have_g_js ) THEN

            a0(i_lo:i_up,2-ii) = xyr2(ylo:yup) ! update upper halo

          ENDIF

          IF ( .NOT. have_g_je ) THEN

            a0(i_lo:i_up,je-1+ii) = xyr1(ylo:yup) ! update lower halo

          ENDIF

        ENDDO


      ELSE              ! use mpitypes

        itype=ns_boundary(ihalo)

        IF ( lnonblock ) THEN    ! use nonblocking communication

          IF ( .NOT. have_g_js ) THEN
            CALL p_isend(itype,a0(i_lo,2),nm,1) ! send to up
            CALL p_irecv(itype,a0(i_lo,j_lo),nm,2) ! recv from up
          ENDIF

          IF ( .NOT. have_g_je ) THEN
            CALL p_isend(itype,a0(i_lo,je-ihalo),np,2) ! send to down
            CALL p_irecv(itype,a0(i_lo,je),np,1) ! recv from down
          ENDIF

          CALL p_wait

        ELSE    ! use blocking communication

          IF ( have_g_js ) THEN

            CALL p_sendrecv(itype,a0(i_lo,je-ihalo),np,a0(i_lo,je),np,3)

          ELSEIF ( have_g_je ) THEN

            CALL p_sendrecv(itype,a0(i_lo,2),nm,a0(i_lo,j_lo),nm,3)

          ELSE

            CALL p_sendrecv(itype,a0(i_lo,2),nm,a0(i_lo,je),np,3)
            CALL p_sendrecv(itype,a0(i_lo,je-ihalo),np,a0(i_lo,j_lo),nm,3)

          ENDIF

        ENDIF

      ENDIF ! lmpitype

    ELSE ! nprocy = 1

      ! nothing needed

    ENDIF ! nprocy = 1


!    lmpitype=lold

#ifdef _PROFILE
    CALL trace_stop ('bounds_exch_2d ew/ns', 3)
    CALL trace_start ('bounds_exch_2d tp', 4)
#endif


    ! make the northern margin cyclic

    IF( lbounds_exch_tp .AND. have_g_js) THEN

      ! define send buffer, only use the ihalo+1 needed lines

      SELECT CASE (ttt)

      CASE('p','p-','u','u+')     !! lines 3 to ihalo+3

        DO ii = 1,ihalo+1

          ylo=(ii-1)*i_size+1
          yup=ii*i_size

          xys1(ylo:yup)      = a0(i_lo:i_up,ii+2)    ! prepare upper send buffer

        ENDDO

      CASE('v','v+','vf','vf+','s','s-')     !! 2 to ihalo+2

        DO ii = 1,ihalo+1

          ylo=(ii-1)*i_size+1
          yup=ii*i_size

          xys1(ylo:yup)      = a0(i_lo:i_up,ii+1)    ! prepare upper send buffer

        ENDDO


      CASE('uu','vv')  ! nothing to be done

      ! jump to the end

      CASE default !! unsupported value for ttt

        CALL STOP_ALL(' bounds_exch : argument unsupported ')

      END SELECT



! ---------------------------------------------------------------------------------

      ne=(nprocx-1)-p_pe                                      ! find the neighbour cpu
      yup=(ihalo+1)*i_size

      IF ( p_pe /= ne ) THEN

        IF ( lnonblock ) THEN

          CALL p_isend(xys1(1:yup),ne,6)                ! exchange with neighbour cpu
          CALL p_irecv(xyr1(1:yup),ne,6)                ! exchange with neighbour cpu

          CALL p_wait

        ELSE

          CALL p_sendrecv(xys1(1:yup),ne,xyr1(1:yup),ne,6)  ! exchange with neighbour cpu

        ENDIF

      ELSE                                                 ! dont do send/recive on the same cpu

        xyr1(1:yup)=xys1(1:yup)

      ENDIF


      DO ii=1,ihalo+1                ! unpack the recive buffer into a help array

        ylo=(ii-1)*i_size+1
        yup=ii*i_size

        buffer(i_lo:i_up,ii,1) = xyr1(ylo:yup) ! update upper halo

      ENDDO

! ---------------------------------------------------------------------------------

      SELECT CASE (ttt)

      CASE('p','p-')                                                     !! p point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='p-')

        DO ii=1,ihalo+1                                                  ! all needed halos
          DO is=1,i_size                                                 ! size from i_lo:i_up

            il=i_lo-1+is                                                 ! left index
            ir=i_up+1-is                                                 ! right index

            a0(il,3-ii) = multip *  buffer(ir,ii,1)                      ! line (-)3 -> 2 ; (-)4 -> 1

          END DO
        END DO

      CASE('v','v+','vf','vf+')                                          !! v point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='v'.or.ttt=='vf')

        DO ii=2,ihalo+1                                                  ! all needed halos
          DO is=1,i_size                                                 ! size from i_lo:i_up

            il=i_lo-1+is                                                 ! left index
            ir=i_up+1-is                                                 ! right index

            a0(il,3-ii) =  multip * buffer(ir,ii,1)                      ! line (-)3 -> 1

          END DO
        END DO


        IF (ttt=='vf'.OR.ttt=='vf+') THEN                                !! v point with/without sign change

          ih=0
          IF ( p_pe <= (nprocx-1)/2 .OR. p_pe == nprocxy ) ih=i_size
          IF ( p_pe == ne ) ih =i_size/2

          DO is=1,ih                                                   ! size from i_lo:i_up

            il=i_lo-1+is                                               ! left index
            ir=i_up+1-is                                               ! right index

            a0(il,2) =  multip * buffer(ir,1,1)                        ! line (-)2 -> 2

          END DO

        ENDIF

      CASE('u','u+')                                                     !! u point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='u')

        DO ii=1,ihalo+1                                                  ! all needed halos
          DO is=1,i_size-1                                               ! size from i_lo:i_up-1

            il=i_lo-1+is                                                 ! left index
            ir=i_up-is                                                   ! right index

            a0(il,3-ii) =  multip * buffer(ir,ii,1)                      ! line (-)3 -> 2 ; (-)4 -> 1

          END DO
        END DO

      CASE('s','s-','sf','sf-')                                          !! psi point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='s-'.OR.ttt=='sf-')

        DO ii=2,ihalo+1                                                  ! all needed halos
          DO is=1,i_size-1                                               ! size from i_lo:i_up-1

            il=i_lo-1+is                                                 ! left index
            ir=i_up-is                                                   ! right index

            a0(il,3-ii) = multip * buffer(ir,ii,1)                       ! line (-)3 -> 1

          END DO
        END DO

        IF (ttt=='sf'.OR.ttt=='sf-') THEN                                !! v point with/without sign change

          ih=0
          IF ( p_pe <= (nprocx-1)/2 .OR. p_pe == nprocxy ) ih=i_size
          IF ( p_pe == ne ) ih =i_size/2

          DO is=1,i_size-1                                             ! size from i_lo:i_up-1

            il=i_lo-1+is                                               ! left index
            ir=i_up-is                                                 ! right index

            a0(il,2) = multip * buffer(ir,1,1)                         ! line (-)2 -> 2

          END DO

        ENDIF

      END SELECT

    ENDIF     ! end of the northern margin treatment

    te = p_time()
    t2d = t2d + te-ts
    n2d = n2d + 1

#ifdef bounds_exch_check

    aa1(:,:) = a0(:,:)
    nn=0
    DO jj=1,je
       DO ii=1,ie
          IF (aa1(ii,jj).NE.aa0(ii,jj)) THEN
            nn=nn+1
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(nn)
    IF ( p_pe .EQ. p_io ) THEN
       IF ( nn .EQ. 0 )   WRITE(0,*) 'attn: pe ',p_pe,' of ',nprocxy,' needless bounds_exch! ',text
    ENDIF
    DEALLOCATE(aa0,aa1)

#endif


    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_2d(a0,text)
      ELSE
        CALL para_check_2d(a0,'bounds_exch_2d')
      ENDIF
    ENDIF

#ifdef _PROFILE
  CALL trace_stop ('bounds_exch_2d tp', 4)
  CALL trace_stop ('bounds_exch_2d total', 2)
#endif

  END SUBROUTINE bounds_exch_halo_2d


  SUBROUTINE bounds_exch_halo_3d(ihalo,ttt,a0,text)

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

    ! Exchanges boundaries of 2D arrays
    ! adopted to support tripolar type grid (hh, 01/2006)
    ! adopted to larger halos (hh, 09/2009)

#ifdef bounds_exch_check
    INTEGER NN,jj
    REAL(wp),ALLOCATABLE :: aa0(:,:,:),aa1(:,:,:)
#endif

    INTEGER :: ihalo ! single halo width
    INTEGER :: i_up  ! upper i-index
    INTEGER :: i_lo  ! lower i-index
    INTEGER :: j_up  ! upper j-index
    INTEGER :: j_lo  ! lower j-index

    INTEGER :: itype
    INTEGER :: xlo,ylo,xup,yup,i_size,j_size

    REAL(wp), INTENT(INOUT) :: a0(:,:,:)
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    INTEGER kk, nm, np, ii, k

    INTEGER  :: is, ir, il, ne, ih                                   ! for northern boundary
    REAL(wp) :: multip

#ifdef _PROFILE
    CALL trace_start ('bounds_exch_3d total', 5)
    CALL trace_start ('bounds_exch_3d ew/ns', 6)
#endif

    ts = p_time()

    kk = UBOUND(a0,3)

!   calculate field bounds incl. halos
    i_up=ie+ihalo-1
    i_lo=2-ihalo
    j_up=je+ihalo-1
    j_lo=2-ihalo
    i_size=i_up-i_lo+1
    j_size=j_up-j_lo+1

#ifdef bounds_exch_check
    ALLOCATE(aa0(i_lo:i_up,j_lo:j_up,kk),aa1(i_lo:i_up,j_lo:j_up,kk))
    aa0(:,:,:) = a0(:,:,:)
#endif


!    lold=lmpitype
!    lmpitype=.false.


    ! x-direction

    ! prepare sendrecv buffers if needed

    IF( ( .NOT. lmpitype) .OR. nprocx == 1 .OR. p_pe == nprocxy ) THEN   !use sendrecv buffers

      DO ii=1,ihalo

        xlo=(ii-1)*je*kk+1
        xup=ii*je*kk

!        WRITE(0,*)p_pe,'xlo',xlo,'xup',xup,'ii',ii,'je',je,'kk',kk
        xys1(xlo:xup) = RESHAPE(a0(1+ii,1:je,1:kk), (/je*kk/))  ! prepare left send buffer
        xys2(xlo:xup) = RESHAPE(a0(ie-1-ihalo+ii,1:je,1:kk), (/je*kk/)) ! prepare right send buffer

      ENDDO

    ENDIF

    ! do mpi communication if needed
    IF ( nprocx > 1 .AND. p_pe < nprocxy ) THEN                        ! mpi communication needed

      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

#ifdef NOLAND

      IF (iproc(p_pe,1).NE.0) THEN ! do not communicate if neighbour p_pe is pure land patch

        IF (iproc(nm,1).EQ.0) nm=MPI_PROC_NULL
        IF (iproc(np,1).EQ.0) np=MPI_PROC_NULL

      ELSE                          !local p_pe is land only

        nm=MPI_PROC_NULL
        np=MPI_PROC_NULL

      ENDIF

#endif

      IF ( .NOT. lmpitype ) THEN
        ! use sendrecv buffers

        xup=ihalo*je*kk

        IF (lnonblock) THEN           ! use nonblocking communication

          CALL p_isend(xys1(1:xup),nm,1) ! send to left
          CALL p_irecv(xyr2(1:xup),nm,2) ! recv from left
          CALL p_isend(xys2(1:xup),np,2) ! send to right
          CALL p_irecv(xyr1(1:xup),np,1) ! recv from right

          CALL p_wait

        ELSE           ! use blocking communication

          CALL p_sendrecv(xys1(1:xup),nm,xyr1(1:xup),np,1) ! send to left, recv from right
          CALL p_sendrecv(xys2(1:xup),np,xyr2(1:xup),nm,2) ! send to right, recv from left

        ENDIF

      ELSE    ! use mpitypes

        IF (kk.EQ.ke) itype=ew_boundary_ke(ihalo)
        IF (kk.EQ.kep) itype=ew_boundary_kep(ihalo)


        IF (lnonblock) THEN           ! use nonblocking communication

          CALL p_isend(itype,a0(2,1,1),nm,1)  ! send to left
          CALL p_irecv(itype,a0(2-ihalo,1,1),nm,2) ! recv from right
          CALL p_isend(itype,a0(ie-ihalo,1,1),np,2) ! send to right
          CALL p_irecv(itype,a0(ie,1,1),np,1)  ! recv from left

          CALL p_wait

        ELSE ! use blocking communication

          CALL p_sendrecv(itype,a0(2,1,1),nm,a0(ie,1,1),np,2)  ! send to left, recv from right
          CALL p_sendrecv(itype,a0(ie-ihalo,1,1),np,a0(2-ihalo,1,1),nm,2) ! send to right, recv from left

        ENDIF

      ENDIF

    ELSE                       ! nprocx == 1 ; no mpi communication but copy needed

      xup=ihalo*je*kk

      xyr1(1:xup)=xys1(1:xup) ! copy left send buffer to right recv buffer
      xyr2(1:xup)=xys2(1:xup) ! copy right send buffer to left recv buffer

    ENDIF

    ! update halos in x direction

    IF ( (.NOT. lmpitype) .OR. nprocx == 1 .OR. p_pe == nprocxy ) THEN   !use sendrecv buffers


      DO ii=1,ihalo

        xlo=(ii-1)*je*kk+1
        xup=ii*je*kk

        IF ( icycli/=0 .OR. .NOT. have_g_is ) THEN ! update left halo

          a0(1-ihalo+ii,1:je,1:kk)=RESHAPE(xyr2(xlo:xup), (/je,kk/))

        ENDIF

        IF ( icycli/=0 .OR. .NOT. have_g_ie ) THEN ! update right halo

          a0(ie-1+ii,1:je,1:kk)=RESHAPE(xyr1(xlo:xup), (/je,kk/))

        ENDIF

      ENDDO



    ENDIF !use sendrecv buffers

!    lold=lmpitype=lold

    ! y-direction

!    lold=lmpitype
!    lmpitype=.false.

    ! Note that there is no action required if nprocy==1

    IF ( nprocy > 1 .AND. p_pe < nprocxy ) THEN

      IF ( .NOT. lmpitype ) THEN          ! use sendrecv buffers

        DO ii=1,ihalo

          ylo=(ii-1)*i_size*kk+1
          yup=ii*i_size*kk

          xys1(ylo:yup) = RESHAPE(a0(i_lo:i_up,ii+1,1:kk), (/i_size*kk/))    ! prepare upper send buffer
          xys2(ylo:yup) = RESHAPE(a0(i_lo:i_up,je-ii,1:kk), (/i_size*kk/))  ! prepare lower send bufffer

        ENDDO

      ENDIF

      ! Get processor numbers of neighbors
      ! For the sake of simplicity  the p_sendrecv call is as
      ! for periodic boundary conditions

      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      IF ( have_g_js ) nm=MPI_PROC_NULL
      IF ( have_g_je ) np=MPI_PROC_NULL

#ifdef NOLAND

      IF (iproc(p_pe,1).NE.0) THEN ! do not communicate if neighbour p_pe is pure land patch

        IF (iproc(nm,1).EQ.0) nm=MPI_PROC_NULL
        IF (iproc(np,1).EQ.0) np=MPI_PROC_NULL

      ELSE                          !local p_pe is land only

        nm=MPI_PROC_NULL
        np=MPI_PROC_NULL

      ENDIF

#endif

      IF ( .NOT. lmpitype ) THEN             !use sendrecv buffers

        yup=ihalo*i_size*kk

        IF ( lnonblock ) THEN            ! use nonblocking communication

          CALL p_isend(xys1(1:yup),nm,1)
          CALL p_irecv(xyr2(1:yup),nm,2)
          CALL p_isend(xys2(1:yup),np,2)
          CALL p_irecv(xyr1(1:yup),np,1)
          CALL p_wait

        ELSE           ! use blocking communication

          CALL p_sendrecv(xys1(1:yup),nm,xyr1(1:yup),np,1)
          CALL p_sendrecv(xys2(1:yup),np,xyr2(1:yup),nm,2)

        ENDIF


        DO ii=1,ihalo

          ylo=(ii-1)*i_size*kk+1
          yup=ii*i_size*kk

          IF(.NOT. have_g_js) THEN

            a0(i_lo:i_up,2-ii,1:kk) = RESHAPE(xyr2(ylo:yup),(/i_size,kk/)) ! update upper halo

          ENDIF

          IF(.NOT. have_g_je) THEN

            a0(i_lo:i_up,je-1+ii,1:kk) = RESHAPE(xyr1(ylo:yup), (/i_size,kk/)) ! update lower halo

          ENDIF

        ENDDO

      ELSE              ! use mpitypes

!!$ the 3d mpidatatypes do not work on our ibm ; as workaround we use 2d mpidatatypes
!!$        IF (kk.EQ.ke) itype=ns_boundary_ke(ihalo)
!!$        IF (kk.EQ.kep) itype=ns_boundary_kep(ihalo)
           itype=ns_boundary(ihalo)


        IF ( lnonblock ) THEN    ! use nonblocking communication

!!$ the 3d mpidatatypes do not work on our ibm ; as workaround we use 2d mpidatatypes
!!$          IF ( .NOT. have_g_js ) THEN
!!$            CALL p_isend(itype,a0(i_lo,2,1),nm,1)        ! send to up
!!$            CALL p_irecv(itype,a0(i_lo,j_lo,1),nm,2)     ! recv from up
!!$          ENDIF
!!$
!!$          IF ( .NOT. have_g_je ) THEN
!!$            CALL p_isend(itype,a0(i_lo,je-ihalo,1),np,2) ! send to down
!!$            CALL p_irecv(itype,a0(i_lo,je,1),np,1)       ! recv from down
!!$          ENDIF
!!$
!!$          CALL p_wait


           IF ( have_g_js ) nm=MPI_PROC_NULL
           IF ( have_g_je ) np=MPI_PROC_NULL

           DO k=1,kk
              CALL p_isend(itype,a0(i_lo,2,k),nm,1)        ! send to up
              CALL p_irecv(itype,a0(i_lo,j_lo,k),nm,2)     ! recv from up
              CALL p_isend(itype,a0(i_lo,je-ihalo,k),np,2) ! send to down
              CALL p_irecv(itype,a0(i_lo,je,k),np,1)       ! recv from down
              CALL p_wait
           ENDDO


        ELSE    ! use blocking communication

!!$ the 3d mpidatatypes do not work on our ibm ; as workaround we use 2d mpidatatypes
!!$          IF ( have_g_js ) THEN
!!$
!!$            CALL p_sendrecv(itype,a0(i_lo,je-ihalo,1),np,a0(i_lo,je,1),np,3)
!!$
!!$          ELSEIF ( have_g_je ) THEN
!!$
!!$            CALL p_sendrecv(itype,a0(i_lo,2,1),nm,a0(i_lo,j_lo,1),nm,3)
!!$
!!$          ELSE
!!$
!!$            CALL p_sendrecv(itype,a0(i_lo,2,1),nm,a0(i_lo,je,1),np,3)
!!$            CALL p_sendrecv(itype,a0(i_lo,je-ihalo,1),np,a0(i_lo,j_lo,1),nm,3)
!!$
!!$          ENDIF

          IF ( have_g_js ) THEN
             DO k=1,kk
                CALL p_sendrecv(itype,a0(i_lo,je-ihalo,k),np,a0(i_lo,je,k),np,3)
             ENDDO
          ELSEIF ( have_g_je ) THEN
             DO k=1,kk
                CALL p_sendrecv(itype,a0(i_lo,2,k),nm,a0(i_lo,j_lo,k),nm,3)
             ENDDO
          ELSE
             DO k=1,kk
                CALL p_sendrecv(itype,a0(i_lo,2,k),nm,a0(i_lo,je,k),np,3)
                CALL p_sendrecv(itype,a0(i_lo,je-ihalo,k),np,a0(i_lo,j_lo,k),nm,3)
             ENDDO
          ENDIF
        ENDIF

      ENDIF ! lmpitype

    ELSE ! nprocy = 1

      ! nothing needed

    ENDIF ! nprocy = 1


!    lmpitype=lold

#ifdef _PROFILE
  CALL trace_stop ('bounds_exch_3d ew/ns', 6)
  CALL trace_start ('bounds_exch_3d tp', 7)
#endif

    ! make the northern margin cyclic

    IF( lbounds_exch_tp .AND. have_g_js) THEN

      ! define send buffer, only use the ihalo+1 needed lines

      SELECT CASE (ttt)

      CASE('p','p-','u','u+')     !! lines 3 to ihalo+3

        DO ii = 1,ihalo+1

          ylo=(ii-1)*i_size*kk+1
          yup=ii*i_size*kk

          xys1(ylo:yup)      = RESHAPE(a0(i_lo:i_up,ii+2,1:kk),(/i_size*kk/))    ! prepare upper send buffer

        ENDDO

      CASE('v','v+','vf','vf+','s','s-')     !! 2 to ihalo+2

        DO ii = 1,ihalo+1

          ylo=(ii-1)*i_size*kk+1
          yup=ii*i_size*kk

          xys1(ylo:yup)      = RESHAPE(a0(i_lo:i_up,ii+1,1:kk), (/i_size*kk/))   ! prepare upper send buffer

        ENDDO


      CASE('uu','vv')

        ! nothing to be done

      CASE default !! unsupported value for ttt

        CALL STOP_ALL(' bounds_exch : argument unsupported ')

      END SELECT



! ---------------------------------------------------------------------------------

      ne=(nprocx-1)-p_pe                                      ! find the neighbour cpu

      yup=(ihalo+1)*i_size*kk

      IF ( p_pe /= ne ) THEN

        IF ( lnonblock ) THEN

          CALL p_isend(xys1(1:yup),ne,6)                ! exchange with neighbour cpu
          CALL p_irecv(xyr1(1:yup),ne,6)                ! exchange with neighbour cpu

          CALL p_wait

        ELSE

          CALL p_sendrecv(xys1(1:yup),ne,xyr1(1:yup),ne,6)  ! exchange with neighbour cpu

        ENDIF

      ELSE                                                 ! dont do send/recive on the same cpu

        xyr1(1:yup)=xys1(1:yup)

      ENDIF


      DO ii=1,ihalo+1

        ylo=(ii-1)*i_size*kk+1
        yup=ii*i_size*kk

        buffer(i_lo:i_up,ii,1:kk) = RESHAPE(xyr1(ylo:yup),(/i_size,kk/)) ! update upper halo

      ENDDO

! ---------------------------------------------------------------------------------

      SELECT CASE (ttt)

      CASE('p','p-')                                                     !! p point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='p-')

        DO ii=1,ihalo+1                                                  ! all needed halos
          DO is=1,i_size                                                 ! size from i_lo:i_up

            il=i_lo-1+is                                                 ! left index
            ir=i_up+1-is                                                 ! right index

            a0(il,3-ii,1:kk) =  multip * buffer(ir,ii,1:kk)              ! line (-)3 -> 2 ; (-)4 -> 1

          END DO
        END DO

      CASE('v','v+','vf','vf+')                                          !! v point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='v'.OR.ttt=='vf')

        DO ii=2,ihalo+1                                                  ! all needed halos
          DO is=1,i_size                                                 ! size from i_lo:i_up

            il=i_lo-1+is                                                 ! left index
            ir=i_up+1-is                                                 ! right index

            a0(il,3-ii,1:kk) = multip * buffer(ir,ii,1:kk)               ! line -3 -> 1

          END DO
        END DO


        IF (ttt=='vf'.OR.ttt=='vf+') THEN                                !! v point with/without sign change

          ih=0
          IF ( p_pe <= (nprocx-1)/2 .OR. p_pe == nprocxy ) ih=i_size
          IF ( p_pe == ne ) ih =i_size/2

          DO is=1,ih                                                     ! size from i_lo:i_up

            il=i_lo-1+is                                                 ! left index
            ir=i_up+1-is                                                 ! right index

            a0(il,2,1:kk) =  multip * buffer(ir,1,1:kk)                  ! line (-)2 -> 2

          END DO

        ENDIF

      CASE('u','u+')                                                     !! u point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='u')

        DO ii=1,ihalo+1                                                  ! all needed halos
          DO is=1,i_size-1                                               ! size from i_lo:i_up-1

            il=i_lo-1+is                                                 ! left index
            ir=i_up-is                                                   ! right index

            a0(il,3-ii,1:kk) =  multip * buffer(ir,ii,1:kk)              ! line (-)3 -> 2 ; (-)4 -> 1

          END DO
        END DO

      CASE('s','s-','sf','sf-')                                          !! psi point with/without sign change

        multip = MERGE(-1.0_wp, 1.0_wp, ttt=='s-'.OR.ttt=='sf-')

        DO ii=2,ihalo+1                                                  ! all needed halos
          DO is=1,i_size-1                                               ! size from i_lo:i_up-1

            il=i_lo-1+is                                                 ! left index
            ir=i_up-is                                                   ! right index

            a0(il,3-ii,1:kk) = multip * buffer(ir,ii,1:kk)               ! line (-)3 -> 1

          END DO
        END DO

        IF (ttt=='sf'.OR.ttt=='sf-') THEN                                !! v point with/without sign change

          ih=0
          IF ( p_pe <= (nprocx-1)/2 .OR. p_pe == nprocxy ) ih=i_size
          IF ( p_pe == ne ) ih =i_size/2

          DO is=1,i_size-1                                               ! size from i_lo:i_up-1

            il=i_lo-1+is                                                 ! left index
            ir=i_up-is                                                   ! right index

            a0(il,2,1:kk) = multip * buffer(ir,1,1:kk)                   ! line (-)2 -> 2

          END DO

        ENDIF



      END SELECT

    ENDIF     ! end of the northern margin treatment


    te = p_time()
    t3d = t3d + te-ts
    n3d = n3d + 1

#ifdef bounds_exch_check

    aa1(:,:) = a0(:,:)
    nn=0
    DO jj=1,je
       DO ii=1,ie
          IF (aa1(ii,jj,1:kk).NE.aa0(ii,jj,1:kk)) THEN
            nn=nn+1
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(nn)
    IF ( p_pe .EQ. p_io ) THEN
       IF ( nn .EQ. 0 )   WRITE(0,*) 'attn: pe ',p_pe,' of ',nprocxy,' needless bounds_exch! ',text
    ENDIF
    DEALLOCATE(aa0,aa1)

#endif


    IF(p_nprocs > nprocxy) THEN
      ! Test mode
      IF(PRESENT(text)) THEN
        CALL para_check_3d(a0,text)
      ELSE
        CALL para_check_3d(a0,'bounds_exch_3d')
      ENDIF
    ENDIF

#ifdef _PROFILE
  CALL trace_stop ('bounds_exch_3d tp', 7)
  CALL trace_stop ('bounds_exch_3d total', 5)
#endif

  END SUBROUTINE bounds_exch_halo_3d


  SUBROUTINE print_stats

    REAL(dp) :: t2, t3
    INTEGER :: n

    IF(p_pe==p_io) THEN
      WRITE(nerr,*) 'Number of 2D boundary exchanges: ',n2d
      WRITE(nerr,*) 'Number of 3D boundary exchanges: ',n3d
    ENDIF

    DO n=0,nprocxy-1
      t2 = t2d
      t3 = t3d
      CALL p_bcast(t2,n)
      CALL p_bcast(t3,n)
      IF(p_pe==p_io) THEN
        WRITE(nerr,'(a,i4,a,2f10.3)') 'PE: ',n, &
             ' Times for 2D/3D boundary exchanges: ',t2,t3
      ENDIF
    ENDDO

  END SUBROUTINE print_stats



      END MODULE MO_BOUNDSEXCH
