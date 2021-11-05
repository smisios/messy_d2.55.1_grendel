MODULE mo_parallel

  USE mo_kind,   ONLY: dp, wp
  USE mo_param1, ONLY: ke, ie, je, ie_g, je_g, icycli
#ifdef MESSY
  USE messy_main_constants_mem, ONLY: nerr
  USE mo_mpi, ONLY: p_pe, p_io, p_nprocs, p_all_comm, p_max, p_min, p_sum, p_bcast, &
                    p_abort, p_gather, p_send, p_recv, p_isend, p_wait, p_barrier, p_sendrecv!, nerr
#else
  USE mo_wncdf,  ONLY: write_ncdf_single_array
  USE mo_mpi
#endif

  IMPLICIT NONE

  !> max possible number of MPI tasks, can be increased if needed
  INTEGER, PARAMETER :: maxproc = 15359
  !> every task replicates hdom_overlap_[xy] rows and columns of each
  !> of its neighbors, resulting in a total overlap of domains of
  !> dom_overlap_[xy]
  INTEGER, PARAMETER :: hdom_overlap_x = 1, dom_overlap_x = 2*hdom_overlap_x
  INTEGER, PARAMETER :: hdom_overlap_y = 1, dom_overlap_y = 2*hdom_overlap_y

  ! Number of subdivisions in x/y direction

  INTEGER :: nprocx, nprocy
  INTEGER :: nprocxy

  ! Our own offset

  INTEGER :: p_ioff, p_joff

  ! Flag if we have the boundaries

  LOGICAL :: have_g_is, have_g_ie, have_g_js, have_g_je

  ! Start of the inner domains in x/y direction
  ! p_lim_x(0) = 2
  ! p_lim_x(i) = start of inner domain i
  ! p_lim_x(nprocx) = ie_g

  INTEGER, PRIVATE :: p_lim_x(0:maxproc), p_lim_y(0:maxproc)

  ! For every processor: number in x/y direction (0-based)

  INTEGER :: p_num_x(0:maxproc), p_num_y(0:maxproc)

  ! Global offsets and sizes for each processor (both for outer domains)

  INTEGER, PRIVATE :: p_ioff_g(0:maxproc), p_joff_g(0:maxproc)
  INTEGER, PRIVATE :: p_size_x(0:maxproc), p_size_y(0:maxproc)


#ifdef NOLAND
  INTEGER, ALLOCATABLE :: iproc(:,:)
#endif

  INTERFACE para_check
    MODULE PROCEDURE para_check_2d
    MODULE PROCEDURE para_check_3d
  END INTERFACE

  INTERFACE scatter
    MODULE PROCEDURE scatter_2d
    MODULE PROCEDURE scatter_3d
!    MODULE PROCEDURE scatter_isend_2d
!    MODULE PROCEDURE scatter_isend_3d
  END INTERFACE

  INTERFACE gather
!    MODULE PROCEDURE gather_2d
    MODULE PROCEDURE gather_mpi_2d
    MODULE PROCEDURE gather_i2d
    MODULE PROCEDURE gather_3d
!    MODULE PROCEDURE gather_ptr_2d
!    MODULE PROCEDURE gather_ptr_3d
  END INTERFACE


  INTERFACE global_sum
    MODULE PROCEDURE global_sum_i
    MODULE PROCEDURE global_sum_r
    MODULE PROCEDURE global_sum_1d
    MODULE PROCEDURE global_sum_2d
    MODULE PROCEDURE global_sum_3d
  END INTERFACE

CONTAINS

  !-----------------------------------------------------------------------

  SUBROUTINE p_deco
    USE MO_COMMO1, only : lbounds_exch_tp
    !   domain decomposition

    INTEGER :: i, nx, ny
    INTEGER :: ii

    IF(p_nprocs > 1) THEN

#ifdef DEBUG
      WRITE(nerr,*) 'Process ',p_pe,' of ',p_nprocs,' is alive'
#endif

      ! set some variables

      IF(p_pe==0) THEN

        ! nprocx and nprocy must be set by the calling process

        IF(nprocx==0 .OR. nprocy==0) THEN
          WRITE(nerr,*) 'ERROR: nprocx or nprocy not set'
          CALL p_abort
        ENDIF

!        IF ( lbounds_exch_tp ) THEN
!           IF (MOD(nprocx,2) /= 0 .AND. nprocx  /= 1)  THEN
!              WRITE(nerr,*) 'ERROR: for the tripolar setup  nprocx should be even (or one) '
!              CALL p_abort
!           ENDIF
!        ENDIF

        nprocxy = nprocx*nprocy

        IF(nprocxy /= p_nprocs .AND. nprocxy /= p_nprocs-1) THEN
          WRITE(nerr,*)'Number of processors = ',p_nprocs
          WRITE(nerr,*)'nprocx = ',nprocx,' nprocy = ',nprocy
          WRITE(nerr,*)'Number of processors doesnt fit!'
          CALL p_abort
        ENDIF

        IF(((ie_g-dom_overlap_x)/nprocx)<3 .OR. ((je_g-dom_overlap_y)/nprocy)<3) THEN
          WRITE(nerr,*)'Decomposed domain gets too small'
          WRITE(nerr,*)'We need at least 3 rows in every direction'
          CALL p_abort
        ENDIF

      ENDIF

      ! broadcast nprocx and nprocy

      CALL p_bcast(nprocx,0)
      CALL p_bcast(nprocy,0)
      nprocxy = nprocx*nprocy

      ! Decomposition - compute domain limits

      if ( lbounds_exch_tp ) then
         DO i=0,(nprocx/2)
            p_lim_x(i) = dom_overlap_x + i*(ie_g-dom_overlap_x)/nprocx
            ii=nprocx-i
            p_lim_x(ii) = (ie_g - (p_lim_x(i)-dom_overlap_x))
         ENDDO
      else
         DO i=0,nprocx
            p_lim_x(i) = dom_overlap_x + i*(ie_g-dom_overlap_x)/nprocx
         ENDDO
      endif

      DO i=0,nprocy
        p_lim_y(i) = dom_overlap_y + i*(je_g-dom_overlap_y)/nprocy
      ENDDO

      ! Set number of processors in x and y direction

      DO i=0,nprocx-1
        p_num_x(i:nprocxy-1:nprocx) = i
      ENDDO
      DO i=0,nprocy-1
        p_num_y(i*nprocx:(i+1)*nprocx-1) = i
      ENDDO

      ! Offsets and sizes

      DO i=0,nprocxy-1
        nx = p_num_x(i)
        ny = p_num_y(i)

        p_ioff_g(i) = p_lim_x(nx) - dom_overlap_x
        p_joff_g(i) = p_lim_y(ny) - dom_overlap_y
        p_size_x(i) = p_lim_x(nx+1) - p_lim_x(nx) + dom_overlap_x
        p_size_y(i) = p_lim_y(ny+1) - p_lim_y(ny) + dom_overlap_y
      ENDDO

      ! Get our own values

      IF(p_pe<nprocxy) THEN
        ie = p_size_x(p_pe)
        je = p_size_y(p_pe)
        p_ioff = p_ioff_g(p_pe)
        p_joff = p_joff_g(p_pe)
        have_g_is = p_num_x(p_pe) == 0
        have_g_ie = p_num_x(p_pe) == nprocx-1
        have_g_js = p_num_y(p_pe) == 0
        have_g_je = p_num_y(p_pe) == nprocy-1
      ELSE
        ie = ie_g
        je = je_g
        p_ioff = 0
        p_joff = 0
        have_g_is = .TRUE.
        have_g_ie = .TRUE.
        have_g_js = .TRUE.
        have_g_je = .TRUE.
      ENDIF

#ifdef DEBUG
      WRITE(nerr,'(a,i2,a,2i4,a,2i4)')'Proc ',p_pe,' offset: ',p_ioff,p_joff, &
           ' Size: ',ie,je
#endif

    ELSE

      nprocx = 1
      nprocy = 1
      nprocxy = 1
      p_lim_x(0) = dom_overlap_x
      p_lim_x(1) = ie_g
      p_lim_y(0) = dom_overlap_y
      p_lim_y(1) = je_g
      p_num_x(0) = 0
      p_num_y(0) = 0
      p_ioff_g(0) = 0
      p_joff_g(0) = 0
      p_size_x(0) = ie_g
      p_size_y(0) = je_g
      ie = ie_g
      je = je_g
      p_ioff = 0
      p_joff = 0
      have_g_is = .TRUE.
      have_g_ie = .TRUE.
      have_g_js = .TRUE.
      have_g_je = .TRUE.

#ifdef DEBUG
      WRITE(nerr,*) 'Running on a single processor'
#endif

    ENDIF

!if (p_pe.eq.0) then
!     WRITE(0,*) 'Am Ende von deco'
!     DO i=0,nprocxy-1
!       write(0,*)'p_num_x(',i,')=', p_num_x(i)
!       write(0,*)'p_size_x(',i,')=', p_size_x(i)
!       write(0,*)'p_ioff_g(',i,')=', p_ioff_g(i)
!     enddo
!     DO i=0,nprocxy-1
!       write(0,*)'p_lim_x(',i,')=', p_lim_x(i)
!     enddo
!     DO i=0,nprocxy-1
!       write(0,*)'p_num_y(',i,')=', p_num_y(i)
!       write(0,*)'p_size_y(',i,')=', p_size_y(i)
!       write(0,*)'p_joff_g(',i,')=', p_joff_g(i)
!     enddo
!     DO i=0,nprocxy-1
!       write(0,*)'p_lim_y(',i,')=', p_lim_y(i)
!     enddo
!ndif

  END SUBROUTINE p_deco

  !-----------------------------------------------------------------------

  SUBROUTINE gather_arr(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

    REAL(wp),    INTENT(IN)  :: arrl(:,:)
    REAL(wp),    INTENT(OUT) :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

!    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
!    REAL(wp), ALLOCATABLE :: aux(:,:)

     CALL gather(arrl,arrg,pe)

!    IF(p_pe/=pe) THEN
!      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
!    ELSE
!      DO n=0,nprocxy-1
!        iis = 1
!        iie = p_size_x(n)
!        jjs = 1
!        jje = p_size_y(n)
!        ALLOCATE(aux(iie,jje))
!        IF (n==pe) THEN
!          aux = arrl
!        ELSE
!          CALL p_recv(aux,n,1111)
!        ENDIF
!        ! Copy only outer limits into arrg
!        IF(p_num_x(n) /= 0        ) iis = iis+1
!        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
!        IF(p_num_y(n) /= 0        ) jjs = jjs+1
!        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
!        ioff = p_ioff_g(n)
!        joff = p_joff_g(n)
!        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
!        DEALLOCATE(aux)
!      ENDDO
!    ENDIF

  END SUBROUTINE gather_arr

  SUBROUTINE gather3_arr(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL(wp),    INTENT(IN)  :: arrl(:,:,:)
    REAL(wp),    INTENT(OUT) :: arrg(:,:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL(wp), ALLOCATABLE :: aux(:,:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje,ke))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,1:ke) = aux(iis:iie,jjs:jje,1:ke)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather3_arr

  !> Gathers all local array parts into a global array on PE pe

  SUBROUTINE gather_mpi_2d(arrl, arrg, pe)
#ifdef _PROFILE
    USE mo_profile,      ONLY: trace_start, trace_stop
#endif

    REAL(wp),    INTENT(IN)  :: arrl(:,:)
    REAL(wp)                 :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    REAL(wp), ALLOCATABLE :: aux1(:,:)
    REAL(wp), ALLOCATABLE :: auxg(:,:)!,aux2(:,:)
    INTEGER :: total_cols
    INTEGER :: i,j
    INTEGER :: n, iis, iie, jjs, jje, ioff, joff!,k
    INTEGER :: ie_max, je_max,offset,jjsl,jjel

!
!       write(0,*)'In my_gatherv before p_gatherv nprocxy=',nprocxy

    ie_max=MAXVAL(p_size_x)
    je_max=MAXVAL(p_size_y)

    ALLOCATE(aux1(ie_max,je_max))
    aux1 = 0.0_wp

!        if (p_pe.eq.pe) then
    total_cols = je_max*nprocxy
!         write(0,*)'total_cols =', total_cols
    ALLOCATE(auxg(ie_max,total_cols))
    auxg = 0.0_wp
!        endif

    iie = p_size_x(p_pe)
    jje = p_size_y(p_pe)

    DO i=1,iie
      DO j=1,jje

        aux1(i,j) = arrl(i,j)

      ENDDO
    ENDDO

!        do i=1,iie
!         do j=1,jje
!
!          write(0,*)'p_pe = ',p_pe,'arrl(',i,',',j,')=', arrl(i,j)
!
!         enddo
!        enddo

!       write(0,*)'In my_gather before p_gather'
!! Case nprocxy=1
    IF (nprocxy.GT.1) THEN
      CALL p_gather(aux1,auxg,pe)
    ELSE
      auxg = aux1
    ENDIF

!       write(0,*)'In my_gather after p_gather'

!        if (p_pe.eq.pe) then
!
!          do i=1, ie_max
!           do j=1,total_cols
!
!            write(0,*)'auxg(',i,',',j,')=',auxg(i,j)
!
!           enddo
!          enddo
!
!        endif

    IF(p_pe.EQ.pe) THEN
!       write(0,*)'In my_gather before do loop nprocxy=',nprocxy
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        offset = (n*je_max) + 1
        jjsl = offset
        jjel = offset + p_size_y(n) - 1

!           ALLOCATE(aux2(iie,jje))
!
!          offset = (n*je_max) + 1
!
!         do i = 1, iie
!          k = 1
!
!          do j = offset, offset+jje-1
!
!           aux2(i,k) = auxg(i,j)
!           k=k+1
!
!          enddo
!         enddo
!
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjsl = jjsl+1
        IF(p_num_y(n) /= nprocy-1 ) jjel = jjel-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
!        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux2(iis:iie,jjs:jje)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = auxg(iis:iie,jjsl:jjel)

!        DEALLOCATE(aux2)
      ENDDO
!       CALL write_ncdf_single_array('dump_arrg.nc', arrg, shape(arrg), 'arrg')
    ENDIF


    DEALLOCATE(aux1)
!        if (p_pe.eq.pe) then
    DEALLOCATE(auxg)
!        endif


  END SUBROUTINE gather_mpi_2d



  SUBROUTINE gather_2d(arrl, arrg, pe)

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif

    REAL(wp),    INTENT(IN)  :: arrl(:,:)
    REAL(wp)                 :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL(wp), ALLOCATABLE :: aux(:,:)



    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
        DEALLOCATE(aux)
      ENDDO
    ENDIF


  END SUBROUTINE gather_2d

  !> Gathers all local array parts into a global array on PE pe
  !> @param arrl local array part to copy
  !> @param arrg global array to collect the copies in
  !> @param pe parallel process identifier to gather at (only this one must
  !> pass a valid arrg argument)
  SUBROUTINE gather_i2d(arrl, arrg, pe)
    INTEGER, INTENT(IN)  :: arrl(:,:)
    INTEGER              :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    INTEGER, ALLOCATABLE :: aux(:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather_i2d

  !> Gathers all local array parts into a global array on PE pe
  SUBROUTINE gather_ptr_2d(arrl, arrg, pe)
    REAL(wp),    INTENT(IN)  :: arrl(:,:)
    REAL(wp),   POINTER      :: arrg(:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff
    REAL(wp), ALLOCATABLE :: aux(:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        ALLOCATE(aux(iie,jje))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje) = aux(iis:iie,jjs:jje)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather_ptr_2d

  SUBROUTINE gather_3d(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL(wp),    INTENT(IN)  :: arrl(:,:,:)
    REAL(wp)                 :: arrg(:,:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: kk, n, iis, iie, jjs, jje, ioff, joff
    REAL(wp), ALLOCATABLE :: aux(:,:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        kk = UBOUND(arrl,3)

        ALLOCATE(aux(iie,jje,kk))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,:) = aux(iis:iie,jjs:jje,:)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather_3d

  SUBROUTINE gather_ptr_3d(arrl, arrg, pe)

    ! Gathers all local array parts into a global array on PE pe

    REAL(wp),    INTENT(IN)  :: arrl(:,:,:)
    REAL(wp), POINTER        :: arrg(:,:,:)
    INTEGER, INTENT(IN)  :: pe

    INTEGER :: kk, n, iis, iie, jjs, jje, ioff, joff
    REAL(wp), ALLOCATABLE :: aux(:,:,:)

    IF(p_pe/=pe) THEN
      IF(p_pe<nprocxy) CALL p_send(arrl,pe,1111)
    ELSE
      DO n=0,nprocxy-1
        iis = 1
        iie = p_size_x(n)
        jjs = 1
        jje = p_size_y(n)
        kk = UBOUND(arrl,3)

        ALLOCATE(aux(iie,jje,kk))
        IF (n==pe) THEN
          aux = arrl
        ELSE
          CALL p_recv(aux,n,1111)
        ENDIF
        ! Copy only outer limits into arrg
        IF(p_num_x(n) /= 0        ) iis = iis+1
        IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
        IF(p_num_y(n) /= 0        ) jjs = jjs+1
        IF(p_num_y(n) /= nprocy-1 ) jje = jje-1
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,:) = aux(iis:iie,jjs:jje,:)
        DEALLOCATE(aux)
      ENDDO
    ENDIF

  END SUBROUTINE gather_ptr_3d

  !-----------------------------------------------------------------------

!!$  SUBROUTINE scatter_arr(arrg, arrl, pe)
!!$
!!$    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
!!$    ! Since this routine is not used in perfomance critical parts
!!$    ! we use the simple broadcast version
!!$
!!$    REAL(wp),    INTENT(INOUT) :: arrg(:,:)
!!$    REAL(wp),    INTENT(OUT)   :: arrl(:,:)
!!$    INTEGER, INTENT(IN)    :: pe
!!$
!!$    CALL p_bcast(arrg,pe)
!!$
!!$    arrl(:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
!!$
!!$  END SUBROUTINE scatter_arr

!!$  SUBROUTINE scatter3_arr(arrg, arrl, pe)
!!$
!!$    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
!!$    ! Since this routine is not used in perfomance critical parts
!!$    ! we use the simple broadcast version
!!$
!!$    REAL(wp),    INTENT(INOUT) :: arrg(:,:,:)
!!$    REAL(wp),    INTENT(OUT)   :: arrl(:,:,:)
!!$    INTEGER, INTENT(IN)    :: pe
!!$
!!$    CALL p_bcast(arrg,pe)
!!$
!!$    arrl(:,:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je,:)
!!$
!!$  END SUBROUTINE scatter3_arr


!!$  SUBROUTINE scatter_old(arrg, arrl, pe)
!!$
!!$    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs
!!$    ! Since this routine is not used in perfomance critical parts
!!$    ! we use the simple broadcast version
!!$
!!$    REAL(wp),    INTENT(INOUT) :: arrg(:,:)
!!$    REAL(wp),    INTENT(OUT)   :: arrl(:,:)
!!$    INTEGER, INTENT(IN)    :: pe
!!$
!!$
!!$    CALL p_bcast(arrg,pe)
!!$
!!$    arrl(:,:) = arrg(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
!!$
!!$  END SUBROUTINE scatter_old

#ifndef __SX__
  SUBROUTINE scatter_isend_2d(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs

    ! This version uses non blocking communication

    TYPE send_buffer
       REAL(wp), ALLOCATABLE :: aux(:,:)
    END TYPE send_buffer

    TYPE(send_buffer), ALLOCATABLE :: send(:)

    REAL(wp)                   :: arrg(:,:)
    REAL(wp),    INTENT(OUT)   :: arrl(:,:)
    INTEGER, INTENT(IN)    :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff

    IF (p_pe == pe ) THEN
       ALLOCATE(send(0:nprocxy-1))
       DO n = 0,nprocxy-1
          iie = p_size_x(n)
          jje = p_size_y(n)
          ALLOCATE(send(n)%aux(iie,jje))
       ENDDO
    ENDIF

    IF (p_pe == pe ) THEN
       DO n=0,nprocxy-1
          iis = 1
          iie = p_size_x(n)
          jjs = 1
          jje = p_size_y(n)
          ioff = p_ioff_g(n)
          joff = p_joff_g(n)

          send(n)%aux(iis:iie,jjs:jje)=arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje)

          IF (n /= pe) THEN
             CALL p_isend(send(n)%aux,n,1112)
          ELSE
             arrl = send(n)%aux
          ENDIF
      ENDDO

! p_waitall from echam is faster

      CALL p_wait
    ELSE
       CALL p_recv(arrl,pe,1112)
    ENDIF

    IF (p_pe == pe ) THEN
       DO n = 0,nprocxy-1
          DEALLOCATE(send(n)%aux)
       ENDDO
       DEALLOCATE(send)
    ENDIF

   END SUBROUTINE scatter_isend_2d
#endif

   SUBROUTINE scatter_2d(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs

    REAL(wp)                   :: arrg(:,:)
    REAL(wp),    INTENT(OUT)   :: arrl(:,:)
    REAL(wp),    ALLOCATABLE   :: aux(:,:)
    INTEGER, INTENT(IN)    :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff

    IF (p_pe == pe ) THEN
       DO n=0,nprocxy-1
          iis = 1
          iie = p_size_x(n)
          jjs = 1
          jje = p_size_y(n)
          ALLOCATE(aux(iie,jje))

          ioff = p_ioff_g(n)
          joff = p_joff_g(n)

!          IF(p_num_x(n) /= 0        ) iis = iis+1
!          IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
!          IF(p_num_y(n) /= 0        ) jjs = jjs+1
!          IF(p_num_y(n) /= nprocy-1 ) jje = jje-1

          aux(iis:iie,jjs:jje)=arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje)

          IF (n /= pe) THEN
             CALL p_send(aux,n,1112)
          ELSE
             arrl = aux
          ENDIF
          DEALLOCATE (AUX)
      ENDDO
    ELSE
       CALL p_recv(arrl,pe,1112)
    ENDIF

   END SUBROUTINE scatter_2d


  SUBROUTINE scatter_3d(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs

    REAL(wp)                   :: arrg(:,:,:)
    REAL(wp),    INTENT(OUT)   :: arrl(:,:,:)
    REAL(wp),    ALLOCATABLE   :: aux(:,:,:)
    INTEGER, INTENT(IN)    :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff, kk

    IF (p_pe == pe ) THEN
       DO n=0,nprocxy-1
          iis = 1
          iie = p_size_x(n)
          jjs = 1
          jje = p_size_y(n)

          kk = UBOUND(arrl,3)

          ALLOCATE(aux(iie,jje,kk))

          ioff = p_ioff_g(n)
          joff = p_joff_g(n)

!          IF(p_num_x(n) /= 0        ) iis = iis+1
!          IF(p_num_x(n) /= nprocx-1 ) iie = iie-1
!          IF(p_num_y(n) /= 0        ) jjs = jjs+1
!          IF(p_num_y(n) /= nprocy-1 ) jje = jje-1

          aux(iis:iie,jjs:jje,:)=arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,:)

          IF (n /= pe) THEN
             CALL p_send(aux,n,1112)
          ELSE
             arrl = aux
          ENDIF
          DEALLOCATE (AUX)
       ENDDO
    ELSE
       CALL p_recv(arrl,pe,1112)
    ENDIF

   END SUBROUTINE scatter_3d

#ifndef __SX__
  SUBROUTINE scatter_isend_3d(arrg, arrl, pe)

    ! Scatters global array (arrg) on PE pe to the local array (arrl) on all PEs

    ! This version uses non blocking communication

    TYPE send_buffer
       REAL(wp), ALLOCATABLE :: aux(:,:,:)
    END TYPE send_buffer

    TYPE(send_buffer), ALLOCATABLE :: send(:)

    REAL(wp)                   :: arrg(:,:,:)
    REAL(wp),    INTENT(OUT)   :: arrl(:,:,:)
    INTEGER, INTENT(IN)    :: pe

    INTEGER :: n, iis, iie, jjs, jje, ioff, joff, kk

    IF (p_pe == pe ) THEN
       ALLOCATE(send(0:nprocxy-1))
       DO n = 0,nprocxy-1
          iie = p_size_x(n)
          jje = p_size_y(n)
          kk = UBOUND(arrl,3)
          ALLOCATE(send(n)%aux(iie,jje,kk))
       ENDDO
    ENDIF

    IF (p_pe == pe ) THEN
       DO n=0,nprocxy-1
          iis = 1
          iie = p_size_x(n)
          jjs = 1
          jje = p_size_y(n)

          ioff = p_ioff_g(n)
          joff = p_joff_g(n)

          send(n)%aux(iis:iie,jjs:jje,:)=arrg(ioff+iis:ioff+iie,joff+jjs:joff+jje,:)

          IF (n /= pe) THEN
             CALL p_isend(send(n)%aux,n,1112)
          ELSE
             arrl = send(n)%aux
          ENDIF
       ENDDO

! p_waitall from echam is faster

      CALL p_wait

    ELSE
       CALL p_recv(arrl,pe,1112)
    ENDIF

    IF (p_pe == pe ) THEN
       DO n = 0,nprocxy-1
          DEALLOCATE(send(n)%aux)
       ENDDO
       DEALLOCATE(send)
    ENDIF

   END SUBROUTINE scatter_isend_3d
#endif

!-----------------------------------------------------------------------


  SUBROUTINE read_slice_old(iunit,arr)

    ! Reads a 2D array and scatters it to all

    INTEGER,INTENT(IN) :: iunit
    REAL(wp), INTENT(OUT) :: arr(:,:)

    REAL(wp) arr_g(ie_g,je_g)

    IF(p_pe==p_io) READ(iunit) arr_g
    CALL p_bcast(arr_g,p_io)

    arr(:,:) = arr_g(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

  END SUBROUTINE read_slice_old

  SUBROUTINE read_slice(iunit,arr_l)

    ! Reads a 2D array and scatters it to all

    INTEGER,INTENT(IN) :: iunit
    REAL(wp), INTENT(OUT) :: arr_l(:,:)

    REAL(wp),ALLOCATABLE :: arr_g(:,:)

    IF(p_pe==p_io) THEN
       ALLOCATE(arr_g(ie_g,je_g))
    ELSE
              ALLOCATE(arr_g(0,0))
    ENDIF

    IF(p_pe==p_io) READ(iunit) arr_g

    CALL scatter(arr_g,arr_l,p_io)

    DEALLOCATE(Arr_g)

  END SUBROUTINE read_slice

  !-----------------------------------------------------------------------

  SUBROUTINE spool_slice(iunit)

    ! Reads first element of a slice only

    INTEGER,INTENT(IN) :: iunit
    REAL(dp)           :: first

    IF(p_pe==p_io) READ(iunit) first


  END SUBROUTINE spool_slice

  !-----------------------------------------------------------------------

  SUBROUTINE write_slice(iunit,arr)

    ! Gathers a 2D array and writes it to iunit

    INTEGER,INTENT(IN) :: iunit
    REAL(wp), INTENT(IN) :: arr(:,:)

    REAL(wp), ALLOCATABLE :: arr_g(:,:)

    IF(p_pe==p_io) then
       ALLOCATE(arr_g(ie_g,je_g))
    else
       ALLOCATE(arr_g(0,0))
    ENDIF

    CALL gather(arr,arr_g,p_io)
    IF(p_pe==p_io) WRITE(iunit) arr_g

    DEALLOCATE(arr_g)

  END SUBROUTINE write_slice

  SUBROUTINE write_slice_sp(iunit,arr)

    ! Gathers a 2D array and writes it to iunit
    USE mo_kind

    INTEGER,INTENT(IN) :: iunit
    REAL(wp), INTENT(IN) :: arr(:,:)

    REAL(wp), ALLOCATABLE :: arr_g(:,:)

    IF(p_pe==p_io) then
       ALLOCATE(arr_g(ie_g,je_g))
    else
       ALLOCATE(arr_g(0,0))
    ENDIF

    CALL gather(arr,arr_g,p_io)
    IF(p_pe==p_io) WRITE(iunit) REAL(arr_g,sp)

    DEALLOCATE(arr_g)

  END SUBROUTINE write_slice_sp

  SUBROUTINE write_extra_sp(date, code, level, iunit,arr)

    ! Gathers a 2D array and writes it to iunit
    USE mo_kind

    INTEGER, INTENT(IN) :: iunit, date, code, level
    REAL(wp), INTENT(IN) :: arr(:,:)

    IF(p_pe==p_io) THEN
       WRITE(iunit) int(date, i4), int(code, i4), int(level, i4), int(ie_g*je_g, i4)
    END IF
    CALL write_slice_sp(iunit,arr)

  END SUBROUTINE write_extra_sp

  SUBROUTINE write_extra(i1,i2,i3,iunit,arr)

    ! Gathers a 2D array and writes it to iunit
    USE mo_kind

    INTEGER,INTENT(IN) :: iunit,i1,i2,i3
    INTEGER(i8)         :: iheader(4)
    REAL(wp), INTENT(IN) :: arr(:,:)

    iheader(1) = INT(i1, i8)
    iheader(2) = INT(i2, i8)
    iheader(3) = INT(i3, i8)
    iheader(4) = INT(ie_g*je_g, i8)

    IF(p_pe==p_io) WRITE(iunit) iheader(1),iheader(2),iheader(3),iheader(4)
    CALL write_slice(iunit,arr)

  END SUBROUTINE write_extra


  !-----------------------------------------------------------------------

  SUBROUTINE global_sum_r(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

    ! Build global sum for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call

    REAL(wp), INTENT(INOUT):: s0
    REAL(wp), INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

    REAL(wp) s(10)
    INTEGER n

    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

    CALL global_sum_1d(s(1:n))

    s0 = s(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = s(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = s(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = s(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = s(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = s(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = s(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = s(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = s(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = s(n); ENDIF

  END SUBROUTINE global_sum_r

  !-----------------------------------------------------------------------

  SUBROUTINE global_sum_i(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

    ! Build global sum for integer scalars
    ! For performance reasons we permit up to 10 arguments in a single call

    INTEGER, INTENT(INOUT):: s0
    INTEGER, INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

    REAL(wp) :: s(10)
    INTEGER :: n

    s(1) = REAL(s0, wp)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = REAL(s1, wp); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = REAL(s2, wp); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = REAL(s3, wp); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = REAL(s4, wp); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = REAL(s5, wp); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = REAL(s6, wp); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = REAL(s7, wp); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = REAL(s8, wp); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = REAL(s9, wp); ENDIF

    CALL global_sum_1d(s(1:n))

    s0 = NINT(s(1))
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = NINT(s(n)); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = NINT(s(n)); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = NINT(s(n)); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = NINT(s(n)); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = NINT(s(n)); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = NINT(s(n)); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = NINT(s(n)); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = NINT(s(n)); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = NINT(s(n)); ENDIF

  END SUBROUTINE global_sum_i

  !-----------------------------------------------------------------------

  SUBROUTINE global_sum_1d(s)

    ! Build global sum for real 1D array

    REAL(wp), INTENT(INOUT):: s(:)

    REAL(wp) :: r(SIZE(s)), q(SIZE(s)), errmax, err
    INTEGER n, i

    ! If we are running in test mode, use result from last PE and
    ! check if the others are approximatly right, else do real summation

    n = SIZE(s)

    IF(p_nprocs>nprocxy) THEN

      ! Test mode
      IF(p_pe>=nprocxy) THEN
        q(:) = s(:) ! Save s
        s(:) = 0._wp! Remove our contribution to global sum
      ENDIF

      ! Get sum on working PEs
      r = p_sum(s)

      ! Check against saved value on test PE
      IF(p_pe==nprocxy) THEN
        errmax = 0._wp
        DO i=1,n
          IF(q(i)/=r(i)) THEN
            err = ABS(q(i)-r(i))/MAX(ABS(q(i)),ABS(r(i)))
            errmax = MAX(errmax,err)
          ENDIF
        ENDDO
        WRITE(nerr,*) 'Global Sum Max Err = ',errmax
        IF (errmax > 1.e-5_wp) CALL p_abort
        s(:) = q(:) ! Restore s
      ENDIF

      ! We use the value from test PE in order to get always
      ! identical results during test
      CALL p_bcast(s,nprocxy)

    ELSE

      r = p_sum(s)
      s(:) = r(:)

    ENDIF

  END SUBROUTINE global_sum_1d

  !-----------------------------------------------------------------------

  SUBROUTINE global_sum_2d(s)

    ! Build global sum for real 2D array

    REAL(wp), INTENT(INOUT):: s(:,:)

    REAL(wp) r(SIZE(s))

    IF(p_nprocs>nprocxy) THEN
      ! Test mode - use global_sum_1d
      r = RESHAPE(s,(/SIZE(s)/))
      CALL global_sum_1d(r)
      s = RESHAPE(r,SHAPE(s))
    ELSE
      r = p_sum(RESHAPE(s,(/SIZE(s)/)))
      s = RESHAPE(r,SHAPE(s))
    ENDIF

  END SUBROUTINE global_sum_2d
  !-----------------------------------------------------------------------

  SUBROUTINE global_sum_3d(s)

    ! Build global sum for real 2D array

    REAL(wp), INTENT(INOUT):: s(:,:,:)

    REAL(wp) r(SIZE(s))

    IF(p_nprocs>nprocxy) THEN
      ! Test mode - use global_sum_1d
      r = RESHAPE(s,(/SIZE(s)/))
      CALL global_sum_1d(r)
      s = RESHAPE(r,SHAPE(s))
    ELSE
      r = p_sum(RESHAPE(s,(/SIZE(s)/)))
      s = RESHAPE(r,SHAPE(s))
    ENDIF

  END SUBROUTINE global_sum_3d
  !-----------------------------------------------------------------------


  SUBROUTINE global_sum_2d_pio(arr,su)

    USE MO_KIND
    USE MO_PARAM1, only : ie_g,je_g
    USE MO_COMMO1, only : lbounds_exch_tp

    INTEGER             ::  i,j,jb
    REAL(wp), INTENT(IN)    ::  arr(:,:)
    REAL(wp), INTENT(INOUT) ::  su

    REAL(wp), ALLOCATABLE       ::  arr_g(:,:)

    if (p_pe==p_io) then
       ALLOCATE(arr_g(ie_g,je_g))
    else
       ALLOCATE(arr_g(0,0))
    endif

!    call gather_arr(arr,arr_g,p_io)
    call gather(arr,arr_g,p_io)

    jb=2
    if ( lbounds_exch_tp ) jb=3

    if ( p_pe == p_io ) then
       su = 0._wp
       do i=2,ie_g-1
          do j=jb,je_g-1
             su=su+arr_g(i,j)
          enddo
       enddo
    endif
    call p_bcast(su,p_io)

    DEALLOCATE(arr_g)

  END SUBROUTINE global_sum_2d_pio



  SUBROUTINE global_mean_2d(arr,su)

    USE MO_KIND
    USE MO_PARAM1, ONLY : ie,je,ie_g,je_g
    USE MO_COMMO1, ONLY : dlxp, dlyp, weto, lbounds_exch_tp

    ! Build global mean for real 2D array

    INTEGER             ::  i,j,jb
    REAL(wp), INTENT(IN)    ::  arr(:,:)
    REAL(wp), INTENT(INOUT) ::  su

    REAL(wp)                ::  as
    REAL(wp)                ::  arr2(ie,je), sarr(ie,je)

    REAL(wp), ALLOCATABLE       ::  arr_g(:,:)
    REAL(wp), ALLOCATABLE       ::  sarr_g(:,:)

    IF (p_pe==p_io) THEN
       ALLOCATE(arr_g(ie_g,je_g),sarr_g(ie_g,je_g))
    ELSE
       ALLOCATE(arr_g(0,0),sarr_g(0,0))
    ENDIF

    DO i=1,ie
       DO j=1,je
          arr2(i,j)=arr(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
          sarr(i,j)=dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
       ENDDO
    ENDDO

    CALL gather(arr2,arr_g,p_io)
    CALL gather(sarr,sarr_g,p_io)

    IF ( p_pe == p_io ) THEN

       jb=2
       if ( lbounds_exch_tp ) jb=3
       su = 0._wp
       as = 0._wp
       DO i=2,ie_g-1
          DO j=jb,je_g-1
             su=su+arr_g(i,j)
             as=as+sarr_g(i,j)
          ENDDO
       ENDDO
       IF (as == 0._wp) THEN
         su = 0._wp
       ELSE
         su = su / as
       ENDIF
    ENDIF

    CALL p_bcast(su,p_io)

    IF (p_pe==p_io) THEN
       DEALLOCATE(arr_g,sarr_g)
    ENDIF

  END SUBROUTINE global_mean_2d
  !-----------------------------------------------------------------------

  SUBROUTINE global_max(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

    ! Build global max for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call

    REAL(wp), INTENT(INOUT):: s0
    REAL(wp), INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

    REAL(wp) s(10), r(10)
    INTEGER n
    s(:) = 0._wp
    r(:) = 0._wp
    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

    r = p_max(s)

    s0 = r(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF

  END SUBROUTINE global_max

  !-----------------------------------------------------------------------

  SUBROUTINE global_min(s0,s1,s2,s3,s4,s5,s6,s7,s8,s9)

    ! Build global min for real scalars
    ! For performance reasons we permit up to 10 arguments in a single call

    REAL(wp), INTENT(INOUT):: s0
    REAL(wp), INTENT(INOUT), OPTIONAL:: s1,s2,s3,s4,s5,s6,s7,s8,s9

    REAL(wp) s(10), r(10)
    INTEGER n
    s(:) = 0._wp
    r(:) = 0._wp

    s(1) = s0
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s(n) = s1; ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s(n) = s2; ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s(n) = s3; ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s(n) = s4; ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s(n) = s5; ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s(n) = s6; ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s(n) = s7; ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s(n) = s8; ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s(n) = s9; ENDIF

    r = p_min(s)

    s0 = r(1)
    n = 1
    IF (PRESENT(s1)) THEN; n = n+1; s1 = r(n); ENDIF
    IF (PRESENT(s2)) THEN; n = n+1; s2 = r(n); ENDIF
    IF (PRESENT(s3)) THEN; n = n+1; s3 = r(n); ENDIF
    IF (PRESENT(s4)) THEN; n = n+1; s4 = r(n); ENDIF
    IF (PRESENT(s5)) THEN; n = n+1; s5 = r(n); ENDIF
    IF (PRESENT(s6)) THEN; n = n+1; s6 = r(n); ENDIF
    IF (PRESENT(s7)) THEN; n = n+1; s7 = r(n); ENDIF
    IF (PRESENT(s8)) THEN; n = n+1; s8 = r(n); ENDIF
    IF (PRESENT(s9)) THEN; n = n+1; s9 = r(n); ENDIF

  END SUBROUTINE global_min

  !-----------------------------------------------------------------------

  REAL(wp) FUNCTION global_array_sum(arr)

    ! Builds the global sum of a 2D array by gathering it on one PE,
    ! calculating th sum on this PE and broadcasting the result
    ! As opposed to calculating the local sum and the using a global_sum call,
    ! this method should give identical results indepentently of the number
    ! of processors, at least when optimization is switched off.

    REAL(wp), INTENT(IN) :: arr(:,:)

    INTEGER i,j
    REAL(wp) arr_g(ie_g,je_g), s, sum

!    CALL gather_arr(arr,arr_g,0)
    CALL gather(arr,arr_g,0)

    IF(p_pe==0) THEN
      sum = 0._wp
      DO j=2,je_g-1
        DO i=2,ie_g-1
          sum = sum + arr_g(i,j)
        ENDDO
      ENDDO
    ENDIF

    CALL p_bcast(sum,0)

    global_array_sum = sum

    IF(p_pe==nprocxy) THEN
      s = 0._wp
      DO j=2,je-1
        DO i=2,ie-1
          s = s + arr(i,j)
        ENDDO
      ENDDO

      IF(s/=sum) THEN
        WRITE(nerr,*) 'global_array_sum: ',s,sum,ABS(s-sum)
        CALL p_abort
      ENDIF
    ENDIF

  END FUNCTION global_array_sum

  !-----------------------------------------------------------------------

  SUBROUTINE para_check_2d(arr,text,lev)

    ! If running in test mode, checks if the results on the PEs running
    ! in parallel are identical with the results on the test PE
    ! This works only if optimization is switched off

    REAL(wp), INTENT(IN) :: arr(:,:)
    CHARACTER (LEN=*), INTENT(IN) :: text
    INTEGER, INTENT(IN), OPTIONAL :: lev

    REAL(wp), ALLOCATABLE :: aux(:,:)
    INTEGER iie, jje, ioff, joff, i, j, n, num

    IF(p_nprocs<=nprocxy) RETURN


    IF(p_pe<nprocxy) THEN
      CALL p_send(arr,nprocxy,1111)
    ELSE
      num = 0
      DO n=0,nprocxy-1
        iie = p_size_x(n)
        jje = p_size_y(n)
        ioff = p_ioff_g(n)
        joff = p_joff_g(n)
        ALLOCATE(aux(iie,jje))
        CALL p_recv(aux,n,1111)
        DO j=1,jje
          DO i=1,iie
            IF(aux(i,j) /= arr(i+ioff,j+joff)) THEN
              IF(num==0) THEN
                IF(PRESENT(lev)) THEN
                  WRITE(nerr,*) 'Consistency Check Error! K=',lev,text
                ELSE
                  WRITE(nerr,*) 'Consistency Check Error! ',text
                ENDIF
              ENDIF
              WRITE(nerr,'(3i4,2e25.16)') n,i,j,aux(i,j),arr(i+ioff,j+joff)
              num = num+1
            ENDIF
          ENDDO
        ENDDO
        DEALLOCATE(aux)
      ENDDO
      IF(num>0) THEN
        WRITE(nerr,*) num,' errors'
        CALL p_abort
      ENDIF
    ENDIF

    CALL p_barrier(p_all_comm)

  END SUBROUTINE para_check_2d

  !-----------------------------------------------------------------------

  SUBROUTINE para_check_3d(arr,text)

    ! interface for 3D arrays

    REAL(dp), INTENT(IN) :: arr(:,:,:)
    CHARACTER (LEN=*), INTENT(IN) :: text

    INTEGER :: k

    DO k=1,UBOUND(arr,3)
      CALL para_check_2d(arr(:,:,k),text,k)
    ENDDO

  END SUBROUTINE para_check_3d

  !-----------------------------------------------------------------------

  SUBROUTINE stop_all(text)

    ! Does an emergency stop by calling p_abort

    CHARACTER (LEN=*), INTENT(IN) :: text

    WRITE(nerr, '(A)') text

    CALL p_abort

  END SUBROUTINE stop_all

  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------

  SUBROUTINE sethalo2(a0,a1,ihalo)

! ao(ie,je), a1(ie+2,je+2)


    REAL(wp), INTENT(INOUT) :: a0(:,:),a1(:,:)
    INTEGER nm, np, ihalo,i,j
    REAL(wp) xr1(je),xr2(je),yr1(ie),yr2(ie)
    REAL(wp) xs1(je),xs2(je),ys1(ie),ys2(ie)

    a1(:,:) = 0._wp

    DO i=1,ie
       DO j=1,je
          a1(i+1,j+1)=a0(i,j)
       ENDDO
    ENDDO

    ! x-direction halo=2

!    write(0,*)'xdir'

    xs1(:) = a0(3,:)
    xs2(:) = a0(ie-2,:)

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
      ! Get processor numbers of neighbors
      nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
      np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

      CALL p_sendrecv(xs1,nm,xr1,np,201)
      CALL p_sendrecv(xs2,np,xr2,nm,202)

    ELSE
      xr1(:)=xs1(:)
      xr2(:)=xs2(:)
    ENDIF

    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
      DO j=1,je
         a1(1,j+1) = xr2(j)
      ENDDO
    ENDIF

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
      DO j=1,je
         a1(ie+2,j+1) = xr1(j)
      ENDDO
    ENDIF

    ! y-direction

!    write(0,*)'ydir'

      ys1(:) = a0(:,3)
      ys2(:) = a0(:,je-2)

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN

      ! Get processor numbers of neighbors
      nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
      np = MOD(p_pe+nprocxy+nprocx,nprocxy)

      CALL p_sendrecv(ys1,nm,yr1,np,301)
      CALL p_sendrecv(ys2,np,yr2,nm,302)

   ELSE
      yr1(:)=ys1(:)
      yr2(:)=ys2(:)
   ENDIF

   IF(icycli/=0 .OR. .NOT. have_g_js) THEN
      DO i=1,ie
         a1(i+1, 1) = yr2(i)
      ENDDO
   ENDIF
   IF(icycli/=0 .OR. .NOT. have_g_je) THEN
      DO i=1,ie
         a1(i+1,je+2) = yr1(i)
      ENDDO
   ENDIF


  END SUBROUTINE sethalo2

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE sethalon(ttt,a0,a1,ihalo,cart_comm)
     USE MO_COMMO1, only : lbounds_exch_tp
#ifndef NOMPI
    INCLUDE 'mpif.h'
#endif
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    INTEGER cart_comm

    REAL(wp), INTENT(inout) :: a0(:,:),a1(:,:)
    INTEGER nm, np, ihalo,i,j,ii,ihm,ne,il,ir
    REAL(wp) xr1(je,ihalo-1),xr2(je,ihalo-1),yr1(ie+2*(ihalo-1),ihalo-1),yr2(ie+2*(ihalo-1),ihalo-1)
    REAL(wp) xs1(je,ihalo-1),xs2(je,ihalo-1),ys1(ie+2*(ihalo-1),ihalo-1),ys2(ie+2*(ihalo-1),ihalo-1)
    REAL(wp) ystp(ie+2*(ihalo-1),ihalo+1),yrtp(ie+2*(ihalo-1),ihalo+1)
#ifdef cartesian_coord
    INTEGER max_dims
    PARAMETER(max_dims=2)
    INTEGER new_rank,left,right, ifehler

    INTEGER dims(max_dims), coords(max_dims), temp_coords(max_dims)
#endif

    ihm=ihalo-1

    a1(:,:) = 0._wp

    DO i=1,ie
       DO j=1,je
          a1(i+ihm,j+ihm)=a0(i,j)
       ENDDO
    ENDDO

    ! x-direction

    !-------- create the send buffers xs1 (left) and xs2 (right) ------------

    DO ii=1,ihm
       xs1(:,ii) = a0(ii+2,:)
       xs2(:,ii) = a0(ie-(ii+1),:)
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
       ! get processor numbers of neighbors
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords

       CALL mpi_cart_shift(cart_comm, 1, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201,cart_comm)
       CALL p_sendrecv(xs2,np,xr2,nm,202,cart_comm)

#else
       nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
       np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

       !      write(0,*)'in sethalon in x direction in mpi_comm_world p_pe=',p_pe, &
       !                 'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201)
       CALL p_sendrecv(xs2,np,xr2,nm,202)

#endif

    ELSE
       xr1(:,:)=xs1(:,:)
       xr2(:,:)=xs2(:,:)
    ENDIF

    !------ update left halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
       DO ii=1,ihm
          DO j=1,je
             a1(ii,j+ihm) = xr2(j,ihalo-ii)
          ENDDO
       ENDDO
    ENDIF

    !------- update right halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
       DO ii=1,ihm
          DO j=1,je
             a1(ie+(ii+ihm),j+ihm) = xr1(j,ii)
          ENDDO
       ENDDO
    ENDIF

    ! y-direction

    !-------- create the send buffers ys1 (upward) and ys2 (downward) ------------

    DO ii=1,ihm
       ys1(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
       ys1(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
       ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo

       ys2(ihalo:ie+ihm,ii) = a0(:,je-(ii+1))                           ! center
       ys2(1:ihm,ii)=a1(1:ihm,je-(ii+1)+ihm)                            ! left halo
       ys2(ie+ihalo:ie+2*ihm,ii)=a1(ie+ihalo:ie+2*ihm,je-(ii+1)+ihm)    ! right halo
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords
       CALL mpi_cart_shift(cart_comm, 0, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301,cart_comm)
       CALL p_sendrecv(ys2,np,yr2,nm,302,cart_comm)
#else
       !      write(0,*)'in sethalon in y direction in mpi_comm_world'
       ! get processor numbers of neighbors
       nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
       np = MOD(p_pe+nprocxy+nprocx,nprocxy)

       !      write(0,*)'in sethalon in y direction in mpi_comm_world p_pe=',p_pe, &
       !                       'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301)
       CALL p_sendrecv(ys2,np,yr2,nm,302)
#endif
    ELSE
       yr1(:,:)=ys1(:,:)
       yr2(:,:)=ys2(:,:)
    ENDIF

    !-------- update north halo-----------------

    IF ( .NOT. have_g_js) THEN
       DO ii=1,ihm
          a1(:,ii) = yr2(:,ihalo-ii)
       ENDDO
       ENDIF

    !-------- update south halo-----------------

    IF ( .NOT. have_g_je) THEN
       DO ii=1,ihm
          a1(:,((je+ihm)+ii)) = yr1(:,ii)
       ENDDO
    ENDIF


    IF( have_g_js .and. lbounds_exch_tp ) THEN

       !-------- create the send buffers ys1 (upward) ------------

!!$       if ( ttt == 'p') then
!!$          DO ii=1,ihm
!!$             ys1(ihalo:ie+ihm,ii) = a0(:,ii+4)                                ! center
!!$             ys1(1:ihm,ii)=a1(1:ihm,ii+4+ihm)                                 ! left halo
!!$             ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+4+ihm) ! right halo
!!$          enddo
!!$       endif
!!$
!!$       if ( ttt == 'u+') then
!!$          DO ii=1,ihm
!!$             ys1(ihalo:ie+ihm,ii) = a0(:,ii+4)                                ! center
!!$             ys1(1:ihm,ii)=a1(1:ihm,ii+4+ihm)                                 ! left halo
!!$             ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+4+ihm) ! right halo
!!$          enddo
!!$       endif
!!$
!!$       if ( ttt == 'v+') then
!!$          DO ii=1,ihm
!!$             ys1(ihalo:ie+ihm,ii) = a0(:,ii+3)                                ! center
!!$             ys1(1:ihm,ii)=a1(1:ihm,ii+3+ihm)                                 ! left halo
!!$             ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+3+ihm) ! right halo
!!$          enddo
!!$       endif

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+1)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+1+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+1+ihm) ! right halo
          ENDDO
       ENDIF


       !-------- find the communication partner ne -----------------------

       ne=(nprocx-1)-p_pe

       IF ( p_pe /= ne ) THEN
!          CALL p_sendrecv(ys1,ne,yr2,ne,303)
          CALL p_sendrecv(ystp,ne,yrtp,ne,303)
       ELSE
!          yr2(:,:)=ys1(:,:)
          yrtp(:,:)=ystp(:,:)
       ENDIF

!       do ii=1,ie+2*ihm
!       print*,'yr2',ii,yr2(ii,1),yr2(ii,2),yr2(ii,3)
!       enddo

       !-------- update north halo-----------------

!!$       if ( ttt == 'p') then
!!$          DO ii=1,ihm
!!$             DO i=1,ie+(2*ihm)
!!$                il=i
!!$                ir=(ie+2*ihm)+1-i
!!$                a1(il,ii) = yr2(ir,ihalo-ii)
!!$             ENDDO
!!$          ENDDO
!!$       endif
!!$
!!$       if ( ttt == 'u+') then
!!$          DO ii=1,ihm
!!$             DO i=1,ie+(2*ihm)-1
!!$                il=i
!!$                ir=(ie+2*ihm)-i
!!$                a1(il,ii) = yr2(ir,ihalo-ii)
!!$             ENDDO
!!$          ENDDO
!!$       endif
!!$
!!$       if ( ttt == 'v+') then
!!$          DO ii=1,ihm
!!$             DO i=1,ie+(2*ihm)
!!$                il=i
!!$                ir=(ie+2*ihm)+1-i
!!$                a1(il,ii) = yr2(ir,ihalo-ii)
!!$             ENDDO
!!$          ENDDO
!!$       endif

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)-1
                il=i
                ir=(ie+2*ihm)-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
          ii=ihalo+1
          DO i=1,ie+(2*ihm)/2
             il=i
             ir=(ie+2*ihm)+1-i
             a1(il,ii) = yrtp(ir,ihalo+2-ii)
          ENDDO
       ENDIF

    ENDIF


  END SUBROUTINE sethalon

  !-----------------------------------------------------------------------

  SUBROUTINE sethalon_new(ttt,a0,a1,ihalo,cart_comm)
     USE MO_COMMO1, only : lbounds_exch_tp
#ifndef NOMPI
    INCLUDE 'mpif.h'
#endif
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    INTEGER, INTENT(IN), OPTIONAL :: cart_comm
    REAL(wp), INTENT(inout) :: a0(:,:),a1(:,:)
    INTEGER nm, np, ihalo,i,j,ii,ihm,ne,il,ir

    REAL(wp) xr1(je,ihalo),xr2(je,ihalo),yr1(ie+2*(ihalo-1),ihalo),yr2(ie+2*(ihalo-1),ihalo)
    REAL(wp) xs1(je,ihalo),xs2(je,ihalo),ys1(ie+2*(ihalo-1),ihalo),ys2(ie+2*(ihalo-1),ihalo)

    REAL(wp) ystp(ie+2*(ihalo-1),ihalo+1),yrtp(ie+2*(ihalo-1),ihalo+1)
#ifdef cartesian_coord
    INTEGER ifehler
    INTEGER max_dims
    PARAMETER(max_dims=2)
    INTEGER new_rank,left,right

    INTEGER dims(max_dims), coords(max_dims), temp_coords(max_dims)
#endif

    ihm=ihalo-1

    a1(:,:) = 0._wp


    !------- update the inner part of the return field (without halos)

    DO i=1,ie
       DO j=1,je
          a1(i+ihm,j+ihm)=a0(i,j)
       ENDDO
    ENDDO

    ! x-direction

    !-------- create the send buffers xs1 (left) and xs2 (right) ------------

    DO ii=1,ihalo
       xs1(:,ii) = a0(1+ii,:)
       xs2(:,ii) = a0(ie-ii,:)
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
       ! get processor numbers of neighbors
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords

       CALL mpi_cart_shift(cart_comm, 1, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201,cart_comm)
       CALL p_sendrecv(xs2,np,xr2,nm,202,cart_comm)

#else
       nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
       np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

       !      write(0,*)'in sethalon in x direction in mpi_comm_world p_pe=',p_pe, &
       !                 'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201)
       CALL p_sendrecv(xs2,np,xr2,nm,202)

#endif

    ELSE
       xr1(:,:)=xs1(:,:)
       xr2(:,:)=xs2(:,:)
    ENDIF


    !------ update left halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
       DO ii=1,ihalo
          DO j=1,je
             a1(ii,j+ihm) = xr2(j,ihalo+1-ii)
          ENDDO
       ENDDO
    ENDIF

    !------- update right halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
       DO ii=1,ihalo
          DO j=1,je
             a1(ie+ihm+ii-1,j+ihm) = xr1(j,ii)
          ENDDO
       ENDDO
    ENDIF

    ! y-direction

    !-------- create the send buffers ys1 (upward) and ys2 (downward) ------------

    DO ii=1,ihalo
       ys1(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
       ys1(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
       ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo

       ys2(ihalo:ie+ihm,ii) = a0(:,je-(ii+1))                           ! center
       ys2(1:ihm,ii)=a1(1:ihm,je-(ii+1)+ihm)                            ! left halo
       ys2(ie+ihalo:ie+2*ihm,ii)=a1(ie+ihalo:ie+2*ihm,je-(ii+1)+ihm)    ! right halo
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords
       CALL mpi_cart_shift(cart_comm, 0, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301,cart_comm)
       CALL p_sendrecv(ys2,np,yr2,nm,302,cart_comm)
#else
       !      write(0,*)'in sethalon in y direction in mpi_comm_world'
       ! get processor numbers of neighbors
       nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
       np = MOD(p_pe+nprocxy+nprocx,nprocxy)

       !      write(0,*)'in sethalon in y direction in mpi_comm_world p_pe=',p_pe, &
       !                       'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301)
       CALL p_sendrecv(ys2,np,yr2,nm,302)
#endif
    ELSE
       yr1(:,:)=ys1(:,:)
       yr2(:,:)=ys2(:,:)
    ENDIF

    !-------- update north halo-----------------

    IF ( .NOT. have_g_js) THEN
       DO ii=1,ihalo
          a1(:,ii) = yr2(:,ihalo+1-ii)
       ENDDO
       ENDIF

    !-------- update south halo-----------------

    IF ( .NOT. have_g_je) THEN
       DO ii=1,ihalo
          a1(:,(je+ihm+ii-1)) = yr1(:,ii)
       ENDDO
    ENDIF


    IF( have_g_js .and. lbounds_exch_tp ) THEN

       !-------- create the send buffers ystp (upward) ------------

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a0(:,ii+1)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+1+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+1+ihm) ! right halo
          ENDDO
       ENDIF


       !-------- find the communication partner ne -----------------------

       ne=(nprocx-1)-p_pe

       IF ( p_pe /= ne ) THEN
          CALL p_sendrecv(ystp,ne,yrtp,ne,303)
       ELSE
          yrtp(:,:)=ystp(:,:)
       ENDIF


       !-------- update north halo-----------------

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)-1
                il=i
                ir=(ie+2*ihm)-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
          ii=ihalo+1
          DO i=1,ie+(2*ihm)/2
             il=i
             ir=(ie+2*ihm)+1-i
             a1(il,ii) = yrtp(ir,ihalo+2-ii)
          ENDDO
       ENDIF

    ENDIF



  END SUBROUTINE sethalon_new



  SUBROUTINE sethalon_new2(ttt,a1,ihalo,cart_comm)
    USE MO_COMMO1, only : lbounds_exch_tp
#ifndef NOMPI
    INCLUDE 'mpif.h'
#endif
    CHARACTER (LEN=*), INTENT(IN) :: ttt
    INTEGER cart_comm
    REAL(wp), INTENT(inout) :: a1(:,:)
    INTEGER nm, np, ihalo,i,j,ii,ihm,ne,il,ir

    REAL(wp) xr1(je,ihalo),xr2(je,ihalo),yr1(ie+2*(ihalo-1),ihalo),yr2(ie+2*(ihalo-1),ihalo)
    REAL(wp) xs1(je,ihalo),xs2(je,ihalo),ys1(ie+2*(ihalo-1),ihalo),ys2(ie+2*(ihalo-1),ihalo)

    REAL(wp) ystp(ie+2*(ihalo-1),ihalo+1),yrtp(ie+2*(ihalo-1),ihalo+1)
#ifdef cartesian_coord
    INTEGER max_dims,ifehler
    PARAMETER(max_dims=2)
    INTEGER new_rank,left,right

    INTEGER dims(max_dims), coords(max_dims), temp_coords(max_dims)
#endif

    ihm=ihalo-1



    !------- update the inner part of the return field (without halos)

!    DO i=1,ie
!       DO j=1,je
!          a1(i+ihm,j+ihm)=a0(i,j)
!       ENDDO
!    ENDDO

    ! x-direction

    !-------- create the send buffers xs1 (left) and xs2 (right) ------------

    DO ii=1,ihalo
       xs1(:,ii) = a1(ihalo+ii,:)
       xs2(:,ii) = a1(ihalo-1+ie-ii,:)
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
       ! get processor numbers of neighbors
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords

       CALL mpi_cart_shift(cart_comm, 1, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in x direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201,cart_comm)
       CALL p_sendrecv(xs2,np,xr2,nm,202,cart_comm)

#else
       nm = MERGE(p_pe-1,p_pe+nprocx-1,p_num_x(p_pe)/=0)
       np = MERGE(p_pe+1,p_pe-nprocx+1,p_num_x(p_pe)/=nprocx-1)

       !      write(0,*)'in sethalon in x direction in mpi_comm_world p_pe=',p_pe, &
       !                 'nm=',nm,'np=',np

       CALL p_sendrecv(xs1,nm,xr1,np,201)
       CALL p_sendrecv(xs2,np,xr2,nm,202)

#endif

    ELSE
       xr1(:,:)=xs1(:,:)
       xr2(:,:)=xs2(:,:)
    ENDIF


    !------ update left halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_is) THEN
       DO ii=1,ihalo
          DO j=1,je
             a1(ii,j+ihm) = xr2(j,ihalo+1-ii)
          ENDDO
       ENDDO
    ENDIF

    !------- update right halo -----------------

    IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
       DO ii=1,ihalo
          DO j=1,je
             a1(ie+ihm+ii-1,j+ihm) = xr1(j,ii)
          ENDDO
       ENDDO
    ENDIF

    ! y-direction

    !-------- create the send buffers ys1 (upward) and ys2 (downward) ------------

    DO ii=1,ihalo
       ys1(ihalo:ie+ihm,ii) = a1(ihalo-1:ie+ihalo-1,ii+2)                                ! center
       ys1(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
       ys1((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo

       ys2(ihalo:ie+ihm,ii) = a1(ihalo-1:ie+ihalo-1,je-(ii+1))                           ! center
       ys2(1:ihm,ii)=a1(1:ihm,je-(ii+1)+ihm)                            ! left halo
       ys2(ie+ihalo:ie+2*ihm,ii)=a1(ie+ihalo:ie+2*ihm,je-(ii+1)+ihm)    ! right halo
    ENDDO

    !-------- find the communication partners -----------------------

    IF (nprocy>1 .AND. p_pe<nprocxy) THEN
#ifdef cartesian_coord
       CALL mpi_comm_rank(cart_comm, new_rank, ifehler)
       CALL mpi_cart_coords(cart_comm,new_rank,max_dims,coords,ifehler)

       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                   'coords=',coords
       CALL mpi_cart_shift(cart_comm, 0, 1, nm, np, ifehler)
       !      write(0,*)' in sethalon in y direction in cartesian for new_rank=',new_rank,&
       !                  'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301,cart_comm)
       CALL p_sendrecv(ys2,np,yr2,nm,302,cart_comm)
#else
       !      write(0,*)'in sethalon in y direction in mpi_comm_world'
       ! get processor numbers of neighbors
       nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
       np = MOD(p_pe+nprocxy+nprocx,nprocxy)

       !      write(0,*)'in sethalon in y direction in mpi_comm_world p_pe=',p_pe, &
       !                       'nm=',nm,'np=',np

       CALL p_sendrecv(ys1,nm,yr1,np,301)
       CALL p_sendrecv(ys2,np,yr2,nm,302)
#endif
    ELSE
       yr1(:,:)=ys1(:,:)
       yr2(:,:)=ys2(:,:)
    ENDIF

    !-------- update north halo-----------------

    IF ( .NOT. have_g_js) THEN
       DO ii=1,ihalo
          a1(:,ii) = yr2(:,ihalo+1-ii)
       ENDDO
       ENDIF

    !-------- update south halo-----------------

    IF ( .NOT. have_g_je) THEN
       DO ii=1,ihalo
          a1(:,(je+ihm+ii-1)) = yr1(:,ii)
       ENDDO
    ENDIF

    IF( have_g_js .and. lbounds_exch_tp ) THEN

       !-------- create the send buffers ystp (upward) ------------

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a1(ihalo-1:ie+ihalo-1,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a1(ihalo-1:ie+ihalo-1,ii+2)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+2+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+2+ihm) ! right halo
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo+1
             ystp(ihalo:ie+ihm,ii) = a1(ihalo-1:ie+ihalo-1,ii+1)                                ! center
             ystp(1:ihm,ii)=a1(1:ihm,ii+1+ihm)                                 ! left halo
             ystp((ie+ihalo):(ie+2*ihm),ii)=a1((ie+ihalo):(ie+2*ihm),ii+1+ihm) ! right halo
          ENDDO
       ENDIF


       !-------- find the communication partner ne -----------------------

       ne=(nprocx-1)-p_pe

       IF ( p_pe /= ne ) THEN
          CALL p_sendrecv(ystp,ne,yrtp,ne,303)
       ELSE
          yrtp(:,:)=ystp(:,:)
       ENDIF


       !-------- update north halo-----------------

       IF ( ttt == 'p') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'u+') THEN
          DO ii=1,ihalo+1
             DO i=1,ie+(2*ihm)-1
                il=i
                ir=(ie+2*ihm)-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
       ENDIF

       IF ( ttt == 'v+') THEN
          DO ii=1,ihalo
             DO i=1,ie+(2*ihm)
                il=i
                ir=(ie+2*ihm)+1-i
                a1(il,ii) = yrtp(ir,ihalo+2-ii)
             ENDDO
          ENDDO
          ii=ihalo+1
          DO i=1,ie+(2*ihm)/2
             il=i
             ir=(ie+2*ihm)+1-i
             a1(il,ii) = yrtp(ir,ihalo+2-ii)
          ENDDO
       ENDIF

    ENDIF



  END SUBROUTINE sethalon_new2

  SUBROUTINE WRITE_MATRIX(A)
    INTEGER :: I,J
    REAL(wp), DIMENSION(:,:) :: A
    WRITE(0,*)
    DO I = LBOUND(A,1), UBOUND(A,1)
      WRITE(0,*) (A(I,J), J = LBOUND(A,2), UBOUND(A,2))
    END DO


  END SUBROUTINE WRITE_MATRIX

  SUBROUTINE WRITE_MATRIX_1d(A)
    INTEGER :: I
    REAL(wp), DIMENSION(:) :: A
    WRITE(0,*)
    DO I = 1, SIZE(A)
      WRITE(0,*) 'A(',I,')=', A(I)
    END DO
  END SUBROUTINE WRITE_MATRIX_1d



END MODULE mo_parallel
