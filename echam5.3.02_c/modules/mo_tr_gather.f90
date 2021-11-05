! enables two-phase gather
#define _USE_MGATHER
! compares results of mgather with general gather
! slow, for debugging only
#undef _CHECK_MGATHER
! use extended set of timers
#undef _USE_GATHER_TIMER
#if defined(NOMPI) || defined(STANDALONE)
#undef _USE_MGATHER
#endif

MODULE mo_tr_gather

  ! Luis Kornblueh, MPIM, March 2010, initial version
  ! Joerg Behrens, DKRZ, August 2010, two phase gather
  ! Luis Kornblueh, MPIM, August 2010, added
  ! Joerg Behrens, DKRZ, August 2010, fixed intent out bug

  USE mo_kind,          ONLY: dp
  USE mo_exception,     ONLY: message, message_text, finish
  USE mo_decomposition, ONLY: ldc => local_decomposition,  &
                              gdc => global_decomposition, &
                              debug_parallel
  USE mo_mpi,           ONLY: p_nprocs, p_pe, p_io,        &
                              p_recv, p_send, p_real_dp ,  &
                              p_all_comm, npes, &
                              p_barrier, p_communicator_b
  USE mo_transpose,     ONLY: reorder !, gen_inv_map ! op_pj_20110205
#ifdef _USE_GATHER_TIMER
  USE mo_control,       ONLY: ltimer
  USE mo_real_timer,    ONLY: new_timer, timer_start, timer_stop
#endif
! op_pj_20140128+
! replaced by INCLUDE 'mpif.h', see below
!!$#ifdef _USE_MGATHER
!!$  ! PRELIMINARY ONLY ---------------------------------------------------------
!!$  USE mpi,              ONLY : MPI_STATUS_SIZE, MPI_SUCCESS
!!$  ! --------------------------------------------------------------------------
!!$#endif
! op_pj_20140128-

#ifdef _USE_MGATHER
#ifndef __PGI
 USE mpi
#endif
#endif

  IMPLICIT NONE

  PRIVATE
  SAVE ! mz_ht_20120123

!!$! op_pj_20140128+
#ifdef __PGI
#ifdef _USE_MGATHER
  INCLUDE 'mpif.h'
#endif
#endif
!!$! op_pj_20140128-

! gather_field routines expanded to work with variable I/O PE
! specifically can be different from p_io and vary between CHANNELs
  INTERFACE gather_field
    MODULE PROCEDURE gather_gp432
    MODULE PROCEDURE gather_gp32
#ifdef _USE_MGATHER
    MODULE PROCEDURE mgather_gp2
   ! mgather_gp3 needs not to be exposed
#else
    MODULE PROCEDURE gather_gp2
    ! gather_gp3 needs not to be exposed
#endif
  END INTERFACE gather_field

  INTERFACE gather_spectral
    MODULE PROCEDURE gather_sp4
  END INTERFACE gather_spectral

  INTERFACE gather_fourier
    MODULE PROCEDURE gather_sa42
  END INTERFACE gather_fourier

  PUBLIC :: gather_field
  PUBLIC :: gather_spectral
  PUBLIC :: gather_fourier
  ! op_pj_20110204+
  PUBLIC :: gen_inv_map ! generate the inverse of the processor mapping
  ! op_pj_20110204-

  INTEGER, PARAMETER :: tag_gather_gp = 301 ! gridpoint
  INTEGER, PARAMETER :: tag_gather_sp = 401 ! spectral
  INTEGER, PARAMETER :: tag_gather_sa = 501 ! symmetric/anti-symmetric Fourier

  INTEGER, PARAMETER :: tag_gp_sub1   = 311 ! channel for first region
  INTEGER, PARAMETER :: tag_gp_sub2   = 312 ! channel for second region
  INTEGER, PARAMETER :: tag_gp_slab   = 320 ! channel for second phase gather

  ! global variables for two phase gather

! global_decomposition(pe2indx(any_pe))%pe == any_pe
  INTEGER, ALLOCATABLE :: pe2indx(:)

  INTEGER              :: slab_comm
  ! array of communicators to store information for each PE, which can potentially be I/O PE
  ! several other arrays of communicator information are set up below (and
  ! filled by define_slab via call in init_tr_gather)
  INTEGER, ALLOCATABLE :: slab_comm_a(:)

  ! two-phase-gather:
  TYPE gpslab_type
     INTEGER :: glats
     INTEGER :: glate
     INTEGER :: nglat
     INTEGER :: size2
     INTEGER :: disp2
     INTEGER :: owner_ip   ! index for gdc (p_all_comm)
     INTEGER :: owner_pe   ! owner-rank from gdc (p_all_comm)
     INTEGER :: owner_b_pe !  '' within slab_comm
     INTEGER :: src_n
     INTEGER, ALLOCATABLE :: src_ip(:) ! gdc process index of msg.-sources
     INTEGER, ALLOCATABLE :: src_ir(:) ! src region index (1 or 2)
  END TYPE gpslab_type

  ! mz_ht_20120123+ GLOBAL SAVE
!!$  TYPE(gpslab_type), ALLOCATABLE, SAVE :: gpslab(:)
!!$  TYPE(gpslab_type), ALLOCATABLE :: gpslab(:)
  TYPE(gpslab_type), ALLOCATABLE :: gpslab_a(:,:) ! for varying I/O PE
  ! mz_ht_20120123-

  INTEGER :: nslabs
  INTEGER :: idest_gpslab(2) ! defines to which slab we send our gp regions
  INTEGER, ALLOCATABLE :: idest_gpslab_a(:,:)  ! for varying I/O PE
  INTEGER :: slab_id         ! gpslab(slab_id)%owner_pe == p_pe, if slab_id > 0
  INTEGER, ALLOCATABLE :: slab_id_a(:)  ! for varying I/O PE
  INTEGER :: max_rmsg_size_2d, max_rmsg_num
  INTEGER, ALLOCATABLE :: max_rmsg_num_a(:)  ! for varying I/O PE

  LOGICAL :: single_stripe_slab ! (calculated) effective value

#ifdef _USE_GATHER_TIMER
  ! timing modes:
  ! effective switch for using timers (calculated):
  LOGICAL :: timings = .FALSE.
  ! extended timings: enables more timers and overrides ltimer
  LOGICAL, PARAMETER :: extended_timings = .TRUE.

  ! measure modes:
  ! in order to measure gather-performance we need to block accumulated
  ! load imbalance with a barrier, disable for faster runs:
  LOGICAL, PARAMETER :: measure_mode = .TRUE.

  ! another barrier before the second phase
  ! disable for faster runs:
  LOGICAL, PARAMETER :: measure_mode2 = .TRUE.

  ! timers:
  INTEGER :: timer_gp2, timer_gp3
  INTEGER :: timer_mgp2, timer_mgp3
  INTEGER :: timer_sub_gp2, timer_slab_gp2
  INTEGER :: timer_sub_gp3, timer_slab_gp3
  INTEGER :: timer_inject2, timer_inject3
  INTEGER :: timer_fill_gp2, timer_mfill_gp2
  INTEGER :: timer_fill_gp3, timer_mfill_gp3
#endif

  ! module state:
  LOGICAL :: require_init = .TRUE.

  ! debug messages yes/no
  LOGICAL, PARAMETER  :: verbose =  .FALSE.

  ! communication modes for second phase:
  INTEGER, PARAMETER :: sel_gatherv = 1
  INTEGER, PARAMETER :: sel_irecv   = 2
  INTEGER, PARAMETER :: sel_recv    = 3

  ! performance switches:
  !
  ! If we don't work with slabdata directly (e.g.. for parallel IO),
  ! then we can delay the insertion of p_io data. If true: do not
  ! send data from p_io to any other task
  LOGICAL, PARAMETER :: retard_p_io_data = .TRUE.
  !
  ! Note: the following is machine/resolution dependend!
  ! In single_stripe_slab mode we have higher parallelism and less sync
  ! in the first phase but more messages in the second phase
  ! blizzard (Power6, p575, IB): for T63L47 and nproca <= 8 single_stripe_slab
  ! is much faster, for nproca > 8, it is slightly slower
  ! Note: for nproca <= s3_limit: use single_stripe_slab
  INTEGER, PARAMETER :: s3_limit = 8
  ! blizzard (Power6, p575, IB): sel_gatherv is a lot faster
  ! than point-to-point messages
  INTEGER, PARAMETER :: sel_slabgather = sel_gatherv

CONTAINS

  SUBROUTINE init_tr_gather()
    ! if one I/O PE (.TRUE.) or if communicator should be set up
    ! for all PE when several I/O PE used (.FALSE.)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::init_tr_gather'
    INTEGER :: na, nb, ip, ia, nslabs

    IF (.NOT. require_init) RETURN
    require_init = .FALSE.

    IF (.NOT. ALLOCATED(pe2indx)) ALLOCATE(pe2indx(0:p_nprocs-1))
    pe2indx = gen_inv_map(gdc)

#ifdef _USE_GATHER_TIMER
    timer_gp2  = new_timer('gather_gp2')
    timer_mgp2 = new_timer('mgather_gp2')
    timer_gp3  = new_timer('gather_gp3')
    timer_mgp3 = new_timer('mgather_gp3')

    IF (extended_timings) THEN

      timings=.TRUE.

      timer_sub_gp2   = new_timer('subgather_gp2')
      timer_slab_gp2  = new_timer('slabgather_gp2')
      timer_inject2   = new_timer('gather_inject2')
      timer_fill_gp2  = new_timer('fill_gp2')
      timer_mfill_gp2 = new_timer('mfill_gp2')

      timer_sub_gp3   = new_timer('subgather_gp3')
      timer_slab_gp3  = new_timer('slabgather_gp3')
      timer_inject3   = new_timer('gather_inject3')
      timer_fill_gp3  = new_timer('fill_gp3')
      timer_mfill_gp3 = new_timer('mfill_gp3')

    ELSE

      timings = ltimer

    ENDIF

    IF (measure_mode .OR. measure_mode2) timings = .TRUE.
#endif

#ifdef _USE_MGATHER

    nb = SIZE(ldc%mapmesh,1) ! nprocb
    na = SIZE(ldc%mapmesh,2) ! nproca

      ! also in define_slab, here for setting up gpslab_a
      IF (na <= s3_limit .AND. nb > 1) THEN
        nslabs             = 2*na
      ELSE
        nslabs             = na
      ENDIF

    ! within mgather we have no intercommunication between debug
    ! and parallel domain yet, so we have to fall back to the
    ! general gather routines in that case
    !
    ! if I/O PE is variable, set up slab_comm as array with information for each
    ! PE during initialisation, then use arrays of all parameters depending on PE
    ! arrays only used outside of define_slab
    IF ( (debug_parallel == -1) .AND. p_nprocs == ldc%d_nprocs ) THEN
       ! slab_comm also used from other points than write_data - try
       ! communicator array as default
       IF (p_nprocs == 1) THEN
          IF (verbose) &
            CALL message(context, 'call define_slab once to fill array slab_comm_a')
          IF (.NOT. ALLOCATED(slab_comm_a)) ALLOCATE(slab_comm_a(0:1))
          IF (.NOT. ALLOCATED(slab_id_a)) ALLOCATE(slab_id_a(0:1))
          IF (.NOT. ALLOCATED(max_rmsg_num_a))  ALLOCATE(max_rmsg_num_a(0:1))
          IF (.NOT. ALLOCATED(idest_gpslab_a))  ALLOCATE(idest_gpslab_a(0:1,2))
          IF (.NOT. ALLOCATED(gpslab_a)) ALLOCATE(gpslab_a(0:1,nslabs))
          ip = 0
          CALL define_slab(ip)
          slab_comm_a(ip) = slab_comm
          slab_id_a(ip) = slab_id
          max_rmsg_num_a(ip) = max_rmsg_num
          idest_gpslab_a(ip,:) = idest_gpslab
       ELSE  ! IF (p_nprocs == 1)
          IF (verbose) &
            CALL message(context, 'call define_slab to fill array slab_comm_a')
          IF (.NOT. ALLOCATED(slab_comm_a)) ALLOCATE(slab_comm_a(0:p_nprocs-1))
          IF (.NOT. ALLOCATED(slab_id_a)) ALLOCATE(slab_id_a(0:p_nprocs-1))
          IF (.NOT. ALLOCATED(max_rmsg_num_a)) &
            ALLOCATE(max_rmsg_num_a(0:p_nprocs-1))
          IF (.NOT. ALLOCATED(idest_gpslab_a)) &
            ALLOCATE(idest_gpslab_a(0:p_nprocs-1,2))
          IF (.NOT. ALLOCATED(gpslab_a)) ALLOCATE(gpslab_a(0:p_nprocs-1,nslabs))
          DO ip=0, p_nprocs-1
            CALL define_slab(ip)
            slab_comm_a(ip) = slab_comm
            slab_id_a(ip) = slab_id
            max_rmsg_num_a(ip) = max_rmsg_num
            idest_gpslab_a(ip,:) = idest_gpslab
          ENDDO
       ENDIF  ! IF (p_nprocs == 1)
    ENDIF

  CONTAINS

    SUBROUTINE define_slab(p_io_c)
      ! I/O PE
      INTEGER, OPTIONAL, INTENT(in) :: p_io_c
      CHARACTER(len=*), PARAMETER :: context = &
           'mo_tr_gather::init_tr_gather::define_slab'
      INTEGER :: d, s, j, k, glats, glate
      INTEGER :: ia, ip, ip1, ip2, iy, k1, k2
      INTEGER ::  my_ip, io_ip, io_set_a, io_set_b1, io_set_b2
      INTEGER :: ypos2ip(2*na), ypos2k(2*na)
      INTEGER :: yocc(ldc%nlat), src_n, ir, n, nio
      INTEGER :: t_b_rank, ierror, ipos, sndiff
      INTEGER :: ranks_slab(1)

     IF (verbose) &
           CALL message(context, 'start define_slabs')

      io_ip     = pe2indx(p_io_c)
      io_set_a  = gdc(io_ip)%set_a
      io_set_b1 = gdc(io_ip)%set_b
      IF (io_set_b1 < 1 .OR.  io_set_b1 > ldc%nprocb) &
           CALL finish(context,'bad io_ip coords')

      IF (na <= s3_limit .AND. nb > 1) THEN
        single_stripe_slab = .TRUE.
        nslabs             = 2*na
        io_set_b2          = MOD(io_set_b1,ldc%nprocb)+1
        sndiff             = 1
      ELSE
        single_stripe_slab = .FALSE.
        nslabs             = na
        io_set_b2          = io_set_b1
        sndiff             = 2
      ENDIF
      IF (verbose) THEN
        WRITE(message_text,'(a,4i6)') &
             'io_ip, io_set_a, io_set_b1, io_set_b2=', &
             io_ip, io_set_a, io_set_b1, io_set_b2
        CALL message(context, message_text)
      ENDIF

      ! help arrays: walk south on mapmesh, note region k and gdc index ip
      DO ia = 1, na
        ! first region:
        ypos2k(ia)  = 1
        ypos2ip(ia) = ldc%mapmesh(io_set_b1,ia)

        ! second region:
        ypos2k(na+ia) = 2
        ypos2ip(na+ia) = ldc%mapmesh(io_set_b2,na-ia+1)
      ENDDO

      my_ip = pe2indx(p_pe)

      ! glats, glate , nglat, size, disp,  owner_pe:
      iy  = 0
      d   = 0
      s   = 0
      nio = 0
      DO ia = 1, nslabs
        IF (single_stripe_slab) THEN
          ip1  = ypos2ip(ia)
          ip2  = ip1
          ipos = ia
          k1   = ypos2k(ia)
          k2   = k1
        ELSE
          ! this is going to be our slab owner
          ip1  = ypos2ip(iy+1)
          ! use domain border IPs as order for collective gather
          ipos = ldc%mapmesh(1,gdc(ip1)%set_a)
          IF (ipos < 1 .OR. ipos > na) &
               CALL finish(context,'ipos < 1 or ipos > na')
          ! second stripe of slab
          ip2 = ypos2ip(iy+2)
          k1  = ypos2k(iy+1)
          k2  = ypos2k(iy+2)
          iy = iy+2
        ENDIF

        gpslab_a(p_io_c,ipos)%glats = gdc(ip1)%glats(k1)
        gpslab_a(p_io_c,ipos)%glate = gdc(ip2)%glate(k2)

      IF (verbose) THEN
          WRITE(message_text,'(a,5i6)') &
               'ia, ip1, ipos, ilats, ilate=', &
               ia, ip1, ipos, gdc(ip1)%glats(k1), gdc(ip2)%glate(k2)
          CALL message(context, message_text, p_io_c=p_io_c)
        ENDIF

        gpslab_a(p_io_c,ipos)%nglat      = gpslab_a(p_io_c,ipos)%glate-gpslab_a(p_io_c,ipos)%glats+1
        gpslab_a(p_io_c,ipos)%owner_ip   = ip1
        gpslab_a(p_io_c,ipos)%owner_pe   = gdc(ip1)%pe
        gpslab_a(p_io_c,ipos)%owner_b_pe = ipos-1

        d = ldc%nlon * (gpslab_a(p_io_c,ipos)%glats-1)
        s = ldc%nlon * gpslab_a(p_io_c,ipos)%nglat

        gpslab_a(p_io_c,ipos)%disp2 = d
        gpslab_a(p_io_c,ipos)%size2 = s

        IF (gpslab_a(p_io_c,ipos)%owner_ip == io_ip) nio = nio+1
      ENDDO

      ! make sure we have included the io task
      IF (nio /= 1) &
           CALL finish(context,'(nio /= 1)')


      ! set slab_id, check if we do not own more than one slab:
      n       = 0
      slab_id = 0
      DO ia = 1, nslabs
        IF (gpslab_a(p_io_c,ia)%owner_pe == p_pe) THEN
          n       = n+1
          slab_id = ia
        ENDIF
      ENDDO
      IF (n > 1) &
           CALL finish(context,'more than one slab per pe')

      IF (single_stripe_slab) THEN
        ! generate special communicator
        CALL gen_my_comm(slab_comm, gpslab_a(p_io_c,:)%owner_pe)
      ELSE
        slab_comm = p_communicator_b
      ENDIF

      ! check if our rank assumptions are correct (if t_b_rank is indeed equal
      ! to p_io_c on sub-communicator):
      IF (slab_id > 0) THEN
        CALL MPI_COMM_RANK(slab_comm, t_b_rank, ierror)
        IF (verbose) THEN
          WRITE(message_text,'(a,2i6)') &
               'slab_id, t_b_rank=', &
               slab_id, t_b_rank
          ! provide debug info only on PE of interest, not via p_io
          CALL message(context, message_text, p_io_c=p_io_c)
        ENDIF
        IF (p_pe == p_io_c) THEN
          CALL check_my_comm(p_io_c, slab_comm, ranks_slab)
          IF (verbose) THEN
            WRITE(message_text,'(a,3i6)') &
                 'p_io, comm_b->p_io, p_io(group)=', &
                 p_io_c, t_b_rank, ranks_slab(1)
            CALL message(context, message_text, p_io_c=p_io_c)
          ENDIF
          IF (ranks_slab(1) /= t_b_rank) &
               CALL finish(context,'p_io /= t_b_rank')
        ENDIF
        IF ((t_b_rank /= slab_id - 1) &
             .OR. (t_b_rank /= gpslab_a(p_io_c,slab_id)%owner_b_pe)) THEN
          ! in this case, we would have to reorder the slab array
          CALL finish(context,'cannot use collective gather due to rank mismatch')
        ENDIF
      ENDIF

      ! check completness:
      yocc = 0
      DO ia = 1, nslabs
        IF (verbose) THEN
          WRITE(message_text,'(a,3i6)') &
               'ia, gpslab_a(ia)%glats, gpslab_a(ia)%glate=', &
               ia, gpslab_a(p_io_c,ia)%glats, gpslab_a(p_io_c,ia)%glate
          CALL message(context, message_text)
        ENDIF
        DO j = gpslab_a(p_io_c,ia)%glats, gpslab_a(p_io_c,ia)%glate
          yocc(j) = yocc(j)+1
          IF (yocc(j) /= 1) &
               CALL finish(context,'yocc(j) /= 1')
        ENDDO
      ENDDO
      IF (SUM(yocc) /= SIZE(yocc)) &
           CALL finish(context,'SUM(yocc) /= SIZE(yocc)')

      ! destination list:
      idest_gpslab = -1
      IF (.NOT. (p_pe == p_io_c .AND. retard_p_io_data)) THEN
        DO k = 1, 2
          glats = gdc(my_ip)%glats(k)
          glate = gdc(my_ip)%glate(k)
          DO ia = 1, nslabs
            IF ((glats >= gpslab_a(p_io_c,ia)%glats) &
                 .AND. (glate <= gpslab_a(p_io_c,ia)%glate)) THEN
              IF (idest_gpslab(k) > 0) &
                   CALL finish(context,'slabs overlap')
              idest_gpslab(k) = ia
            ENDIF
          ENDDO
          IF (idest_gpslab(k) < 0) &
               CALL finish(context,'slabs not complete')
          IF (verbose) THEN
            WRITE(message_text,'(a,3i6)') &
                 'my_ip, k, idest_gpslab=', &
                 my_ip, k, idest_gpslab(k)
            CALL message(context, message_text)
          ENDIF
        ENDDO
      ENDIF

      ! source lists (only needed at owner):
      DO ia = 1, nslabs
        IF (gpslab_a(p_io_c,ia)%owner_pe == p_pe) THEN
          gpslab_a(p_io_c,ia)%src_n = sndiff*nb
          ALLOCATE(gpslab_a(p_io_c,ia)%src_ip(sndiff*nb))
          ALLOCATE(gpslab_a(p_io_c,ia)%src_ir(sndiff*nb))
          glats = gpslab_a(p_io_c,ia)%glats
          glate = gpslab_a(p_io_c,ia)%glate
          src_n = 0
          DO ip = 1, na*nb
            DO ir = 1, 2
              IF (glats <= gdc(ip)%glats(ir) &
                   .AND. gdc(ip)%glate(ir) <= glate) THEN
                IF (.NOT. (ip == io_ip .AND. retard_p_io_data)) THEN
                  src_n = src_n+1
                  IF (src_n > sndiff*nb) &
                       CALL finish(context,'gpslab_a(ia)%src_ip too small')
                  gpslab_a(p_io_c,ia)%src_ip(src_n) = ip
                  gpslab_a(p_io_c,ia)%src_ir(src_n) = ir
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          IF ( .NOT. (retard_p_io_data)) THEN
            WRITE(0,*) 'src_n, sndiff*nb=',src_n, sndiff*nb
            IF (src_n /= sndiff*nb) &
                 CALL finish(context,'gpslab_a(ia)%src_ip: size mismatch')
          END IF
          gpslab_a(p_io_c,ia)%src_n = src_n
          IF (verbose) THEN
            DO ip = 1, SIZE(gpslab_a(p_io_c,ia)%src_ip)
              WRITE(message_text,'(a,3i6)') &
                   'ia,my_ip,src_ip=',ia,my_ip,gpslab_a(p_io_c,ia)%src_ip(ip)
              CALL message(context, message_text)
              WRITE(message_text,'(a,3i6)') &
                   'ia,my_ip,src_ir=',ia,my_ip,gpslab_a(p_io_c,ia)%src_ir(ip)
              CALL message(context, message_text)
            ENDDO
          ENDIF
        ELSE
          gpslab_a(p_io_c,ia)%src_n = 0
        ENDIF
      ENDDO

      ! check if each slab owner's private data overlabs with his slab:
      IF (slab_id > 0) THEN
        IF (.NOT. (p_pe == p_io_c .AND. retard_p_io_data)) THEN
          glats = gpslab_a(p_io_c,slab_id)%glats
          glate = gpslab_a(p_io_c,slab_id)%glate
          DO ir = 1, 2
            IF (glats <= ldc%glats(ir) .AND. ldc%glate(ir) <= glate) n = n+1
          ENDDO
          IF (n == 0) &
               CALL finish(context, &
                 'no selfoverlap between slab and private data of slab_owner')
        ENDIF
      ENDIF

      ! calculate max msg-size, msg-num
      max_rmsg_size_2d = 0
      DO ip = 1, ldc%d_nprocs
        DO k = 1, 2
          n = gdc(ip)%nglon*gdc(ip)%nglh(k)
          max_rmsg_size_2d = MAX(max_rmsg_size_2d, n)
        ENDDO
      ENDDO
      IF (slab_id > 0) THEN
        max_rmsg_num = gpslab_a(p_io_c,slab_id)%src_n
      ELSE
        max_rmsg_num = 0
      ENDIF
      IF (verbose) THEN
        WRITE(message_text,'(a,2i6)') &
             'max_msg_size_2d, max_msg_num=', &
             max_rmsg_size_2d, max_rmsg_num
        CALL message(context,'end define_slabs')
      ENDIF

    END SUBROUTINE define_slab

    SUBROUTINE gen_my_comm(new_comm, ranks)
      INTEGER, INTENT(out) :: new_comm
      INTEGER, INTENT(in) :: ranks(:)

      INTEGER :: all_group, ierr, new_group

      CALL MPI_COMM_GROUP(p_all_comm, all_group, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::gen_my_comm: bad case (1)')

      CALL MPI_GROUP_INCL(all_group, SIZE(ranks), ranks, new_group, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::gen_my_comm: bad case (2)')

      CALL MPI_COMM_CREATE(p_all_comm, new_group, new_comm, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::gen_my_comm: bad case (3)')
    END SUBROUTINE gen_my_comm

#endif

  END SUBROUTINE init_tr_gather
  !-----------------------------------------------------------------------------
#ifdef _CHECK_MGATHER
  SUBROUTINE cmp2(a, b)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::cmp2'
    REAL(dp), INTENT(in) :: a(:,:)
    REAL(dp), INTENT(in) :: b(:,:)

    INTEGER :: i1, i2

    ! only for debugging
    IF (SIZE(a,1) /= SIZE(b,1)) CALL finish(context,'size error 1')
    IF (SIZE(a,2) /= SIZE(b,2)) CALL finish(context,'size error 2')

    DO i2 = 1, SIZE(a,2)
      DO i1 = 1, SIZE(a,1)
        IF (a(i1,i2) /= b(i1,i2) ) THEN
          WRITE(message_text,'(a,4i6)') &
               'cmp2: i1,i2, a, b=', i1,i2, a(i1,i2), b(i1,i2)
          CALL message(context, message_text)
          CALL finish(context,'ERROR')
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE cmp2
  !-----------------------------------------------------------------------------
  SUBROUTINE cmp3(a, b)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::cmp3'
    REAL(dp), INTENT(in) :: a(:,:,:)
    REAL(dp), INTENT(in) :: b(:,:,:)

    INTEGER :: i1, i2, i3

    ! only for debugging
    IF (SIZE(a,1) /= SIZE(b,1)) CALL finish(context,'size error 1')
    IF (SIZE(a,2) /= SIZE(b,2)) CALL finish(context,'size error 2')

    DO i3 = 1, SIZE(a,3)
      DO i2 = 1, SIZE(a,2)
        DO i1 = 1, SIZE(a,1)
          IF (a(i1,i2,i3) /= b(i1,i2,i3) ) THEN
            WRITE(message_text,'(a,5i6)') &
                 'cmp3: i1,i2,i3, a, b=', &
                 i1, i2, i3, a(i1,i2,i3), b(i1,i2,i3)
            CALL message(context, message_text)
            CALL finish(context,'ERROR')
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE cmp3
#endif
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_gp432(global, global_dims, local, source, lall_write, p_io_c)
    REAL(dp), POINTER                :: global(:,:,:,:)
    INTEGER,           INTENT(in)    :: global_dims(:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c

    REAL(dp), POINTER :: global3d(:,:,:)
    REAL(dp), POINTER :: global3d_c(:,:,:) ! mz_ht_20120123

    INTEGER :: i, nt
    INTEGER :: dimsize4
    INTEGER :: dimsize3 ! mz_ht_20120123
    INTEGER :: p_io_up

    dimsize4 = global_dims(4)

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    NULLIFY(global3d)
    NULLIFY(global3d_c) ! mz_ht_20120123

    IF (dimsize4 == 1) THEN
      IF (p_pe == p_io_up) global3d => global(:,:,:,1)
      CALL gather_gp32 (global3d, global_dims, local(:,:,:,1), source, &
           p_io_c=p_io_c)
    ELSE
      nt = global_dims(3)
! mz_ht_20120123+
      IF (p_pe == p_io_up) &
           ALLOCATE(global3d_c(global_dims(1),global_dims(2),global_dims(4)))
! mz_ht_20120123-
      DO i = 1, nt
! mz_ht_20120123+
!!$        IF (p_pe == p_io) global3d => global(:,:,i,:)
!!$        CALL gather_gp32 (global3d, global_dims, local(:,:,i,:), source)
         ! op_pj_20120125+
         ! NOTE: This copy is required, because MPI_GATHERV operates with
         !       start adresses and offsets, which require consecutive
         !       data layouts in memory. This prerequisite is violated,
         !       however, if rank 3 (instead of 4) is cut out.
         ! op_pj_20120125-
         IF (p_pe == p_io_up) global3d_c(:,:,:) = global(:,:,i,:)
         dimsize3 = global_dims(4)
         CALL gather_gp32 (global3d_c, global_dims, local(:,:,i,:), source, &
                    dimsize3, p_io_c=p_io_c)
         IF (p_pe == p_io_up) global(:,:,i,:) = global3d_c(:,:,:)
! mz_ht_20120123-
      END DO
! mz_ht_20120123+
      IF (ASSOCIATED(global3d_c)) DEALLOCATE(global3d_c)
! mz_ht_20120123-
    ENDIF

  END SUBROUTINE gather_gp432
  !-----------------------------------------------------------------------------
! mz_ht_20120123+
!!$  SUBROUTINE gather_gp32(global, global_dims, local, source, lall_write)
  SUBROUTINE gather_gp32(global, global_dims, local, source, pdimsize3, lall_write, p_io_c)
! mz_ht_20120123-
    REAL(dp), POINTER                :: global(:,:,:)
    INTEGER,           INTENT(in)    :: global_dims(:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    INTEGER, OPTIONAL, INTENT(in)    :: pdimsize3  ! mz_ht_20120123
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c
    REAL(dp), POINTER :: global2d(:,:)

    INTEGER :: dimsize3
    INTEGER :: p_io_up

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    dimsize3 = global_dims(3)
    IF (PRESENT(pdimsize3)) dimsize3 = pdimsize3 ! mz_ht_20120123

    IF (dimsize3 == 1) THEN
      NULLIFY(global2d)
      IF (p_pe == p_io_up) global2d => global(:,:,1)
#ifdef _USE_MGATHER
      CALL mgather_gp2(global2d, local(:,:,1), source, p_io_c=p_io_c)
#else
      CALL gather_gp2(global2d, local(:,:,1), source, p_io_c=p_io_c)
#endif
    ELSE
#ifdef _USE_MGATHER
      CALL mgather_gp3(global, local(:,:,:), source, p_io_c=p_io_c)
#else
      CALL gather_gp3(global, local(:,:,:), source, p_io_c=p_io_c)
#endif
    ENDIF

  END SUBROUTINE gather_gp32
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_gp2(global, local, source, lall_write, p_io_c)
#ifdef STANDALONE
    USE mo_jsbach_transpose, ONLY: jsb_gather_gp2
#endif
    REAL(dp), POINTER                :: global(:,:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c

    REAL(dp), POINTER :: sendbuf(:,:)
    REAL(dp), POINTER :: recvbuf(:,:)

    INTEGER :: i, pe, src, nlon
    LOGICAL :: lreg
    INTEGER :: p_io_up

#ifdef STANDALONE
    CALL jsb_gather_gp2(global, local(:,:), source)
    RETURN
#endif

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    IF (require_init) THEN
      CALL init_tr_gather()
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (measure_mode) CALL p_barrier
    IF (timings) CALL timer_start(timer_gp2)
#endif

    src = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    nlon = ldc%nlon
    lreg = ldc%lreg

    IF (lreg) THEN
      sendbuf => local
    ELSE
      ALLOCATE(sendbuf(ldc%nglon,ldc%nglat))
      CALL reorder(sendbuf,local)
    ENDIF

    IF (p_pe /= p_io_up) THEN
      CALL p_send(sendbuf(:,:), p_io_up, tag_gather_gp)
    ELSE
      DO i = 1, p_nprocs
        pe    = gdc(i)%pe
        ALLOCATE(recvbuf(gdc(i)%nglon,gdc(i)%nglat))
        IF (pe /= p_pe) THEN
          CALL p_recv(recvbuf(:,:), pe, tag_gather_gp)
        ELSE
          recvbuf(:,:) = sendbuf(:,:)
        ENDIF
        IF (src == -1 .OR. (src == 0 .AND. pe == p_io_up) &
             .OR. (src == 1 .AND. pe /= p_io_up)) THEN
#ifdef _USE_GATHER_TIMER
          IF (extended_timings) CALL timer_start(timer_fill_gp2)
#endif
          ! unpack first segment
          global(gdc(i)%glons(1):gdc(i)%glone(1),gdc(i)%glats(1):gdc(i)%glate(1)) &
               = recvbuf(:,:gdc(i)%nglh(1))
          ! unpack second segment
          IF (gdc(i)%nglh(2) > 0) THEN
            IF (gdc(i)%glone(2) > gdc(i)%glons(2)) THEN
              global(gdc(i)%glons(2):gdc(i)%glone(2),gdc(i)%glats(2):gdc(i)%glate(2)) &
                   = recvbuf(:,gdc(i)%nglat/2+1:)
            ELSE
              ! unpacking second segment, split in longitudes
              global(gdc(i)%glons(2):nlon,gdc(i)%glats(2):gdc(i)%glate(2)) &
                   = recvbuf(:nlon-gdc(i)%glons(2)+1,gdc(i)%nglat/2+1:)
              global(:gdc(i)%glone(2),gdc(i)%glats(2):gdc(i)%glate(2))     &
                   = recvbuf(gdc(i)%nglon-gdc(i)%glone(2)+1:,gdc(i)%nglat/2+1:)
            ENDIF
          ENDIF
#ifdef _USE_GATHER_TIMER
          IF (extended_timings) CALL timer_stop(timer_fill_gp2)
#endif
        ENDIF
        DEALLOCATE(recvbuf)
      ENDDO
    ENDIF

    IF (lreg) THEN
      NULLIFY (sendbuf)
    ELSE
      DEALLOCATE (sendbuf)
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_stop(timer_gp2)
#endif

  END SUBROUTINE gather_gp2
  !-----------------------------------------------------------------------------
#ifdef _USE_MGATHER
  SUBROUTINE subgather_gp2(local, slab, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::subgather_gp2'
    ! blocked nglon x nglat => unblocked nglon x nglat
    REAL(dp),           INTENT(in)  :: local(ldc%nglon,ldc%nglat)
    ! (nlon, latsection)
    REAL(dp), OPTIONAL, INTENT(out) :: slab(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: p_io_c

    INTEGER, PARAMETER :: tag(2) = (/ tag_gp_sub1, tag_gp_sub2 /)
    INTEGER :: nreq
    INTEGER, ALLOCATABLE :: ireq(:), istat(:, :)

    INTEGER :: my_ip, ilat1(2),pn(2), ir, islab, islab_ip, islab_pe
    INTEGER :: ierr, isrc, nlon, ip, msize, latshift
    INTEGER :: p_io_up

    REAL(dp), ALLOCATABLE:: rbuf(:,:)

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    ALLOCATE(ireq(2+max_rmsg_num_a(p_io_up)))
    ALLOCATE(istat(MPI_STATUS_SIZE, 2+max_rmsg_num_a(p_io_up)))
    ALLOCATE(rbuf(max_rmsg_size_2d, max_rmsg_num_a(p_io_up)))

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_sub_gp2)
#endif

    nlon = ldc%nlon

    my_ip  = pe2indx(p_pe)

    ! first region:
    ilat1(1) = 1
    pn(1)    = ldc%nglon*ldc%nglh(1)

    ! second region:
    ilat1(2) = ldc%nglh(1)+1
    pn(2)    = ldc%nglon*ldc%nglh(2)

    IF (slab_id_a(p_io_up) > 0) THEN
      latshift = gpslab_a(p_io_up,slab_id_a(p_io_up))%glats-1
    ELSE
      latshift = 0 ! not used
    ENDIF

    nreq = 0
    DO ir = 1, 2
      islab = idest_gpslab_a(p_io_up,ir)
      IF (islab < 1) CYCLE
      islab_ip = gpslab_a(p_io_up,islab)%owner_ip
      IF (my_ip /= islab_ip) THEN
        islab_pe = gpslab_a(p_io_up,islab)%owner_pe
        nreq = nreq +1
        CALL MPI_ISEND(local(1,ilat1(ir)), pn(ir), p_real_dp, islab_pe, tag(ir), &
                       p_all_comm, ireq(nreq), ierr)
      ELSE
        ! "send" to self
        CALL fill_gp2(my_ip, ir, latshift, local(1,ilat1(ir)), slab)
      ENDIF
    ENDDO

    IF (slab_id_a(p_io_up) > 0) THEN
      DO isrc = 1, gpslab_a(p_io_up,slab_id_a(p_io_up))%src_n
        ip = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ip(isrc)
        ir = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ir(isrc)
        msize = gdc(ip)%nglon*gdc(ip)%nglh(ir)
        IF (my_ip /= ip) THEN
          nreq = nreq +1
          CALL MPI_IRECV(rbuf(1,isrc), msize, p_real_dp, gdc(ip)%pe, tag(ir), &
                         p_all_comm, ireq(nreq), ierr)
        ENDIF
      ENDDO
    ENDIF

    CALL MPI_WAITALL(nreq, ireq, istat, ierr)

    IF (slab_id_a(p_io_up) > 0) THEN
      islab_ip = gpslab_a(p_io_up,slab_id_a(p_io_up))%owner_ip
      DO isrc = 1, gpslab_a(p_io_up,slab_id_a(p_io_up))%src_n
        ip = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ip(isrc)
        ir = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ir(isrc)
        IF (my_ip /= ip) CALL fill_gp2(ip, ir, latshift, rbuf(1,isrc), slab)
      ENDDO
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_sub_gp2)
#endif

    DEALLOCATE(ireq)
    DEALLOCATE(istat)
    DEALLOCATE(rbuf)

  END SUBROUTINE subgather_gp2
  !-----------------------------------------------------------------------------
  SUBROUTINE fill_gp2(ip, ir, glat0, buf, slab)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::fill_gp2'
    INTEGER,  INTENT(in)    :: ip, ir, glat0
    REAL(dp), INTENT(in)    :: buf(gdc(ip)%nglon,gdc(ip)%nglh(ir))
    ! (nlon, latsection)
    REAL(dp), INTENT(inout) :: slab(:,:)

    INTEGER :: nlon, glons, glone, glats, glate, nglon

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_mfill_gp2)
#endif

    nlon = ldc%nlon

    glons = gdc(ip)%glons(ir)
    glone = gdc(ip)%glone(ir)
    glats = gdc(ip)%glats(ir)-glat0
    glate = gdc(ip)%glate(ir)-glat0
    nglon = gdc(ip)%nglon
    IF (ir == 1) THEN
      ! unpack first segment
      slab(glons:glone, glats:glate) = buf
    ELSE
      ! unpack second segment
      IF (glone > glons) THEN
        slab(glons:glone, glats:glate) = buf
      ELSE
        ! unpacking second segment, split in longitudes
        slab(glons:nlon, glats:glate) = buf(:nlon-glons+1 ,:)
        slab(:glone    , glats:glate) = buf(nglon-glone+1:,:)
      ENDIF
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_mfill_gp2)
#endif

  END SUBROUTINE fill_gp2
  !-----------------------------------------------------------------------------
  SUBROUTINE slabgather_gp2(slab, global, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::slabgather_gp2'
    ! global is only accessed at p_io, we do a complete overwrite
    ! this does not accept a nonassociated pointer
    REAL(dp), OPTIONAL, INTENT(out) :: global(ldc%nlon,ldc%nlat)
    REAL(dp),           INTENT(in)  :: slab(:,:)
    INTEGER,  OPTIONAL, INTENT(in)  :: p_io_c

    INTEGER  :: ierr, is
    INTEGER  :: ireq(ldc%nproca), istat(MPI_STATUS_SIZE, ldc%nproca), nreq
    REAL(dp) :: gtmp(1,1)
    INTEGER  :: p_io_up, ranks_slab(1)

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_slab_gp2)
#endif

    ! we have three methods for the second phase: collective gather,
    ! blocking recv and nonblocking recv

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    IF (slab_id_a(p_io_up) > 0) THEN

      SELECT CASE (sel_slabgather)

      CASE (sel_gatherv)

        IF (p_pe == p_io_up) THEN
          CALL check_my_comm(p_io_up, slab_comm_a(p_io_up), ranks_slab)
          CALL MPI_GATHERV(slab, gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, global, &
                           gpslab_a(p_io_up,:)%size2, gpslab_a(p_io_up,:)%disp2,               &
                           p_real_dp, ranks_slab(1), slab_comm_a(p_io_up), ierr)
        ELSE
          ! the compiler doesn't know that global is not used by
          ! mpi_gatherv and requires a valid address, so we use gtmp here
          CALL check_my_comm(p_io_up, slab_comm_a(p_io_up), ranks_slab)
          CALL MPI_GATHERV(slab, gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, gtmp, &
                           gpslab_a(p_io_up,:)%size2, gpslab_a(p_io_up,:)%disp2,             &
                           p_real_dp, ranks_slab(1), slab_comm_a(p_io_up), ierr)
        ENDIF

      CASE (sel_irecv)

        IF (p_pe /= p_io_up) THEN
          CALL MPI_SEND(slab, gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, p_io_up, tag_gp_slab, &
                        p_all_comm, ierr)
        ELSE
          ! p_pe == p_io:
          nreq = 0
          DO is = 1, nslabs
            IF (gpslab_a(p_io_up,is)%owner_pe == p_io_up) THEN
               global(:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate) = slab
            ELSE
              nreq = nreq+1
              IF (.NOT. PRESENT(global)) &
                   CALL finish(context,'internal error')
              CALL MPI_IRECV(global(:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate), &
                             gpslab_a(p_io_up,is)%size2, p_real_dp,                 &
                             gpslab_a(p_io_up,is)%owner_pe, tag_gp_slab,            &
                             p_all_comm, ireq(nreq), ierr)
            ENDIF
          ENDDO
          CALL MPI_WAITALL(nreq, ireq, istat, ierr)
        ENDIF

      CASE (sel_recv)

        IF (p_pe /= p_io_up) THEN
          CALL p_send(slab, p_io_up, tag_gp_slab)
        ELSE
          ! p_pe == p_io:
          DO is = 1, nslabs
            IF (gpslab_a(p_io_up,is)%owner_pe == p_io_up) THEN
               global(:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate) = slab
            ELSE
              CALL p_recv(global(:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate), &
                          gpslab_a(p_io_up,is)%owner_pe,tag_gp_slab)
            ENDIF
          ENDDO
        ENDIF

      END SELECT
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_slab_gp2)
#endif

  END SUBROUTINE slabgather_gp2
  !-----------------------------------------------------------------------------
  SUBROUTINE mgather_gp2(global, local, source, lall_write, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::mgather_gp2'
    REAL(dp), POINTER                :: global(:,:)
    REAL(dp),          INTENT(in)    :: local(:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write  ! not used yet (?)
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c

    REAL(dp), ALLOCATABLE :: slab(:,:)
    INTEGER :: src
#ifdef _CHECK_MGATHER
    REAL(dp), POINTER :: ref_global(:,:)
#endif
    INTEGER :: p_io_up

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    IF (require_init) THEN
      CALL init_tr_gather()
    ENDIF
#ifdef _USE_GATHER_TIMER
    IF (measure_mode) CALL p_barrier
#endif


    src = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    ! use general gather for debugging:
    IF ((src /= -1) .OR. p_nprocs /= ldc%d_nprocs) THEN
      CALL gather_gp2(global, local, source, lall_write, p_io_c=p_io_c)
      RETURN
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_start(timer_mgp2)
#endif

    ! here we limit ourselves to the intra domain case:
    ! (src == -1) .AND. (p_nprocs == ldc%d_nprocs)
    ! so no debug support yet
#ifdef NOMPI
    CALL reorder(global,local)
#else
    IF (slab_id_a(p_io_up) > 0) THEN
      !IF (ALLOCATED(slab)) DEALLOCATE(slab)
      ALLOCATE(slab(ldc%nlon,gpslab_a(p_io_up,slab_id_a(p_io_up))%nglat))
      CALL subgather_gp2(local, slab, p_io_c=p_io_c)
    ELSE
      CALL subgather_gp2(local, p_io_c=p_io_c)
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (measure_mode2) CALL p_barrier
#endif

    IF (slab_id_a(p_io_up) > 0) THEN
      IF (p_pe == p_io_up) THEN
        CALL slabgather_gp2(slab, global, p_io_c=p_io_c)
        ! now we have to insert our delayed own data
        IF (retard_p_io_data) CALL my_late_injection(local)
      ELSE
        CALL slabgather_gp2(slab, p_io_c=p_io_c)
      ENDIF
    ENDIF
#endif

#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_stop(timer_mgp2)
#endif

#ifdef _CHECK_MGATHER
    IF (p_pe == p_io_up) THEN
      ALLOCATE(ref_global(ldc%nlon,ldc%nlat))
    ELSE
      NULLIFY(ref_global)
    ENDIF
    CALL gather_gp2(ref_global,local,source, p_io_c=p_io_c)
    IF (p_pe == p_io_up) THEN
      CALL message(context,'check')
      CALL cmp2(ref_global, global)
      DEALLOCATE(ref_global)
    ENDIF
#endif

  CONTAINS

    SUBROUTINE my_late_injection(buf)
      CHARACTER(len=*), PARAMETER :: context = &
           'mo_tr_gather::mgather_gp2::my_late_injection'
      REAL(dp), INTENT(in)  :: buf(ldc%nglon,ldc%nglat)

      INTEGER :: my_ip, p1, p2

#ifdef _USE_GATHER_TIMER
      IF (extended_timings) CALL timer_start(timer_inject2)
#endif

      my_ip = pe2indx(p_pe)
      p1    = 1
      p2    = ldc%nglh(1)+1

      CALL fill_gp2(my_ip, 1, 0, buf(1,p1), global)
      CALL fill_gp2(my_ip, 2, 0, buf(1,p2), global)

#ifdef _USE_GATHER_TIMER
      IF (extended_timings) CALL timer_stop(timer_inject2)
#endif

    END SUBROUTINE my_late_injection

  END SUBROUTINE mgather_gp2
  !-----------------------------------------------------------------------------
  SUBROUTINE subgather_gp3(nk, local, slab, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::subgather_gp3'
    INTEGER,            INTENT(in)  :: nk
    REAL(dp),           INTENT(in)  :: local(ldc%nglon,nk,ldc%nglat)
    REAL(dp), OPTIONAL, INTENT(out) :: slab(:,:,:)  ! (nlon, nk, latsection)
    INTEGER, OPTIONAL, INTENT(in) :: p_io_c

    INTEGER, PARAMETER :: tag(2) = (/ tag_gp_sub1, tag_gp_sub2 /)
    INTEGER :: nreq
    INTEGER, ALLOCATABLE :: ireq(:), istat(:,:)

    INTEGER :: my_ip, ilat1(2),pn(2), ir, islab, islab_ip, islab_pe
    INTEGER :: ierr, isrc, nlon, ip, msize, latshift
    INTEGER :: p_io_up

    REAL(dp), ALLOCATABLE:: rbuf(:,:)

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    ALLOCATE(ireq(2+max_rmsg_num_a(p_io_up)))
    ALLOCATE(istat(MPI_STATUS_SIZE, 2+max_rmsg_num_a(p_io_up)))
    ALLOCATE(rbuf(nk*max_rmsg_size_2d, max_rmsg_num_a(p_io_up)))

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_sub_gp3)
#endif

    nlon = ldc%nlon

    my_ip = pe2indx(p_pe)

    ! first region:
    ilat1(1) = 1
    pn(1)    = ldc%nglon*nk*ldc%nglh(1)

    ! second region:
    ilat1(2) = ldc%nglh(1)+1
    pn(2)    = ldc%nglon*nk*ldc%nglh(2)

    IF (slab_id_a(p_io_up) > 0) THEN
      latshift = gpslab_a(p_io_up,slab_id_a(p_io_up))%glats-1
    ELSE
      latshift = 0 ! not used
    ENDIF

    nreq = 0
    DO ir = 1, 2
      islab = idest_gpslab_a(p_io_up,ir)
      IF (islab < 1) CYCLE
      islab_ip = gpslab_a(p_io_up,islab)%owner_ip
      IF (my_ip /= islab_ip) THEN
        islab_pe = gpslab_a(p_io_up,islab)%owner_pe
        nreq = nreq +1
        CALL MPI_ISEND(local(1,1,ilat1(ir)), pn(ir), p_real_dp, islab_pe, tag(ir), &
             &         p_all_comm, ireq(nreq), ierr)
      ELSE
        CALL fill_gp3(my_ip, ir, nk, latshift, local(1,1,ilat1(ir)), slab)
      ENDIF
    ENDDO

    IF (slab_id_a(p_io_up) > 0) THEN
      DO isrc = 1, gpslab_a(p_io_up,slab_id_a(p_io_up))%src_n
        ip    = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ip(isrc)
        ir    = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ir(isrc)
        msize = gdc(ip)%nglon*nk*gdc(ip)%nglh(ir)
        IF (my_ip /= ip) THEN
          nreq = nreq +1
          CALL MPI_IRECV(rbuf(1,isrc), msize, p_real_dp, gdc(ip)%pe, tag(ir), &
               &         p_all_comm, ireq(nreq), ierr)
        ENDIF
      ENDDO
    ENDIF

    CALL MPI_WAITALL(nreq, ireq, istat, ierr)

    IF (slab_id_a(p_io_up) > 0) THEN
      islab_ip = gpslab_a(p_io_up,slab_id_a(p_io_up))%owner_ip
      DO isrc = 1, gpslab_a(p_io_up,slab_id_a(p_io_up))%src_n
        ip = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ip(isrc)
        ir = gpslab_a(p_io_up,slab_id_a(p_io_up))%src_ir(isrc)
        IF (my_ip /= ip) CALL fill_gp3(ip, ir, nk, latshift, rbuf(1,isrc), slab)
      ENDDO
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_sub_gp3)
#endif

    DEALLOCATE(ireq)
    DEALLOCATE(istat)
    DEALLOCATE(rbuf)

  END SUBROUTINE subgather_gp3
  !-----------------------------------------------------------------------------
  SUBROUTINE fill_gp3(ip, ir, nk, glat0, buf, slab)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::fill_gp3'
    INTEGER,  INTENT(in)    :: ip, ir, nk, glat0
    REAL(dp), INTENT(in)    :: buf(gdc(ip)%nglon,nk,gdc(ip)%nglh(ir))
    REAL(dp), INTENT(inout) :: slab(:,:,:)  ! (nlon, nk, latsection)

    INTEGER:: nlon, glons, glone, glats, glate, nglon

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_mfill_gp3)
#endif

    nlon = ldc%nlon

    glons = gdc(ip)%glons(ir)
    glone = gdc(ip)%glone(ir)
    glats = gdc(ip)%glats(ir)-glat0
    glate = gdc(ip)%glate(ir)-glat0
    nglon = gdc(ip)%nglon

    IF (ir == 1) THEN
      ! unpack first segment
      slab(glons:glone, :,glats:glate) = buf
    ELSE
      ! unpack second segment
      IF (glone > glons) THEN
        slab(glons:glone, :,glats:glate) = buf
      ELSE
        ! unpacking second segment, split in longitudes
        slab(glons:nlon, :, glats:glate) = buf(:nlon-glons+1 , :, :)
        slab(:glone    , :, glats:glate) = buf(nglon-glone+1:, :, :)
      ENDIF
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_mfill_gp3)
#endif

  END SUBROUTINE fill_gp3
  !-----------------------------------------------------------------------------
  SUBROUTINE slabgather_gp3(nk,slab,global, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::slabgather_gp3'
    INTEGER,            INTENT(in)  :: nk
    ! global is only accessed at p_io, we do a complete overwrite
    REAL(dp), OPTIONAL, INTENT(out) :: global(ldc%nlon,nk,ldc%nlat)
    REAL(dp),           INTENT(in)  :: slab(:,:,:)
    INTEGER, OPTIONAL, INTENT(in)  :: p_io_c

    ! mz_ht_20111128+
!!$ INTEGER  :: ierr, is
    INTEGER  :: ierr = 0
    INTEGER  :: is
    ! mz_ht_20111128-
    INTEGER  :: ireq(ldc%nproca), istat(MPI_STATUS_SIZE, ldc%nproca),nreq
    REAL(dp) :: gtmp(1,1) = 0._dp ! mz_hz_20120125 = 0._dp added
    INTEGER  :: p_io_up, ranks_slab(1)

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_start(timer_slab_gp3)
#endif

    ! we have three methods for the second phase: collective gather,
    ! blocking recv and nonblocking recv
    IF (slab_id_a(p_io_up) > 0) THEN

      SELECT CASE (sel_slabgather)

      CASE (sel_gatherv)

        IF (p_pe == p_io_up) THEN
          CALL check_my_comm(p_io_up, slab_comm_a(p_io_up), ranks_slab)
          CALL MPI_GATHERV(slab, nk*gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, global, &
                           nk*gpslab_a(p_io_up,:)%size2, nk*gpslab_a(p_io_up,:)%disp2, &
                           p_real_dp, ranks_slab(1), slab_comm_a(p_io_up), ierr)
        ELSE
          ! the compiler doesn't know that global is not used by mpi_gatherv
          ! and requires a valid address, so we use gtmp here
          CALL check_my_comm(p_io_up, slab_comm_a(p_io_up), ranks_slab)
          CALL MPI_GATHERV(slab, nk*gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, gtmp, &
                           nk*gpslab_a(p_io_up,:)%size2, nk*gpslab_a(p_io_up,:)%disp2, &
                           p_real_dp, ranks_slab(1), slab_comm_a(p_io_up), ierr)
        ENDIF

      CASE (sel_irecv)

        IF (p_pe /= p_io_up) THEN
          CALL MPI_SEND(slab, nk*gpslab_a(p_io_up,slab_id_a(p_io_up))%size2, p_real_dp, p_io_up, tag_gp_slab, &
                        p_all_comm, ierr)
        ELSE
          ! p_pe == p_io_up:
          nreq = 0
          DO is = 1, nslabs
            IF (gpslab_a(p_io_up,is)%owner_pe == p_io_up) THEN
              global(:,:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate)=slab
            ELSE
              nreq = nreq+1
              IF (.NOT. PRESENT(global)) CALL finish(context,'internal error')
              CALL MPI_IRECV(global(:,:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate), &
                             nk*gpslab_a(p_io_up,is)%size2, p_real_dp,                &
                             gpslab_a(p_io_up,is)%owner_pe, tag_gp_slab,              &
                             p_all_comm, ireq(nreq), ierr)
            ENDIF
          ENDDO
          CALL MPI_WAITALL(nreq, ireq, istat, ierr)
        ENDIF

      CASE (sel_recv)

        IF (p_pe /= p_io_up) THEN
          CALL p_send(slab, p_io_up, tag_gp_slab)
        ELSE
          ! p_pe == p_io:
          DO is = 1, nslabs
            IF (gpslab_a(p_io_up,is)%owner_pe == p_io_up) THEN
              global(:,:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate) = slab
            ELSE
              CALL p_recv(global(:,:,gpslab_a(p_io_up,is)%glats:gpslab_a(p_io_up,is)%glate), gpslab_a(p_io_up,is)%owner_pe, tag_gp_slab)
            ENDIF
          ENDDO
        ENDIF

      END SELECT

    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (extended_timings) CALL timer_stop(timer_slab_gp3)
#endif

  END SUBROUTINE slabgather_gp3
  !-----------------------------------------------------------------------------
  SUBROUTINE mgather_gp3(global, local, source, lall_write, p_io_c)
    CHARACTER(len=*), PARAMETER :: context = 'mo_tr_gather::mgather_gp3'
    REAL(dp), POINTER                :: global(:,:,:) ! (nlon,:,nlat)
    REAL(dp),          INTENT(in)    :: local(:,:,:)  ! blocked (nglon,nk,nglat)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write  ! not used yet (?)
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c

    REAL(dp), ALLOCATABLE :: slab(:,:,:), reg_local(:,:,:)
#ifdef _CHECK_MGATHER
    REAL(dp), POINTER :: ref_global(:,:,:)
#endif
    INTEGER :: nk, src, p_io_up

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

    IF (require_init) THEN
      CALL init_tr_gather()
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (measure_mode) CALL p_barrier
#endif

    src = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    ! use general gather for debugging:
    IF ( (src /= -1) .OR. p_nprocs /= ldc%d_nprocs ) THEN
      CALL gather_gp3(global, local, source, lall_write,p_io_c=p_io_c)
      RETURN
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_start(timer_mgp3)
#endif

#ifdef NOMPI
    CALL reorder(global,local)
#else
    ! here we limit ourselves to the intra domain case:
    ! (src == -1) .AND. (p_nprocs == ldc%d_nprocs)
    ! no debug support yet
    nk = SIZE(local,2)

    IF (slab_id_a(p_io_up) > 0) THEN
         !IF (ALLOCATED(slab)) DEALLOCATE(slab)
         ALLOCATE(slab(ldc%nlon,nk,gpslab_a(p_io_up,slab_id_a(p_io_up))%nglat))
    ENDIF
    IF (ldc%lreg) THEN
      IF (slab_id_a(p_io_up) > 0) THEN
        CALL subgather_gp3(nk, local, slab, p_io_c=p_io_c)
      ELSE
        CALL subgather_gp3(nk, local, p_io_c=p_io_c)
      ENDIF
    ELSE
      ALLOCATE(reg_local(ldc%nglon,nk,ldc%nglat))
      CALL reorder(reg_local,local)
      IF (slab_id_a(p_io_up) > 0) THEN
        CALL subgather_gp3(nk, reg_local, slab, p_io_c=p_io_c)
      ELSE
        CALL subgather_gp3(nk, reg_local, p_io_c=p_io_c)
      ENDIF
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (measure_mode2) CALL p_barrier
#endif

    IF (slab_id_a(p_io_up) > 0) THEN
      IF (p_pe == p_io_up) THEN
        IF (SIZE(local,2) /= SIZE(global,2) ) &
             CALL finish(context,'size mismatch')
        CALL slabgather_gp3(nk, slab, global, p_io_c=p_io_c)
        IF (retard_p_io_data) THEN
          IF (ldc%lreg) THEN
            CALL my_late_injection(local)
          ELSE
            CALL my_late_injection(reg_local)
          ENDIF
        ENDIF
      ELSE
        CALL slabgather_gp3(nk, slab, p_io_c=p_io_c)
      ENDIF
    ENDIF

    IF (ALLOCATED(reg_local)) DEALLOCATE(reg_local)
#endif
#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_stop(timer_mgp3)
#endif
#ifdef _CHECK_MGATHER
    IF (p_pe == p_io_up) THEN
      ALLOCATE(ref_global(ldc%nlon,nk,ldc%nlat))
    ELSE
      NULLIFY(ref_global)
    ENDIF
    CALL gather_gp3(ref_global, local, source, p_io_c=p_io_c)
    IF (p_pe == p_io_up) THEN
      CALL message(context,'check')
      CALL cmp3(ref_global, global)
      DEALLOCATE(ref_global)
    ENDIF
#endif

  CONTAINS

    SUBROUTINE my_late_injection(buf)
    CHARACTER(len=*), PARAMETER :: context = &
         'mo_tr_gather::mgather_gp3::my_late_injection'
      REAL(dp), INTENT(in) :: buf(ldc%nglon,nk,ldc%nglat)

      INTEGER :: my_ip, p1, p2

#ifdef _USE_GATHER_TIMER
      IF (extended_timings) CALL timer_start(timer_inject3)
#endif

      my_ip = pe2indx(p_pe)
      p1    = 1
      p2    = ldc%nglh(1)+1
      CALL fill_gp3(my_ip, 1, nk, 0, buf(1,1,p1), global)
      CALL fill_gp3(my_ip, 2, nk, 0, buf(1,1,p2), global)

#ifdef _USE_GATHER_TIMER
      IF (extended_timings) CALL timer_stop(timer_inject3)
#endif

    END SUBROUTINE my_late_injection

  END SUBROUTINE mgather_gp3
#endif
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_gp3(global, local, source, lall_write, p_io_c)
#ifdef STANDALONE
    USE mo_jsbach_transpose, ONLY: jsb_gather_gp3
#endif
    REAL(dp), POINTER                :: global(:,:,:)
    REAL(dp), TARGET,  INTENT(in)    :: local(:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write
    INTEGER, OPTIONAL, INTENT(in)    :: p_io_c

    REAL(dp), POINTER :: sendbuf(:,:,:)
    REAL(dp), POINTER :: recvbuf(:,:,:)

    INTEGER :: i, src, pe, nlon, nk, p_io_up
    LOGICAL :: lreg

    p_io_up = p_io
    IF (PRESENT(p_io_c)) p_io_up = p_io_c

#ifdef STANDALONE
    CALL jsb_gather_gp3(global, local(:,:,:), source)
    RETURN
#endif

    IF (require_init) THEN
      CALL init_tr_gather()
    ENDIF
#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_start(timer_gp3)
#endif
    !
    src = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    !
    nlon = ldc%nlon
    lreg = ldc%lreg
    nk   = SIZE(local,2)

    !
    IF (lreg) THEN
      sendbuf => local
    ELSE
      ALLOCATE(sendbuf(ldc%nglon,nk,ldc%nglat))
      CALL reorder(sendbuf,local)
    ENDIF

    !
    IF (p_pe /= p_io_up) THEN
      CALL p_send(sendbuf(:,:,:), p_io_up, tag_gather_gp)
    ELSE
      DO i = 1, p_nprocs
        pe    = gdc(i)%pe
        ALLOCATE(recvbuf(gdc(i)%nglon,nk,gdc(i)%nglat))
        IF (pe /= p_pe) THEN
          CALL p_recv(recvbuf(:,:,:), pe, tag_gather_gp)
        ELSE
          recvbuf(:,:,:) = sendbuf(:,:,:)
        ENDIF
        IF (src == -1 .OR. (src == 0 .AND. pe == p_io_up) &
             .OR. (src == 1 .AND. pe /= p_io_up)) THEN
#ifdef _USE_GATHER_TIMER
          IF (extended_timings) CALL timer_start(timer_fill_gp3)
#endif
          ! unpack first segment
          global(gdc(i)%glons(1):gdc(i)%glone(1),:,gdc(i)%glats(1):gdc(i)%glate(1)) &
               = recvbuf(:,:,:gdc(i)%nglh(1))
          ! unpack second segment
          IF (gdc(i)%nglh(2) > 0) THEN
            IF (gdc(i)%glone(2) > gdc(i)%glons(2)) THEN
              global(gdc(i)%glons(2):gdc(i)%glone(2),:,gdc(i)%glats(2):gdc(i)%glate(2)) &
                   = recvbuf(:,:,gdc(i)%nglat/2+1:)
            ELSE
              ! unpacking second segment, split in longitudes
              global(gdc(i)%glons(2):nlon,:,gdc(i)%glats(2):gdc(i)%glate(2)) &
                   = recvbuf(:nlon-gdc(i)%glons(2)+1,:,gdc(i)%nglat/2+1:)
              global(:gdc(i)%glone(2),:,gdc(i)%glats(2):gdc(i)%glate(2))     &
                   = recvbuf(gdc(i)%nglon-gdc(i)%glone(2)+1:,:,gdc(i)%nglat/2+1:)
            ENDIF
          ENDIF
#ifdef _USE_GATHER_TIMER
          IF (extended_timings) CALL timer_stop(timer_fill_gp3)
#endif
        ENDIF
        DEALLOCATE(recvbuf)
      ENDDO
    ENDIF

    IF (lreg) THEN
      NULLIFY (sendbuf)
    ELSE
      DEALLOCATE (sendbuf)
    ENDIF

#ifdef _USE_GATHER_TIMER
    IF (timings) CALL timer_stop(timer_gp3)
#endif

  END SUBROUTINE gather_gp3
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sa42(global, global_dims, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:,:,:)
    INTEGER,           INTENT(in)    :: global_dims(:)
    REAL(dp),          INTENT(in)    :: local(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: global2d(:,:)

    INTEGER :: dimsize

    dimsize = global_dims(3) * global_dims(4)

    IF (dimsize == 1) THEN
      IF (p_pe == p_io) global2d => global(:,:,1,1)
      CALL gather_sa2 (global2d, local(:,:,1,1), source)
    ELSE
      CALL gather_sa4 (global, local, source)
    ENDIF

  END SUBROUTINE gather_sa42
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sa4 (global, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:,:,:)
    REAL(dp),          INTENT(in)    :: local(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: recvbuf(:,:,:,:)

    INTEGER :: i, im, src, pe, mp1, llevs, lleve, nlm, nhgl, ke, nk

    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    !
    IF (p_pe /= p_io) THEN
      IF (SIZE(local,1) > 0) THEN
        CALL p_send(local(:,:,:,:), p_io, tag_gather_sa)
      END IF
    ELSE
      DO i = 1, p_nprocs
        llevs  = gdc(i)%llevs
        lleve  = gdc(i)%lleve
        nlm    = gdc(i)%nlm
        nhgl   = gdc(i)%nlat/2
        ke = MIN(lleve,SIZE(global,1))
        nk = ke - llevs + 1
        IF (nk < 1) CYCLE

        pe = gdc(i)%pe
        ALLOCATE (recvbuf(nk,2,nlm,nhgl))
        IF (pe /= p_pe) THEN
          CALL p_recv(recvbuf, pe, tag_gather_sa)
        ELSE
          recvbuf(:,:,:,:) = local(:,:,:,:)
        ENDIF

        IF (src == -1 .OR. (src == 0 .AND. pe == p_io) &
                      .OR. (src == 1 .AND. pe /= p_io)) THEN
          DO im = 1, nlm
            mp1 = gdc(i)%lm(im)+1
            global(llevs:ke,:,mp1,:) = recvbuf(:,:,im,:)
          END DO
        ENDIF
        DEALLOCATE (recvbuf)
      END DO
    END IF

  END SUBROUTINE gather_sa4
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sa2 (global, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:)
    REAL(dp),          INTENT(in)    :: local(:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: recvbuf(:,:)

    INTEGER :: i, indx, pe, src, llevs, lleve, nlnm0, nhgl, ke, nk

    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    !
    IF (p_pe /= p_io) THEN
      indx = -1
      search_index: DO i = 1, SIZE(gdc)
        IF(gdc(i)%pe == p_pe) THEN
          indx = i
          EXIT search_index
        ENDIF
      END DO search_index
      nlnm0  = gdc(indx)%nlnm0
      IF (SIZE(local,1) > 0 .AND. nlnm0 > 0) THEN
        CALL p_send(local(:,:), p_io, tag_gather_sa)
      END IF
    ELSE
      DO i = 1, p_nprocs

        llevs  = gdc(i)%llevs
        lleve  = gdc(i)%lleve
        nlnm0  = gdc(i)%nlnm0
        nhgl   = gdc(i)%nlat/2
        ke     = MIN(lleve,SIZE(global,1))
        nk     = ke-llevs+1
        IF (nk < 1 .OR. nlnm0 < 1) CYCLE

        pe = gdc(i)%pe
        ALLOCATE (recvbuf(nk,nhgl))
        IF (pe /= p_pe) THEN
          CALL p_recv(recvbuf, pe, tag_gather_sa)
        ELSE
          recvbuf(:,:) = local(:,:)
        ENDIF
        IF (src == -1 .OR. (src == 0 .AND. pe == p_io) &
                      .OR. (src == 1 .AND. pe /= p_io)) THEN
          global(llevs:ke,:) = recvbuf(:,:)
        ENDIF
        DEALLOCATE(recvbuf)
      END DO
    END IF

  END SUBROUTINE gather_sa2
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sp4(global, global_dims, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:,:,:)
    INTEGER,           INTENT(in)    :: global_dims(:)
    REAL(dp),          INTENT(in)    :: local(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: global3d(:,:,:)
    REAL(dp), POINTER :: global2d(:,:)

    INTEGER :: dimsize3, dimsize4

    dimsize3 = global_dims(3)
    dimsize4 = global_dims(4)

    IF (dimsize4 /= 1) &
         CALL finish('gather_sp4','dimsize4 /= 1')

    IF (dimsize3 == 1) THEN
      NULLIFY (global2d)
      IF (p_pe == p_io) global2d => global(:,:,1,1)
      CALL gather_sp0 (global2d, local(:,:,1,1), source)
    ELSE
      NULLIFY (global3d)
      IF (p_pe == p_io) global3d => global(:,:,:,1)
      CALL gather_sp3 (global3d, local(:,:,:,1), source)
    ENDIF

  END SUBROUTINE gather_sp4
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sp3 (global, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:,:)
    REAL(dp),          INTENT(in)    :: local(:,:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: recvbuf(:,:,:)

    INTEGER :: i, im, src, pe, nlev, mp1, mp, np, snsp, mpgl

    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    !
    IF (p_pe /= p_io) THEN
      IF (SIZE(local) > 0) THEN
        CALL p_send (local(:,:,:), p_io, tag_gather_sp)
      ENDIF
    ELSE
      DO i = 1, p_nprocs
        snsp = gdc(i)%snsp
        IF (snsp < 1) CYCLE
        nlev = SIZE(global,1)
        pe = gdc(i)%pe
        ALLOCATE(recvbuf(nlev,2,snsp))
        IF (pe /= p_pe) THEN
          CALL p_recv(recvbuf(:,:,:), pe, tag_gather_sp)
        ELSE
          recvbuf(:,:,:) = local(:,:,:)
        ENDIF
        IF (src == -1 .OR. (src == 0 .AND. pe == p_io) &
                     .OR. (src == 1 .AND. pe /= p_io)) THEN
          ! unpack
          mp = 0
          DO im = 1, gdc(i)%nsm
            mp1  = gdc(i)%sm(im)+1
            np   = gdc(i)%snnp(im)
            mpgl = gdc(i)%nmp(mp1)+gdc(i)%snn0(im)
            global(:,:,mpgl+1:mpgl+np) = &
                 recvbuf(:,:,mp+1:mp+np)
            mp = mp+np
          END DO
          IF (mp /= snsp) THEN
            WRITE(message_text,*) 'gather_sp: PE', pe, ',mp/=snsp:', mp, snsp
            CALL message('', message_text)
            CALL finish('gather_sp','mp/=snsp')
          ENDIF
        ENDIF
        DEALLOCATE (recvbuf)
      END DO
    END IF
  END SUBROUTINE gather_sp3
  !-----------------------------------------------------------------------------
  SUBROUTINE gather_sp0 (global, local, source, lall_write)
    REAL(dp), POINTER                :: global(:,:)
    REAL(dp),          INTENT(in)    :: local(:,:)
    INTEGER, OPTIONAL, INTENT(in)    :: source
    LOGICAL, OPTIONAL, INTENT(in)    :: lall_write

    REAL(dp), POINTER :: recvbuf (:,:)

    INTEGER :: i, src, pe, nsnm0, snn0, nlev

    src   = debug_parallel
    IF (debug_parallel >= 0 .AND. PRESENT(source)) src = source

    !
    IF (p_pe /= p_io) THEN
      IF (SIZE(local) > 0) THEN
        CALL p_send (local(:,:), p_io, tag_gather_sp)
      ENDIF
    ELSE
      DO i = 1, p_nprocs
        nsnm0 = gdc(i)%nsnm0
        IF (nsnm0 == 0) CYCLE
        snn0  = gdc(i)%snn0(1)
        nlev  = SIZE(global,1)
        pe = gdc(i)%pe
        ALLOCATE (recvbuf (nlev, nsnm0))
        IF (pe /= p_pe) THEN
          CALL p_recv( recvbuf, pe, tag_gather_sp)
        ELSE
          recvbuf(:,:) = local(:,:)
        ENDIF
        IF (src ==-1 .OR. (src == 0 .AND. pe == p_io) &
                     .OR. (src == 1 .AND. pe /= p_io)) THEN
          ! unpack
          global(:,1+snn0:nsnm0+snn0) = recvbuf(:,:)
        ENDIF
        DEALLOCATE (recvbuf)
      END DO
    END IF

  END SUBROUTINE gather_sp0

    SUBROUTINE check_my_comm(p_io_up, group_comm, sub_group_rank)
      ! SUBROUTINE to calculate rank of p_io_up is in group_comm subcommunicator
      ! used to check that subcommunicator is correctly constructed

      INTEGER, INTENT(in)  :: p_io_up
      INTEGER, INTENT(in)  :: group_comm
      INTEGER, INTENT(out) :: sub_group_rank(1)

      INTEGER :: ierr, all_group, sub_group

      ! make group from p_all_comm
      CALL MPI_COMM_GROUP(p_all_comm, all_group, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::check_my_comm: bad case (1)')
      ! make group from sub communicator
      CALL MPI_COMM_GROUP(group_comm, sub_group, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::check_my_comm: bad case (2)')
      ! compare rank of p_io in p_all_comm with rank in slab_comm
      CALL MPI_GROUP_TRANSLATE_RANKS(all_group, 1, (/ p_io_up /), sub_group,  &
           sub_group_rank, ierr)
      IF (ierr /= MPI_SUCCESS) &
           CALL finish('mo_tr_gather::check_my_comm: bad case (3)')
      IF (sub_group_rank(1) == MPI_UNDEFINED) &
           CALL finish('mo_tr_gather::check_my_comm: MPI_UNDEFINED')
    END SUBROUTINE check_my_comm

! op_pj_20110204+
!==============================================================================
  FUNCTION gen_inv_map (gl_dc)
    USE mo_decomposition, ONLY: pe_decomposed

    TYPE (pe_decomposed), INTENT(in) :: gl_dc(:)

    INTEGER :: gen_inv_map(0:SIZE(gl_dc)-1)

    INTEGER :: i

    gen_inv_map = -1

    DO i = 1, SIZE(gl_dc)
      gen_inv_map(gl_dc(i)% pe) = i
    END DO

    IF(ANY(gen_inv_map==-1))THEN
      WRITE (message_text,*) 'index not found in decomposition table'
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  required:',i
      CALL message ('', TRIM(message_text))
      WRITE (message_text,*) '  found  :',gl_dc% pe
      CALL message ('', TRIM(message_text))
      CALL finish('mo_transpose:gen_inv_map', &
           'index not found in decomposition table')
    END IF
  END FUNCTION gen_inv_map
!==============================================================================
! op_pj_20110204-

END MODULE mo_tr_gather
