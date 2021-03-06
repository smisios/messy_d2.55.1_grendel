! ##########################################################################
#ifndef ICON
! ##########################################################################

MODULE mo_mpi

  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve proper output.

#ifdef _OPENMP
#if (! defined __PGI)
  USE omp_lib,   ONLY: OMP_GET_MAX_THREADS, OMP_SET_NUM_THREADS, &
                       OMP_SET_DYNAMIC
#endif
#ifdef NAG
  USE f90_unix_env, ONLY: getenv
#endif
#endif

#ifdef MESSYMMD
  USE MMD_handle_communicator, ONLY: MMD_get_model_communicator &
       , MMD_Print_Error_Message, MMD_STATUS_OK, MMD_FreeMem_communicator
#endif
#ifndef MESSY
  USE mo_kind
  USE mo_doctor, ONLY: nerr
! mz_ab_20100222+
#else
  USE messy_main_constants_mem, ONLY: dp, sp, i4, i8, nerr !, iouerr
#endif
! mz_ab_20100222-

#ifndef NOMPI
#ifndef __PGI
  USE mpi
#endif
#endif

  IMPLICIT NONE

  PRIVATE                          ! all declarations are private

#ifndef NOMPI
#ifdef __PGI
  INCLUDE  'mpif.h'
#endif
#endif

#if (defined _OPENMP) && (defined __PGI)
  INTEGER, EXTERNAL :: OMP_GET_MAX_THREADS
  EXTERNAL :: OMP_SET_NUM_THREADS
  EXTERNAL :: OMP_SET_DYNAMIC
#endif

  ! subroutines defined, overloaded depending on argument type
  PUBLIC :: ask_omp_threads
  PUBLIC :: p_start, p_stop, p_abort
  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier
  PUBLIC :: p_isend, p_irecv, p_wait, p_wait_any
  PUBLIC :: p_gather, p_max, p_min, p_sum, p_global_sum, p_field_sum
  PUBLIC :: nerr
#ifdef ONESIDED
  PUBLIC :: p_win_create, p_win_fence, p_win_free, p_put
#endif
  PUBLIC :: p_set_communicator
  PUBLIC :: p_probe
  PUBLIC :: p_allgather
#ifdef MESSY
  PUBLIC :: p_time
  PUBLIC :: p_lor
#ifndef NOMPI
  PUBLIC :: p_error, p_real
#endif
  PUBLIC :: npes
#if defined (MPIOM_2000)
  INTEGER, PUBLIC, PARAMETER :: ihalo_max=10
  PUBLIC :: init_MPI_datatypes
  INTEGER, public :: ns_boundary(ihalo_max), ew_boundary(ihalo_max)
  INTEGER, public :: ns_boundary_ke(ihalo_max), ew_boundary_ke(ihalo_max)
  INTEGER, public :: ns_boundary_kep(ihalo_max), ew_boundary_kep(ihalo_max)
  INTEGER :: mmm_mpi_datatype_dp, op_mmm_mpi_reduce_dp
  TYPE min_mean_max_dp
    SEQUENCE
    REAL(dp) :: min, mean, max
    INTEGER :: pe_min, rcount, pe_max
  END TYPE min_mean_max_dp
  PUBLIC :: generate_mpi_struct_type, create_mpi_op
#endif
#endif

#ifndef NOMPI
  PUBLIC :: MPI_INTEGER, MPI_STATUS_SIZE, MPI_SUCCESS, MPI_ANY_SOURCE
!#ifdef HAVE_LIBPNETCDF ! mz_pj_20061107
#if defined(HAVE_LIBPNETCDF) || defined (MESSY)
  PUBLIC :: MPI_INFO_NULL, MPI_OFFSET_KIND
#endif
  PUBLIC :: MPI_SUM, MPI_MIN, MPI_MAX, MPI_LOR
#endif

  ! real data type matching real type of MPI implementation

  PUBLIC :: p_real_dp

  ! logical switches

  PUBLIC :: p_parallel, p_parallel_io

  ! PE identifier

  PUBLIC :: p_pe, p_io, p_nprocs

  ! communicator

  PUBLIC :: p_all_comm
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d

  ! old fashioned method (MPI-1)

!!$#ifndef NOMPI
!!$  INCLUDE 'mpif.h'
!!$#endif

  ! general run time information

#ifndef NOMPI
  INTEGER :: version, subversion   ! MPI version
#endif

  ! MPI call inherent variables

#ifndef NOMPI
  INTEGER :: p_error                     ! MPI error number

  INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV
#endif

  ! public parallel run information

  LOGICAL, SAVE :: p_parallel, p_parallel_io ! mz_pj_20050610 SAVE
  ! um_ak_20120601+ SAVE added
  INTEGER, SAVE  :: p_pe     ! this is the PE number of this task
  INTEGER, SAVE  :: p_io     ! PE number of PE handling IO
  INTEGER, SAVE  :: p_nprocs ! number of available PEs (processors)
  ! um_ak_20120601-

  ! communicator sets

  INTEGER :: p_all_comm       ! replaces MPI_COMM_WORLD in one application
  INTEGER :: p_communicator_a ! for Set A
  INTEGER :: p_communicator_b ! for Set B
  INTEGER :: p_communicator_d ! for debug node

  ! non blocking calls

#ifndef NOMPI
#ifndef _SX
  INTEGER, ALLOCATABLE, SAVE :: p_request(:) ! request values for non blocking calls
#else
  INTEGER, POINTER    , SAVE :: p_request(:) ! request values for non blocking calls
#endif
  INTEGER :: p_irequest ! the first p_irequest values of p_request are in use
  INTEGER :: p_mrequest ! actual size of p_request
  INTEGER, PARAMETER :: p_request_alloc_size = 4096
#endif

  ! module intrinsic names

  INTEGER :: mype                  ! this is the PE number of this task
#ifndef NOMPI
  INTEGER :: iope                  ! PE able to do IO
#endif
  INTEGER :: npes                  ! number of available PEs

  INTEGER :: nbcast                ! counter for broadcasts for debugging

  ! MPI transfer types

  INTEGER :: p_real_dp
#ifndef NOMPI
  INTEGER :: p_real_sp
  INTEGER :: p_int_i4
  INTEGER :: p_int_i8

  ! native types

  INTEGER :: p_int       ! maybe switched by compiler options therefor reset
  INTEGER :: p_real      ! maybe switched by compiler options therefor reset

  INTEGER :: p_bool
  INTEGER :: p_char

  INTEGER :: p_ig, p_rg
  INTEGER :: p_i4, p_i8
  INTEGER :: p_sp, p_dp
#endif

  ! for checking out integer and real variables separat. KIND values usually
  ! overlap and give the byte size or are defined as sequence separate for
  ! both groups.

  INTEGER, PARAMETER :: real_type    = 1
  INTEGER, PARAMETER :: integer_type = 2

  ! define generic interfaces to allow proper compiling
  ! with picky compilers like NAG f95 for clean argument checking and
  ! shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real_5d
#if defined(MESSY) && defined(MPIOM_2000)
! mz_ap_20100922
     MODULE PROCEDURE p_isend_mpitype_2d
#endif
  END INTERFACE

#ifdef ONESIDED
  INTERFACE p_win_create
     MODULE PROCEDURE p_win_create_real_2d
     MODULE PROCEDURE p_win_create_real_3d
  END INTERFACE

  INTERFACE p_win_fence
     MODULE PROCEDURE p_win_fence_int
     MODULE PROCEDURE p_win_fence_int_1d
  END INTERFACE

  INTERFACE p_win_free
     MODULE PROCEDURE p_win_free_int
     MODULE PROCEDURE p_win_free_int_1d
  END INTERFACE

  INTERFACE p_put
     MODULE PROCEDURE p_put_real_2d
     MODULE PROCEDURE p_put_real_3d
  END INTERFACE
#endif

  INTERFACE p_recv
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_real_4d
#if defined(MESSY) && defined(MPIOM_2000)
     MODULE PROCEDURE p_irecv_mpitype_2d
#endif
     ! ju_ch_20110429+
     MODULE PROCEDURE p_irecv_int_2d
     MODULE PROCEDURE p_irecv_int_3d
     ! ju_ch_20110429-
  END INTERFACE

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
#if defined(MESSY) && defined(MPIOM_2000)
     MODULE PROCEDURE p_sendrecv_mpitype_2d
#endif
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
!!$     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_char_1d
     ! um_ak_20130702+
     MODULE PROCEDURE p_bcast_real_sp_0d
     MODULE PROCEDURE p_bcast_real_sp_1d
     MODULE PROCEDURE p_bcast_real_sp_2d
     MODULE PROCEDURE p_bcast_real_sp_3d
     MODULE PROCEDURE p_bcast_real_sp_4d
     MODULE PROCEDURE p_bcast_int_i4_1d
     MODULE PROCEDURE p_bcast_int_i8_1d
     ! um_ak_20130702-
  END INTERFACE

  INTERFACE p_probe
     MODULE PROCEDURE p_probe_real
     MODULE PROCEDURE p_probe_int
     MODULE PROCEDURE p_probe_bool
     MODULE PROCEDURE p_probe_char
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_1d2d
#if defined(MESSY) && defined(MPIOM_2000)
     MODULE PROCEDURE p_gather_real_2d2d
#endif
  END INTERFACE

  ! mz_bk_20120725+
  INTERFACE p_allgather
     MODULE PROCEDURE p_allgather_real_2d2d
     MODULE PROCEDURE p_allgather_real_3d2d ! op_pj_20160811
  END INTERFACE
  ! mz_bk_20120725-

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_0d
     MODULE PROCEDURE p_sum_0d_int     ! ju_ch_20110429
     MODULE PROCEDURE p_sum_1d
#ifdef MESSY
     MODULE PROCEDURE p_sum_2d
     MODULE PROCEDURE p_sum_2d_int     ! ju_nt_20160122
     MODULE PROCEDURE p_sum_3d
#endif
  END INTERFACE

  INTERFACE p_lor
     MODULE PROCEDURE p_lor_0d
  END INTERFACE p_lor

  INTERFACE p_global_sum
     MODULE PROCEDURE p_global_sum_1d
  END INTERFACE

  INTERFACE p_field_sum
     MODULE PROCEDURE p_field_sum_1d
  END INTERFACE

CONTAINS

  SUBROUTINE p_start(model_name)

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    USE mod_prism_proto, ONLY: prism_ok
#endif
#endif

#ifndef MESSY
    USE mo_util_string, ONLY: toupper
! mz_ab_20100222+
!!$#else
!!$    USE messy_main_tools, ONLY: toupper=>ucase
#endif
! mz_ab_20100222-

    CHARACTER(len=*), INTENT(in), OPTIONAL :: model_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds

    INTEGER :: io_size, integer_io_size, integer_byte_size

    INTEGER      :: iig = 0
    INTEGER (i4) :: ii4 = 0_i4
    INTEGER (i8) :: ii8 = 0_i8

    REAL         :: rrg = 0.0
    REAL (sp)    :: rsp = 0.0_sp
    REAL (dp)    :: rdp = 0.0_dp

    ! temporary array to distibute the determined MPI types

    INTEGER :: p_send(7)

    ! variables used for determing the I/O PE

    LOGICAL :: liope
    INTEGER, ALLOCATABLE :: iope_table(:)
    CHARACTER(len=132) :: io_pe_message

    CHARACTER(len=132) :: yname

    ! variables used for determing the OpenMP threads
    ! suitable as well for coupled models

#if (defined _OPENMP) && (! defined NAG)
    CHARACTER(len=32) :: env_name
    CHARACTER(len=32) :: thread_num
    INTEGER :: env_threads, threads
#if defined(__GFORTRAN__)
    INTRINSIC :: getenv
#else
#ifndef LF
    EXTERNAL :: getenv
#endif
#endif
#endif

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)
    INTEGER :: prism_model_number
    CHARACTER(len=132) :: prism_model_name

    EXTERNAL :: prism_abort_proto
    EXTERNAL :: prism_init_comp_proto
    EXTERNAL :: prism_get_localcomm_proto
#endif
#endif

    ! loop index

    INTEGER :: jp

    ! Executable statements:

    IF (PRESENT(model_name)) THEN
      yname = TRIM(model_name)
    ELSE
#if defined(MPIOM_2000) || defined(MPIOM_13B)
      yname = 'mpi-om'
#else
      yname = 'echam5'
#endif
    END IF

    nbcast = 0

    io_pe_message(:) = ' '

    ! start MPI

#ifndef NOMPI
    CALL MPI_INIT (p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_INIT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#endif

    ! create communicator for this process alone before
    ! potentially joining MPI2

#ifndef NOMPI
#if defined (__prism) && defined (use_comm_MPI1)

    prism_model_name = TRIM(yname)

    CALL prism_init_comp_proto (prism_model_number, TRIM(prism_model_name), &
         p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) ' prism_init_comp_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort1')
    ENDIF

    CALL prism_get_localcomm_proto(p_all_comm, p_error)

    IF (p_error /= prism_ok) THEN
      WRITE (nerr,*) ' prism_get_localcomm_proto failed'
      CALL prism_abort_proto(prism_model_number, TRIM(yname),'abort2')
    ENDIF

#else

#ifdef MESSYMMD
    CALL MMD_get_model_communicator (p_all_comm, p_error)
     if(p_error /= MMD_STATUS_OK)   then
        WRITE (nerr,'(a)') ' MMD_get_model_communicator failed.'
        call MMD_Print_Error_Message ( nerr, p_error)
        CALL p_abort
     end if
#else
    CALL MPI_COMM_DUP (MPI_COMM_WORLD, p_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_DUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF
#endif
#endif
#endif

    ! get local PE identification

#ifndef NOMPI
    CALL MPI_COMM_RANK (p_all_comm, mype, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    ELSE
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a)') ' PE ', mype, ' started.'
#endif
    END IF
#else
    mype = 0
#endif

    ! get number of available PEs

#ifndef NOMPI
    CALL MPI_COMM_SIZE (p_all_comm, npes, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF
#else
    npes = 1
#endif

    ! for non blocking calls

#ifndef NOMPI
    p_mrequest = p_request_alloc_size
    ALLOCATE(p_request(p_mrequest))
    p_irequest = 0
#endif

    ! look for a dedicated IO PE

#ifndef NOMPI
! mz_kk_20081107+
!kk If all MPI processes can do IO, assign pe 0 as p_io, else use
!   SPECIAL_MPI_IO
#ifdef SPECIAL_MPI_IO
! mz_kk_20081107-
#ifndef OPENMPI
    CALL MPI_ATTR_GET (p_all_comm, MPI_IO, iope, liope, p_error)
#else
! op_pj_20111129+
!!$    CALL MPI_ATTR_GET (MPI_COMM_WORLD, MPI_IO, iope, liope, p_error)
    CALL MPI_COMM_GET_ATTR(p_all_comm, MPI_IO, iope, liope, p_error)
! op_pj_20111129-
#ifdef MESSYMMD
    if(p_error /= MMD_STATUS_OK)   then
       WRITE (nerr,'(a)') ' MMD does not work in case of OPENMPI'
       CALL p_abort
    end if
#endif
#endif

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_ATTR_GET failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (iope == MPI_ANY_SOURCE) THEN

       ! all nodes can do IO

       IF (mype == 0) THEN
          WRITE (io_pe_message,'(a)') &
               '  All nodes can do I/O, select PE 0 for I/O handling.'
       END IF
       p_io = 0

    ELSE

       ALLOCATE (iope_table(0:npes-1))

       IF (liope) THEN
          iope_table(mype) = iope
          CALL MPI_GATHER (iope_table(mype), 1, MPI_INTEGER, &
               iope_table, 1, MPI_INTEGER,       &
               0, p_all_comm, p_error)

          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GATHER failed.'
             WRITE (nerr,'(a,i4)') ' Error =  ', p_error
             CALL p_abort
          END IF

          IF (mype == 0) THEN
             ! Now select the first given iope from table as IO PE.
             WRITE (io_pe_message,'(a,i3,a)') &
                  '  Select PE ', iope_table(0), ' for I/O handling.'
             p_io = iope_table(0)
          END IF
          CALL MPI_BCAST (p_io, 1, MPI_INTEGER, 0, p_all_comm, p_error)

       ELSE
          ! if no dedicated IO PE is given, use PE 0
          p_io = 0

          WRITE (io_pe_message,'(a)') &
               '  No dedicated I/O PE, select PE 0 for I/O handling.'
       END IF

       DEALLOCATE (iope_table)

    END IF
#else
    p_io = 0
#endif
! mz_kk_20081107+
#else
    p_io = 0
#endif
! mz_kk_20081107-

    ! Information ...

    IF (mype == 0) THEN
       WRITE (nerr,'(/,a,a,a)') ' ', &
            TRIM(yname), ' MPI interface runtime information:'
    END IF

     IF (npes < 2) THEN
       p_parallel = .FALSE.
       p_parallel_io = .TRUE.   ! can always do I/O
       IF (mype == 0) THEN
          WRITE (nerr,'(a)') '  Single processor run.'
       END IF
       p_pe = 0
       p_nprocs = 1
    ELSE
       p_parallel = .TRUE.
       IF (mype == p_io) THEN
          p_parallel_io = .TRUE.
       ELSE
          p_parallel_io = .FALSE.
       END IF
       IF (mype == 0) THEN
          WRITE (nerr,'(a,i4,a)') '  Run on ', npes, ' processors.'
       END IF
       p_pe = mype
       p_nprocs = npes
    END IF

#ifdef _OPENMP

    ! Expect that PE 0 did got the information of OMP_NUM_THREADS.
    ! That might be wrong in the coupled case when the model is
    ! started via MPI dynamic process creation. So we have to check
    ! the environment variable too.

    IF (mype == 0) THEN

#ifndef MESSY
       env_name = toupper(TRIM(yname)) // '_THREADS'
#else
! op_pj_20160802+
!!$    env_name = TRIM(yname)
!!$    CALL toupper(env_name)
!!$    env_name = TRIM(env_name)//'_THREADS'
       env_name = 'OMP_NUM_THREADS'
! op_pj_20160802-
#endif

#ifdef NAG
      CALL getenv(TRIM(env_name), thread_num, errno=p_error)
      IF (p_error /= 0) THEN
#else
      CALL getenv(TRIM(env_name), thread_num)
      IF (thread_num /= ' ') THEN
#endif
        READ(thread_num,*) env_threads
      ELSE
        WRITE (nerr,'(a)') ' Number of OpenMP threads unknown!'
        WRITE (nerr,'(a,a,a,/,a,a,a)') &
             ' Environment variable ', TRIM(env_name), ' either not set,', &
             ' or not available to ', TRIM(yname), ' root PE.'
        CALL p_abort
      ENDIF
      threads = env_threads
   ENDIF

#ifndef NOMPI
    ! Make number of threads from environment available to all model PEs

    CALL MPI_BCAST (threads, 1, MPI_INTEGER, 0, p_all_comm, p_error)
#endif

    ! Inform on OpenMP thread usage


    CALL OMP_SET_DYNAMIC(.FALSE.)

    CALL OMP_SET_NUM_THREADS(threads)

    threads = OMP_GET_MAX_THREADS()

    IF (mype == 0) THEN

       ! write out thread usage of master PE

       WRITE (nerr,'(a,i4,a,i4,a)') &
            '  PE ', mype, ': using ', &
            threads, &
            ' OpenMP threads.'

#ifndef NOMPI
       ! recevie and write the remaining MPI PEs thread number

       DO jp = 1, npes-1
          CALL MPI_RECV(threads, 1, MPI_INTEGER, jp, jp, &
               p_all_comm, p_status, p_error)
          WRITE (nerr,'(a,i4,a,i4,a)') &
               '  PE ', jp, ': using ', &
               threads, &
               ' OpenMP threads.'
       ENDDO
#endif
    ELSE
#ifndef NOMPI
       ! send threads per PE to master PE for write out

       CALL MPI_SEND (threads, 1, MPI_INTEGER, 0, mype, &
               p_all_comm, p_error)
#endif
    ENDIF


#endif

    ! inform on I/O PE situation

    IF (mype == 0) THEN
       WRITE (nerr, '(a)') io_pe_message(1:LEN_TRIM(io_pe_message))
    END IF

#ifndef NOMPI

    IF (p_parallel) THEN

       ! lets check the available MPI version

       CALL MPI_GET_VERSION (version, subversion, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_GET_VERSION failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) THEN
          WRITE (nerr,'(a,i1,a1,i1)') &
               '  Used MPI version: ', version, '.', subversion
       END IF

       ! due to a possible circular dependency with mo_machine and other
       ! modules, we determine here locally the I/O size of the different
       ! kind types (assume 8 bit/byte. This is than used for determing
       ! the right MPI send/receive type parameters.

       ! first get the native INTEGER size

       integer_byte_size = BIT_SIZE(iig)/8

       ! and inquire for the I/O size (is independent of byte or word
       ! values ...)

       INQUIRE (iolength=io_size) iig
       integer_io_size = io_size
       p_ig = io_size/integer_io_size*integer_byte_size

       ! and the native REAL size

       INQUIRE (iolength=io_size) rrg
       p_rg = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte INTEGER
       ! (might be 8 byte both, or only the 4 byte available ...

       INQUIRE (iolength=io_size) ii4
       p_i4 = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) ii8
       p_i8 = io_size/integer_io_size*integer_byte_size

       ! find now the size of usual 4 byte and 8 byte REAL
       ! (might be 8 byte both)

       INQUIRE (iolength=io_size) rsp
       p_sp = io_size/integer_io_size*integer_byte_size
       INQUIRE (iolength=io_size) rdp
       p_dp = io_size/integer_io_size*integer_byte_size

       ! testing this variables

       p_int_i4  = p_type (i4, integer_type)
       p_int_i8  = p_type (i8, integer_type)
       p_real_sp = p_type (sp, real_type)
       p_real_dp = p_type (dp, real_type)

       IF (mype == 0) THEN

          IF (p_ig == p_i4) THEN
             p_int = p_int_i4
          ELSE IF (p_ig == p_i8) THEN
             p_int = p_int_i8
          END IF

          IF (p_rg == p_sp) THEN
             p_real = p_real_sp
          ELSE IF (p_rg == p_dp) THEN
             p_real = p_real_dp
          END IF

       END IF

       p_send(1) = p_int
       p_send(2) = p_int_i4
       p_send(3) = p_int_i8
       p_send(4) = p_real
       p_send(5) = p_real_sp
       p_send(6) = p_real_dp

       CALL MPI_BCAST (p_send, 6, MPI_INTEGER, 0, p_all_comm, p_error)

#ifdef DEBUG
       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a)') ' MPI_BCAST for send/receive types failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
       END IF
#endif

       p_int      = p_send(1)
       p_int_i4   = p_send(2)
       p_int_i8   = p_send(3)
       p_real     = p_send(4)
       p_real_sp  = p_send(5)
       p_real_dp  = p_send(6)

       ! set logical and character types to native types

       p_bool = MPI_LOGICAL
       p_char = MPI_CHARACTER

       IF (mype == 0) THEN
          WRITE (nerr,'(/)')
          IF (p_real_sp == MPI_REAL) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_REAL for KIND sp'
          ELSE IF (p_real_sp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND sp'
          END IF

          IF (p_real_dp == MPI_DOUBLE_PRECISION) THEN
             WRITE (nerr,'(a)') &
                  ' Selected type: MPI_DOUBLE_PRECISION for KIND dp'
          END IF

          IF (p_int_i4 == MPI_INTEGER4) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER4 for KIND i4'
          ELSE IF (p_int_i4 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i4'
          END IF

          IF (p_int_i8 == MPI_INTEGER8) THEN
             WRITE (nerr,'(a)') ' Selected type: MPI_INTEGER8 for KIND i8'
          END IF
          WRITE (nerr,'(/)')
       END IF
    END IF

#endif



#ifdef DEBUG
    WRITE (nerr,'(a)')    ' Transfer types:'
    WRITE (nerr,'(a,i4)') '  INTEGER generic:', p_int
    WRITE (nerr,'(a,i4)') '  INTEGER 4 byte :', p_int_i4
    WRITE (nerr,'(a,i4)') '  INTEGER 8 byte :', p_int_i8
    WRITE (nerr,'(a,i4)') '  REAL generic   :', p_real
    WRITE (nerr,'(a,i4)') '  REAL single    :', p_real_sp
    WRITE (nerr,'(a,i4)') '  REAL double    :', p_real_dp
#endif
#ifdef MESSY
    CALL messy_setup(0)
#endif

  END SUBROUTINE p_start

  SUBROUTINE p_stop

    ! finish MPI and clean up all PEs

#ifndef NOMPI
    ! to prevent abort due to unfinished communication
    CALL p_barrier(p_all_comm)

#if defined(MESSY) && defined(ECHAM5)
    CALL messy_finalize ! ub_ak_20181030
#endif
    CALL MPI_FINALIZE (p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
    p_parallel = .FALSE.
    DEALLOCATE(p_request)
#ifdef MESSYMMD
    CALL MMD_FreeMem_Communicator
#endif
#endif

  END SUBROUTINE p_stop

  SUBROUTINE p_abort

    ! this routine should be used instead of abort, util_abort() or STOP
    ! in all routines for proper clean up of all PEs

#if defined(ECHAM5)
    EXTERNAL util_exit
#endif

#ifndef NOMPI
#ifndef MESSY
    CALL MPI_ABORT (MPI_COMM_WORLD, 0, p_error)
#else
    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)  ! mz_pj_20030313
#endif

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF
#else
#if defined(ECHAM5)
    CALL util_exit(1)
#else
#ifndef NOMPI
    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)
#else
    STOP
#endif
#endif
! mz_ab_20100504-
#endif

  END SUBROUTINE p_abort

  FUNCTION p_type (kind_type, var_type) RESULT (p_message_type)

!!$    USE mo_kind  ! op_pj_20100503: globally used

    INTEGER              :: p_message_type
    INTEGER, INTENT(in)  :: kind_type, var_type

    IF (var_type == integer_type) THEN
       IF (kind_type == i8) THEN
          p_message_type = check_type_i8 ()
       ELSE IF (kind_type == i4) THEN
          p_message_type = check_type_i4 ()
       END IF
    ELSE IF (var_type == real_type) THEN
       IF (kind_type == dp) THEN
          p_message_type = check_type_dp ()
       ELSE IF (kind_type == sp) THEN
          p_message_type = check_type_sp ()
       END IF
    END IF

  END FUNCTION p_type

  FUNCTION check_type_i4 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    INTEGER (i4) :: buf_int_i4(1,8)
    INTEGER (i4) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_int_i4(1,1) = a

       CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER4, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER4, 1, 2, &
            p_all_comm, p_status, p_error)
       b = buf_int_i4(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)
          CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 3, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error)
          c =  buf_int_i4(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER4 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER4 anyway'
          END IF

          p_type = MPI_INTEGER4

       ELSE
          CALL MPI_SEND (MPI_INTEGER8, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)
          CALL MPI_SEND (buf_int_i4(1,1), 1, MPI_INTEGER8, 1, 5, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_int_i4(1,5), 1, MPI_INTEGER8, 1, 6, &
               p_all_comm, p_status, p_error)
          b = buf_int_i4(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_int_i4(1,1), p_i4, MPI_BYTE, 1, 7, &
                  p_all_comm, p_error)
             CALL MPI_RECV (buf_int_i4(1,5), p_i4, MPI_BYTE, 1, 8, &
                  p_all_comm, p_status, p_error)
             c =  buf_int_i4(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a)') &
                     ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_INTEGER8 anyway'
             END IF

             p_type = MPI_INTEGER8

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER4, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_int_i4(1,3) = buf_int_i4(1,7)
       CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER4, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            p_all_comm, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 3, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 4, &
               p_all_comm, p_error)
       ELSE IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), 1, MPI_INTEGER8, 0, 5, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), 1, MPI_INTEGER8, 0, 6, &
               p_all_comm, p_error)
       END IF

       IF (datatype == MPI_INTEGER8) THEN
          CALL MPI_RECV (buf_int_i4(1,7), p_i4, MPI_BYTE, 0, 7, &
               p_all_comm, p_status, p_error)
          buf_int_i4(1,3) = buf_int_i4(1,7)
          CALL MPI_SEND (buf_int_i4(1,3), p_i4, MPI_BYTE, 0, 8, &
               p_all_comm, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i4

  FUNCTION check_type_i8 () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER (i8) :: buf_int_i8(1,8)
    INTEGER (i8) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_int_i8(1,1) = a

       CALL MPI_SEND (buf_int_i8(1,1), 1, MPI_INTEGER8, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_int_i8(1,5), 1, MPI_INTEGER8, 1, 2, &
            p_all_comm, p_status, p_error)
       b = buf_int_i8(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_int_i8(1,1), p_i8, MPI_BYTE, 1, 3, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_int_i8(1,5), p_i8, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error)
          c =  buf_int_i8(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_INTEGER8 and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_INTEGER8 anyway'
          END IF

          p_type = MPI_INTEGER8

       ELSE
          WRITE (nerr,'(a)') ' MPI_INTEGER8 not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_int_i8(1,7), 1, MPI_INTEGER8, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), 1, MPI_INTEGER8, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (buf_int_i8(1,7), p_i8, MPI_BYTE, 0, 3, &
            p_all_comm, p_status, p_error)
       buf_int_i8(1,3) = buf_int_i8(1,7)
       CALL MPI_SEND (buf_int_i8(1,3), p_i8, MPI_BYTE, 0, 4, &
            p_all_comm, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_i8

  FUNCTION check_type_sp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    INTEGER :: datatype

    REAL (sp) :: buf_real_sp(1,8)
    REAL (sp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_real_sp(1,1) = a

       CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_REAL, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_REAL, 1, 2, &
            p_all_comm, p_status, p_error)
       b = buf_real_sp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (MPI_BYTE, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)
          CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 3, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error)
          c =  buf_real_sp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') ' Warning: MPI_REAL and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_REAL anyway'
          END IF

          p_type = MPI_REAL

       ELSE
          CALL MPI_SEND (MPI_DOUBLE_PRECISION, 1, MPI_INTEGER, 1, 0, &
               p_all_comm, p_error)
          CALL MPI_SEND (buf_real_sp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 5, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_real_sp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 6, &
               p_all_comm, p_status, p_error)
          b = buf_real_sp(1,5)

          IF (a == b) THEN
             CALL MPI_SEND (buf_real_sp(1,1), p_sp, MPI_BYTE, 1, 7, &
                  p_all_comm, p_error)
             CALL MPI_RECV (buf_real_sp(1,5), p_sp, MPI_BYTE, 1, 8, &
                  p_all_comm, p_status, p_error)
             c =  buf_real_sp(1,5)

             IF (a /= c) THEN
                WRITE (nerr,'(a,a)') &
                     ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE ', &
                     'not equivalent'
                WRITE (nerr,'(a)') &
                     ' Using MPI_DOUBLE_PRECISION anyway'
             END IF

             p_type = MPI_DOUBLE_PRECISION

          END IF
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_REAL, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_real_sp(1,3) = buf_real_sp(1,7)
       CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_REAL, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (datatype, 1, MPI_INTEGER, 0, 0, &
            p_all_comm, p_status, p_error)

       IF (datatype == MPI_BYTE) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 3, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 4, &
               p_all_comm, p_error)
       ELSE IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 5, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 6, &
               p_all_comm, p_error)
       END IF

       IF (datatype == MPI_DOUBLE_PRECISION) THEN
          CALL MPI_RECV (buf_real_sp(1,7), p_sp, MPI_BYTE, 0, 7, &
               p_all_comm, p_status, p_error)
          buf_real_sp(1,3) = buf_real_sp(1,7)
          CALL MPI_SEND (buf_real_sp(1,3), p_sp, MPI_BYTE, 0, 8, &
               p_all_comm, p_error)
       END IF
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_sp

  FUNCTION check_type_dp () RESULT (p_type)

    INTEGER :: p_type

#ifndef NOMPI
    REAL (dp) :: buf_real_dp(1,8)
    REAL (dp) :: a, b, c

    p_type = MPI_DATATYPE_NULL

    IF (mype == 0) THEN

       a = HUGE(a)
       buf_real_dp(1,1) = a

       CALL MPI_SEND (buf_real_dp(1,1), 1, MPI_DOUBLE_PRECISION, 1, 1, &
            p_all_comm, p_error)
       CALL MPI_RECV (buf_real_dp(1,5), 1, MPI_DOUBLE_PRECISION, 1, 2, &
            p_all_comm, p_status, p_error)
       b = buf_real_dp(1,5)

       IF (a == b) THEN
          CALL MPI_SEND (buf_real_dp(1,1), p_dp, MPI_BYTE, 1, 3, &
               p_all_comm, p_error)
          CALL MPI_RECV (buf_real_dp(1,5), p_dp, MPI_BYTE, 1, 4, &
               p_all_comm, p_status, p_error)
          c =  buf_real_dp(1,5)

          IF (a /= c) THEN
             WRITE (nerr,'(a)') &
                  ' Warning: MPI_DOUBLE_PRECISION and MPI_BYTE not equivalent'
             WRITE (nerr,'(a)') ' Using MPI_DOUBLE_PRECISION anyway'
          END IF

          p_type = MPI_DOUBLE_PRECISION

       ELSE
          WRITE (nerr,'(a)') ' MPI_DOUBLE_PRECISION not available.'
       END IF
    END IF

    IF (mype == 1) THEN
       CALL MPI_RECV (buf_real_dp(1,7), 1, MPI_DOUBLE_PRECISION, 0, 1, &
            p_all_comm, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), 1, MPI_DOUBLE_PRECISION, 0, 2, &
            p_all_comm, p_error)

       CALL MPI_RECV (buf_real_dp(1,7), p_dp, MPI_BYTE, 0, 3, &
            p_all_comm, p_status, p_error)
       buf_real_dp(1,3) = buf_real_dp(1,7)
       CALL MPI_SEND (buf_real_dp(1,3), p_dp, MPI_BYTE, 0, 4, &
            p_all_comm, p_error)
    END IF
#else
    p_type = 0
#endif

  END FUNCTION check_type_dp

  ! communicator set up

  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

#ifndef NOMPI
    INTEGER :: all_debug_pes(SIZE(mapmesh))

    INTEGER :: group_world, group_a, group_b, group_d
    INTEGER :: p_communicator_tmp

    INTEGER :: n, members

    INTEGER :: ranks(1) = 0

    ! first set global group

    CALL MPI_COMM_GROUP (p_all_comm, group_world, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_GROUP failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    ! communicator is p_all_comm

    IF (debug_parallel >= 0 ) THEN

       CALL MPI_GROUP_INCL (group_world, 1, ranks, group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype == 0) p_communicator_d = p_communicator_tmp

       DO n = 1, SIZE(mapmesh)
          all_debug_pes(n) = n
       END DO

       CALL MPI_GROUP_INCL (group_world, SIZE(mapmesh), all_debug_pes, &
            group_d, p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_d, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       IF (mype /= 0) p_communicator_d = p_communicator_tmp

    ELSE
       p_communicator_d = p_all_comm
    END IF

    DO n = 0, nproca-1
       members = nprocb
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(:,n), group_a, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_a, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_a = p_communicator_tmp

    END DO

    ! create groups for set Bs

    DO n = 0, nprocb-1
       members = nproca
       CALL MPI_GROUP_INCL (group_world, members, mapmesh(n,:), group_b, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_GROUP_INCL failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF

       CALL MPI_COMM_CREATE (p_all_comm, group_b, p_communicator_tmp, &
            p_error)

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_COMM_CREATE failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', p_error
          CALL p_abort
       END IF
       IF(p_communicator_tmp/=MPI_COMM_NULL) &
         p_communicator_b = p_communicator_tmp

    END DO

    CALL MPI_BARRIER (p_all_comm, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' PE: ', mype, ' MPI_BARRIER failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', p_error
       CALL p_abort
    END IF

    IF (debug_parallel >= 0 .AND. mype == 0) THEN
      p_communicator_a = p_communicator_d
      p_communicator_b = p_communicator_d
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,3i8)') &
         'p_set_communicator on PE ', mype, ': ', &
         p_communicator_d, &
         p_communicator_a, &
         p_communicator_b
#endif
#endif
  END SUBROUTINE p_set_communicator

!=========================================================================

  ! send implementation

  SUBROUTINE p_send_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real

  SUBROUTINE p_send_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_1d

  SUBROUTINE p_send_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_char

! non-blocking sends
#ifndef _SX
  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest))
      tmp(:) = p_request(:)
      DEALLOCATE(p_request)
      ALLOCATE(p_request(p_mrequest+p_request_alloc_size))
      p_request(1:p_mrequest) = tmp(:)
      p_mrequest = p_mrequest+p_request_alloc_size
      DEALLOCATE(tmp)
    ENDIF
#endif

  END SUBROUTINE p_inc_request
#else
  SUBROUTINE p_inc_request
    INTEGER, POINTER     :: tmp(:)

#ifndef NOMPI
    p_irequest = p_irequest + 1
    IF (p_irequest > p_mrequest) THEN
      ALLOCATE(tmp(p_mrequest+p_request_alloc_size))
      tmp(:p_mrequest) = p_request(:)
      DEALLOCATE(p_request)
      p_request => tmp
      p_mrequest = p_mrequest+p_request_alloc_size
    ENDIF
#endif

  END SUBROUTINE p_inc_request
#endif

  SUBROUTINE p_isend_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real

  SUBROUTINE p_isend_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_isend_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_real_dp, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_int, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, SIZE(buffer), p_bool, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(inout) :: buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, LEN(buffer), p_char, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', mype, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_char

#ifdef ONESIDED
  ! some MPI2 things

  SUBROUTINE p_win_create_real_2d(buffer,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(out) :: p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_WIN_CREATE(buffer,SIZE(buffer)*p_rg,p_rg,MPI_INFO_NULL,p_comm,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_WIN_CREATE on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
   END SUBROUTINE p_win_create_real_2d

  SUBROUTINE p_win_create_real_3d(buffer,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(out) :: p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_WIN_CREATE(buffer,SIZE(buffer)*p_rg,p_rg,MPI_INFO_NULL,p_comm,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_WIN_CREATE on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif

  END SUBROUTINE p_win_create_real_3d

  SUBROUTINE p_win_fence_int(p_win)

    INTEGER, INTENT(in) :: p_win

#ifndef NOMPI

    CALL MPI_WIN_FENCE(0,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FENCE on ', mype, ' for win', p_win, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
  END SUBROUTINE p_win_fence_int

  SUBROUTINE p_win_fence_int_1d(p_win)

    INTEGER, INTENT(in) :: p_win(:)
    INTEGER :: i

#ifndef NOMPI

    do i=1,SIZE(p_win)
    CALL MPI_WIN_FENCE(0,p_win(i),p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FENCE on ', mype, ' for win', p_win(i), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
    end do

#endif
  END SUBROUTINE p_win_fence_int_1d

  SUBROUTINE p_win_free_int(p_win)

    INTEGER, INTENT(in) :: p_win

#ifndef NOMPI

    CALL MPI_WIN_FREE(p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FREE on ', mype, ' for win', p_win, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
  END SUBROUTINE p_win_free_int

  SUBROUTINE p_win_free_int_1d(p_win)

    INTEGER, INTENT(in) :: p_win(:)
    INTEGER :: i

#ifndef NOMPI

    do i=1,SIZE(p_win)
    CALL MPI_WIN_FREE(p_win(i),p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_WIN_FREE on ', mype, ' for win', p_win(i), ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
    end do

#endif
  END SUBROUTINE p_win_free_int_1d

  SUBROUTINE p_put_real_2d(buffer,p_destination,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_PUT(buffer,SIZE(buffer),p_real,p_destination,0,SIZE(buffer),p_real,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a)') ' MPI_PUT on ', mype,'to',p_destination',&
                               'through win', p_win',' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif

#endif
   END SUBROUTINE p_put_real_2d

  SUBROUTINE p_put_real_3d(buffer,p_destination,p_win,comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_win
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_PUT(buffer,SIZE(buffer),p_real,p_destination,0,SIZE(buffer),p_real,p_win,p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a)') ' MPI_PUT on ', mype,'to',p_destination',&
                               'through win', p_win',' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

   END SUBROUTINE p_put_real_3d
#endif /* not use_comm_MPI1 */
  ! recv implementation

  SUBROUTINE p_recv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real

  SUBROUTINE p_recv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
    buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = 0
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_int, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = .FALSE.
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = .FALSE.
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_1d

  SUBROUTINE p_recv_bool_2d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = .FALSE.
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = .FALSE.
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

#ifdef MESSY
      buffer = .FALSE.
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), p_bool, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(out) :: buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, LEN(buffer), p_char, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_char

! non-blocking receives

  SUBROUTINE p_irecv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, 1, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real

  SUBROUTINE p_irecv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_irecv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_real_dp, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_4d

! ju_ch_20110429+
  SUBROUTINE p_irecv_int_2d (buffer, p_source, p_tag, p_count, comm)

    INTEGER,  INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_int_i4, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_int_i4, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_2d

  SUBROUTINE p_irecv_int_3d (buffer, p_source, p_tag, p_count, comm)

    INTEGER,  INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_int_i4, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, SIZE(buffer), p_int_i4, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', mype, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_int_3d
! ju_ch_20110429-

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, SIZE(sendbuf), p_real_dp, p_dest,   p_tag, &
                        recvbuf, SIZE(recvbuf), p_real_dp, p_source, p_tag, &
                        p_comm, p_status, p_error)

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i4,a,i6,a)') ' MPI_SENDRECV by ', mype, &
            ' to ', p_dest, ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_sendrecv_real_4d

  ! bcast implementation

  SUBROUTINE p_bcast_real (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real

  SUBROUTINE p_bcast_real_1d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_1d

  SUBROUTINE p_bcast_real_2d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_2d

  SUBROUTINE p_bcast_real_3d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (buffer, p_source, comm)

    REAL (dp), INTENT(inout) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_dp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_int_i4 (buffer, p_source, comm)

    INTEGER (i4), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i4, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (buffer, p_source, comm)

    INTEGER (i8), INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8

!!$  SUBROUTINE p_bcast_int_1d (buffer, p_source, comm)
!!$
!!$    INTEGER, INTENT(inout) :: buffer(:)
!!$    INTEGER, INTENT(in)    :: p_source
!!$    INTEGER, OPTIONAL, INTENT(in) :: comm
!!$
!!$#ifndef NOMPI
!!$    INTEGER :: p_comm
!!$
!!$    IF (PRESENT(comm)) THEN
!!$       p_comm = comm
!!$    ELSE
!!$       p_comm = p_all_comm
!!$    ENDIF
!!$
!!$#ifdef DEBUG
!!$    nbcast = nbcast+1
!!$#endif
!!$
!!$    IF (npes == 1) THEN
!!$       RETURN
!!$    ELSE
!!$       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
!!$            p_comm, p_error)
!!$    ENDIF
!!$
!!$#ifdef DEBUG
!!$    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
!!$            ' with broadcast number ', nbcast, ' succesfull.'
!!$
!!$     IF (p_error /= MPI_SUCCESS) THEN
!!$       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
!!$            ' failed.'
!!$       WRITE (nerr,'(a,i4)') ' Error = ', p_error
!!$       CALL p_abort
!!$    END IF
!!$#endif
!!$#endif
!!$
!!$  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_2d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_4d


  SUBROUTINE p_bcast_bool (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, 1, p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (buffer, p_source, comm)

    LOGICAL, INTENT(inout) :: buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_bool, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (buffer, p_source, comm)

    CHARACTER(len=*),  INTENT(inout) :: buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, LEN(buffer), p_char, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_char

#ifdef MESSY
  FUNCTION p_time() RESULT(time_in_seconds)

    REAL(dp) :: time_in_seconds

#ifndef NOMPI
    DOUBLE PRECISION :: t

    t = MPI_Wtime()
#else
    REAL :: t

    CALL cpu_time(t)
#endif

    time_in_seconds = REAL(t,dp)

  END FUNCTION p_time

  FUNCTION p_lor_0d (zfield, comm) RESULT (p_lor)

    LOGICAL,          INTENT(in)  :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    LOGICAL                       :: p_lor

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_lor, 1, p_bool, &
            MPI_LOR, p_comm, p_error)
    ELSE
       p_lor = zfield
    END IF
#else
    p_lor = zfield
#endif

  END FUNCTION p_lor_0d
#endif

  SUBROUTINE p_bcast_char_1d (buffer, p_source, comm)

    CHARACTER (*), INTENT(inout) :: buffer(:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                      :: lexlength,i,flength

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       lexlength=LEN(buffer(1))
       flength=SIZE(buffer)
       lexlength=lexlength*flength
       CALL MPI_BCAST (buffer, lexlength, p_char, p_source, p_comm, p_error)
#ifdef DEBUG
       WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

       IF (p_error /= MPI_SUCCESS) THEN
          WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
               ' failed.'
          WRITE (nerr,'(a,i4)') ' Error = ', p_error
          CALL p_abort
       END IF
#endif
    ENDIF
#endif

  END SUBROUTINE p_bcast_char_1d

  ! um_ak_20130702+
  SUBROUTINE p_bcast_real_sp_0d (buffer, p_source, comm)

    REAL(sp), INTENT(INOUT)          :: buffer
    INTEGER,  INTENT(IN)             :: p_source
    INTEGER,  INTENT(IN),   OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
!!$    CALL MPI_BCAST (buffer, 1, p_real, p_source, &
       CALL MPI_BCAST (buffer, 1, p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_sp_0d

  SUBROUTINE p_bcast_real_sp_1d (buffer, p_source, comm)

    REAL(sp), INTENT(INOUT)          :: buffer(:)
    INTEGER,  INTENT(IN)             :: p_source
    INTEGER,  INTENT(IN),   OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
!!$    CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_sp_1d

  SUBROUTINE p_bcast_real_sp_2d (buffer, p_source, comm)

    REAL(sp), INTENT(INOUT)          :: buffer(:,:)
    INTEGER,  INTENT(IN)             :: p_source
    INTEGER,  INTENT(IN), OPTIONAL   :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
!!$    CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_sp_2d

  SUBROUTINE p_bcast_real_sp_3d (buffer, p_source, comm)

    REAL(sp), INTENT(INOUT)          :: buffer(:,:,:)
    INTEGER,  INTENT(IN)             :: p_source
    INTEGER,  INTENT(IN),   OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
!!$    CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_sp_3d

  SUBROUTINE p_bcast_real_sp_4d (buffer, p_source, comm)

    REAL(sp), INTENT(INOUT)          :: buffer(:,:,:,:)
    INTEGER,  INTENT(IN)             :: p_source
    INTEGER,  INTENT(IN),   OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
!!$    CALL MPI_BCAST (buffer, SIZE(buffer), p_real, p_source, &
       CALL MPI_BCAST (buffer, SIZE(buffer), p_real_sp, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_real_sp_4d

  SUBROUTINE p_bcast_int_i4_1d (buffer, p_source, comm)

    INTEGER(I4), INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int_i4, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i4_1d

  SUBROUTINE p_bcast_int_i8_1d (buffer, p_source, comm)

    INTEGER(I8), INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

#ifdef DEBUG
    nbcast = nbcast+1
#endif

    IF (npes == 1) THEN
       RETURN
    ELSE
       CALL MPI_BCAST (buffer, SIZE(buffer), p_int_i8, p_source, &
            p_comm, p_error)
    ENDIF

#ifdef DEBUG
    WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' with broadcast number ', nbcast, ' succesfull.'

     IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_bcast_int_i8_1d
  ! um_ak_20130702-

  ! probe implementation

  SUBROUTINE p_probe_real (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    REAL (dp), INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i = 1, p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_real_dp, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO

#endif

  END SUBROUTINE p_probe_real

  SUBROUTINE p_probe_int (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_int, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_int

  SUBROUTINE p_probe_bool (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    LOGICAL,   INTENT(in)  :: buffer
    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_bool, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_bool

  SUBROUTINE p_probe_char (buffer, p_tagcount, p_tagtable, p_source, &
&                          p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in)  :: buffer
    INTEGER,           INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,           INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm, i
    LOGICAL :: flag

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    p_tag = -1
    DO WHILE (p_tag == -1)
       DO i=1,p_tagcount
          CALL MPI_IPROBE (MPI_ANY_SOURCE, p_tagtable(i), p_comm, &
&                          flag, p_status, p_error)
#ifdef DEBUG
          IF (p_error /= MPI_SUCCESS) THEN
             WRITE (nerr,'(a,i4,a,i4,a)') ' MPI_IPROBE on ', mype, &
                  ' for tag ', p_tagtable(i), ' failed.'
             WRITE (nerr,'(a,i4)') ' Error = ', p_error
             CALL p_abort
          END IF
#endif
          IF (flag) THEN
             p_source = p_status(MPI_SOURCE)
             p_tag = p_status(MPI_TAG)
             CALL MPI_GET_COUNT(p_status, p_char, p_count, p_error)
#ifdef DEBUG
             IF (p_error /= MPI_SUCCESS) THEN
                WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_GET_COUNT on ', mype, &
                     ' for tag ', p_tag, ' from ' , p_source, ' failed.'
                WRITE (nerr,'(a,i4)') ' Error = ', p_error
                CALL p_abort
             END IF
#endif
             EXIT
          ELSE
             p_tag = -1
          END IF
       END DO
    END DO


#endif

  END SUBROUTINE p_probe_char

  SUBROUTINE p_wait
#ifndef NOMPI
    INTEGER :: i

    DO i = 1, p_irequest
        CALL MPI_WAIT(p_request(i), p_status, p_error)
    END DO
    p_irequest = 0
#endif
  END SUBROUTINE p_wait

  SUBROUTINE p_wait_any(return_pe)

    INTEGER, INTENT(out) :: return_pe
#ifndef NOMPI
    INTEGER :: i

    CALL MPI_WAITANY(p_irequest, p_request, i, p_status, p_error)
    IF (i == MPI_UNDEFINED) THEN
      p_irequest = 0
      return_pe = -1
    ELSE
      return_pe = p_status(MPI_SOURCE)
    ENDIF
#else
    return_pe = 0
#endif

  END SUBROUTINE p_wait_any

  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm
#ifndef NOMPI
    INTEGER :: com
! mz_kk_20081107+
!!$    com = MPI_COMM_WORLD; IF(PRESENT(comm)) com = comm
    ! With the introduction of MMD is p_all_comm only part of MPI_COMM_WORLD
    com = p_all_comm; IF(PRESENT(comm)) com = comm
! mz_kk_20081107-
    CALL MPI_BARRIER (com, p_error)

!#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', mype, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
!#endif
#endif

  END SUBROUTINE p_barrier

  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp)                      :: p_sum
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d

! ju_ch_20110429+
  FUNCTION p_sum_0d_int (zfield, comm) RESULT (p_sum)

    INTEGER                       :: p_sum
    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_int_i4, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d_int
! ju_ch_20110429-

  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d

#ifdef MESSY
  FUNCTION p_sum_2d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d

  FUNCTION p_sum_2d_int (zfield, comm) RESULT (p_sum)

    INTEGER,           INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_int_i4, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d_int

  FUNCTION p_sum_3d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2),SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_3d
#endif

  FUNCTION p_global_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum
    REAL(dp)                      :: pe_sums(SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_REDUCE (zfield, pe_sums, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_io, p_comm, p_error)
       p_sum = SUM(pe_sums)
    ELSE
       p_sum = SUM(zfield)
    END IF
#else
    p_sum = SUM(zfield)
#endif

  END FUNCTION p_global_sum_1d

  FUNCTION p_field_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_REDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
            MPI_SUM, p_io, p_comm, p_error)
        IF (.NOT. p_parallel_io) p_sum = 0.0_dp
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_field_sum_1d

  FUNCTION p_max_0d (zfield, comm) RESULT (p_max)

    REAL(dp)                      :: p_max
    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d

  FUNCTION p_max_1d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_2d

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, SIZE(zfield), p_real_dp, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_3d

  FUNCTION p_min_0d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, 1, p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_0d

  FUNCTION p_min_1d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_min, SIZE(zfield), p_real_dp, &
            MPI_MIN, p_comm, p_error)
    ELSE
       p_min = zfield
    END IF
#else
    p_min = zfield
#endif

  END FUNCTION p_min_3d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real_dp, &
                     recvbuf, SIZE(recvbuf), p_real_dp, &
                     p_dest, p_comm, p_error)
#else
     recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif
   END SUBROUTINE p_gather_real_1d2d

   ! mz_bk_20120725+
   SUBROUTINE p_allgather_real_2d2d (sendbuf, recvbuf, comm)

    REAL(dp),          INTENT(inout) :: sendbuf(:,:), recvbuf(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_ALLGATHER(sendbuf, SIZE(sendbuf), p_real, &
                     recvbuf, SIZE(sendbuf), p_real, &
                     p_comm, p_error)
#else
     recvbuf(:,:) = sendbuf(:,:)
#endif
   END SUBROUTINE p_allgather_real_2d2d
   ! mz_bk_20120725-

   ! op_pj_20160408+
   SUBROUTINE p_allgather_real_3d2d(recvbuf, sendbuf, comm)

     IMPLICIT NONE

     ! I/O
     REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: recvbuf
     REAL(DP), DIMENSION(:,:),   INTENT(IN)  :: sendbuf
     INTEGER,  OPTIONAL,         INTENT(IN)  :: comm

     ! LOCAL
#ifndef NOMPI
     INTEGER :: p_error
     INTEGER :: p_comm

     IF (PRESENT(comm)) THEN
        p_comm = comm
     ELSE
        p_comm = p_all_comm
     ENDIF

     CALL MPI_ALLGATHER(sendbuf, SIZE(sendbuf)  , p_real_dp, &
          recvbuf, SIZE(recvbuf)/SIZE(recvbuf,3), p_real_dp, &
          p_comm, p_error)
!!$  CALL MPI_GATHER(sendbuf, SIZE(sendbuf),                 p_real_dp, &
!!$                  recvbuf, SIZE(recvbuf)/SIZE(recvbuf,3), p_real_dp, &
!!$                  p_io, p_comm, p_error)
!!$  CALL MPI_BCAST(recvbuf, SIZE(recvbuf), p_real_dp, p_io, p_comm, p_error)
#else
     recvbuf(:,:,LBOUND(recvbuf,3)) = sendbuf(:,:)
#endif

   END SUBROUTINE p_allgather_real_3d2d
   ! op_pj_20160408-

   SUBROUTINE ask_omp_threads

#ifdef _OPENMP
   INTEGER :: threads,jp

    threads = OMP_GET_MAX_THREADS()

    IF (mype == 0) THEN

       ! write out thread usage of master PE

       WRITE (nerr,'(a,i4,a,i4,a)') &
            'MPIOM:  PE ', mype, ': using ', &
            threads, &
            ' OpenMP threads.'
#ifndef NOMPI
       ! recevie and write the remaining MPI PEs thread number

       DO jp = 1, npes-1
          CALL MPI_RECV(threads, 1, MPI_INTEGER, jp, jp, &
               p_all_comm, p_status, p_error)
          WRITE (nerr,'(a,i4,a,i4,a)') &
               'MPIOM:  PE ', jp, ': using ', &
               threads, &
               ' OpenMP threads.'
       ENDDO

    ELSE

       ! send threads per PE to master PE for write out

       CALL MPI_SEND (threads, 1, MPI_INTEGER, 0, mype, &
               p_all_comm, p_error)
#endif

    ENDIF

#endif

   END SUBROUTINE ask_omp_threads

#if defined(MESSY) && defined(MPIOM_2000)
!mz_ap_20100922+

  SUBROUTINE p_isend_mpitype_2d (p_mpitype, buffer, p_destination, p_tag, p_count, comm)
    INTEGER,   INTENT(in) :: p_mpitype
    REAL(dp), INTENT(inout) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_ISEND (buffer, p_count, p_mpitype, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_ISEND (buffer, 1, p_mpitype, p_destination, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#endif

  END SUBROUTINE p_isend_mpitype_2d

  SUBROUTINE p_irecv_mpitype_2d (p_mpitype, buffer, p_source, p_tag, p_count, comm)
    INTEGER,   INTENT(in) :: p_mpitype
    REAL(dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL p_inc_request
       CALL MPI_IRECV (buffer, p_count, p_mpitype, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    ELSE
       CALL p_inc_request
       CALL MPI_IRECV (buffer, 1, p_mpitype, p_source, p_tag, &
            p_comm, p_request(p_irequest), p_error)
    END IF

#endif

  END SUBROUTINE p_irecv_mpitype_2d

  SUBROUTINE p_sendrecv_mpitype_2d (p_mpitype,sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)
    INTEGER , INTENT(in)           :: p_mpitype
!    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    REAL(dp), INTENT(in)           :: sendbuf
    INTEGER,  INTENT(in)           :: p_dest
!    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    REAL(dp), INTENT(out)          :: recvbuf
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    CALL MPI_SENDRECV (sendbuf, 1, p_mpitype, p_dest,   p_tag, &
                        recvbuf, 1, p_mpitype, p_source, p_tag, &
                        p_comm, p_status, p_error)

#endif

  END SUBROUTINE p_sendrecv_mpitype_2d

SUBROUTINE p_gather_real_2d2d (sendbuf, recvbuf, p_dest, comm)

    REAL(dp),          INTENT(inout) :: sendbuf(:,:), recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

     CALL MPI_GATHER(sendbuf, SIZE(sendbuf), p_real, &
                     recvbuf, SIZE(sendbuf), p_real, &
                     p_dest, p_comm, p_error)

   END SUBROUTINE p_gather_real_2d2d


   SUBROUTINE init_MPI_datatypes(ie,je,ke)
#ifndef NOMPI
!      USE MO_PARAM1, ONLY : IE,JE,KE
!      use fields
!      use MPI_handles

     IMPLICIT NONE

     INTEGER, INTENT(IN):: ie,je,ke
     INTEGER :: ihalo


     ! do for all halos up to ihalo_max
     DO ihalo=1,ihalo_max

!     section 1 - North/South boundaries
!--------------------------------------------------------------------------

      ! 2d case
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,                &
            ie+2*(ihalo-1), je+2*(ihalo-1), 1,  &
            ie+2*(ihalo-1), ihalo,               &
                                     .FALSE., .FALSE.,                    &
            NS_BOUNDARY(ihalo) )

      ! 3d case (ke level at cell center)
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,                &
            ie+2*(ihalo-1), je+2*(ihalo-1), ke,  &
            ie+2*(ihalo-1), ihalo,               &
                                     .TRUE., .FALSE.,                    &
            NS_BOUNDARY_KE(ihalo) )

      ! 3d case (ke+1 level at cell interfaces)
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,                &
            ie+2*(ihalo-1), je+2*(ihalo-1), ke+1,  &
            ie+2*(ihalo-1), ihalo,               &
                                     .TRUE., .FALSE.,                    &
            NS_BOUNDARY_KEP(ihalo) )

!     section 2 - East/West boundaries
!--------------------------------------------------------------------------

      ! 2d case
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,              &
            ie+2*(ihalo-1),je+2*(ihalo-1), 1,  &
            ihalo, je,              &
                                     .FALSE., .TRUE.,                    &
            EW_boundary(ihalo) )

      ! 3d case (ke level at cell center)
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,              &
            ie+2*(ihalo-1),je+2*(ihalo-1), ke,  &
            ihalo, je,              &
                                     .TRUE., .TRUE.,                    &
            EW_boundary_KE(ihalo) )

      ! 3d case (ke+1 level at cell interfaces)
       CALL generate_MPI_subarray_type(MPI_DOUBLE_PRECISION,              &
            ie+2*(ihalo-1),je+2*(ihalo-1), ke+1,  &
            ihalo, je,              &
                                     .TRUE., .TRUE.,                    &
            EW_boundary_KEP(ihalo) )

     ENDDO


     ! reduction specific types
     !---------------------------------------------------------------------
     CALL generate_MPI_struct_type(2, (/3, 3/), (/p_real_dp, p_int/), &
          mmm_mpi_datatype_dp)
     CALL create_MPI_op(mmm_mpi_reduce_dp, .TRUE., op_mmm_mpi_reduce_dp)

#endif

   END SUBROUTINE init_MPI_datatypes

#ifndef NOMPI
!--------------------------------------------------------------------------
      subroutine generate_MPI_subarray_type(old_type, NX, NY, NZ,       &
                                           MX, MY, dim3, row_align,    &
                                           new_type)

!      use MPI_handles

      implicit none

      integer       :: NX,NY,NZ     ! dimensions of the large array
      integer       :: MX,MY        ! dimensions of the sub-array
      integer       :: old_type     ! old MPI data type
      integer       :: new_type     ! new MPI data type
      logical       :: dim3         ! is the sub-array two
                                    !    or three dimensional ?
      logical       :: row_align    ! .TRUE. if two adjacent elements
                                    !    are aligned along rows

!     This abstract MPI data type generator generates a data type
!     that describs a sub array
!
!     B( X1:X1+MX-1, Y1:Y1+MY-1, Z1 )          - if dim3 .eq. .FALSE.
!
!          or
!
!     B( X1:X1+MX-1, Y1:Y1+MY-1, Z0:Z0+NZ-1 )  - if dim3 .eq. .TRUE.
!
!       of a given 3 dimensional array
!
!     A( X0:X0+NX-1, Y0:Y0+NY-1, Z0:Z0+NZ-1 ) of type old_type.
!
!     Please note that definition below does not depend on the variables
!     X0, Y0, Z0, X1, Y1, Z1. These variables
!     can be arbitrary, as long as the following conditions are met.
!
!     X0 .le. X1,     X1+MX-1 .le. X0+NX-1
!     Y0 .le. Y1,     Y1+MY-1 .le. Y0+NY-1
!     Z0 .le. Z1,     Z1      .le. Z0+NZ-1
!     (in short - if B is a true sub-array of A
!
!     The logical variable row_align indicates where two "adjacent"
!     sub arrays are located.
!
!     If row_alingn .eq. .TRUE. then two adjacent sub-arrays are
!     aligned along rows, i.e. if B starts at (X1,Y1,Z1), the next
!     sub-array starts at (X1+MX,Y1,Z1).
!     If row_alingn .eq. .FALSE. then two adjacent sub-arrays are
!     aligned along columns, i.e. if B starts at (X1,Y1,Z1), the next
!     sub-array starts at (X1,Y1+MY,Z1).
!
!     Local variables

      integer (kind=MPI_ADDRESS_KIND)     :: sizeof_old_type, lb
      integer (kind=MPI_ADDRESS_KIND)     :: old_extent, new_extent
      integer       :: MPI_type_1dim      ! 1-dimensional building block
      integer       :: MPI_type_2dim      ! 2-dimensional building block
      integer       :: MPI_type_3dim      ! 3-dimensional building block
      integer       :: ierror

!     1. Define 1-dimensional building block
!--------------------------------------------------------------------------

      call MPI_Type_contiguous(MX, old_type, MPI_type_1dim, ierror)
      if ( ierror .ne. MPI_SUCCESS ) THEN
         call p_abort
      end if

!     2. Define 2-dimensional building block
!--------------------------------------------------------------------------

      call MPI_Type_get_extent(old_type, lb, sizeof_old_type, ierror)
      if ( ierror .ne. MPI_SUCCESS ) THEN
         call p_abort
      end if

      call MPI_Type_create_hvector(MY, 1, NX*sizeof_old_type,          &
                                  MPI_type_1dim,                       &
                                  MPI_type_2dim, ierror)
      if ( ierror .ne. MPI_SUCCESS ) THEN
         call p_abort
      end if

!     3. Define 3-dimensional building block
!--------------------------------------------------------------------------

      call MPI_Type_create_hvector(NZ, 1, NX*NY*sizeof_old_type,       &
                                  MPI_type_2dim,                      &
                                  MPI_type_3dim, ierror)
      if ( ierror .ne. MPI_SUCCESS ) THEN
         call p_abort
      end if

!     4. Resize the extent of the data type to allow proper alignment
!--------------------------------------------------------------------------

      if ( dim3 ) then
         call MPI_Type_get_extent(MPI_type_3dim, lb, old_extent, ierror)
         if ( ierror .ne. MPI_SUCCESS ) THEN
            call p_abort
         end if

         if ( row_align ) then
            new_extent = MX * sizeof_old_type
         else
            new_extent = MY * NX * sizeof_old_type
         endif

         call MPI_Type_create_resized(MPI_type_3dim, lb, new_extent,  &
                                     new_type, ierror)
         if ( ierror .ne. MPI_SUCCESS ) THEN
            call p_abort
         end if
      else
         call MPI_Type_get_extent(MPI_type_2dim, lb, old_extent, ierror)
         if ( ierror .ne. MPI_SUCCESS ) THEN
            call p_abort
         end if

         if ( row_align ) then
            new_extent = MX * sizeof_old_type
         else
            new_extent = MY * NX * sizeof_old_type
         endif

         call MPI_Type_create_resized(MPI_type_2dim, lb, new_extent,  &
                                     new_type, ierror)
         if ( ierror .ne. MPI_SUCCESS ) THEN
            call p_abort
         end if
      endif

!     5. Commit the MPI data type
!--------------------------------------------------------------------------

      call MPI_Type_commit(new_type, ierror)
!      write(iouerr,*) 'MPI_Type_commit',new_type
      if ( ierror .ne. MPI_SUCCESS ) THEN
         call p_abort
      end if


      end subroutine generate_MPI_subarray_type

      SUBROUTINE generate_MPI_struct_type(nblk, nelem, component_type, newtype)

        INTEGER, INTENT(IN)    :: nblk
        INTEGER, INTENT(IN)    :: nelem(0:nblk-1), component_type(0:nblk-1)
        INTEGER, INTENT(INOUT) :: newtype

        INTEGER             :: ierr, i
        INTEGER(KIND = MPI_ADDRESS_KIND) :: noffset(0:nblk-1), extent ,lb


        noffset(0) = 0
        DO i=0, nblk-2
          CALL MPI_TYPE_GET_EXTENT(component_type(i), lb, extent, ierr)
          noffset(i+1) = noffset(i) + nelem(i) * extent
        ENDDO

        CALL MPI_TYPE_CREATE_STRUCT(nblk, nelem, noffset, component_type, newtype, ierr)
        IF (ierr /= MPI_SUCCESS) THEN
          CALL p_abort
        END IF

        CALL MPI_TYPE_COMMIT(newtype, ierr)
        IF (ierr /= MPI_SUCCESS) THEN
          CALL p_abort
        END IF

      END SUBROUTINE generate_MPI_struct_type

      SUBROUTINE create_MPI_op(func, lcommute, newop)

        EXTERNAL func
        LOGICAL, INTENT(IN)  :: lcommute
        INTEGER, INTENT(OUT) :: newop
        INTEGER              :: ierr


        CALL MPI_OP_CREATE(func, lcommute, newop, ierr)
        IF (ierr /= MPI_SUCCESS) THEN
          CALL p_abort
        ENDIF

      END SUBROUTINE create_MPI_op

#endif

   SUBROUTINE mmm_mpi_reduce_dp(a, b, n, dtype)
     INTEGER, INTENT(in) :: n, dtype
     TYPE(min_mean_max_dp), INTENT(in) :: a(n)
     TYPE(min_mean_max_dp), INTENT(inout) :: b(n)

     INTEGER :: i, rpe, count_total

     DO i = 1, n
       rpe = MERGE(a(i)%pe_min, b(i)%pe_min, a(i)%min < b(i)%min)
       b(i)%min = MIN(a(i)%min, b(i)%min)
       b(i)%pe_min = rpe

       count_total = a(i)%rcount + b(i)%rcount
       b(i)%mean = (a(i)%mean * REAL(a(i)%rcount, dp) &
            + b(i)%mean * REAL(b(i)%rcount, dp)) &
            / REAL(count_total, dp)
       b(i)%rcount = count_total

       rpe = MERGE(a(i)%pe_max, b(i)%pe_max, a(i)%max >= b(i)%max)
       b(i)%max = MAX(a(i)%max, b(i)%max)
       b(i)%pe_max = rpe
     END DO
   END SUBROUTINE mmm_mpi_reduce_dp

!mz_ap_20100922-
#endif

END MODULE mo_mpi

! ##########################################################################
#else
!#include "mo_mpi_icon.inc"
#include "../messy/bmil/mo_mpi_icon.inc"
#endif
! ##########################################################################
