#define _DIAGONALEX_
! **************************************************************************
MODULE messy_main_mpi_bi
! **************************************************************************

  USE messy_main_constants_mem, ONLY: dp, sp, i4, i8, iouerr

#if defined(ECHAM5) || defined(MBM_CLAMS) || defined (MBM_MPIOM)

  USE mo_mpi,           ONLY: p_parallel_io               &  ! LOGICAL
                            , p_io, p_pe, p_nprocs        &  ! INTEGER
                            , p_bcast                     &  ! SUBROUTINES
                            , p_send, p_recv, p_sendrecv  &  ! SUBROUTINES
                            , p_all_comm                  &  ! INTEGER
                            , p_abort, p_parallel         &
                            , p_sum                       &  ! FUNCTION
                            , p_set_communicator          &  ! mz_ab_20100307+
                            , p_barrier                   &  ! mz_ab_20100518
                            , npes                        &  ! um_ak_20130502
                            , p_isend, p_irecv, p_wait    &  ! mz_ht_20140801
                            , p_max, p_lor

#ifndef NOMPI
  USE mo_mpi,           ONLY: p_error, p_real          ! um_ak_20130502
#endif
#endif

! ### ADDITINAL REQUESTS ...

#if defined(ECHAM5) || defined(MBM_MPIOM)
#ifdef MPIOM_2000
  USE mo_mpi,           ONLY: init_mpi_datatypes             !mz_ap_20101126
#endif
#endif

#if defined(ECHAM5)
  ! NON-MESSY; TO BE REPLACED BY error_bi and info_bi from BLATHER
  USE mo_exception,     ONLY: finish, message                ! SUBROUTINES
  ! DECOMPOSITION
  USE mo_decomposition,    ONLY: dcg => global_decomposition       &
                               , dcl => local_decomposition

  ! TRANSFORMATION
  USE mo_transpose,        ONLY: gather_gp, scatter_gp, indx, reorder &
                               , gather_sa, scatter_sa  &
                               , gather_sp, scatter_sp
  USE mo_tr_gather,        ONLY: gather_field
  USE mo_tr_allgather,     ONLY: allgather_field
#endif

#if defined(MBM_MPIOM)
  USE mo_mpi,              ONLY: p_start, p_stop
#endif

#if defined(ECHAM5) || defined(MBM_MPIOM)
  USE messy_mpiom_mem_e5,  ONLY:  gather_arr, scatter_arr &
                               ,  allgather_arr
#endif

#if defined(MBM_CLAMS)
  USE mo_mpi,           ONLY: p_isend, p_irecv            &
                            , p_start, p_stop, p_wait     &
                            , p_barrier
#endif

#if defined(COSMO)
  USE data_parallel,      ONLY: imp_reals, imp_grib,  imp_integers     &
                              , imp_byte, imp_character, imp_logical   &
                              , icomm_world, p_all_comm => icomm_world &
                              , p_nprocs => nproc, p_pe => my_world_id &
                              , nprocx, nprocy, isubpos, icomm_compute &
                              , num_compute, nboundlines, my_cart_id   &
                              , icomm_cart, ldatatypes, ncomm_type     &
                              , my_cart_neigh, sendbuf, isendbuflen    &
                              , icomm_cart, nboundlines, num_compute   &
                              , lasync_io
  USE parallel_utilities, ONLY: distribute_field, gather_field         &
                              , gather_values, global_values, ij_local
  USE environment,        ONLY: exchg_boundaries
#ifndef NOMPI
#ifndef __PGI
  USE mpi
#endif
#endif
#endif

#if defined(CESM1) && defined(SPMD)
  USE spmd_utils,         ONLY:  p_nprocs=> npes,          &
       imp_integers=>mpi_integer,imp_integers8=> mpi_integer8, &
       imp_character=>mpi_character,   &
       imp_logical=>mpi_logical, &
       imp_reals=> mpi_real8, imp_reals4=>mpi_real4,          &
       p_pe => iam, iam, p_all_comm=>mpicom, icomm_world=>mpicom, &
       p_real_dp=>mpi_real8, MPI_STATUS_SIZE
  USE mo_mpi, ONLY: MPI_SUM, MPI_MAX, MPI_LOR
#endif

  ! op_bk_20140116+
#if defined(ICON)
  USE mo_mpi,           ONLY: my_process_is_stdio              &  ! LOGICAL FUNC
       &                    , p_io=>process_mpi_stdio_id       &  ! INTEGER
       &                    , get_my_mpi_all_id                &  ! p_pe
       &                    , get_my_mpi_all_comm_size         &  ! p_nprocs
!!$    &                    , get_my_mpi_work_id               &  ! p_pe
       &                    , get_my_mpi_work_comm_size        &  ! p_nprocs
!!$    &                    , num_work_procs                   &  ! p_nprocs
       &                    , p_work0=>process_mpi_all_workroot_id &
       &                    , p_bcast_icon=>p_bcast            &  ! SUBROUTINES
       &                    , p_send, p_recv, p_sendrecv       &  ! SUBROUTINES
       &                    , p_isend, p_irecv                 &
!!$       &                    , p_all_comm=>process_mpi_all_comm &  ! INTEGER
       &                    , p_all_comm=>p_comm_work          &  ! INTEGER
       &                    , p_abort=>abort_mpi               &
       &                    , my_process_is_mpi_all_parallel   &  ! p_parallel
!!$    &                    , p_sum                            &  ! SUBROUTINE
!!$    &                    , my_process_is_mpi_all_parallel   &  ! p_parallel
!!$    &                    , my_process_is_stdio               ! p_parallel_io
       &                    , my_process_is_work               &
       &                    , my_process_is_io                 &
       &                    , p_comm_work                      &
       &                    , p_comm_work_io                   &
       &                    , p_comm_work_2_io                 &
       &                    , p_wait                           &
       &                    , p_max, p_lor
  USE mo_exception,     ONLY: finish, message                ! SUBROUTINES
#ifndef NOMPI
#if !defined (__SUNPRO_F95)
  USE mpi
#endif
#endif
#endif
! ... defined(ICON)
  ! op_bk_20130828-


! ######################################################################
  IMPLICIT NONE
  PUBLIC
  SAVE
! ######################################################################

  ! op_pj_20121001+
  INTERFACE p_allgather
     MODULE PROCEDURE p_allgather_2d1d
     MODULE PROCEDURE p_allgather_3d2d ! op_pj_20160408
     MODULE PROCEDURE p_allgather_4d3d
  END INTERFACE
  INTERFACE p_scatter
     MODULE PROCEDURE p_scatter_4d3d
  END INTERFACE
  ! op_pj_20121001-

#if ! defined(ECHAM5)
   TYPE decomp
      LOGICAL :: lreg = .TRUE.
   END TYPE decomp
   TYPE(decomp), POINTER :: dcg
   TYPE(decomp)          :: dcl
   INTERFACE reorder
      MODULE PROCEDURE reorder2
      MODULE PROCEDURE reorder3
      MODULE PROCEDURE reorder4
   END INTERFACE
#endif

#if defined(ECHAM5)
   PUBLIC :: dcl  ! um_ak_20120601
   PUBLIC :: dcg  ! um_ak_20130618
   PUBLIC :: p_pe ! um_ak_20120601
#endif

#if defined(COSMO)
#ifndef NOMPI
#ifdef __PGI
   INCLUDE 'mpif.h'
#endif
#endif
   INTEGER(I4), PARAMETER :: p_io = 0     ! CPU writing LOG-File
   INTEGER, PARAMETER :: iexchg_MPI_type_len    = 200
            ! length of global vector iexchg_MPI_type used in environment.f90

   INTEGER :: iexchg_MPI_types(iexchg_MPI_type_len) = MPI_DATATYPE_NULL
                  ! List of MPI data types used in exchg_datatypes
                  ! routine. If set to MPI_DATATYPE_NULL, the
                  ! corresponding entry has not been properly set up.
   INTEGER ::  iexchg_counts(iexchg_MPI_type_len)
                  ! List of counts of these data types
                  ! used in exchg_datatypes routine.
                  ! Has meaningful value only if the corresponding
                  ! vector element of iexchg_MPI_types is not
                  ! MPI_DATATYPE_NULL.
   LOGICAL :: p_parallel_io = .FALSE. ! TRUE if p_pe = p_io
   LOGICAL :: p_parallel    = .TRUE.  ! parallel environment

   INTERFACE global_fields
      MODULE PROCEDURE      &
           global_1d_intfield,   &
           global_1d_realfield,  &
           global_2d_intfield,   &
           global_2d_realfield,  &
           global_3d_intfield,   &
           global_3d_realfield
   END INTERFACE

   PUBLIC :: p_pe
   PUBLIC :: exchg_trac_boundaries
   PUBLIC :: exchg_trac_boundaries2
!!$EXTERNAL :: MPI_BCAST
#endif

#if defined(CESM1)
   PUBLIC :: finish, message
#if defined(SPMD)
   LOGICAL :: p_parallel_io = .FALSE. ! TRUE if p_pe = p_io
   LOGICAL :: p_parallel    = .TRUE.  ! parallel environment
   INTEGER(I4), PARAMETER :: p_io = 0 ! CPU writing LOG-File
#endif
#endif

! op_bk_20170223+
#if defined(ICON)
   INTEGER, PUBLIC :: p_pe, p_nprocs
   LOGICAL, PUBLIC :: p_parallel, p_parallel_io
   INTEGER, PARAMETER :: iexchg_MPI_type_len    = 200
   ! length of global vector iexchg_MPI_type used in environment.f90

   INTEGER :: iexchg_MPI_types(iexchg_MPI_type_len) = MPI_DATATYPE_NULL
                  ! List of MPI data types used in exchg_datatypes
                  ! routine. If set to MPI_DATATYPE_NULL, the
                  ! corresponding entry has not been properly set up.
   INTEGER :: iexchg_counts(iexchg_MPI_type_len)
                  ! List of counts of these data types
                  ! used in exchg_datatypes routine.
                  ! Has meaningful value only if the corresponding
                  ! vector element of iexchg_MPI_types is not
                  ! MPI_DATATYPE_NULL.
   INTEGER :: imp_reals = MPI_DOUBLE_PRECISION
   INTEGER :: imp_integers = MPI_INTEGER

   INTEGER, PUBLIC :: p_io0

   PUBLIC :: p_io
   PUBLIC :: p_work0
   PUBLIC :: p_abort
   PUBLIC :: p_all_comm
   PUBLIC :: p_comm_work
   PUBLIC :: p_comm_work_io
   PUBLIC :: p_comm_work_2_io

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_0d_int
     MODULE PROCEDURE p_sum_0d
     MODULE PROCEDURE p_sum_1d
     MODULE PROCEDURE p_sum_2d
     MODULE PROCEDURE p_sum_2d_int
     MODULE PROCEDURE p_sum_3d
     MODULE PROCEDURE p_sum_4d
  END INTERFACE
  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_real_single
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_real_1d_single
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_real_2d_single
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_real_5d
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_int_i8_1d
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_char_1d
  END INTERFACE

!!$  PUBLIC :: p_allgather
!!$  PUBLIC :: p_scatter
  PUBLIC :: p_sum
  PUBLIC :: p_bcast
  ! ub_ak_20180116+
  !PUBLIC :: messy_mpi_initialize
  PUBLIC :: messy_mpi_setup
  ! ub_ak_20180116-
#endif
! op_bk_20170223-

#if defined(BLANK)
  INTEGER, PARAMETER :: MPI_SUM = -1           ! MPI dummy (mpif.h or USE mpi)
#endif

!#if defined(CESM1) && defined(SPMD)
!#define NOMPI
!#endif
!#ifdef NOMPI
#if defined(BLANK) || (defined(CESM1) && !defined(SPMD)) || defined(VERTICO)
  ! DUMMY FOR (NON-)PARALLEL ENVIRONMENT
  INTEGER(I4), PARAMETER :: p_io = 0           ! CPU writing LOG-File
  INTEGER, PARAMETER :: p_pe         = 0       ! CPU number
  INTEGER, PARAMETER :: p_nprocs     = 1       ! number of parallel CPUs
  INTEGER, PARAMETER :: p_all_comm = 1         ! communicator dummy
! op_pj_20120618+
!!$  INTEGER(I4), PARAMETER :: icomm_world   = 1  ! communicator dummy
  INTEGER, PARAMETER :: icomm_world   = 1  ! communicator dummy
! op_pj_20120618-
  INTEGER, PARAMETER :: imp_integers = 1       ! MPI dummy
  INTEGER, PARAMETER :: imp_reals    = 2       ! MPI dummy
  INTEGER, PARAMETER :: p_real_dp    = 2       ! MPI dummy
  INTEGER, PARAMETER :: imp_logical  = 3       ! MPI dummy
  INTEGER, PARAMETER :: imp_character = 4      ! MPI dummy
  LOGICAL            :: p_parallel_io = .TRUE. ! TRUE if p_pe = p_io
  LOGICAL            :: p_parallel    = .FALSE.  ! parallel environment
!... defined(BLANK) || (defined(CESM1) && !defined(SPMD))  || defined(VERTICO)
#endif
!#endif

#if (defined COSMO) || defined(BLANK) || defined(CESM1) || defined(VERTICO)
  INTEGER :: itag = 88

  INTERFACE gather_gp
     MODULE PROCEDURE gather_gp432 ! gather gridp. field (nlon,nlev,ntrac,nlat)
                                   !                  or (nlon,nlev,nlat,1)
                                   !                  or (nlon,nlat)
     MODULE PROCEDURE gather_gp32  ! gather gridp. field (nlon,nlev,nlat)
                                   !                 or  (nlon,nlat,1)
     MODULE PROCEDURE gather_gp2   ! gather only m=0 wave number (nlon,nlat)
  END INTERFACE

  INTERFACE scatter_gp
     MODULE PROCEDURE scatter_gp432! scatter gridp. field (nlon,nlev,ntrac,nlat)
                                   !                   or (nlon,nlev,nlat,1)
     MODULE PROCEDURE scatter_gp32 ! scatter gridp. field (nlon,nlev,nlat)
                                   !                   or (nlon,nlat,1)
     MODULE PROCEDURE scatter_gp2  ! scatter gridp. field (nlon,nlat)
     ! um_ak_20121022+
     MODULE PROCEDURE scatter_gp432_sp   ! scatter gridp. field
                                         !                (nlon,nlev,ntrac,nlat)
                                         !             or (nlon,nlev,nlat,1)
     MODULE PROCEDURE scatter_gp32_sp    ! scatter gridp. field (nlon,nlev,nlat)
                                         !             or (nlon,nlat,1)
     MODULE PROCEDURE scatter_gp2_sp     ! scatter gridp. field (nlon,nlat)
     ! um_ak_20121022-
#ifdef CESM1
     MODULE PROCEDURE scatter_gp1  ! scatter field (ncol)
#endif
  END INTERFACE

  INTERFACE p_bcast
     MODULE PROCEDURE              &
          distribute_kind8,        &
          distribute_kind4,        &
          distribute_oneinteger,   &
          distribute_dp,           &
          distribute_sp,           &
          distribute_onedouble,    &
          distribute_onesingle,    &
          distribute_logical,      &
          distribute_onelogical,   &
          distribute_character,    &
          distribute_onecharacter, &
          bcast_4d,                &
          bcast_4d_sp,             &
          bcast_3d,                &  ! mz_ab_20140401
          bcast_2d                    ! mz_ab_20140401
  END INTERFACE

  ! op_pj_20120925+
  INTERFACE p_sum
     MODULE PROCEDURE p_sum_0d
     MODULE PROCEDURE p_sum_1d
     MODULE PROCEDURE p_sum_2d
  END INTERFACE p_sum
  ! op_pj_20120925-
  INTERFACE p_lor
     MODULE PROCEDURE p_lor_0d
  END INTERFACE p_lor
  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
  END INTERFACE p_max
#endif
! ... (defined COSMO) || defined(BLANK) || defined(CESM1) || defined(VERTICO)

#if defined(COSMO) || defined(BLANK) || defined(CESM1)
  INTERFACE p_send
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE
  INTERFACE p_recv
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE
!... defined(COSMO) || defined(BLANK) || defined(CESM1)
#endif

! ######################################################################
 CONTAINS
! ######################################################################

! ############################################################################
! ############################################################################

! ############################################################################
! ############################################################################

! ############################################################################
! ############################################################################
#if ! defined(ECHAM5)
!==============================================================================
  SUBROUTINE reorder2 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:)

#if (defined CRAY) || (defined sun) || (defined NAG) || (defined __SX__)
    CALL util_reshape(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2))
#else
    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))
#endif

  END SUBROUTINE reorder2
!------------------------------------------------------------------------------
  SUBROUTINE reorder3 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:)
    INTEGER :: k

    DO k=1,SIZE(x,2)
      y(:,k,:) = RESHAPE (x(:,k,:),(/SIZE(y,1),SIZE(y,3)/),(/0._dp/))
    END DO

  END SUBROUTINE reorder3
!------------------------------------------------------------------------------
  SUBROUTINE reorder4 (y,x)
    REAL(dp) ,INTENT(out) :: y (:,:,:,:)
    REAL(dp) ,INTENT(in)  :: x (:,:,:,:)
    INTEGER :: k, l

    DO l=1,SIZE(x,3)
      DO k=1,SIZE(x,2)
        y(:,k,l,:) = RESHAPE (x(:,k,l,:),(/SIZE(y,1),SIZE(y,4)/),(/0._dp/))
      END DO
    END DO

  END SUBROUTINE reorder4
!==============================================================================

!... ! defined(ECHAM5)
#endif
! ############################################################################
! ############################################################################

! ############################################################################
! ############################################################################
#if (defined COSMO) || defined(BLANK) || defined(CESM1) || defined(VERTICO) || defined(MESSYDWARF)
!------------------------------------------------------------------------------
  SUBROUTINE main_mpi_setup
    !
    ! Author: Astrid Kerkweg, Uni-Mainz, Mar 2008
    !
    ! This subroutines initializes some variables needed in MESSy
    !
    IF (p_nprocs > 1 ) THEN
       p_parallel = .TRUE.
    ELSE
       p_parallel = .FALSE.
    ENDIF

    p_parallel_io = (p_pe == p_io)

    ALLOCATE(dcg) ! mz_ab_20131203

  END SUBROUTINE main_mpi_setup
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE p_abort(pstr,qstr)
!
! Author: Astrid Kerkweg, Uni-Mainz, Mar 2008
!
! This subroutines calls the COSMO model_abort routine
!
#ifdef COSMO
  USE environment,      ONLY: model_abort
#endif
#ifdef CESM1
 use abortutils,   only: endrun
#endif
  IMPLICIT NONE
  INTRINSIC :: PRESENT, TRIM

  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: pstr
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: qstr

  CHARACTER(LEN=256) :: lpstr
  CHARACTER(LEN=256) :: lqstr

  IF (PRESENT(pstr) ) THEN
     lpstr=pstr
  ELSE
     lpstr='messy_main_mpi_bi'
  ENDIF
  IF (PRESENT(qstr) ) THEN
     lqstr=qstr
  ELSE
     lqstr='forced exit'
  ENDIF

#ifdef COSMO
  CALL model_abort(p_pe,42,TRIM(lpstr),TRIM(lqstr))
#endif
#ifdef CESM1
  CALL endrun(lpstr//','//lqstr)
#endif
END SUBROUTINE p_abort
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE bcast_4d (buffer, sender)

  IMPLICIT NONE
  INTRINSIC :: INT, SIZE

  REAL(dp), INTENT(INOUT) :: buffer(:,:,:,:)
  INTEGER,  INTENT(IN)    :: sender

  ! LOCAL
  INTEGER(i4) :: isender
  INTEGER     :: i,j,k

  isender=INT(sender,i4)

  DO i=1,SIZE(buffer,4)
     DO j=1,SIZE(buffer,3)
        DO k=1,SIZE(buffer,2)
           CALL p_bcast(buffer(:,k,j,i), isender)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE bcast_4d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE bcast_4d_sp (buffer, sender)

  IMPLICIT NONE
  INTRINSIC :: INT, SIZE

  REAL,     INTENT(INOUT) :: buffer(:,:,:,:)
  INTEGER,  INTENT(IN)    :: sender

  ! LOCAL
  INTEGER(i4) :: isender
  INTEGER     :: i,j,k

  isender=INT(sender,i4)

  DO i=1,SIZE(buffer,4)
     DO j=1,SIZE(buffer,3)
        DO k=1,SIZE(buffer,2)
           CALL p_bcast(buffer(:,k,j,i), isender)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE bcast_4d_sp
!------------------------------------------------------------------------------

! mz_ab_20140401+
!------------------------------------------------------------------------------
SUBROUTINE bcast_3d (buffer, sender)

  IMPLICIT NONE
  INTRINSIC :: INT, SIZE

  REAL(dp), INTENT(INOUT) :: buffer(:,:,:)
  INTEGER,  INTENT(IN)    :: sender

  ! LOCAL
  INTEGER(i4) :: isender
  INTEGER     :: i,j,k

  isender=INT(sender,i4)

  DO j=1,SIZE(buffer,3)
     DO k=1,SIZE(buffer,2)
        CALL p_bcast(buffer(:,k,j), isender)
     ENDDO
  ENDDO

END SUBROUTINE bcast_3d
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE bcast_2d (buffer, sender)

  IMPLICIT NONE
  INTRINSIC :: INT, SIZE

  REAL(dp), INTENT(INOUT) :: buffer(:,:)
  INTEGER,  INTENT(IN)    :: sender

  ! LOCAL
  INTEGER(i4) :: isender
  INTEGER     :: i,j,k

  isender=INT(sender,i4)

  DO k=1,SIZE(buffer,2)
     CALL p_bcast(buffer(:,k), isender)
  ENDDO

END SUBROUTINE bcast_2d
!------------------------------------------------------------------------------
! mz_ab_20140401-

!------------------------------------------------------------------------------
!
! SUBROUTINE p_bcast(buffer, isender, icomm, ibufferlen, idatatype,
!                                ierrorcode)
!
!------------------------------------------------------------------------------
!
! Description:
!  p_bcast is a generic name for several subroutines that distribute
!  values from one processor to all others. Depending on the type of the
!  first argument, the appropriate procedure is chosen.
!
! Method:
!  With the MPI_BCAST command the buffer is distributed to all other
!  processors.
!
!  This routines are -in principle- copies of the COSMO routines
!  distribute_values (parallel_utilities.f90) but they arguments that are
!  not given by the p_bcast calls in the MESSy submodel interface layer are
!  made OPTIONAL parameter here.
!
!------------------------------------------------------------------------------

! Following are the different subroutines

!------------------------------------------------------------------------------

!==============================================================================

!+ Subroutine for array of kind=8 integers

SUBROUTINE distribute_kind8 (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
INTEGER (KIND=i8), INTENT(INOUT) :: buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender   ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):

INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4)  :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length
INTEGER (KIND=i4) :: p_idatatype   ! local data type
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,              &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_kind8

!==============================================================================

!==============================================================================

!+ Subroutine for array of kind=4 integers

SUBROUTINE distribute_kind4 (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

INTEGER (KIND=i4), INTENT(INOUT) :: buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender   ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                         &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4)  :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length
INTEGER (KIND=i4) :: p_idatatype        ! local data type
INTEGER (KIND=i4) :: p_icomm    ! local communicator
!- End of header
!------------------------------------------------------------------------------
IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,                 &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_kind4

!==============================================================================

!==============================================================================

!+ Subroutine for one model integers

SUBROUTINE distribute_oneinteger(buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
INTEGER (KIND=i4),  INTENT(INOUT) :: buffer  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)     :: isender ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::                          &
  ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_idatatype   ! local data type
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_integers
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST ( buffer, 1, p_idatatype, isender, p_icomm &
                 , p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, 1, p_idatatype, isender, p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_oneinteger

!==============================================================================

!==============================================================================

!+ Subroutine for array of doubles

SUBROUTINE distribute_dp (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=dp),  INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)  :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen, & ! length of the buffer
  idatatype,  & ! type of buffer
  icomm         ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode    ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_ibufferlen ! local buffer length
INTEGER (KIND=i4) :: p_idatatype  ! local data type
INTEGER (KIND=i4) :: p_icomm      ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif
IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_dp

!==============================================================================

!==============================================================================

!+ Subroutine for one double

SUBROUTINE distribute_onedouble (buffer, isender, ibufferlen, icomm,  &
                             idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
REAL (KIND=dp),  INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)  :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length
INTEGER (KIND=i4) :: p_idatatype        ! local data type
INTEGER (KIND=i4) :: p_icomm    ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = 1
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
 CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onedouble

!==============================================================================

!==============================================================================

!+ Subroutine for array of singles

SUBROUTINE distribute_sp (buffer, isender, ibufferlen  &
                        , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=sp),    INTENT(INOUT) ::  buffer(:) ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode     ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length
INTEGER (KIND=i4) :: p_idatatype   ! local data type
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_sp

!==============================================================================

!==============================================================================

!+ Subroutine for one single

SUBROUTINE distribute_onesingle (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)
  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
REAL (KIND=sp),    INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode   ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode   ! local error code
INTEGER (KIND=i4) :: p_ibufferlen   ! local buffer length
INTEGER (KIND=i4) :: p_idatatype    ! local data type
INTEGER (KIND=i4) :: p_icomm        ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = 1
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_reals
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onesingle

!==============================================================================

!==============================================================================

!+ Subroutine for array of default logicals

SUBROUTINE distribute_logical  (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
LOGICAL,           INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode  ! local error code
INTEGER (KIND=i4) :: p_ibufferlen  ! local buffer length
INTEGER (KIND=i4) :: p_idatatype   ! local data type
INTEGER (KIND=i4) :: p_icomm       ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_logical
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
  CALL MPI_BCAST (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, p_ibufferlen, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
  ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_logical

!==============================================================================

!==============================================================================

!+ Subroutine for one default logical

SUBROUTINE distribute_onelogical (buffer, isender, ibufferlen  &
                                , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments
! Array arguments with intent(inout):
LOGICAL,           INTENT(INOUT) :: buffer     ! buffer to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL ::  ierrorcode   ! error code

!------------------------------------------------------------------------------

! Local Scalars
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_idatatype  ! local data type
INTEGER (KIND=i4) :: p_icomm      ! local communicator
!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_logical
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF
#ifdef COSMO
 CALL MPI_BCAST (buffer, 1, p_idatatype, isender, p_icomm &
               , p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, 1, p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
p_ierrorcode = 0
#endif

IF (PRESENT(ierrorcode)) THEN
   ierrorcode = p_ierrorcode
ENDIF

END SUBROUTINE distribute_onelogical

!==============================================================================

!==============================================================================

!+ Subroutine for array of characters

SUBROUTINE distribute_character (buffer, isender, ibufferlen  &
                               , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: CHAR, ICHAR, PRESENT, SIZE
#ifdef COSMO
  EXTERNAL  :: MPI_COMM_RANK
#endif

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments


! Array arguments with intent(inout):
! op_pj_20110313+
!!$CHARACTER (LEN=100), INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted
CHARACTER (LEN=*), INTENT(INOUT) :: buffer(:)  ! buffer to be broadcasted
! op_pj_20110313-

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)      :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors

! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode  ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER (KIND=i4) :: p_ierrorcode       ! local error code
INTEGER (KIND=i4) :: p_ibufferlen       ! local buffer length
INTEGER (KIND=i4) :: p_icomm            ! local communicator
INTEGER (KIND=i4) ::  my_comm_id, i, j

INTEGER    :: intbuf(300)      ! Standard integer ! mz_ab_20130515 for exec checksum
!- End of header
!------------------------------------------------------------------------------


 IF (PRESENT(ibufferlen)) THEN
   p_ibufferlen = ibufferlen
ELSE
   p_ibufferlen = SIZE(buffer)
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF
#if defined(COSMO) || (defined(CESM1) && defined(SPMD))
  CALL MPI_COMM_RANK(p_icomm, my_comm_id, p_ierrorcode)
#else
  my_comm_id=0
#endif
  DO i=1,p_ibufferlen
    IF (my_comm_id == isender) THEN
       ! um_ak_20110721+
       !DO j=1,100
       DO j=1,LEN(buffer(i))
       ! um_ak_20110721-
        intbuf(j) = ICHAR ( buffer(i)(j:j) )
      ENDDO
    ENDIF

#ifdef COSMO
    CALL MPI_BCAST (intbuf, 300, imp_integers, isender  &
                  , p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (intbuf, 300, imp_integers, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
    p_ierrorcode = 0
#endif

    IF (my_comm_id /= isender ) THEN
      !um_ak_20110721+
      !DO j=1,100
      DO j=1,LEN(buffer(i))
      !um_ak_20110721-
        buffer(i)(j:j) = CHAR (intbuf(j) )
      ENDDO
    ENDIF
  ENDDO

! and this would be the normal way
! CALL MPI_BCAST (buffer, ibufferlen, MPI_CHARACTER, isender,   &
!                 icomm, implcode)

  IF (PRESENT(ierrorcode)) THEN
    ierrorcode = p_ierrorcode
  ENDIF

END SUBROUTINE distribute_character

!==============================================================================

!==============================================================================

!+ Subroutine for one word of characters

SUBROUTINE distribute_onecharacter (buffer, isender, ibufferlen  &
                              , icomm, idatatype, ierrorcode)

  IMPLICIT NONE

  INTRINSIC :: LEN, PRESENT

!------------------------------------------------------------------------------
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!
! Subroutine / Function arguments

! Array arguments with intent(inout):
CHARACTER (LEN=*), INTENT(INOUT) ::  buffer     ! character to be broadcasted

! Scalar arguments with intent(in):
INTEGER (KIND=i4), INTENT(IN)    :: isender    ! sending processor

INTEGER (KIND=i4), INTENT(IN),  OPTIONAL ::                         &
  ibufferlen,         & ! length of the buffer
  idatatype,          & ! type of buffer
  icomm                 ! involved group of processors


! Scalar arguments with intent(out):
INTEGER (KIND=i4), INTENT(OUT), OPTIONAL :: ierrorcode            ! error code

!------------------------------------------------------------------------------

! Local Scalars and Arrays
INTEGER (KIND=i4) :: p_ierrorcode ! local error code
INTEGER (KIND=i4) :: p_idatatype  ! local data type
INTEGER (KIND=i4) :: p_icomm      ! local communicator

!- End of header
!------------------------------------------------------------------------------

IF (PRESENT(idatatype)) THEN
   p_idatatype = idatatype
ELSE
   p_idatatype = imp_character
ENDIF

IF (PRESENT(icomm)) THEN
   p_icomm = icomm
ELSE
   p_icomm = icomm_world
ENDIF

#ifdef COSMO
 CALL MPI_BCAST (buffer, len(buffer), p_idatatype, isender,              &
                 p_icomm, p_ierrorcode)
#elif (defined(CESM1) && defined(SPMD))
  CALL mpibcast (buffer, len(buffer), p_idatatype, isender,   &
                  p_icomm)
  p_ierrorcode = 0
#else
 p_ierrorcode = 0
#endif

 IF (PRESENT(ierrorcode)) THEN
    ierrorcode = p_ierrorcode
 ENDIF

END SUBROUTINE distribute_onecharacter

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

! op_pj_20120925+
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
       CALL MPI_ALLREDUCE (zfield, p_lor, 1, imp_logical, &
            MPI_LOR, p_comm, p_error)
    ELSE
       p_lor = zfield
    END IF
#else
    p_lor = zfield
#endif

  END FUNCTION p_lor_0d

  FUNCTION p_max_0d (zfield, comm) RESULT (p_max)

    REAL(DP),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(DP)                      :: p_max

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_max, 1, imp_reals, &
            MPI_MAX, p_comm, p_error)
    ELSE
       p_max = zfield
    END IF
#else
    p_max = zfield
#endif

  END FUNCTION p_max_0d

  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d

  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d

  FUNCTION p_sum_2d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d
! op_pj_20120925-

!==============================================================================
  ! um_ak_20140415+
  !SUBROUTINE gather_gp432 (gl, lc, gl_dc, source, lg32)
  SUBROUTINE gather_gp432 (gl, lc, gl_dc, source, lg32, lsum)
  ! um_ak_20140415-

  !
  ! gather global grid point field from pe-s (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !

    IMPLICIT NONE
    INTRINSIC :: SIZE

    ! INPUT
    REAL(dp),             POINTER     :: gl   (:,:,:,:) ! global field
    REAL(dp), TARGET,     INTENT(in)  :: lc   (:,:,:,:) ! local  field
    TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc          ! global decomposition
    INTEGER, OPTIONAL    ,INTENT(in)  :: source         ! source to gather from
    !                                                 ! -1=all;0=p_io;1=not p_io
    LOGICAL, OPTIONAL,    INTENT(in)  :: lg32           ! um_ak_20110603
    LOGICAL, OPTIONAL,    INTENT(IN)  :: lsum           ! um_ak_20140415
    ! LOCAL
    INTEGER(i4) :: size4 ! size of 4th dimension
    INTEGER(i4) :: size3 ! size of 3rd dimension ! mz_ab_20140401
    !
    INTEGER :: i
    REAL(dp), POINTER :: gl3(:,:,:)
    !
    ! call 2D gather routine if 4th dimension size is 1
    ! else loop over 3th index
    !
#ifdef COSMO
    IF (p_parallel_io) size4 = SIZE(gl,4)
    CALL p_bcast (size4, p_io)
    DO i=1,size4
       NULLIFY (gl3)
       gl3 => gl(:,:,:,i)
       CALL gather_gp32 (gl3, lc(:,:,:,i), gl_dc, &
            source=source, lg32=lg32, lsum=lsum)! um_ak_20140415 lsum added
    END DO
#elif defined(CESM1)
    IF (p_parallel_io) size4 = SIZE(gl,4)
    CALL p_bcast (size4, p_io)
    IF (size4 == 1) THEN
       NULLIFY (gl3)
       gl3 => gl(:,:,:,1)
       CALL gather_gp32 (gl3, lc(:,:,:,1), gl_dc, &
            source=source, lg32=lg32)
    ELSE
       IF (p_parallel_io) size3 = SIZE(gl,3)
       CALL p_bcast (size3, p_io)
       DO i=1,size3
          NULLIFY (gl3)
          gl3 => gl(:,:,i,:)
          CALL gather_gp32 (gl3, lc(:,:,i,:), gl_dc, &
               source=source, lg32=lg32)
       END DO
    END IF
#else
    gl(:,:,:,:) = lc(:,:,:,:)
#endif
    !
  END SUBROUTINE gather_gp432
!------------------------------------------------------------------------------
  ! um_ak_20140415+
  !SUBROUTINE gather_gp32 (gl, lc, gl_dc,source, lg32)
  SUBROUTINE gather_gp32 (gl, lc, gl_dc,source, lg32, lsum)
  ! um_ak_20140415-
  !
  ! gather global grid point field from pe-s (nlon,nlev,nlat) or (nlon,nlat,1)
  !
    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL(dp)             ,POINTER     :: gl   (:,:,:) ! global field
    REAL(dp)    , TARGET ,INTENT(in)  :: lc   (:,:,:) ! local  field
    TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc        ! global decomposition
    INTEGER, OPTIONAL    ,INTENT(in)  :: source       ! source to gather from
    !                                                 ! -1=all;0=p_io;1=not p_io
    LOGICAL, OPTIONAL    , INTENT(in) :: lg32         ! um_ak_20110603
    LOGICAL, OPTIONAL    , INTENT(in) :: lsum         ! um_ak_20140415

    ! LOCAL
    INTEGER(i4) :: isize ! size of 3rd dimension
    INTEGER     :: i
    REAL(dp),POINTER :: gl2(:,:)
    LOGICAL     :: l3
    !
    ! call 2D gather routine if 3rd dimension size is 1
    ! else call 3D gather routine
    !
#ifdef COSMO
!!$    l3=.TRUE.
!!$    IF (PRESENT(lg32)) THEN
!!$       IF (lg32) l3=.FALSE.
!!$    ENDIF
!!$
!!$    IF (l3) THEN
       IF (p_parallel_io) isize = SIZE(gl,3)
       CALL p_bcast(isize, p_io)
       NULLIFY (gl2)
       DO i=1,isize
          !IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,:,i)
          CALL gather_gp2 (gl2, lc(:,:,i),lsum=lsum)! um_ak_20140415 lsum added
       ENDDO
!!$    ELSE
!!$       IF (p_parallel_io) isize = SIZE(gl,2)
!!$       CALL p_bcast(isize, p_io)
!!$       NULLIFY (gl2)
!!$       DO i=1,isize
!!$          !IF (p_pe == p_io) gl2 => gl(:,:,i)
!!$          gl2 => gl(:,i,:)
!!$          CALL gather_gp2 (gl2, lc(:,i,:),lsum=lsum)! um_ak_20140415 lsum added
!!$       ENDDO
!!$
!!$    ENDIF

#elif defined(CESM1)
    l3=.FALSE.
    IF (PRESENT(lg32)) THEN
       IF (lg32) l3=.TRUE.
    ENDIF

    IF (l3) THEN
       IF (p_parallel_io) isize = SIZE(gl,3)
       CALL p_bcast(isize, p_io)
       NULLIFY (gl2)
       DO i=1,isize
          gl2 => gl(:,:,i)
          CALL gather_gp2 (gl2, lc(:,:,i))
       ENDDO
    ELSE
       IF (p_parallel_io) isize = SIZE(gl,2)
       CALL p_bcast(isize, p_io)
       NULLIFY (gl2)
       DO i=1,isize
          !IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,i,:)
          CALL gather_gp2 (gl2, lc(:,i,:))
       ENDDO
    ENDIF
#else
    gl(:,:,:) = lc(:,:,:)
#endif
  END SUBROUTINE gather_gp32
!------------------------------------------------------------------------------
  ! um_ak_20140415+
  !SUBROUTINE gather_gp2 (gl, lc)
  SUBROUTINE gather_gp2 (gl, lc, lsum)
  ! um_ak_20140415-
#ifdef CESM1
    use phys_grid, only: gather_chunk_to_field
#endif

    IMPLICIT NONE
    INTRINSIC :: SIZE

    REAL(dp),POINTER               :: gl   (:,:) ! global field
    REAL(dp), TARGET  ,INTENT(in)  :: lc   (:,:) ! local  field
    LOGICAL, OPTIONAL , INTENT(in) :: lsum   ! um_ak_20140415 lsum added
    !
    INTEGER(i4) :: l_size(2), g_size(2)
    INTEGER     :: ierror = 0
    !
#ifdef COSMO
    IF (p_pe == p_io) THEN
       g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
    END IF
    CALL p_bcast (g_size, p_io)
    l_size = (/ SIZE(lc,1), SIZE(lc,2) /)

    CALL gather_field(lc,l_size(1),l_size(2),gl,g_size(1),g_size(2) &
         ,p_io, ierror, lsum) ! um_ak_20140415 lsum added

    IF (ierror /= 0) THEN
       CALL p_abort('gather_gp2','messy_main_mpi_bi')
       RETURN
    ENDIF
#elif  defined(CESM1)
    call gather_chunk_to_field(1,1,1,SIZE(gl,1),lc,gl)
#else
    gl(:,:) = lc(:,:)
#endif
  END SUBROUTINE gather_gp2
!==============================================================================


!==============================================================================
  SUBROUTINE scatter_gp432 (gl, lc, gl_dc, lg32)
  !
  ! scatter global grid point field from pe-s (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
  ! HOMMESE:    (ncol,nlev,ntrac,1) -> (nproma,nlev,ntrac,ngpblks)
  !             (ncol,nlev,1,1)     -> (nproma,nlev,ngpblks,1)
  !             (ncol,1,1,1)        -> (nproma,ngpblks,1,1)

    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL(dp), POINTER                 :: gl   (:,:,:,:) ! global field
    REAL(dp), TARGET ,INTENT(out)     :: lc   (:,:,:,:) ! local  field
    TYPE(DECOMP),OPTIONAL,INTENT(in)  :: gl_dc          ! global decomposition
    LOGICAL, OPTIONAL,    INTENT(in)  :: lg32           ! um_ak_20110603
    ! LOCAL
    INTEGER(i4)      :: size4 ! size of 4th dimension
    INTEGER(i4)      :: size3 ! size of 3rd dimension ! mz_ab_20140401
    INTEGER          :: i
    REAL(dp),POINTER :: gl3(:,:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3th index
    !
#if defined(COSMO)
    IF (p_parallel_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    NULLIFY(gl3)
    DO i=1,size4
!       IF (p_pe == p_io) gl3 => gl(:,:,:,i)
       gl3 => gl(:,:,:,i)
       CALL scatter_gp32 (gl3, lc(:,:,:,i), gl_dc, lg32=lg32)! um_ak_20110603 lg32 added)
    END DO
#elif defined(CESM1)
#ifndef HOMMESE
    IF (p_parallel_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    IF (size4 == 1) THEN
       NULLIFY (gl3)
       IF (p_parallel_io) gl3 => gl(:,:,:,1)
       CALL scatter_gp32 (gl3, lc(:,:,:,1), gl_dc, lg32=lg32)
    ELSE
       IF (p_parallel_io) size3 = SIZE(gl,3)
       CALL p_bcast (size3, p_io)
       DO i=1,size3
          NULLIFY (gl3)
          IF (p_parallel_io) gl3 => gl(:,:,i,:)
          CALL scatter_gp32 (gl3, lc(:,:,i,:), gl_dc, lg32=lg32)
       END DO
    END IF
#else
!!$    IF (p_parallel_io) size3 = SIZE(gl,3)
!!$    CALL p_bcast (size3, p_io)
!!$    DO i=1,size3
!!$       NULLIFY (gl3)
!!$       IF (p_parallel_io) gl3 => gl(:,:,i,:)
!!$       CALL scatter_gp32 (gl3, lc(:,:,i,:), gl_dc, lg32=lg32)
!!$    END DO
    IF (p_parallel_io) size3 = SIZE(gl,3) !
    CALL p_bcast (size3, p_io)
    IF (size3 == 1) THEN
       IF (p_parallel_io) gl3 => gl(:,:,1,:) !
       CALL scatter_gp32 (gl3, lc(:,:,:,1), gl_dc, lg32=lg32)
    ELSE
       DO i=1,size3
          NULLIFY (gl3)
          IF (p_parallel_io) gl3 => gl(:,:,i,:)
          CALL scatter_gp32 (gl3, lc(:,:,i,:), gl_dc, lg32=lg32)
       END DO
    ENDIF
#endif
#else
    lc(:,:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3),1:SIZE(lc,4))
#endif
    !
  END SUBROUTINE scatter_gp432
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp32 (gl, lc, gl_dc, lg32)

#ifdef CESM1
    use phys_grid, only: scatter_field_to_chunk
#endif
    IMPLICIT NONE

    INTRINSIC :: SIZE
  !
  ! send global grid point field to pe-s (nlon,nlev,nlat) or (nlon,nlat,1)
  !

  REAL(dp), POINTER                 :: gl   (:,:,:) ! global field
  REAL(dp), TARGET ,INTENT(out)     :: lc   (:,:,:) ! local  field
  TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc        ! global decomposition
  LOGICAL, OPTIONAL    , INTENT(in) :: lg32         ! um_ak_20110603

  ! LOCAL
  REAL(dp),POINTER :: gl2(:,:)
  INTEGER(i4)      :: isize    ! size of 3rd dimension
  INTEGER(I4)      :: g_size(3)
  INTEGER          :: i
  LOGICAL          :: l3
  !
  ! call 2D scatter routine if 3rd dimension size is 1
  ! else call 3D scatter routine
  !
#ifdef COSMO
!!$    l3=.TRUE.
!!$    IF (PRESENT(lg32)) THEN
!!$       IF (lg32) l3=.FALSE.
!!$    ENDIF
!!$
!!$    IF (l3) THEN
       IF (p_pe == p_io) isize = (SIZE(gl,3))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,:,i)
          CALL scatter_gp2 (gl2, lc(:,:,i), gl_dc)
       ENDDO
!!$    ELSE
!!$       IF (p_pe == p_io) isize = (SIZE(gl,2))
!!$       CALL p_bcast (isize, p_io)
!!$       DO i=1,isize
!!$          NULLIFY(gl2)
!!$          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
!!$          gl2 => gl(:,i,:)
!!$          CALL scatter_gp2 (gl2, lc(:,i,:), gl_dc)
!!$       ENDDO
!!$    ENDIF
#elif defined(CESM1)
    l3=.FALSE.
    IF (PRESENT(lg32)) THEN
       IF (lg32) l3=.TRUE.
    ENDIF

    IF (l3) THEN
       IF (p_pe == p_io) isize = (SIZE(gl,3))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          IF (p_parallel_io) gl2 => gl(:,:,i)
          CALL scatter_gp2 (gl2, lc(:,:,i), gl_dc)
       ENDDO
    ELSE
!!$       IF (p_pe == p_io) isize = (SIZE(gl,2))
!!$       CALL p_bcast (isize, p_io)
!!$       IF (isize == 1) THEN
!!$          NULLIFY(gl2)
!!$          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
!!$          IF (p_parallel_io) gl2 => gl(:,1,:)
!!$          CALL scatter_gp2 (gl2, lc(:,1,:), gl_dc)
!!$       ELSE
          IF (p_pe == p_io) THEN
             g_size = (/ SIZE(gl,1), SIZE(gl,2), SIZE(gl,3) /)
          END IF
          CALL p_bcast (g_size, p_io)
          IF (p_pe /= p_io) THEN
             ALLOCATE(gl(g_size(1),g_size(2),g_size(3)))
          END IF
          call scatter_field_to_chunk(1, g_size(2), 1, g_size(1), gl, lc)
!!$       ENDIF
    ENDIF
#else
    lc(:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3))
#endif
    !
  END SUBROUTINE scatter_gp32
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp2 (gl, lc, gl_dc, ijsender)
  !
  ! send global 2D grid point field to local pe-s (nlon,nlat)
  !
#ifdef CESM1
    use phys_grid, only: scatter_field_to_chunk
#endif
    IMPLICIT NONE
    INTRINSIC :: SIZE

  REAL(dp)                 ,POINTER     :: gl   (:,:) ! global field
  REAL(dp)         ,TARGET ,INTENT(OUT) :: lc   (:,:) ! local  field
  TYPE(decomp),   OPTIONAL ,INTENT(IN)  :: gl_dc !global decomposition
  INTEGER,        OPTIONAL, INTENT(IN)  :: ijsender
    !
  INTEGER(I4) :: l_size(2), g_size(2)
  INTEGER(I4) :: ierror = 0
  INTEGER(I4) :: isender
  !
#ifdef COSMO
  IF (PRESENT(ijsender)) THEN
     isender = INT(ijsender,I4)
  ELSE
     isender = 0_I4
  END IF

  IF (p_pe == p_io) THEN
     g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
  END IF
  CALL p_bcast (g_size, p_io)
  l_size = (/ SIZE(lc,1), SIZE(lc,2) /)

  CALL distribute_field(gl,g_size(1),g_size(2),lc,l_size(1),l_size(2) &
       , isender, ierror )

  IF (ierror /= 0) THEN
     CALL p_abort('scatter_gp2','messy_main_mpi_bi')
     RETURN
  ENDIF
#elif defined(CESM1)
  IF (p_pe == p_io) THEN
     g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
  END IF
  CALL p_bcast (g_size, p_io)
  IF (p_pe /= p_io) THEN
    ALLOCATE(gl(g_size(1),g_size(2)))
  END IF
  call scatter_field_to_chunk(1,1,1, g_size(1),gl,lc)
#else
    lc(:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2))
!... COSMO ... elif
#endif
  END SUBROUTINE scatter_gp2
!------------------------------------------------------------------------------
#ifdef CESM1
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp1 (gl, lc, gl_dc, ijsender)
  !
  ! send global 2D grid point field to local pe-s (nlon,nlat)
  !
#ifdef CESM1
    use phys_grid, only: scatter_field_to_chunk
#endif
    IMPLICIT NONE
    INTRINSIC :: SIZE

  REAL(dp)                 ,POINTER     :: gl   (:) ! global field
  REAL(dp)         ,TARGET ,INTENT(OUT) :: lc   (:,:) ! local  field
  TYPE(decomp),   OPTIONAL ,INTENT(IN)  :: gl_dc !global decomposition
  INTEGER,        OPTIONAL, INTENT(IN)  :: ijsender
    !
  INTEGER(I4) :: l_size(2), g_size(1)
  INTEGER(I4) :: ierror = 0
  INTEGER(I4) :: isender
  !
  IF (p_pe == p_io) THEN
     g_size = (/ SIZE(gl,1) /)
  END IF
  CALL p_bcast (g_size, p_io)
  IF (p_pe /= p_io) THEN
    ALLOCATE(gl(g_size(1)))
  END IF
  call scatter_field_to_chunk(1,1,1, g_size(1),gl,lc)
  END SUBROUTINE scatter_gp1
!------------------------------------------------------------------------------
#endif

!==============================================================================
! um_ak_20121022+
!==============================================================================
  SUBROUTINE scatter_gp432_sp (gl, lc, gl_dc, lg32)
  !
  ! scatter global grid point field from pe-s (nlon,nlev,ntrac,nlat)
  !                                       or (nlon,nlev,nlat ,1   )
  !                                       or (nlon,nlat,1    ,1   )
  !
    IMPLICIT NONE

    INTRINSIC :: SIZE

    REAL, POINTER                 :: gl   (:,:,:,:) ! global field
    REAL, TARGET ,INTENT(out)     :: lc   (:,:,:,:) ! local  field
    TYPE(DECOMP),OPTIONAL,INTENT(in)  :: gl_dc          ! global decomposition
    LOGICAL, OPTIONAL,    INTENT(in)  :: lg32           ! um_ak_20110603
    ! LOCAL
    INTEGER(i4)      :: size4 ! size of 4th dimension
    INTEGER          :: i
    REAL,POINTER :: gl3(:,:,:)
    !
    ! call 3D scatter routine if 4th dimension size is 1
    ! else loop over 3th index
    !
#ifdef COSMO
    IF (p_parallel_io) size4 = (SIZE(gl,4))
    CALL p_bcast (size4, p_io)
    NULLIFY(gl3)
    DO i=1,size4
!       IF (p_pe == p_io) gl3 => gl(:,:,:,i)
       gl3 => gl(:,:,:,i)
       CALL scatter_gp32_sp (gl3, lc(:,:,:,i), gl_dc, lg32=lg32)! um_ak_20110603 lg32 added)
    END DO
#else
    lc(:,:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3),1:SIZE(lc,4))
#endif
    !
  END SUBROUTINE scatter_gp432_sp
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp32_sp (gl, lc, gl_dc, lg32)

    IMPLICIT NONE

    INTRINSIC :: SIZE
  !
  ! send global grid point field to pe-s (nlon,nlev,nlat) or (nlon,nlat,1)
  !

  REAL, POINTER                 :: gl   (:,:,:) ! global field
  REAL, TARGET ,INTENT(out)     :: lc   (:,:,:) ! local  field
  TYPE(decomp),OPTIONAL,INTENT(in)  :: gl_dc        ! global decomposition
  LOGICAL, OPTIONAL    , INTENT(in) :: lg32         ! um_ak_20110603

  ! LOCAL
  REAL,POINTER :: gl2(:,:)
  INTEGER(i4)      :: isize    ! size of 3rd dimension
  INTEGER          :: i
  LOGICAL          :: l3
  !
  ! call 2D scatter routine if 3rd dimension size is 1
  ! else call 3D scatter routine
  !
#ifdef COSMO
    l3=.TRUE.
    IF (PRESENT(lg32)) THEN
       IF (lg32) l3=.FALSE.
    ENDIF

    IF (l3) THEN
       IF (p_pe == p_io) isize = (SIZE(gl,3))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,:,i)
          CALL scatter_gp2_sp (gl2, lc(:,:,i), gl_dc)
       ENDDO
    ELSE
       IF (p_pe == p_io) isize = (SIZE(gl,2))
       CALL p_bcast (isize, p_io)
       DO i=1,isize
          NULLIFY(gl2)
          !       IF (p_pe == p_io) gl2 => gl(:,:,i)
          gl2 => gl(:,i,:)
          CALL scatter_gp2_sp (gl2, lc(:,i,:), gl_dc)
       ENDDO
    ENDIF
#else
    lc(:,:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2),1:SIZE(lc,3))
#endif
    !
  END SUBROUTINE scatter_gp32_sp
!------------------------------------------------------------------------------
  SUBROUTINE scatter_gp2_sp (gl, lc, gl_dc, ijsender)
  !
  ! send global 2D grid point field to local pe-s (nlon,nlat)
  !
    IMPLICIT NONE
    INTRINSIC :: SIZE

  REAL                 ,POINTER     :: gl   (:,:) ! global field
  REAL         ,TARGET ,INTENT(OUT) :: lc   (:,:) ! local  field
  TYPE(decomp),   OPTIONAL ,INTENT(IN)  :: gl_dc !global decomposition
  INTEGER,        OPTIONAL, INTENT(IN)  :: ijsender
    !
  INTEGER(I4) :: l_size(2), g_size(2)
  INTEGER(I4) :: ierror = 0
  INTEGER(I4) :: isender
  REAL(dp), ALLOCATABLE :: gl_dp(:,:)
  REAL(dp), ALLOCATABLE :: lc_dp(:,:)
  !
#ifdef COSMO
  IF (PRESENT(ijsender)) THEN
     isender = INT(ijsender,I4)
  ELSE
     isender = 0_I4
  END IF

  IF (p_pe == p_io) THEN
     g_size = (/ SIZE(gl,1), SIZE(gl,2) /)
  END IF
  CALL p_bcast (g_size, p_io)
  l_size = (/ SIZE(lc,1), SIZE(lc,2) /)
  ALLOCATE(gl_dp(SIZE(gl,1),SIZE(gl,2)))
  ALLOCATE(lc_dp(SIZE(lc,1),SIZE(lc,2)))
  gl_dp = REAL(gl,dp)

  CALL distribute_field(gl_dp,g_size(1),g_size(2),lc_dp,l_size(1),l_size(2) &
       , isender, ierror )

  IF (ierror /= 0) THEN
     CALL p_abort('scatter_gp2','messy_main_mpi_bi')
     RETURN
  ENDIF
  lc = REAL(lc_dp)
  DEALLOCATE(lc_dp,gl_dp)
#else
    lc(:,:) = gl(1:SIZE(lc,1),1:SIZE(lc,2))
!... COSMO
#endif
  END SUBROUTINE scatter_gp2_sp
!------------------------------------------------------------------------------

!!$!==============================================================================
!!$  SUBROUTINE reorder (y,x)
!!$
!!$    IMPLICIT NONE
!!$    INTRINSIC :: RESHAPE, SIZE
!!$
!!$    REAL(dp) ,INTENT(out) :: y (:,:)
!!$    REAL(dp) ,INTENT(in)  :: x (:,:)
!!$
!!$#if (defined CRAY) || (defined sun) || (defined NAG) || (defined __SX__) || defined(__PGI)
!!$    CALL util_reshape(y, x, SIZE(y,1)*SIZE(y,2), SIZE(x,1)*SIZE(x,2))
!!$#else
!!$    y = RESHAPE (x,(/SIZE(y,1),SIZE(y,2)/),(/0._dp/))
!!$#endif
!!$
!!$  END SUBROUTINE reorder
!!$!------------------------------------------------------------------------------

#ifdef I2CINC
! ----------------------------------------------------------------------
  SUBROUTINE switch_par_utilities(flag)

    ! SWITCH BETWEEN COSMO AND INT2COSMO DIMENSIONS OF FIELDS
    ! FOR parallel_utilities

    ! INT2COSMO
    USE data_grid_lm,         ONLY: istartpar_i2c   => istartpar    &
                                  , iendpar_i2c     => iendpar      &
                                  , jstartpar_i2c   => jstartpar    &
                                  , jendpar_i2c     => jendpar      &
                                  , ie2lm_tot, je2lm_tot, kelm_tot  &
                                  , kelm, ie2lm_max, je2lm_max      &
                                  , ie2lm, je2lm
    USE data_grid_in,         ONLY: istartpar_in, iendpar_in        &
                                  , jstartpar_in, jendpar_in        &
                                  , ie_in,     je_in,     ke_in     &
                                  , ie_in_tot, je_in_tot, ke_in_tot &
                                  , ie_in_max, je_in_max
    USE data_int2lm_parallel, ONLY: nboundlines_in                  &
                                  , ie_in_tot_red, je_in_tot_red    &
                                  , isubpos_in_red                  &
                                  , isubpos_in, isubpos_coarse

    USE data_int2lm_parallel, ONLY: nboundlines_i2c => nboundlines &
                                  , isubpos_i2c     => isubpos

    ! COSMO
    USE data_modelconfig, ONLY: istartpar_c4   => istartpar &
                              , iendpar_c4     => iendpar   &
                              , jstartpar_c4   => jstartpar &
                              , jendpar_c4     => jendpar   &
                              , ke, ke_tot, ie, je, ie_tot  &
                              , je_tot, ie_max, je_max
    ! COSMO/INT2COSMO
    USE data_parallel,    ONLY: nboundlines_c4 => nboundlines  &
                              , isubpos_c4     => isubpos      &
                              , icomm_cart, imp_integers       &
                              , nprocx, nprocy, nproc, nprocio &
#ifdef COSMOv509
                              , lcompute_pe                    &
#endif
                              , my_cart_id

    USE parallel_utilities, ONLY: init_par_utilities
    USE messy_main_constants_mem, ONLY: iouerr

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: flag

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='switch_par_utilities'

    SELECT CASE (flag)
    CASE(1)
       ! SWITCH FROM COSMO TO INT2COSMO
       CALL init_par_utilities (ie2lm, je2lm, kelm     &
            , ie2lm_tot, je2lm_tot, kelm_tot           &
            , ie2lm_max, je2lm_max                     &
            , istartpar_i2c, iendpar_i2c               &
            , jstartpar_i2c, jendpar_i2c               &
            , nproc, nprocx, nprocy, nprocio           &
            , isubpos_i2c, nboundlines_i2c, icomm_cart &
            , my_cart_id, imp_reals, imp_integers      &
#ifdef COSMOv509
            , lcompute_pe                              &
#endif
            )
    CASE(2)
       ! SWITCH FROM INT2COSMO TO COSMO
       CALL init_par_utilities (ie, je, ke, ie_tot, je_tot, ke_tot     &
            , ie_max, je_max, istartpar_c4, iendpar_c4                 &
            , jstartpar_c4, jendpar_c4, nproc, nprocx, nprocy, nprocio &
            , isubpos_c4, nboundlines_c4, icomm_cart, my_cart_id       &
            , imp_reals, imp_integers &
#ifdef COSMOv509
            , lcompute_pe                                              &
#endif
            )
    CASE(3)
       ! SWITCH FROM INT2COSMO TO INT2COSMO-IN
       CALL init_par_utilities (ie_in, je_in, ke_in                    &
            , ie_in_tot_red, je_in_tot_red, ke_in_tot, ie_in_max, je_in_max    &
            , istartpar_in, iendpar_in, jstartpar_in, jendpar_in       &
            , nproc, nprocx, nprocy, nprocio, isubpos_in_red           &
            , nboundlines_in, icomm_cart, my_cart_id                   &
            , imp_reals, imp_integers &
#ifdef COSMOv509
            , lcompute_pe                                              &
#endif
            )
    CASE(4)
       ! SWITCH FROM INT2COSMO TO COSMO-OUT
       CALL init_par_utilities (ie_in, je_in, ke_in                    &
            , ie_in_tot_red , je_in_tot_red , ke_in_tot, ie_in_max, je_in_max&
            , istartpar_in, iendpar_in, jstartpar_in, jendpar_in       &
            , nproc, nprocx, nprocy, nprocio, isubpos_in_red           &
            , nboundlines_in, icomm_cart, my_cart_id                   &
            , imp_reals, imp_integers &
#ifdef COSMOv509
            , lcompute_pe                                              &
#endif
            )

    CASE(5) ! 5 only required for gathering of fland_in_tot
       ! SWITCH FROM INT2COSMO TO INT2COSMO-IN
       CALL init_par_utilities (ie_in, je_in, ke_in                    &
            , ie_in_tot, je_in_tot, ke_in_tot, ie_in_max, je_in_max    &
            , istartpar_in, iendpar_in, jstartpar_in, jendpar_in       &
            , nproc, nprocx, nprocy, nprocio, isubpos_in               &
            , nboundlines_in, icomm_cart, my_cart_id                   &
            , imp_reals, imp_integers &
#ifdef COSMOv509
            , lcompute_pe                                              &
#endif
            )

    CASE DEFAULT
       ! NO FURTHER SWITCHING OPTIONS
       write(iouerr,*) ' THIS SWITCHING OPTION IS NOT DEFINED! '
    END SELECT

  END SUBROUTINE switch_par_utilities

! -------------------------------------------------------------------------
!... I2CINC
#endif

#if !defined(BLANK) && !defined(CESM1) && ! defined(VERTICO)
! -------------------------------------------------------------------------
  SUBROUTINE global_1d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:)
    INTEGER, INTENT(OUT) :: outfield(:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF
    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_1d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_2d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:,:)
    INTEGER, INTENT(OUT) :: outfield(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_2d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_3d_intfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)  :: dims
    INTEGER, INTENT(IN)  :: infield(:,:,:)
    INTEGER, INTENT(OUT) :: outfield(:,:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_integers, noper, icomm_cart, ierror)

  END SUBROUTINE global_3d_intfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_1d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:)
    REAL(dp), INTENT(OUT) :: outfield(:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)

  END SUBROUTINE global_1d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_2d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:,:)
    REAL(dp), INTENT(OUT) :: outfield(:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)

  END SUBROUTINE global_2d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE global_3d_realfield(infield, outfield, dims, tag, yerrmsg, ierror)

    IMPLICIT NONE

    ! parameter
    INTEGER, INTENT(IN)   :: dims
    REAL(dp), INTENT(IN)  :: infield(:,:,:)
    REAL(dp), INTENT(OUT) :: outfield(:,:,:)
    CHARACTER(LEN=3), INTENT(IN) :: tag

    INTEGER, INTENT(OUT)          :: ierror
    CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg
    ! LOCAL
    INTEGER :: noper, dim

    ierror = 0
    yerrmsg = '                                        '

    IF (tag == 'SUM') THEN
       noper = MPI_SUM
    ELSEIF (tag == 'MAX') THEN
       noper = MPI_MAX
    ELSEIF (tag == 'MIN') THEN
       noper = MPI_MIN
    ELSE
       ierror  = 1
       yerrmsg = 'no valid operation type in global_fields'
       RETURN
    ENDIF

    CALL MPI_ALLREDUCE                                       &
      (infield, outfield, dims, imp_reals, noper, icomm_cart, ierror)


  END SUBROUTINE global_3d_realfield
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE exchange_boundaries(kdims, field, ierror, yerrmsg )

    USE data_modelconfig, ONLY: jstartpar, jendpar, ie, je
    USE data_runcontrol,  ONLY: lperi_x, lperi_y, l2dim, nnow

    IMPLICIT NONE

    INTEGER, INTENT(IN)                       :: kdims(24)
    REAL(dp), INTENT(INOUT), DIMENSION(:,:,:) :: field
    INTEGER, INTENT(OUT)                      :: ierror
    CHARACTER(LEN=*), INTENT(OUT)             :: yerrmsg

    itag = itag+1

    CALL exchg_boundaries                                            &
             (nnow+39, sendbuf, isendbuflen, imp_reals, icomm_cart   &
             , num_compute, ie, je, kdims, jstartpar, jendpar        &
             , nboundlines, nboundlines, my_cart_neigh               &
             , lperi_x, lperi_y, l2dim                               &
             , itag, ldatatypes, ncomm_type, ierror, yerrmsg         &
             , field(:,:,:))


  END SUBROUTINE exchange_boundaries
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------

  SUBROUTINE exchg_trac_boundaries2 ( idim, jdim, kdim, trac_dim,  &
                 jstartpar, jendpar, nlines, nboundlines,         &
                 neighbors, ntag, ierror, yerrmsg, trac_field)

    !----------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine performs the boundary exchange for the full tracer field.
    !
    !----------------------------------------------------------------------------

#ifndef NOMPI
#ifndef __PGI
  USE mpi
#endif
#endif

    IMPLICIT NONE

#ifndef NOMPI
#ifdef __PGI
  INCLUDE 'mpif.h'
#endif
#endif

    ! Subroutine arguments
    INTEGER, INTENT(IN)         ::    &
         idim, jdim,         & ! horizontal dimensions of the fields
         kdim,               & ! vertical dimensions  of the tracer field
         trac_dim,           & ! number dimension of the tracer field
         jstartpar,          & ! start index in j-direction
         jendpar,            & ! end index in j-direction
         nlines,             & ! number of lines that have to be exchanged
                               ! (<= nboundlines)
         nboundlines,        & ! number of overlapping boundary lines
         neighbors(4),       & ! process-ids of the neighbors in the grid
         ntag                  ! tag of the message

    INTEGER, INTENT (OUT) :: ierror    ! error status variable

    CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg ! for MPI error message

    REAL (dp), DIMENSION(:,:,:,:), POINTER :: trac_field

    ! LOCAL
    INTEGER :: izmplcode            ! for MPI error code
    INTEGER :: i,j,k, jt, ind, datasize
    INTEGER :: ireq(4)
    INTEGER :: status (MPI_STATUS_SIZE)
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: sbuffer
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: rbuffer

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '

  ! check whether nlines <= nboundlines
  IF (nlines > nboundlines) THEN
    ierror  = 9011
    yerrmsg     = ' *** nlines > nboundlines *** '
    RETURN
  ENDIF

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
  ! Horizontal exchange

  datasize = nlines * (jendpar - jstartpar + 1) * kdim * trac_dim

  ALLOCATE (sbuffer(datasize))
  ALLOCATE (rbuffer(datasize))

  ! left neighbor is present
  IF (neighbors(1) /= -1) THEN
     ind = 1
     DO jt = 1, trac_dim
        DO i = nboundlines+1, nboundlines + nlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO

     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(1), ntag, icomm_cart, ireq(1), ierror)
     IF (ierror /= 0) RETURN

  END IF

  ! right neighbor exists
  IF (neighbors(3) /= -1) THEN
     ! 2.) send buffer
     ind = 1
     DO jt = 1, trac_dim
        DO i = idim - nboundlines - nlines + 1, idim-nboundlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(3), ntag+1, icomm_cart, ireq(3), ierror)
     IF (ierror /= 0) RETURN

     ! 1.) receive buffer of right neighbor
     rbuffer(:) = 0._dp

      CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(3), ntag, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = idim - nboundlines + 1,idim - nboundlines + nlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO

  END IF

  ! left neighbors exists recv buffer
  IF (neighbors(1) /= -1) THEN

     rbuffer(:) = 0._dp

     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(1), ntag+1, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = nboundlines + 1 - nlines, nboundlines
           DO j = jstartpar, jendpar
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  END IF

  IF (neighbors(1) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(1), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 1'
        RETURN
     ENDIF
  ENDIF
  IF (neighbors(3) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(3), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 3'
        RETURN
     ENDIF
  ENDIF



  DEALLOCATE (rbuffer, sbuffer)

  ! *************************************************
  ! EXCHANGE with northern and southern neighbors
  ! *************************************************

  ! DIMENSION buffer
  datasize =  nlines * idim * kdim * trac_dim
  ALLOCATE (sbuffer(datasize))
  ALLOCATE (rbuffer(datasize))

  ! upper neighbor is present
  IF (neighbors(2) /= -1) THEN
      ! send buffer
     ind = 1
     DO jt = 1, trac_dim
        DO i = 1 , idim
           DO j = jdim - nboundlines - nlines + 1, jdim - nboundlines
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(2), ntag+2, icomm_cart, ireq(2), ierror)
     IF (ierror /= 0) RETURN

  ENDIF

  ! lower neighbor is present
  IF (neighbors(4) /= -1) THEN
      ! send buffer to lower neighbor
     ind = 1
     DO jt = 1, trac_dim
        DO i = 1 , idim
           DO j = nboundlines + 1, nboundlines + nlines
              DO k= 1, kdim
                 sbuffer(ind) = trac_field(i,j,jt,k)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
     CALL MPI_ISEND (sbuffer, datasize,  imp_reals, &
          neighbors(4), ntag+3, icomm_cart, ireq(4), ierror)
     IF (ierror /= 0) RETURN

     rbuffer(:) = 0._dp

     ! receive buffer of lower neighbor
     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(4), ntag+2, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = 1, idim
           DO j = nboundlines-nlines+1, nboundlines
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  ENDIF

  ! upper neighbor is present
  IF (neighbors(2) /= -1) THEN
     ! receive buffer of upper neighbor
     rbuffer(:) = 0._dp

     CALL MPI_RECV(rbuffer, datasize,  imp_reals &
          , neighbors(2), ntag+3, icomm_cart, MPI_STATUS_IGNORE, ierror)
     IF (ierror /= 0) RETURN

     ind = 1
     DO jt = 1, trac_dim
        DO i = 1, idim
!           DO j = jdim -nboundlines - nlines +1 , jdim - nboundlines
           DO j = jdim -nboundlines + 1 , jdim - nboundlines + nlines
              DO k= 1, kdim
                 trac_field(i,j,jt,k) = rbuffer(ind)
                 ind =ind +1
              END DO
           END DO
        END DO
     END DO
  ENDIF

  IF (neighbors(2) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(2), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 2'
        RETURN
     ENDIF
  ENDIF
  IF (neighbors(4) /= -1) THEN
     ! wait for the completion of the last send to neighbors(1)
     CALL MPI_WAIT (ireq(4), status, izmplcode)
     IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_WAIT 4'
        RETURN
     ENDIF
  ENDIF
   DEALLOCATE (rbuffer, sbuffer)


END SUBROUTINE exchg_trac_boundaries2

!==============================================================================
!==============================================================================
!+ This subroutine performs the data exchange between boundaries
!------------------------------------------------------------------------------

SUBROUTINE exchg_trac_boundaries                                             &
               ( icase, imp_type, icomm, idim, jdim,   &
                 kdim, trac_dim, jstartpar, jendpar, nlines, nboundlines,    &
                 neighbors, ntag, lmpi_types, ntype, ierror, yerrmsg,        &
                 trac_field)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine performs the boundary exchange of up to 20 variables. Only
!   one variable has to occur, the others are optional.
!
! Method:
!   At the moment there are 3 different MPI-communications implemented:
!     1) immediate send, blocking receive and wait
!     2) immediate receive, blocking send and wait
!     3) MPI_SendRecv
!   Also there is the choice of an explicit buffering (putbuf, getbuf) or
!   implicit buffering (MPI-Datatypes) of the data to be send.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    icase,              & ! tag for exchange scenario
!    isendbuflen,        & ! length of sendbuffer
    imp_type,           & ! determines the REAL type used
    icomm,              & ! communicator for virtual cartesian topology
    idim, jdim,         & ! horizontal dimensions of the fields
    kdim,               & ! array for the vertical dimensions of var1..var20
    trac_dim,           & ! number of tracers
    jstartpar,          & ! start index in j-direction
    jendpar,            & ! end index in j-direction
    nlines,             & ! number of lines that have to be exchanged
                          ! (<= nboundlines)
    nboundlines,        & ! number of overlapping boundary lines
    neighbors(4),       & ! process-ids of the neighbors in the grid
    ntag,               & ! tag of the message
    ntype                 ! indicates how the communication should be done

  LOGICAL, INTENT(IN) :: lmpi_types ! whether implicit (with MPI-Datatypes)
                                    ! or explicit
                                    ! (putbuf, getbuf) buffering of data is used

  INTEGER,         INTENT (OUT) :: ierror      ! error status variable

  CHARACTER(LEN=*), INTENT(OUT) :: yerrmsg     ! for MPI error message

  REAL(dp),       INTENT (INOUT)      ::    &
!    sendbuf (isendbuflen, 8),& ! send buffer
    trac_field(idim,jdim,trac_dim,kdim)


  ! LOCAL
  INTEGER ::       &
    ! the following numbers are for filling/emptying the buffers for each
    ! neighbor
    izlo_lr, izup_lr, jzlo_lr, jzup_lr,     & ! left , receive
    izlo_rr, izup_rr, jzlo_rr, jzup_rr,     & ! right, receive
    izlo_ur, izup_ur, jzlo_ur, jzup_ur,     & ! upper, receive
    izlo_dr, izup_dr, jzlo_dr, jzup_dr,     & ! down , receive
    izlo_ls, izup_ls, jzlo_ls, jzup_ls,     & ! left , send
    izlo_rs, izup_rs, jzlo_rs, jzup_rs,     & ! right, send
    izlo_us, izup_us, jzlo_us, jzup_us,     & ! upper, send
    izlo_ds, izup_ds, jzlo_ds, jzup_ds,     & ! down , send
    nzcount_ls, nzcount_rs,     & ! counting the values
    nzcount_us, nzcount_ds,     & ! counting the values
    nzcount_lr, nzcount_rr,     & ! counting the values
    nzcount_ur, nzcount_dr,     & ! counting the values
    nzrequest(MPI_STATUS_SIZE), & ! for MPI-receive
    nzstatus (MPI_STATUS_SIZE), & ! for MPI-WAIT
    ncount, type_handle,        & ! return values from setup_data_type
    MPI_neighbors(4), i,        & ! same as neighbors, if neighbor exists
    ilocalreq(4)                  ! the local requests for the ISEND and IRECV

  REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: sendbuf
  INTEGER  :: isendbuflen
  INTEGER  :: izmplcode                   ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  ierror     = 0
  izmplcode  = 0
  yerrmsg    = '    '

  isendbuflen = MAX(idim,jdim)* kdim*trac_dim*nlines
  ALLOCATE (sendbuf(isendbuflen, 8))
 ! check whether nlines <= nboundlines
  IF (nlines > nboundlines) THEN
    ierror  = 9011
    yerrmsg     = ' *** nlines > nboundlines *** '
    RETURN
  ENDIF

  ! Determine the start- and end-indices for routines putbuf, getbuf
  izlo_ls = nboundlines + 1
  izup_ls = nboundlines + nlines
  jzlo_ls = jstartpar
  jzup_ls = jendpar

  izlo_lr = nboundlines + 1 - nlines
  izup_lr = nboundlines
  jzlo_lr = jstartpar
  jzup_lr = jendpar

  izlo_us = 1
  izup_us = idim
  jzlo_us = jdim - nboundlines - nlines + 1
  jzup_us = jdim - nboundlines

  izlo_ur = 1
  izup_ur = idim
  jzlo_ur = jdim - nboundlines + 1
  jzup_ur = jdim - nboundlines + nlines

  izlo_rs = idim - nboundlines - nlines + 1
  izup_rs = idim - nboundlines
  jzlo_rs = jstartpar
  jzup_rs = jendpar

  izlo_rr = idim - nboundlines + 1
  izup_rr = idim - nboundlines + nlines
  jzlo_rr = jstartpar
  jzup_rr = jendpar

  izlo_ds = 1
  izup_ds = idim
  jzlo_ds = nboundlines + 1
  jzup_ds = nboundlines + nlines

  izlo_dr = 1
  izup_dr = idim
  jzlo_dr = nboundlines + 1 - nlines
  jzup_dr = nboundlines

  nzcount_lr = 0
  nzcount_rr = 0
  nzcount_ur = 0
  nzcount_dr = 0
  nzcount_ls = 0
  nzcount_rs = 0
  nzcount_us = 0
  nzcount_ds = 0

  ! Fix list of neighbors (use MPI_PROC_NULL rather than -1 to indicate
  ! missing neighbor).
  DO i= 1, 4
    IF ( neighbors(i) /= -1 ) THEN
      MPI_neighbors(i) = neighbors(i)
    ELSE
      MPI_neighbors(i) = MPI_PROC_NULL
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- Section 2: Determine the necessary datatypes
!------------------------------------------------------------------------------

  IF (lmpi_types) THEN
    ! Exchange with left and right neighbor
    IF ( iexchg_MPI_types(2*icase-1) == MPI_DATATYPE_NULL ) THEN
      IF ( MPI_neighbors(1) /= MPI_PROC_NULL ) THEN
        CALL setup_data_type( trac_field, trac_dim,                        &
           idim, jdim, kdim, izlo_ls, izup_ls, jzlo_ls, jzup_ls,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase-1) = type_handle
        iexchg_counts   (2*icase-1) = ncount
      ELSE
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_rr, izup_rr, jzlo_rr, jzup_rr,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase-1) = type_handle
        iexchg_counts   (2*icase-1) = ncount
      ENDIF
    ENDIF

    ! Exchange with upper and lower neighbor
    IF ( iexchg_MPI_types(2*icase) == MPI_DATATYPE_NULL ) THEN
      IF ( MPI_neighbors(2) /= MPI_PROC_NULL ) THEN
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_us, izup_us, jzlo_us, jzup_us,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase) = type_handle
        iexchg_counts   (2*icase) = ncount
      ELSE
        CALL setup_data_type(trac_field, trac_dim,                         &
           idim, jdim, kdim, izlo_dr, izup_dr, jzlo_dr, jzup_dr,           &
           ierror, yerrmsg, imp_type, ncount, type_handle )
        iexchg_MPI_types(2*icase) = type_handle
        iexchg_counts   (2*icase) = ncount
      ENDIF
    ENDIF
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Exchange with immediate Send and blocking Recv
!------------------------------------------------------------------------------

  IF   (ntype == 1) THEN

    IF (lmpi_types) THEN

      !------------------------------------------------------------------------
      !- Section 3.1: exchange with left and right neighbors using datatypes
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        CALL MPI_ISEND ( trac_field(izlo_ls,jzlo_ls,1,1) &
                       , iexchg_counts(2*icase-1), &
                         iexchg_MPI_types(2*icase-1), MPI_neighbors(1),      &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present
        ! send the data
        CALL MPI_ISEND ( trac_field(izlo_rs,jzlo_rs,1,1) &
                       , iexchg_counts(2*icase-1), &
                         iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present
        ! receive the data
        CALL MPI_RECV (trac_field(izlo_rr,jzlo_rr,1,1) &
                     , iexchg_counts(2*icase-1),   &
                       iexchg_MPI_types(2*icase-1), MPI_neighbors(3),        &
                       ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data
        CALL MPI_RECV ( trac_field(izlo_lr,jzlo_lr,1,1) &
                      , iexchg_counts(2*icase-1),  &
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),       &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.2: exchange with upper and lower neighbors using datatypes
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! send the data
        CALL MPI_ISEND (trac_field(izlo_us,jzlo_us,1,1) &
                      , iexchg_counts(2*icase),    &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! send the data
        CALL MPI_ISEND (trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),    &
                        iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        CALL MPI_RECV (trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),    &
                       iexchg_MPI_types(2*icase), MPI_neighbors(4),         &
                       ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        CALL MPI_RECV (trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),     &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),         &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ELSE

      !------------------------------------------------------------------------
      !- Section 3.3: exchange with left and right neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present

        ! determine start- and end-indices for routine putbuf
        nzcount_ls = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),   &
                         ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! right neighbor is present

        ! determine start- and end-indices for routine putbuf
        nzcount_rs = 0
        CALL putbuf ( trac_field, trac_dim,&
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),   &
                         ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! receive the data
        CALL MPI_RECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! left neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! wait for the completion of the last send to neighbors(1)
        ! to safely reuse sendbuf(1,1)
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! wait for the completion of the last send to neighbors(3)
        ! to safely reuse sendbuf(1,3)
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
          IF (izmplcode /= 0) THEN
            ierror  = izmplcode
            yerrmsg = 'MPI_WAIT'
            RETURN
          ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 3.4: exchange with lower and upper neigh. using explict buff.
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),   &
                         ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present

        nzcount_ds = 0
        CALL putbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        ! send the data
        CALL MPI_ISEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),   &
                         ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_ISEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! lower neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

      IF (neighbors(2) /= -1) THEN
        ! upper neighbor is present
        ! receive the data

        CALL MPI_RECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),   &
                        ntag, icomm, nzrequest, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_RECV'
          RETURN
        ENDIF

        CALL getbuf (  trac_field, trac_dim,                                 &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

      IF (neighbors(2) /= -1) THEN
        ! wait for the completion of the last send to neighbors(2)
        ! to safely reuse sendbuf(1,2)
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        ! wait for the completion of the last send to neighbors(4)
        ! to safely reuse sendbuf(1,4)
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ENDIF

!------------------------------------------------------------------------------
!- Section 4: Exchange with immediate Recv and blocking Send
!------------------------------------------------------------------------------

  ELSEIF (ntype == 2) THEN

    IF (lmpi_types) THEN

      !------------------------------------------------------------------------
      !- Section 4.1: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_rr,jzlo_rr,1,1) &
                     , iexchg_counts(2*icase-1), &
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(3),      &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_lr,jzlo_lr,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        ! now send the data to the left neighbor
        CALL MPI_SEND ( trac_field(izlo_ls,jzlo_ls,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        ! now send the data to the right neighbor
        CALL MPI_SEND ( trac_field(izlo_rs,jzlo_rs,1,1), iexchg_counts(2*icase-1),&
                        iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.2: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),  &
                        iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV (trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),  &
                        iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_SEND (trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),  &
                       iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                       ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_SEND (trac_field(izlo_us,jzlo_us,1,1), iexchg_counts(2*icase),  &
                       iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                       ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF
      ENDIF

    ELSE

      !------------------------------------------------------------------------
      !- Section 4.3: exchange with left and right neighbors
      !------------------------------------------------------------------------

      IF (neighbors(3) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,7), isendbuflen, imp_type, neighbors(3),  &
                        ntag, icomm, ilocalreq(3), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,5), isendbuflen, imp_type, neighbors(1),  &
                        ntag, icomm, ilocalreq(1), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(1) /= -1) THEN
        nzcount_ls = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )

        CALL MPI_SEND ( sendbuf(1,1), nzcount_ls, imp_type, neighbors(1),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        nzcount_rs = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )

        CALL MPI_SEND ( sendbuf(1,3), nzcount_rs, imp_type, neighbors(3),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(3) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(3), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      IF (neighbors(1) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(1), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !------------------------------------------------------------------------
      !- Section 4.4: exchange with upper and lower neighbors
      !------------------------------------------------------------------------

      IF (neighbors(2) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,6), isendbuflen, imp_type, neighbors(2),  &
                        ntag, icomm, ilocalreq(2), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        CALL MPI_IRECV ( sendbuf(1,8), isendbuflen, imp_type, neighbors(4),  &
                        ntag, icomm, ilocalreq(4), izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_IRECV'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(4) /= -1) THEN
        nzcount_ds = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )

        CALL MPI_SEND ( sendbuf(1,4), nzcount_ds, imp_type, neighbors(4),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        nzcount_us = 0
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )

        CALL MPI_SEND ( sendbuf(1,2), nzcount_us, imp_type, neighbors(2),    &
                        ntag, icomm, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_SEND'
          RETURN
        ENDIF
      ENDIF

      IF (neighbors(2) /= -1) THEN
        CALL MPI_WAIT (ilocalreq(2), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF


      IF (neighbors(4) /= -1) THEN
        ! Now wait until the data have arrived
        CALL MPI_WAIT (ilocalreq(4), nzstatus, izmplcode)
        IF (izmplcode /= 0) THEN
          ierror  = izmplcode
          yerrmsg = 'MPI_WAIT'
          RETURN
        ENDIF

        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

    ENDIF

!------------------------------------------------------------------------------
!- Section 5: Exchange with SendRecv
!------------------------------------------------------------------------------

  ELSEIF (ntype == 3) THEN

    IF (lmpi_types) THEN

      !--------------------------------------------------------------------------
      !- Section 5.1: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_ls,jzlo_ls,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag,                                              &
                          trac_field(izlo_rr,jzlo_rr,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                          ntag,                                              &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-1'
        RETURN
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.2: Receive data from the left and send to the right neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_rs,jzlo_rs,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(3),     &
                          ntag+1,                                            &
                          trac_field(izlo_lr,jzlo_lr,1,1), iexchg_counts(2*icase-1),&
                          iexchg_MPI_types(2*icase-1), MPI_neighbors(1),     &
                          ntag+1,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-2'
        RETURN
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.3: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_us,jzlo_us,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                          ntag+2,                                            &
                          trac_field(izlo_dr,jzlo_dr,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                          ntag+2,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-3'
        RETURN
      ENDIF


      !--------------------------------------------------------------------------
      !- Section 5.4: Receive data from the upper and send to the lower neighbor
      !--------------------------------------------------------------------------

      CALL MPI_SENDRECV ( trac_field(izlo_ds,jzlo_ds,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(4),       &
                          ntag+3,                                            &
                          trac_field(izlo_ur,jzlo_ur,1,1), iexchg_counts(2*icase),  &
                          iexchg_MPI_types(2*icase), MPI_neighbors(2),       &
                          ntag+3,                                            &
                          icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV-4'
        RETURN
      ENDIF

    ELSE

      !--------------------------------------------------------------------------
      !- Section 5.5: Send data to the left and receive from the right neighbor
      !--------------------------------------------------------------------------

      nzcount_ls = 0
      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ls, izup_ls, jzlo_ls, jzup_ls, nzcount_ls, 1 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,1), nzcount_ls,  imp_type, MPI_neighbors(1), ntag,    &
             sendbuf(1,7), isendbuflen, imp_type, MPI_neighbors(3), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rr, izup_rr, jzlo_rr, jzup_rr, nzcount_rr, 7 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.6: Send data to the right and receive from the left neighbor
      !--------------------------------------------------------------------------

      nzcount_rs = 0
      IF (MPI_neighbors(3) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_rs, izup_rs, jzlo_rs, jzup_rs, nzcount_rs, 3 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,3), nzcount_rs,  imp_type, MPI_neighbors(3), ntag,    &
             sendbuf(1,5), isendbuflen, imp_type, MPI_neighbors(1), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(1) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_lr, izup_lr, jzlo_lr, jzup_lr, nzcount_lr, 5 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.7: Send data to the upper and receive from the lower neighbor
      !--------------------------------------------------------------------------

      nzcount_us = 0
      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_us, izup_us, jzlo_us, jzup_us, nzcount_us, 2 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,2), nzcount_us,  imp_type, MPI_neighbors(2), ntag,    &
             sendbuf(1,8), isendbuflen, imp_type, MPI_neighbors(4), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_dr, izup_dr, jzlo_dr, jzup_dr, nzcount_dr, 8 )
      ENDIF

      !--------------------------------------------------------------------------
      !- Section 5.8: Send data to the lower and receive from the upper neighbor
      !--------------------------------------------------------------------------

      nzcount_ds = 0
      IF (MPI_neighbors(4) /= MPI_PROC_NULL) THEN
        CALL putbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ds, izup_ds, jzlo_ds, jzup_ds, nzcount_ds, 4 )
      ENDIF

      CALL MPI_SENDRECV                                                      &
           ( sendbuf(1,4), nzcount_ds,  imp_type, MPI_neighbors(4), ntag,    &
             sendbuf(1,6), isendbuflen, imp_type, MPI_neighbors(2), ntag,    &
             icomm, nzstatus, izmplcode)
      IF (izmplcode /= 0) THEN
        ierror  = izmplcode
        yerrmsg = 'MPI_SENDRECV'
        RETURN
      ENDIF

      IF (MPI_neighbors(2) /= MPI_PROC_NULL) THEN
        CALL getbuf ( trac_field, trac_dim,                                  &
                      sendbuf, isendbuflen, idim, jdim, kdim,                &
                      izlo_ur, izup_ur, jzlo_ur, jzup_ur, nzcount_ur, 6 )
      ENDIF

    ENDIF

  ENDIF

  DEALLOCATE (sendbuf)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE exchg_trac_boundaries

!==============================================================================
!==============================================================================
!+ This subroutine puts all necessary values into sendbuf
!------------------------------------------------------------------------------

SUBROUTINE putbuf  (trac_field, trac_dim,                    &
                    sendbuf, isendbuflen, idim, jdim, kdim,  &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine puts the necessary values from the present variables
!   (determined by ilo, iup, jlo, jup) into sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim, &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of the sendbuf

  INTEGER, INTENT (INOUT)      ::    &
    ncount                          ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field (idim, jdim, trac_dim, kdim)   ! first field that has to occur
! Local variables

  INTEGER :: i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 2: Put data into the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O''Neill to
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt =1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              sendbuf (nzc,nentry) = trac_field(i,j,jt,k)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE putbuf

!==============================================================================
!==============================================================================
!+ This subroutine gets all necessary values from sendbuf
!------------------------------------------------------------------------------

SUBROUTINE getbuf  (trac_field, trac_dim, &
                    sendbuf, isendbuflen, idim, jdim, kdim,                 &
                    ilo, iup, jlo, jup, ncount, nentry )

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine gets the necessary values for the present variables
!   (determined by ilo, iup, jlo, jup) from sendbuf.
!
! Method:
!   Check which variables are present.
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    isendbuflen,                  & ! length of sendbuffer
    idim, jdim, kdim,         & ! dimensions of the fields
    trac_dim,     &
    ilo, iup, jlo, jup,           & ! start- and end-indices
    nentry                          ! specifies the row of sendbuf to be used

  INTEGER, INTENT (INOUT)  :: ncount      ! counts the variables

  REAL (dp), INTENT (INOUT)            ::    &
    sendbuf (isendbuflen, 8),     & ! send buffer
    trac_field(idim, jdim, trac_dim,kdim)   ! first field that has to occur

! Local variables

  INTEGER ::   i, j, jt, k, nzc

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!- Section 2: Get data from the buffer
!------------------------------------------------------------------------------

  ! use nzc as a local counter  (based on a work from Mike O''Neill to
  ! improve vectorization of putbuf and getbuf)
  nzc = ncount

  ! first variable that has to be present
  DO k = 1, kdim
     DO jt = 1, trac_dim
        DO j = jlo, jup
           DO i = ilo, iup
              nzc = nzc + 1
              trac_Field(i,j,jt,k) = sendbuf (nzc,nentry)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! put nzc to global counter
  ncount = nzc

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE getbuf

!==============================================================================
!==============================================================================
!+ Calls MPI_BARRIER
!------------------------------------------------------------------------------

SUBROUTINE comm_barrier (icomm, ierror, yerrmsg)

!------------------------------------------------------------------------------
!
! Description:
!
! Method:
!   MPI-routine MPI_BARRIER
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)  :: icomm            ! communicator to be used

  INTEGER, INTENT (OUT) :: ierror           ! error-status variable

  CHARACTER (LEN=*), INTENT (OUT) :: yerrmsg  ! for MPI error message

! Local variables

  INTEGER ::  izmplcode                     ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

ierror    = 0
yerrmsg   = '   '
izmplcode = 0

CALL MPI_BARRIER (icomm, izmplcode)

IF (izmplcode /= 0) THEN
  ierror = izmplcode
  yerrmsg = 'MPI_BARRIER'
ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE comm_barrier

!==============================================================================
!==============================================================================
!+ This subroutine defines and allocates MPI data types for exchg_boundaries
!------------------------------------------------------------------------------

SUBROUTINE setup_data_type(trac_field, trac_dim,                             &
       idim, jdim, kdim, ilo, iup, jlo, jup,                                 &
       ierror, yerrmsg, imp_type, ncount, type_handle)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates and commits a data type consisting of a
!   certain subgrid of up to 14 variables. The subgrid is specified
!   via  idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup . Only
!   one variable has to occur, the other are optional.
!
!   The values of "var1(ilo, jlo, klo), ncount, datatype"
!   should constitute the right entries for the dummy variables
!   "BUF, COUNT, DATATYPE" in a call to
!   MPI_SEND/MPI_RECV. This should describe the subgrid for all
!   up to 14 variables.
!
!   As a consequence, the same data type can only be used, if the same
!   variables (same position in memory !) are used.
!
!  Author: C. Pospiech, IBM
!
!------------------------------------------------------------------------------

! Subroutine arguments
  INTEGER, INTENT (IN)         ::    &
    imp_type,                     & ! determines the REAL type used
    idim, jdim,                   & ! horizontal dimensions of the fields
    kdim,                     & ! vertical dimensions of var01..var20
    trac_dim, &
    ilo, iup, jlo, jup              ! start- and end-indices

  INTEGER, INTENT (OUT)      ::    &
    ncount,                       & ! how many copies of type_handle
    type_handle                     ! handle for MPI data type

  INTEGER,           INTENT (OUT) ::  ierror       ! error status variable

  CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg       ! for MPI error message

  REAL(dp), INTENT (INOUT) ::  trac_field(idim,jdim,trac_dim, kdim)

! Local variables

  INTEGER ::               &
    nzc,                   &
    sect1d, sect2d,        &   ! Variables to hold intermediate
    sect3d, sect4d,        &
    meta_vect,             &   ! MPI data types
    meta_addr,             &   ! Vector of addresses of the varxx
    meta_disp,             &   ! Displacements of varxx in memory
    meta_blklen,           &   ! some intermediate variable
    num_meta_entries,      &   ! how many varxx are present
    disp(2), blocklen(2),  &   ! Variables needed to define
    vartype(2),            &   ! MPI data type of a certain extent
    sizeofreal                 ! size of data type in byte

  INTEGER :: izmplcode                   ! for MPI error code

!------------------------------------------------------------------------------
!- End of header -
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Section 1: Initializations
!------------------------------------------------------------------------------

  sect3d = MPI_DATATYPE_NULL
  sect4d = MPI_DATATYPE_NULL

!------------------------------------------------------------------------------
!- Section 2: Set up of MPI data types *** subarrays
!------------------------------------------------------------------------------

  ! set up 1-dimensional section
  nzc = iup - ilo + 1
  CALL MPI_TYPE_CONTIGUOUS  (nzc, imp_type, sect1d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_CONTIGUOUS TRAC'
    RETURN
  ENDIF

  ! set up 2-dimensional section
  nzc = jup - jlo +1
  CALL MPI_TYPE_EXTENT(imp_type, sizeofreal, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_EXTENT TRAC'
    RETURN
  ENDIF
  CALL MPI_TYPE_HVECTOR    (nzc, 1, idim*sizeofreal,         &
                            sect1d, sect2d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-2 TRAC'
    RETURN
  ENDIF

!US: this must be done for every optional entry
  ! set up 3-dimensional section
  CALL MPI_TYPE_HVECTOR    (trac_dim, 1, idim*jdim*sizeofreal,    &
                              sect2d, sect3d, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-3 TRAC'
    RETURN
  ENDIF
  CALL MPI_TYPE_HVECTOR    (kdim, 1, idim*jdim*sizeofreal,    &
                              sect3d, type_handle, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_HVECTOR-4 TRAC'
    RETURN
  ENDIF

!------------------------------------------------------------------------------
!- Section 3: Set up of MPI data types *** meta structure from all varxx
!------------------------------------------------------------------------------

  meta_addr = 0 ! initialize

  CALL MPI_ADDRESS (trac_field, meta_addr, izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_ADDRESS-01 TRAC'
    RETURN
  ENDIF

!!$  meta_disp(:)   = meta_addr(:) - meta_addr(1)
!!$  meta_blklen(:) = 1
!!$  CALL MPI_TYPE_STRUCT  (num_meta_entries, meta_blklen, meta_disp, sect3d,  &
!!$                         meta_vect, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

!------------------------------------------------------------------------------
!- Section 4: Reset extent of this new data type by defining upper bound
!------------------------------------------------------------------------------

!US Bin mir nicht sicher, wozu das gut sein soll
!!$  blocklen(:)   = 1
!!$  disp(1)     = 0
!!$  disp(2)     = sizeofreal*(iup - ilo + 2)
!!$  vartype(1)  = meta_vect
!!$  vartype(2)  = MPI_UB
!!$  CALL MPI_TYPE_STRUCT     (2, blocklen, disp, vartype, type_handle, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

!------------------------------------------------------------------------------
!- Section 5: Commit the data type
!------------------------------------------------------------------------------

  CALL MPI_TYPE_COMMIT       (type_handle,izmplcode)
  IF (izmplcode /= 0) THEN
    ierror  = izmplcode
    yerrmsg = 'MPI_TYPE_COMMIT'
    RETURN
  ENDIF
  ncount = 1

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE setup_data_type

!==============================================================================
!... !defined(BLANK) && !defined(CESM1) && ! defined(VERTICO)
#endif
! ############################################################################

! ############################################################################
#if defined(CESM1)

  subroutine finish(name, text)

    ! use abortutils, only: endrun

    CHARACTER(*) :: name
    CHARACTER(*), OPTIONAL :: text

    write(*,*) name
    if (present(text)) write(*,*) text


    ! call endrun
    stop

  end subroutine finish

  SUBROUTINE message (name, text)
    CHARACTER (*) :: name
    CHARACTER (*), OPTIONAL :: text

    write(*,*) name
    if (present(text)) write(*,*) text

  END SUBROUTINE MESSAGE

!... CESM1
#endif

!... (defined COSMO) || defined(BLANK) || defined(CESM1) || defined(VERTICO)
#endif
! ############################################################################
! ############################################################################

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
! op_pj_20121001+
!------------------------------------------------------------------------------
SUBROUTINE p_allgather_2d1d(recvbuf, sendbuf, comm)

  USE messy_main_constants_mem, ONLY: dp
#if defined(ECHAM5) || defined(ICON)
  USE mo_mpi, ONLY: p_real_dp                      ! op_pj_20121001
#endif
#if defined(COSMO)
  USE data_parallel,    ONLY: p_real_dp => imp_reals
#endif
#if defined(CESM1)
  USE spmd_utils,       ONLY: p_real_dp => mpi_real8
#endif
#if defined(MBM_CLAMS)
  USE mo_mpi, ONLY: p_real_dp                      ! ju_nt_20160304
#endif
#if defined(MBM_MPIOM)
  USE mo_mpi, ONLY: p_real_dp                      ! mz_ap_20161103
#endif

  IMPLICIT NONE

  ! I/O
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: recvbuf
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: sendbuf
  INTEGER,  OPTIONAL,       INTENT(IN)  :: comm

  ! LOCAL
#ifndef NOMPI
  INTEGER :: p_error
  INTEGER :: p_comm

  IF (PRESENT(comm)) THEN
     p_comm = comm
  ELSE
     p_comm = p_all_comm
  ENDIF

!!$  CALL MPI_ALLGATHER(sendbuf, SIZE(sendbuf),   p_real_dp, &
!!$                     recvbuf, SIZE(recvbuf,1), p_real_dp, &
!!$                     p_comm, p_error)
  CALL MPI_GATHER(sendbuf, SIZE(sendbuf),   p_real_dp, &
                  recvbuf, SIZE(recvbuf,1), p_real_dp, &
                  p_io, p_comm, p_error)
  CALL MPI_BCAST(recvbuf, SIZE(recvbuf), p_real_dp, p_io, p_comm, p_error)
#else
  recvbuf(:,LBOUND(recvbuf,2)) = sendbuf(:)
#endif

END SUBROUTINE p_allgather_2d1d
!------------------------------------------------------------------------------
! op_pj_20121001-
! op_pj_20121027+
!------------------------------------------------------------------------------
SUBROUTINE p_allgather_4d3d(recvbuf, sendbuf, comm)

  USE messy_main_constants_mem, ONLY: dp
#if defined(ECHAM5) || defined(ICON)
  USE mo_mpi, ONLY: p_real_dp                      ! op_pj_20121001
#endif
#if defined(COSMO)
  USE data_parallel,    ONLY: p_real_dp => imp_reals
#endif
#if defined(MBM_CLAMS)
  USE mo_mpi, ONLY: p_real_dp                      ! ju_nt_20160304
#endif
#if defined(MBM_MPIOM)
  USE mo_mpi, ONLY: p_real_dp                      ! ap_mz_20161103
#endif

  IMPLICIT NONE

  ! I/O
  REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: recvbuf
  REAL(DP), DIMENSION(:,:,:),   INTENT(IN)  :: sendbuf
  INTEGER,  OPTIONAL,           INTENT(IN)  :: comm

  ! LOCAL
#ifndef NOMPI
  INTEGER :: p_error
  INTEGER :: p_comm

  IF (PRESENT(comm)) THEN
     p_comm = comm
  ELSE
     p_comm = p_all_comm
  ENDIF

  CALL MPI_ALLGATHER(sendbuf, SIZE(sendbuf)                , p_real_dp, &
                     recvbuf, SIZE(recvbuf)/SIZE(recvbuf,4), p_real_dp, &
                     p_comm, p_error)
!!$  CALL MPI_GATHER(sendbuf, SIZE(sendbuf),                 p_real_dp, &
!!$                  recvbuf, SIZE(recvbuf)/SIZE(recvbuf,4), p_real_dp, &
!!$                  p_io, p_comm, p_error)
!!$  CALL MPI_BCAST(recvbuf, SIZE(recvbuf), p_real_dp, p_io, p_comm, p_error)
#else
  recvbuf(:,:,:,LBOUND(recvbuf,4)) = sendbuf(:,:,:)
#endif

END SUBROUTINE p_allgather_4d3d
!------------------------------------------------------------------------------
! op_pj_20121027-
! op_pj_20160408+
!------------------------------------------------------------------------------
SUBROUTINE p_allgather_3d2d(recvbuf, sendbuf, comm)

  USE messy_main_constants_mem, ONLY: dp
#if defined(ECHAM5) || defined(ICON)
  USE mo_mpi, ONLY: p_real_dp
#endif
#if defined(COSMO)
  USE data_parallel,    ONLY: p_real_dp => imp_reals
#endif
#if defined(MBM_CLAMS)
  USE mo_mpi, ONLY: p_real_dp                      ! ju_nt_20160304
#endif
#if defined(MBM_MPIOM)
  USE mo_mpi, ONLY: p_real_dp                      ! mz_ap_20161103
#endif

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

  CALL MPI_ALLGATHER(sendbuf, SIZE(sendbuf)                , p_real_dp, &
                     recvbuf, SIZE(recvbuf)/SIZE(recvbuf,3), p_real_dp, &
                     p_comm, p_error)
!!$  CALL MPI_GATHER(sendbuf, SIZE(sendbuf),                 p_real_dp, &
!!$                  recvbuf, SIZE(recvbuf)/SIZE(recvbuf,3), p_real_dp, &
!!$                  p_io, p_comm, p_error)
!!$  CALL MPI_BCAST(recvbuf, SIZE(recvbuf), p_real_dp, p_io, p_comm, p_error)
#else
  recvbuf(:,:,LBOUND(recvbuf,3)) = sendbuf(:,:)
#endif

END SUBROUTINE p_allgather_3d2d
!------------------------------------------------------------------------------
! op_pj_20160408-
! op_pj_20121027+
!------------------------------------------------------------------------------
SUBROUTINE p_scatter_4d3d(sendbuf, recvbuf, comm)

  USE messy_main_constants_mem, ONLY: dp
#if defined(ECHAM5) || defined(ICON)
  USE mo_mpi, ONLY: p_real_dp                      ! op_pj_20121001
#endif
#if defined(COSMO)
  USE data_parallel,    ONLY: p_real_dp => imp_reals
#endif
#if defined(MBM_CLAMS)
  USE mo_mpi, ONLY: p_real_dp                      ! ju_nt_20160304
#endif
#if defined(MBM_MPIOM)
  USE mo_mpi, ONLY: p_real_dp                      ! mz_ap_20161103
#endif

  IMPLICIT NONE
#ifdef NOMPI
  INTRINSIC :: LBOUND
#endif

  ! I/O
  REAL(DP), DIMENSION(:,:,:,:), INTENT(IN)  :: sendbuf
  REAL(DP), DIMENSION(:,:,:),   INTENT(OUT) :: recvbuf
  INTEGER,  OPTIONAL,           INTENT(IN)  :: comm

  ! LOCAL
#ifndef NOMPI
  INTEGER :: p_error
  INTEGER :: p_comm

  IF (PRESENT(comm)) THEN
     p_comm = comm
  ELSE
     p_comm = p_all_comm
  ENDIF

! op_pj_20130606+
!!$  CALL MPI_Scatter(sendbuf, SIZE(sendbuf)/SIZE(sendbuf,4), p_real_dp, &
  CALL MPI_Scatter(sendbuf, SIZE(recvbuf), p_real_dp, &
! op_pj_20130606-
                   recvbuf, SIZE(recvbuf), p_real_dp, p_io, &
                   p_comm, p_error)
#else
  recvbuf(:,:,:) = sendbuf(:,:,:,LBOUND(sendbuf,4))
#endif

END SUBROUTINE p_scatter_4d3d
!------------------------------------------------------------------------------
! op_pj_20121027-

!------------------------------------------------------------------------------
SUBROUTINE get_node_ids(status, comm, id)

#ifndef NOMPI
#ifndef __PGI
  USE mpi
#endif
#endif
  USE messy_main_tools, ONLY: unique

  IMPLICIT NONE
  INTRINSIC :: TRIM, ADJUSTL

#ifndef NOMPI
#ifdef __PGI
   INCLUDE 'mpif.h'
#endif
#endif

  ! I/O
  INTEGER,                            INTENT(OUT) :: status
  INTEGER,                            INTENT(IN)  :: comm ! communiator
  INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: id  ! node-IDs

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'get_node_ids'
  INTEGER :: npes
  INTEGER :: p_error
  INTEGER :: my_pe
  INTEGER :: len, i, j
#ifndef NOMPI
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), DIMENSION(:), ALLOCATABLE :: names
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), DIMENSION(:), ALLOCATABLE :: uniqnames
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME)  :: local_name
#else
  CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: names
  CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: uniqnames
#endif

#ifndef NOMPI
  CALL MPI_COMM_SIZE(comm, npes, p_error)
  IF (p_error /= MPI_SUCCESS) THEN
     WRITE(iouerr,*) 'ERROR: MPI_COMM_SIZE failed (',substr,')'
     status = 1
     RETURN
  END IF
  CALL MPI_COMM_RANK(comm, my_pe, p_error)
  IF (p_error /= MPI_SUCCESS) THEN
     IF (p_error == MPI_ERR_COMM) THEN
        ! current taks is not a valid rank in comm
        my_pe = -1
     ELSE
        WRITE(iouerr,*) 'ERROR: MPI_COMM_RANK failed (',substr,')'
        status = 1
        RETURN
     END IF
  END IF
#else
  npes = 1
  my_pe = 0
  status = 0
#endif

  !write(*,*) substr,': PE = ',my_pe, ' out of ', npes-1

  ALLOCATE(id(0:npes-1))
  ALLOCATE(names(0:npes-1))

#ifndef NOMPI
  CALL MPI_GET_PROCESSOR_NAME(local_name, len, p_error )
  IF (p_error /= MPI_SUCCESS) THEN
     WRITE(iouerr,*) 'ERROR: MPI_GET_PROCESSOR_NAME failed (',substr,')'
     status = 1
     RETURN
  END IF
#else
  names(0) = 'node01'
  status = 0
#endif

  !write(*,*) substr,': ',my_pe, names(my_pe)

#ifndef NOMPI
  IF (my_pe >= 0) THEN
     CALL MPI_ALLGATHER(local_name, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                        names(:),     MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                        comm, p_error)
     IF (p_error /= MPI_SUCCESS) THEN
        WRITE(iouerr,*) 'ERROR: MPI_ALLGATHER failed (',substr,')'
        status = 1
        RETURN
     END IF
     names(my_pe) = local_name
  END if
#else
  ! nothing to do
#endif

  ! create a list of unique node names
  CALL unique(names, uniqnames)

  DO i=0, npes-1
     DO j=1, SIZE(uniqnames)
        IF ( TRIM(ADJUSTL(names(i))) == TRIM(ADJUSTL(uniqnames(j))) ) THEN
           id(i) = j-1
           EXIT
        END IF
     END DO
  END DO

  ! DEBUG OUTPUT
  IF (my_pe == 0) THEN
     DO i=0, npes-1
        WRITE(*,*) substr,': PE = ',i ,'; node = ', TRIM(ADJUSTL(names(i))) &
             , '; id = ', id(i)
     END DO
  END IF

  ! FREE MEMORY
  DEALLOCATE(names, uniqnames)

  status = 0

END SUBROUTINE get_node_ids
!------------------------------------------------------------------------------

#if defined(COSMO) || defined(BLANK) || defined(CESM1)
! um_ak_20130822+

  ! send implementation

  SUBROUTINE p_send_real (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, 1, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real

  ! ---------------------------------------------------------

  SUBROUTINE p_send_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real_1d

  ! ---------------------------------------------------------

  SUBROUTINE p_send_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real_2d

  ! ---------------------------------------------------------

  SUBROUTINE p_send_real_3d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real_3d

  ! ---------------------------------------------------------

  SUBROUTINE p_send_real_4d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real_4d

  ! ---------------------------------------------------------

  SUBROUTINE p_send_real_5d (buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    ELSE
       CALL MPI_SEND (buffer, SIZE(buffer), imp_reals, p_destination, p_tag, &
            p_comm, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', p_pe, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_send_real_5d

  ! ---------------------------------------------------------
  ! recv implementation

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, 1, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real_1d

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
    buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real_2d

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real_3d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real_3d

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real_4d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real_4d

  ! ---------------------------------------------------------

  SUBROUTINE p_recv_real_5d (buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error
    INTEGER :: p_status(MPI_STATUS_SIZE)   ! standard information of MPI_RECV

#ifdef MESSY
      buffer = 0._dp
#endif
    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), imp_reals, p_source, p_tag, &
            p_comm, p_status, p_error)
    END IF

#ifdef DEBUG
    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', p_pe, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', p_error
       CALL p_abort
    END IF
#endif
!... MPI
#endif

  END SUBROUTINE p_recv_real_5d
! um_ak_20130822-
!... defined(COSMO) || defined(BLANK) || defined(CESM1)
#endif

!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================
#ifdef ICON

  ! ---------------------------------------------------------------------------
  ! ub_ak_20180116+
  !SUBROUTINE messy_mpi_initialize
  SUBROUTINE messy_mpi_setup
  ! ub_ak_20180116-
    !
    ! This subroutines initializes some variables needed in MESSy
    !
    p_parallel = my_process_is_mpi_all_parallel()
    p_parallel_io = my_process_is_stdio()
    ! op_bk_20170126+
    p_pe = get_my_mpi_all_id()
    p_nprocs = get_my_mpi_all_comm_size()
!!$  p_pe = get_my_mpi_work_id()
!!$  p_nprocs = get_my_mpi_work_comm_size()
!!$  p_nprocs = num_work_procs
    p_io0 = p_work0
!!$  p_parallel_io = (p_pe == p_io0)
    ! op_bk_20170126-
    ! ub_ak_20180116+
  END SUBROUTINE messy_mpi_setup
  !END SUBROUTINE messy_mpi_initialize
    ! ub_ak_20180116-
  ! ---------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE putbuf  (trac_field, trac_dim,     &
       sendbuf, isendbuflen, idim, jdim, kdim,  &
       ilo, iup, jlo, jup, ncount, nentry )

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine puts the necessary values from the present variables
    !   (determined by ilo, iup, jlo, jup) into sendbuf.
    !
    ! Method:
    !   Check which variables are present.
    !
    !--------------------------------------------------------------------------

    ! Subroutine arguments
    INTEGER, INTENT (IN)         ::    &
         isendbuflen,                  & ! length of sendbuffer
         idim, jdim, kdim,         & ! dimensions of the fields
         trac_dim, &
         ilo, iup, jlo, jup,           & ! start- and end-indices
         nentry                          ! specifies the row of the sendbuf

    INTEGER, INTENT (INOUT)      ::    &
         ncount                          ! counts the variables

    REAL (dp), INTENT (INOUT)            ::    &
         sendbuf (isendbuflen, 8),     & ! send buffer
         trac_field (idim, jdim, trac_dim, kdim)   ! first field that has to occur
    ! Local variables

    INTEGER :: i, j, jt, k, nzc

    !--------------------------------------------------------------------------
    !- End of header -
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !- Section 1: Initializations
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !- Section 2: Put data into the buffer
    !--------------------------------------------------------------------------

    ! use nzc as a local counter  (based on a work from Mike O'Neill to
    ! improve vectorization of putbuf and getbuf)
    nzc = ncount

    ! first variable that has to be present
    DO k = 1, kdim
       DO jt =1, trac_dim
          DO j = jlo, jup
             DO i = ilo, iup
                nzc = nzc + 1
                sendbuf (nzc,nentry) = trac_field(i,j,jt,k)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! put nzc to global counter
    ncount = nzc

    !--------------------------------------------------------------------------
    !- End of the Subroutine
    !--------------------------------------------------------------------------
  END SUBROUTINE putbuf
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE getbuf  (trac_field, trac_dim,    &
       sendbuf, isendbuflen, idim, jdim, kdim, &
       ilo, iup, jlo, jup, ncount, nentry )

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine gets the necessary values for the present variables
    !   (determined by ilo, iup, jlo, jup) from sendbuf.
    !
    ! Method:
    !   Check which variables are present.
    !
    !--------------------------------------------------------------------------

    ! Subroutine arguments
    INTEGER, INTENT (IN)         ::    &
         isendbuflen,                  & ! length of sendbuffer
         idim, jdim, kdim,             & ! dimensions of the fields
         trac_dim,                     &
         ilo, iup, jlo, jup,           & ! start- and end-indices
         nentry                       ! specifies the row of sendbuf to be used

    INTEGER, INTENT (INOUT)  :: ncount      ! counts the variables

    REAL (dp), INTENT (INOUT)            ::    &
         sendbuf (isendbuflen, 8),     & ! send buffer
         trac_field(idim, jdim, trac_dim,kdim)   ! first field that has to occur

    ! Local variables

    INTEGER ::   i, j, jt, k, nzc

    !--------------------------------------------------------------------------
    !- End of header -
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !- Section 1: Initializations
    !--------------------------------------------------------------------------


    !--------------------------------------------------------------------------
    !- Section 2: Get data from the buffer
    !--------------------------------------------------------------------------

    ! use nzc as a local counter  (based on a work from Mike O'Neill to
    ! improve vectorization of putbuf and getbuf)
    nzc = ncount

    ! first variable that has to be present
    DO k = 1, kdim
       DO jt = 1, trac_dim
          DO j = jlo, jup
             DO i = ilo, iup
                nzc = nzc + 1
                trac_Field(i,j,jt,k) = sendbuf (nzc,nentry)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! put nzc to global counter
    ncount = nzc

    !--------------------------------------------------------------------------
    !- End of the Subroutine
    !--------------------------------------------------------------------------
  END SUBROUTINE getbuf
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE comm_barrier (icomm, ierror, yerrmsg)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !
    ! Method:
    !   MPI-routine MPI_BARRIER
    !
    !--------------------------------------------------------------------------

    ! Subroutine arguments
    INTEGER, INTENT (IN)  :: icomm            ! communicator to be used

    INTEGER, INTENT (OUT) :: ierror           ! error-status variable

    CHARACTER (LEN=*), INTENT (OUT) :: yerrmsg  ! for MPI error message

    ! Local variables

    INTEGER ::  izmplcode                     ! for MPI error code

    !--------------------------------------------------------------------------
    !- End of header -
    !--------------------------------------------------------------------------

    ierror    = 0
    yerrmsg   = '   '
    izmplcode = 0

    CALL MPI_BARRIER (icomm, izmplcode)

    IF (izmplcode /= 0) THEN
       ierror = izmplcode
       yerrmsg = 'MPI_BARRIER'
    ENDIF

    !--------------------------------------------------------------------------
    !- End of the Subroutine
    !--------------------------------------------------------------------------

  END SUBROUTINE comm_barrier
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE setup_data_type(trac_field, trac_dim,     &
       idim, jdim, kdim, ilo, iup, jlo, jup,           &
       ierror, yerrmsg, imp_type, ncount, type_handle)

    !--------------------------------------------------------------------------
    !
    ! Description:
    !   This subroutine allocates and commits a data type consisting of a
    !   certain subgrid of up to 14 variables. The subgrid is specified
    !   via  idim, jdim, kdim, ilo, iup, jlo, jup, klo, kup . Only
    !   one variable has to occur, the other are optional.
    !
    !   The values of "var1(ilo, jlo, klo), ncount, datatype"
    !   should constitute the right entries for the dummy variables
    !   "BUF, COUNT, DATATYPE" in a call to
    !   MPI_SEND/MPI_RECV. This should describe the subgrid for all
    !   up to 14 variables.
    !
    !   As a consequence, the same data type can only be used, if the same
    !   variables (same position in memory !) are used.
    !
    !  Author: C. Pospiech, IBM
    !
    !--------------------------------------------------------------------------

    ! Subroutine arguments
    INTEGER, INTENT (IN)         ::    &
         imp_type,                     & ! determines the REAL type used
         idim, jdim,                   & ! horizontal dimensions of the fields
         kdim,                         & ! vertical dimensions of var01..var20
         trac_dim,                     &
         ilo, iup, jlo, jup              ! start- and end-indices

    INTEGER, INTENT (OUT)      ::      &
         ncount,                       & ! how many copies of type_handle
         type_handle                     ! handle for MPI data type

    INTEGER,           INTENT (OUT) ::  ierror       ! error status variable

    CHARACTER (LEN=*), INTENT(OUT)  :: yerrmsg       ! for MPI error message

    REAL(dp), INTENT (INOUT) ::  trac_field(idim,jdim,trac_dim, kdim)

    ! Local variables

    INTEGER ::               &
         nzc,                   &
         sect1d, sect2d,        &   ! Variables to hold intermediate
         sect3d, sect4d,        &
         meta_vect,             &   ! MPI data types
         meta_addr,             &   ! Vector of addresses of the varxx
         meta_disp,             &   ! Displacements of varxx in memory
         meta_blklen,           &   ! some intermediate variable
         num_meta_entries,      &   ! how many varxx are present
         disp(2), blocklen(2),  &   ! Variables needed to define
         vartype(2),            &   ! MPI data type of a certain extent
         sizeofreal                 ! size of data type in byte

    INTEGER :: izmplcode                   ! for MPI error code

    !--------------------------------------------------------------------------
    !- End of header -
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !- Section 1: Initializations
    !--------------------------------------------------------------------------

    sect3d = MPI_DATATYPE_NULL
    sect4d = MPI_DATATYPE_NULL

    !--------------------------------------------------------------------------
    !- Section 2: Set up of MPI data types *** subarrays
    !--------------------------------------------------------------------------

    ! set up 1-dimensional section
    nzc = iup - ilo + 1
    CALL MPI_TYPE_CONTIGUOUS  (nzc, imp_type, sect1d, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_CONTIGUOUS TRAC'
       RETURN
    ENDIF

    ! set up 2-dimensional section
    nzc = jup - jlo +1
    CALL MPI_TYPE_EXTENT(imp_type, sizeofreal, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_EXTENT TRAC'
       RETURN
    ENDIF
    CALL MPI_TYPE_HVECTOR    (nzc, 1, idim*sizeofreal,         &
         sect1d, sect2d, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_HVECTOR-2 TRAC'
       RETURN
    ENDIF

    !US: this must be done for every optional entry
    ! set up 3-dimensional section
    CALL MPI_TYPE_HVECTOR    (trac_dim, 1, idim*jdim*sizeofreal,    &
         sect2d, sect3d, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_HVECTOR-3 TRAC'
       RETURN
    ENDIF
    CALL MPI_TYPE_HVECTOR    (kdim, 1, idim*jdim*sizeofreal,    &
         sect3d, type_handle, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_HVECTOR-4 TRAC'
       RETURN
    ENDIF

    !--------------------------------------------------------------------------
    !- Section 3: Set up of MPI data types *** meta structure from all varxx
    !--------------------------------------------------------------------------

    meta_addr = 0 ! initialize

    CALL MPI_ADDRESS (trac_field, meta_addr, izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_ADDRESS-01 TRAC'
       RETURN
    ENDIF

!!$  meta_disp(:)   = meta_addr(:) - meta_addr(1)
!!$  meta_blklen(:) = 1
!!$  CALL MPI_TYPE_STRUCT  (num_meta_entries, meta_blklen, meta_disp, sect3d,  &
!!$                         meta_vect, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

    !--------------------------------------------------------------------------
    !- Section 4: Reset extent of this new data type by defining upper bound
    !--------------------------------------------------------------------------

    !US Bin mir nicht sicher, wozu das gut sein soll
!!$  blocklen(:)   = 1
!!$  disp(1)     = 0
!!$  disp(2)     = sizeofreal*(iup - ilo + 2)
!!$  vartype(1)  = meta_vect
!!$  vartype(2)  = MPI_UB
!!$  CALL MPI_TYPE_STRUCT     (2, blocklen, disp, vartype, type_handle, izmplcode)
!!$  IF (izmplcode /= 0) THEN
!!$    ierror  = izmplcode
!!$    yerrmsg = 'MPI_TYPE_STRUCT'
!!$    RETURN
!!$  ENDIF

    !--------------------------------------------------------------------------
    !- Section 5: Commit the data type
    !--------------------------------------------------------------------------

    CALL MPI_TYPE_COMMIT       (type_handle,izmplcode)
    IF (izmplcode /= 0) THEN
       ierror  = izmplcode
       yerrmsg = 'MPI_TYPE_COMMIT'
       RETURN
    ENDIF
    ncount = 1

    !--------------------------------------------------------------------------
    !- End of the Subroutine
    !--------------------------------------------------------------------------

  END SUBROUTINE setup_data_type
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_0d_int (zfield, comm) RESULT (p_sum)

    INTEGER                       :: p_sum
    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, 1, p_int_i4, &
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, imp_integers, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d_int
  !----------------------------------------------------------------------------

  ! op_bk_20141211+
  !----------------------------------------------------------------------------
  FUNCTION p_sum_2d_int (zfield, comm) RESULT (p_sum)

    INTEGER,           INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error  ! op_pj_20120925

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_integers, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d_int
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_0d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, 1, imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_0d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_1d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_2d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_2d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_3d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2),SIZE(zfield,3))

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_3d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  FUNCTION p_sum_4d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2),SIZE(zfield,3),SIZE(zfield,4))

#ifndef NOMPI
!#ifdef MPI
    INTEGER :: p_comm
    INTEGER :: p_error

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_all_comm
    ENDIF

    IF (p_parallel) THEN
!!$       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), p_real_dp, &
       CALL MPI_ALLREDUCE (zfield, p_sum, SIZE(zfield), imp_reals, &
            MPI_SUM, p_comm, p_error)
    ELSE
       p_sum = zfield
    END IF
#else
    p_sum = zfield
#endif

  END FUNCTION p_sum_4d
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE p_bcast_real(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real

  SUBROUTINE p_bcast_real_single(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (sp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_single

  SUBROUTINE p_bcast_real_1d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_1d

  SUBROUTINE p_bcast_real_1d_single(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (sp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_1d_single

  SUBROUTINE p_bcast_real_2d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_2d

  SUBROUTINE p_bcast_real_2d_single(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (sp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_2d_single

  SUBROUTINE p_bcast_real_3d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_real_5d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_real_5d

  SUBROUTINE p_bcast_int_i4(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER (i4), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER (i8), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER,   INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_i8_1d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER(i8), INTENT(inout) :: t_buffer(:)
    INTEGER,     INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_i8_1d

  SUBROUTINE p_bcast_int_2d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER,   INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER,   INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    INTEGER,   INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_int_4d

  SUBROUTINE p_bcast_bool(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    LOGICAL,   INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    LOGICAL,   INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    LOGICAL,   INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    LOGICAL,   INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    LOGICAL,   INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    CHARACTER(len=*),   INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_char

  SUBROUTINE p_bcast_char_1d(t_buffer, p_source, comm)
#ifndef NOMPI
    USE mo_mpi, ONLY:  p_comm_work
#endif
    CHARACTER(*),   INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = p_comm_work
    ENDIF

    CALL p_bcast_icon(t_buffer, p_source, p_comm)
#endif
  END SUBROUTINE p_bcast_char_1d
  !----------------------------------------------------------------------------
  ! op_bk_20141211-

#endif
!==============================================================================
!==============================================================================
!==============================================================================
!==============================================================================

! **************************************************************************
END MODULE messy_main_mpi_bi
! **************************************************************************
