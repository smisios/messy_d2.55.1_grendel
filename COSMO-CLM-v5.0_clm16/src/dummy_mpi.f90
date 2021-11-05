! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Dummy interfaces for MPI routines
!-------------------------------------------------------------------------------
!
! Description:
!   This file provides dummy interfaces for the MPI-calls that can be
!   used on sequential platforms where MPI is not available. In this
!   way the compiler will not complain about "unsatisfied external
!   references".
!   Some routines, however, really have to set the variables according
!   to the sequential mode of the program. These routines are:
!    
!     - mpi_comm_size
!     - mpi_comm_rank
!     - mpi_comm_group
!     - mpi_error_string
!     - mpi_abort
!
!   All other routines set the error code to -9999, indicating that they
!   should not be called in sequential mode. 
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! !VERSION!  !DATE!     Ulrich Schaettler
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================

!-------------------------------------------------------------------------------

!*******************************************************************************

!-------------------------------------------------------------------------------
SUBROUTINE MPI_TEST(REQUEST, FLAG, STATUS, IERROR)

  INCLUDE "mpif.h"

  INTEGER, INTENT(in)   ::  REQUEST
  INTEGER, INTENT(out)  ::  STATUS(MPI_STATUS_SIZE), IERROR
  LOGICAL, INTENT(out)   :: FLAG

  flag = .true.
  ierror = 0
  status(:) = 0

END SUBROUTINE MPI_TEST

!-------------------------------------------------------------------------------

SUBROUTINE MPI_GET_ADDRESS (var01, meta_addr, izmplcode)

  REAL    :: var01 (:,:,:)
  INTEGER :: meta_addr
  INTEGER, INTENT(out) :: izmplcode

  izmplcode = 0

END SUBROUTINE MPI_GET_ADDRESS

!-------------------------------------------------------------------------------

SUBROUTINE MPI_SIZEOF (xtype, isize, izmplcode)

  INTEGER :: xtype
  INTEGER :: isize
  INTEGER, INTENT(out) :: izmplcode

  izmplcode = -1

END SUBROUTINE MPI_SIZEOF

!-------------------------------------------------------------------------------

SUBROUTINE MPI_TYPE_MATCH_SIZE  (MPI_TYPECLASS, isize, xtype, izmplcode)

  INTEGER :: MPI_TYPECLASS
  INTEGER :: isize
  INTEGER :: xtype
  INTEGER, INTENT(out) :: izmplcode

  izmplcode = -1

END SUBROUTINE MPI_TYPE_MATCH_SIZE

!-------------------------------------------------------------------------------

SUBROUTINE MPI_TYPE_CREATE_STRUCT  (num_meta_entries, meta_blklen, meta_disp,   &
                                sect3d, meta_vect, izmplcode)

  INTEGER :: num_meta_entries
  INTEGER :: meta_blklen(:)
  INTEGER :: meta_disp(:)
  INTEGER :: sect3d(:)
  INTEGER :: meta_vect
  INTEGER, INTENT(out) :: izmplcode

  izmplcode = 0

END SUBROUTINE MPI_TYPE_CREATE_STRUCT

!-------------------------------------------------------------------------------

SUBROUTINE mpi_error_string  (ierrorcode, ystring, iresultlen, ierror)
  INTEGER, INTENT (IN)  :: ierrorcode
  INTEGER, INTENT (OUT) :: iresultlen, ierror
  CHARACTER (LEN=*), INTENT (OUT)        :: ystring

! Begin SUBROUTINE mpi_error_string
  ystring    = ' *** Only dummy mpi available *** '
  iresultlen = 0
  ierror     = -9999
END SUBROUTINE mpi_error_string

!-------------------------------------------------------------------------------

SUBROUTINE mpi_abort  (icomm, ierrorcode, ierror)
  INTEGER, INTENT (IN)  :: icomm, ierrorcode, ierror

! Begin SUBROUTINE mpi_abort
  IF (ierrorcode /= 0) THEN
    STOP 1
  ELSE
    STOP
  ENDIF

END SUBROUTINE mpi_abort

!-------------------------------------------------------------------------------

SUBROUTINE mpi_init  (ierror)
  INTEGER, INTENT (OUT) :: ierror

! Begin SUBROUTINE mpi_init
  ierror = 0
END SUBROUTINE mpi_init

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_size  (icomm, isize, ierror)
  INTEGER, INTENT (IN)  :: icomm
  INTEGER, INTENT (OUT) :: isize, ierror

! Begin SUBROUTINE mpi_comm_size
  isize  = 1
  ierror = 0
END SUBROUTINE mpi_comm_size

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_rank  (icomm, irank, ierror)
  INTEGER, INTENT (IN)  :: icomm
  INTEGER, INTENT (OUT) :: irank, ierror

! Begin SUBROUTINE mpi_comm_rank
  irank  = 0
  ierror = 0
END SUBROUTINE mpi_comm_rank

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_group  (icomm, igroup, ierror)
  INTEGER, INTENT (IN)  :: icomm
  INTEGER, INTENT (OUT) :: igroup, ierror

! Begin SUBROUTINE mpi_comm_group
  igroup = 0
  ierror = 0
END SUBROUTINE mpi_comm_group

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_split (icomm, icolor, ikey, inewcomm, ierror)
  INTEGER, INTENT (IN)  :: icomm, icolor, ikey
  INTEGER, INTENT (OUT) :: inewcomm, ierror

! Begin SUBROUTINE mpi_comm_split
  inewcomm   = -1
  ierror     = -9999
END SUBROUTINE mpi_comm_split

!-------------------------------------------------------------------------------

SUBROUTINE MPI_COMM_GET_ATTR(COMM, COMM_KEYVAL, ATTRIBUTE_VAL, FLAG, IERROR)

  INCLUDE "mpif.h"

  INTEGER, INTENT(in)  ::  COMM, COMM_KEYVAL
  INTEGER, INTENT(out) ::  IERROR
  INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(out) :: ATTRIBUTE_VAL
  LOGICAL, intent(out) ::  FLAG

  flag = .TRUE.
  IERROR = 0
  ATTRIBUTE_VAL = 0

END SUBROUTINE mpi_comm_get_attr

!-------------------------------------------------------------------------------

SUBROUTINE mpi_bcast (realbuffer, icount, idatatype, iroot, icomm, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (INOUT) :: realbuffer (icount)
  INTEGER, INTENT (IN)    :: idatatype, iroot, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_bcast
  realbuffer(:) = 0.0
  ierror        = -9999
END SUBROUTINE mpi_bcast

!-------------------------------------------------------------------------------

SUBROUTINE mpi_cart_create  (icomm_old, ndims, idimens, lperiods,           &
                             lreorder, icomm_cart, ierror)
  INTEGER, INTENT (IN)  :: icomm_old, ndims
  INTEGER, INTENT (IN)  :: idimens(ndims)
  LOGICAL, INTENT (IN)                   :: lperiods(ndims), lreorder
  INTEGER, INTENT (OUT) :: icomm_cart, ierror

! Begin SUBROUTINE mpi_cart_create
  icomm_cart = 0
  ierror     = -9999
END SUBROUTINE mpi_cart_create

!-------------------------------------------------------------------------------

SUBROUTINE mpi_cart_coords (icomm, irank, ndims, icoords, ierror)
  INTEGER, INTENT (IN)  :: icomm, irank, ndims
  INTEGER, INTENT (OUT) :: icoords(*), ierror

! Begin SUBROUTINE mpi_cart_coords
  icoords(1:ndims) = 0
  ierror     = -9999
END SUBROUTINE mpi_cart_coords

!-------------------------------------------------------------------------------

SUBROUTINE mpi_cart_rank  (icomm, icoords, irank, ierror)
  INTEGER, INTENT (IN)  :: icomm, icoords(*)
  INTEGER, INTENT (OUT) :: irank, ierror

! Begin SUBROUTINE mpi_cart_rank
  irank  = 0
  ierror = -9999
END SUBROUTINE mpi_cart_rank

!-------------------------------------------------------------------------------

SUBROUTINE mpi_group_incl  (igroup, n, iranks, inewgroup, ierror)
  INTEGER, INTENT (IN)  :: igroup, n
  INTEGER, INTENT (IN)  :: iranks(n)
  INTEGER, INTENT (OUT) :: inewgroup, ierror

! Begin SUBROUTINE mpi_group_incl
  inewgroup = 0
  ierror    = -9999
END SUBROUTINE mpi_group_incl

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_create  (icomm, igroup, inewcomm, ierror)
  INTEGER, INTENT (IN)  :: icomm, igroup
  INTEGER, INTENT (OUT) :: inewcomm, ierror

! Begin SUBROUTINE mpi_comm_create
  inewcomm = 0
  ierror   = -9999
END SUBROUTINE mpi_comm_create

!-------------------------------------------------------------------------------

SUBROUTINE mpi_gather (rsendbuf, isendcount, isendtype,                       &
                       rrecvbuf, irecvcount, irecvtype, iroot, icomm, ierror)
  INTEGER, INTENT (IN)    :: isendcount, irecvcount
  REAL,    INTENT (IN)    :: rsendbuf (isendcount)
  REAL,    INTENT (OUT)   :: rrecvbuf (irecvcount)
  INTEGER, INTENT (IN)    :: isendtype, irecvtype, iroot, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_gather
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_gather

!-------------------------------------------------------------------------------

SUBROUTINE mpi_gatherv (rsendbuf, isendcount, isendtype,                      &
                        rrecvbuf, irecvcounts, idispls, irecvtype,            &
                        iroot, icomm, ierror)
  INTEGER, INTENT (IN)    :: isendcount, irecvcounts(*), idispls(*)
  REAL,    INTENT (IN)    :: rsendbuf (isendcount)
  REAL,    INTENT (OUT)   :: rrecvbuf (isendcount)
  INTEGER, INTENT (IN)    :: isendtype, irecvtype, iroot, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_gatherv
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_gatherv

!-------------------------------------------------------------------------------
SUBROUTINE mpi_wait  (irequest, istatus, ierror)
  INTEGER, INTENT (INOUT)  :: irequest
  INTEGER, INTENT (OUT) :: istatus(*), ierror

! Begin SUBROUTINE mpi_wait
  istatus(1) = 0
  ierror     = -9999
END SUBROUTINE mpi_wait

!-------------------------------------------------------------------------------

SUBROUTINE mpi_waitall  (icount, irequest, istatus, ierror)
  INTEGER, INTENT (IN)     :: icount
  INTEGER, INTENT (INOUT)  :: irequest(*)
  INTEGER, INTENT (OUT)    :: istatus (*), ierror

! Begin SUBROUTINE mpi_waitall
  irequest(1:icount) = 0
  istatus (1:icount) = 0
  ierror      = -9999
END SUBROUTINE mpi_waitall

!-------------------------------------------------------------------------------

SUBROUTINE mpi_get_count  (istatus, idatatype, icount, ierror)
  INTEGER, INTENT (IN)  :: istatus(*), idatatype
  INTEGER, INTENT (OUT) :: icount, ierror

! Begin SUBROUTINE mpi_get_count
  icount = 0
  ierror = -9999
END SUBROUTINE mpi_get_count

!-------------------------------------------------------------------------------

SUBROUTINE mpi_send  (realbuf, icount, idatatype, idest, itag, icomm, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (IN)    :: realbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, idest, itag, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_send
  ierror = -9999
END SUBROUTINE mpi_send

!-------------------------------------------------------------------------------

SUBROUTINE mpi_isend  (realbuf, icount, idatatype, idest, itag, icomm,        &
                       irequest, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (IN)    :: realbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, idest, itag, icomm
  INTEGER, INTENT (OUT)   :: irequest, ierror

! Begin SUBROUTINE mpi_isend
  irequest    = 0
  ierror      = -9999
END SUBROUTINE mpi_isend

!-------------------------------------------------------------------------------

SUBROUTINE mpi_recv  (realbuf, icount, idatatype, isource, itag, icomm,      &
                      istatus, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (OUT)   :: realbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, isource, itag, icomm
  INTEGER, INTENT (OUT)   :: istatus(*), ierror

! Begin SUBROUTINE mpi_recv
  realbuf(1:icount) = 0.0
  istatus(1:icount) = 0
  ierror     = -9999
END SUBROUTINE mpi_recv

!-------------------------------------------------------------------------------

SUBROUTINE mpi_irecv  (realbuf, icount, idatatype, isource, itag, icomm,     &
                     irequest, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (OUT)   :: realbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, isource, itag, icomm
  INTEGER, INTENT (OUT)   :: irequest , ierror

! Begin SUBROUTINE mpi_irecv
  irequest   = 0
  realbuf(:) = 0.0
  ierror     = -9999
END SUBROUTINE mpi_irecv

!-------------------------------------------------------------------------------

SUBROUTINE mpi_barrier  (icomm, ierror)
  INTEGER, INTENT (IN)  :: icomm
  INTEGER, INTENT (OUT) :: ierror

! Begin SUBROUTINE mpi_barrier
  ierror = -9999
END SUBROUTINE mpi_barrier

!-------------------------------------------------------------------------------

SUBROUTINE mpi_reduce  (rsendbuf, rrecvbuf, icount, idatatype, iop, iroot,   &
                        icomm, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (IN)    :: rsendbuf (icount)
  REAL,    INTENT (OUT)   :: rrecvbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, iop, iroot, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_reduce
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_reduce

!-------------------------------------------------------------------------------

SUBROUTINE mpi_scatter    (rsendbuf, iscount, isdatatype,         &
                           rrecvbuf, ircount, irdatatype,         &
                           irecv, icomm, ierror)

  INTEGER, INTENT (IN)    :: iscount, ircount
  REAL,    INTENT (IN)    :: rsendbuf (iscount)
  REAL,    INTENT (OUT)   :: rrecvbuf (ircount)
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, irecv, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_scatter  
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_scatter  

!-------------------------------------------------------------------------------

SUBROUTINE mpi_scatterv   (rsendbuf, iscounts, idispls, isdatatype,         &
                           rrecvbuf, ircount ,          irdatatype,         &
                           irecv, icomm, ierror)

  INTEGER, INTENT (IN)    :: iscounts(*), ircount, idispls(*)
  REAL,    INTENT (IN)    :: rsendbuf (*)
  REAL,    INTENT (OUT)   :: rrecvbuf (ircount)
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, irecv, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_scatterv
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_scatterv

!-------------------------------------------------------------------------------

SUBROUTINE mpi_allgather  (rsendbuf, iscount, isdatatype,            &
                           rrecvbuf, ircount, irdatatype,            &
                           icomm, ierror)  

  INTEGER, INTENT (IN)    :: iscount, ircount
  REAL,    INTENT (IN)    :: rsendbuf (iscount)
  REAL,    INTENT (OUT)   :: rrecvbuf (ircount)
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_allgather
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_allgather

!-------------------------------------------------------------------------------

SUBROUTINE mpi_allgatherv (rsendbuf, iscount, isdatatype,                     &
                           rrecvbuf, ircount, iidispl, irdatatype,            &
                           icomm, ierror)  

  INTEGER, INTENT (IN)    :: iscount, ircount, iidispl
  REAL,    INTENT (IN)    :: rsendbuf (iscount)
  REAL,    INTENT (OUT)   :: rrecvbuf (ircount)
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_allgatherv
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_allgatherv

!-------------------------------------------------------------------------------

SUBROUTINE mpi_alltoall  (rsendbuf, iscounts, isdatatype,           &
                           rrecvbuf, ircounts, irdatatype,           &
                           icomm, ierror) 

  INTEGER, INTENT (IN)    :: iscounts(*), ircounts(*)
  REAL,    INTENT (IN)    :: rsendbuf (*)
  REAL,    INTENT (OUT)   :: rrecvbuf (iscounts(1))
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_alltoall
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_alltoall

!-------------------------------------------------------------------------------

SUBROUTINE mpi_alltoallv  (rsendbuf, iscounts, isdispl, isdatatype,           &
                           rrecvbuf, ircounts, irdispl, irdatatype,           &
                           icomm, ierror)  

  INTEGER, INTENT (IN)    :: iscounts(*), ircounts(*), isdispl(*), irdispl(*)
  REAL,    INTENT (IN)    :: rsendbuf (*)
  REAL,    INTENT (OUT)   :: rrecvbuf (iscounts(1))
  INTEGER, INTENT (IN)    :: isdatatype, irdatatype, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_alltoallv 
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_alltoallv 

!-------------------------------------------------------------------------------

SUBROUTINE mpi_allreduce  (rsendbuf, rrecvbuf, icount, idatatype, iop,       &
                           icomm, ierror)
  INTEGER, INTENT (IN)    :: icount
  REAL,    INTENT (IN)    :: rsendbuf (icount)
  REAL,    INTENT (OUT)   :: rrecvbuf (icount)
  INTEGER, INTENT (IN)    :: idatatype, iop, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_allreduce
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_allreduce

!-------------------------------------------------------------------------------

SUBROUTINE mpi_finalize  (ierror)
  INTEGER, INTENT (OUT) :: ierror

! Begin SUBROUTINE mpi_finalize
  ierror = -9999
END SUBROUTINE mpi_finalize

!-------------------------------------------------------------------------------

SUBROUTINE mpi_iprobe (isource, itag, icomm, lflag, istatus, ierror)

  INTEGER, INTENT (IN)  :: isource, itag, icomm
  INTEGER, INTENT (OUT) :: istatus, ierror
  LOGICAL, INTENT (OUT) :: lflag

! Begin SUBROUTINE mpi_iprobe
  ierror  = -9999
  istatus = -1
  lflag   = .FALSE.
END SUBROUTINE mpi_iprobe

!-------------------------------------------------------------------------------

SUBROUTINE mpi_intercomm_create (ilocal_comm, ilocal_leader, ipeer_comm,      &
                                 iremote_leader, itag, inewintercomm, ierror)
  INTEGER, INTENT (IN)  :: ilocal_comm, ilocal_leader, ipeer_comm,            &
                           iremote_leader, itag
  INTEGER, INTENT (OUT) :: inewintercomm, ierror

! Begin SUBROUTINE mpi_intercomm_create
  ierror        = -9999
  inewintercomm = -1
END SUBROUTINE mpi_intercomm_create

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_compare (icomm1, icomm2, iresult, ierror)

  INTEGER, INTENT (IN)  :: icomm1, icomm2
  INTEGER, INTENT (OUT) :: iresult, ierror

INCLUDE "mpif.h"

! Begin SUBROUTINE mpi_comm_compare
  ierror  = -9999
  iresult = MPI_CONGRUENT
END SUBROUTINE mpi_comm_compare

!-------------------------------------------------------------------------------

SUBROUTINE mpi_comm_dup (icommold, icommnew, ierror)

  INTEGER, INTENT (IN)  :: icommold
  INTEGER, INTENT (OUT) :: icommnew, ierror

! Begin SUBROUTINE mpi_comm_dup
  ierror   = -9999
  icommnew = -1
END SUBROUTINE mpi_comm_dup

!-------------------------------------------------------------------------------

SUBROUTINE mpi_group_translate_ranks (iori_group, num_compute, isource_ranks, iclose_group, &
                                 iclose_comm_compute_ranks, ierror)

  INTEGER, INTENT (IN)    :: iori_group, num_compute, isource_ranks(:), iclose_group, &
                             iclose_comm_compute_ranks(:)
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_group_translate_ranks
  ierror      = -9999
END SUBROUTINE mpi_group_translate_ranks

!-------------------------------------------------------------------------------

SUBROUTINE mpi_testsome (imax_request, irecv_request, ioutcount, ioutindex, ioutstatus, ierror)

  INTEGER, INTENT (IN)    :: imax_request, irecv_request(:), ioutcount, ioutindex(:), &
                             ioutstatus(:,:)
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_testsome
  ierror      = -9999
END SUBROUTINE mpi_testsome

!-------------------------------------------------------------------------------

SUBROUTINE mpi_probe  (isource, itag, icomm, ierror)

  INTEGER, INTENT (IN)    :: isource, itag, icomm
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_probe
  ierror      = -9999
END SUBROUTINE mpi_probe

!-------------------------------------------------------------------------------

SUBROUTINE mpi_sendrecv  (rsendbuf, iscount, isdatatype, idest,   istag,     &
                          rrecvbuf, ircount, irdatatype, isource, irtag,     &
                          icomm, ierror)  

  INTEGER, INTENT (IN)    :: iscount, isdatatype, idest,   istag,           &
                             ircount, irdatatype, isource, irtag, icomm
  REAL,    INTENT (IN)    :: rsendbuf (iscount)
  REAL,    INTENT (OUT)   :: rrecvbuf (ircount)
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_sendrecv
  rrecvbuf(:) = 0.0
  ierror      = -9999
END SUBROUTINE mpi_sendrecv

!-------------------------------------------------------------------------------

SUBROUTINE mpi_type_contiguous (icount, ioldtype, inewtype, ierror)

  INTEGER, INTENT (IN)    :: icount, ioldtype
  INTEGER, INTENT (OUT)   :: inewtype, ierror

! Begin SUBROUTINE mpi_type_contiguous
  inewtype    = -1
  ierror      = -9999
END SUBROUTINE mpi_type_contiguous

!-------------------------------------------------------------------------------

SUBROUTINE mpi_type_hvector (icount, ILEN, istride, ioldtype, inewtype, ierror)

  INTEGER, INTENT (IN)    :: icount, ILEN, istride, ioldtype
  INTEGER, INTENT (OUT)   :: inewtype, ierror

! Begin SUBROUTINE mpi_type_hvector
  inewtype    = -1
  ierror      = -9999
END SUBROUTINE mpi_type_hvector

!-------------------------------------------------------------------------------

SUBROUTINE mpi_type_struct (icount, iarrlen,iarrdsp,iarrtype, inewtype, ierror)

  INTEGER, INTENT (IN)    :: icount, iarrlen,iarrdsp,iarrtype
  INTEGER, INTENT (OUT)   :: inewtype, ierror

! Begin SUBROUTINE mpi_type_struct
  inewtype    = -1
  ierror      = -9999
END SUBROUTINE mpi_type_struct

!-------------------------------------------------------------------------------

SUBROUTINE mpi_type_extent (idatatype, iextent, ierror)

  INTEGER, INTENT (IN)    :: idatatype
  INTEGER, INTENT (OUT)   :: iextent, ierror

! Begin SUBROUTINE mpi_type_extent
  iextent     = -1
  ierror      = -9999
END SUBROUTINE mpi_type_extent

!-------------------------------------------------------------------------------

SUBROUTINE mpi_type_commit (idatatype, ierror)

  INTEGER, INTENT (INOUT) :: idatatype
  INTEGER, INTENT (OUT)   :: ierror

! Begin SUBROUTINE mpi_type_commit
  idatatype   = -1
  ierror      = -9999
END SUBROUTINE mpi_type_commit

!-------------------------------------------------------------------------------

SUBROUTINE mpi_address (rtype, iaddress, ierror)

  REAL,    INTENT (IN)    :: rtype(*)
  INTEGER, INTENT (OUT)   :: iaddress, ierror

! Begin SUBROUTINE mpi_address
  iaddress    = -1
  ierror      = -9999
END SUBROUTINE mpi_address

!-------------------------------------------------------------------------------

FUNCTION MPI_NULL_COPY_FN (invec, inoutvec,len, TYPE)
  INTEGER   :: len, TYPE
  INTEGER   :: invec(len), inoutvec(len)
  MPI_NULL_COPY_FN = 0
END FUNCTION MPI_NULL_COPY_FN

FUNCTION MPI_NULL_DELETE_FN (invec, inoutvec,len, TYPE)
  INTEGER   :: len, TYPE
  INTEGER   :: invec(len), inoutvec(len)
  MPI_NULL_DELETE_FN = 0
END FUNCTION MPI_NULL_DELETE_FN

FUNCTION MPI_DUP_FN (invec, inoutvec,len, TYPE)
  INTEGER   :: len, TYPE
  INTEGER   :: invec(len), inoutvec(len)
  MPI_DUP_FN = 0
END FUNCTION MPI_DUP_FN

DOUBLE PRECISION FUNCTION MPI_WTICK ()
  MPI_WTICK = 0.0
END FUNCTION MPI_WTICK

DOUBLE PRECISION FUNCTION MPI_WTIME ()
  MPI_WTIME = 0.0
END FUNCTION MPI_WTIME

!*******************************************************************************

! end dummy_mpi
