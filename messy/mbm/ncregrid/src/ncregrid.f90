! ******************************************************************
! ------------------------------------------------------------------
PROGRAM NCREGRID
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE mo_f2kcli                    ! command line interface
  USE messy_ncregrid_control       ! regridding interface

  IMPLICIT NONE

  INTRINSIC :: TRIM

#ifdef PNCREGRID
! mz_pj_20080317+
!!$  USE mpi       ! mz_kk_20080218
#include <mpif.h>
! mz_pj_20080317-
#endif
 
  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

#ifdef PNCREGRID
  ! mz_kk_20080218+
  INTEGER         :: ierr, root
  TYPE(t_mpi_def) :: my_mpi

  ! MPI Setup
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, my_mpi%rank,  ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, my_mpi%nproc, ierr)
  my_mpi%comm = MPI_COMM_WORLD
  ! mz_kk_20080218-
#endif

  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name

  IF (NARG > 1) THEN 
     WRITE(*,*) 'Too many arguments !'
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  IF (NARG == 0) THEN 
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1,CMD)  

  ! REGRIDDING ...
  RG_CTRL = RG_PROC
  RG_NML  = NML_NEXT
  DO ! endless DO loop (must be terminated with EXIT)
#ifdef PNCREGRID
     ! mz_kk_20080218+
     CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS, TRIM(CMD) &
          , my_mpi=my_mpi, root_pe=root)
     ! mz_kk_20080218-
#else
     CALL REGRID_CONTROL(RG_CTRL, RG_NML, RG_STATUS, TRIM(CMD))
#endif
     IF (RG_STATUS == RGSTAT_STOP) EXIT ! leave endless DO loop
  END DO
  ! END OF REGRIDDING

#ifdef PNCREGRID
  CALL MPI_Finalize(ierr)  ! mz_kk_20080218
#endif

CONTAINS

  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'NCREGRID Version ',NCREGRIDVERS
    WRITE(*,*) 'Author: Patrick Joeckel, MPICH, June 2002'
#ifdef PNCREGRID
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'MPI parallelisation (number of variables)'
    WRITE(*,*) 'Author: Klaus Ketelsen, MPICH, Dec 2007'
#endif
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '-----------------------------------------'
  END SUBROUTINE USAGE

END PROGRAM NCREGRID
! ------------------------------------------------------------------
