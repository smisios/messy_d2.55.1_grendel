MODULE MESSY_MAIN_GRID_NETCDF_BI

IMPLICIT NONE
PRIVATE

PUBLIC :: p_bcast_ncvar
PUBLIC :: p_bcast_ncatt
PUBLIC :: p_bcast_ncdim
PUBLIC :: p_bcast_narray

CONTAINS

! ----------------------------------------------------------
SUBROUTINE P_BCAST_NCVAR(var, proc)
  
  USE MESSY_MAIN_MPI_BI,      ONLY: P_BCAST, p_pe

  USE MESSY_MAIN_GRID_NETCDF, ONLY: t_ncvar, INIT_NCDIM, INIT_NCATT &
                                  , INIT_NCVAR ! op_pj_20180807

  IMPLICIT NONE

  ! I/O
  TYPE(t_ncvar), INTENT(INOUT) :: var
  INTEGER      , INTENT(IN)    :: proc
  ! LOCAL
  INTEGER :: i ! loop variable

  IF (p_pe /= proc) CALL INIT_NCVAR(var) ! op_pj_20180807

  CALL P_BCAST(var%name,  proc)
  CALL P_BCAST(var%ID,    proc)
  CALL P_BCAST(var%XTYPE, proc)
  CALL P_BCAST(var%ndims, proc)

  IF (p_pe /= proc) ALLOCATE(var%dim(var%ndims))
  DO i=1,var%ndims
     IF (p_pe /= proc) CALL INIT_NCDIM(var%dim(i))
     CALL P_BCAST_NCDIM(var%dim(i),  proc)
  END DO

  CALL P_BCAST(var%uid,   proc)
  CALL P_BCAST(var%ustep, proc)
  CALL P_BCAST(var%natts, proc)

  IF (p_pe /= proc) ALLOCATE(var%att(var%natts))
  DO i=1,var%natts
     IF (p_pe /= proc) CALL INIT_NCATT(var%att(i))
     CALL P_BCAST_NCATT(var%att(i), proc)
  END DO

  CALL P_BCAST_NARRAY(var%dat, proc)

END SUBROUTINE P_BCAST_NCVAR
! ----------------------------------------------------------

! ----------------------------------------------------------
SUBROUTINE P_BCAST_NCDIM(dim, proc)
  
  USE MESSY_MAIN_MPI_BI,      ONLY: P_BCAST

  USE MESSY_MAIN_GRID_NETCDF, ONLY: t_ncdim
  
  IMPLICIT NONE

  ! I/O
  TYPE(t_ncdim), INTENT(INOUT) :: dim
  INTEGER      , INTENT(IN)    :: proc

  CALL P_BCAST(dim%name,  proc)
  CALL P_BCAST(dim%ID,    proc)
  CALL P_BCAST(dim%len,   proc)
  CALL P_BCAST(dim%fuid,  proc)
  CALL P_BCAST(dim%varid, proc)

END SUBROUTINE P_BCAST_NCDIM
! ----------------------------------------------------------

! ----------------------------------------------------------
SUBROUTINE P_BCAST_NCATT(att, proc)
  
  USE MESSY_MAIN_MPI_BI,      ONLY: P_BCAST

  USE MESSY_MAIN_GRID_NETCDF, ONLY: t_ncatt
  
  IMPLICIT NONE

  ! I/O
  TYPE (t_ncatt), INTENT(INOUT) :: att
  INTEGER       , INTENT(IN)    :: proc

  CALL P_BCAST(att%name,  proc)
  CALL P_BCAST(att%ID,    proc)
  CALL P_BCAST(att%xtype, proc)
  CALL P_BCAST(att%varid, proc)
  CALL P_BCAST(att%len,   proc)
  CALL P_BCAST_NARRAY(att%dat, proc)

END SUBROUTINE P_BCAST_NCATT
! ----------------------------------------------------------

! ----------------------------------------------------------
SUBROUTINE P_BCAST_NARRAY(array, proc)
  
  USE MESSY_MAIN_MPI_BI,      ONLY: P_BCAST, p_pe
  USE MESSY_MAIN_GRID_NETCDF, ONLY: t_narray, VTYPE_UNDEF, VTYPE_INT &
                                  , VTYPE_CHAR &
                                  , VTYPE_DOUBLE, VTYPE_REAL, VTYPE_BYTE &
                                  , INIT_NARRAY
 
  IMPLICIT NONE
  INTRINSIC :: PRODUCT

  ! I/O
  TYPE(t_narray), INTENT(INOUT) :: array
  INTEGER       , INTENT(IN)    :: proc
  ! LOCAL
  INTEGER :: i ! loop variable
  INTEGER :: s ! size of array  ! op_pj_20180807

  IF (p_pe /= proc) CALL INIT_NARRAY(array)

  CALL P_BCAST(array%type,  proc)
  CALL P_BCAST(array%n,     proc)
  IF (p_pe /= proc) ALLOCATE(array%dim(array%n))
  DO i=1, array%n
     CALL P_BCAST(array%dim(i), proc)
  END DO

  IF (array%n == 0) RETURN  ! op_pj_20180807
  s = PRODUCT(array%dim(:)) ! op_pj_20180807

  SELECT CASE (array%type)
  CASE(VTYPE_UNDEF)
     ! NOTHING TO DO
  CASE(VTYPE_INT)
     IF (p_pe /= proc) ALLOCATE(array%vi(s))
     CALL P_BCAST(array%vi(:), proc)
  CASE(VTYPE_REAL)
     IF (p_pe /= proc) ALLOCATE(array%vr(s))
     CALL P_BCAST(array%vr(:), proc)
  CASE(VTYPE_DOUBLE)
     IF (p_pe /= proc) ALLOCATE(array%vd(s))
     CALL P_BCAST(array%vd(:), proc)
  CASE(VTYPE_BYTE)
     IF (p_pe /= proc) ALLOCATE(array%vb(s))
     CALL P_BCAST(array%vb(:), proc)
  CASE(VTYPE_CHAR)
     IF (p_pe /= proc) ALLOCATE(array%vc(s))
     CALL P_BCAST(array%vc(:), proc)
  END SELECT

END SUBROUTINE P_BCAST_NARRAY
! ----------------------------------------------------------

END MODULE MESSY_MAIN_GRID_NETCDF_BI
