      SUBROUTINE READ_NETCDF_VAR(ncid,desc,arr,klev)
!**************************************************************************
!
! Reads a variable from a NETCDF file and distributes it to all PEs
!
! The NETCDF File is only accessed by p_io
!
!**************************************************************************

      USE mo_param1
      USE mo_parallel

      implicit none

      INTEGER ncid, klev
      CHARACTER (LEN=*) desc
      REAL(wp) arr(ie,je,klev)

      REAL(wp),ALLOCATABLE :: arr_g(:,:,:)
      CHARACTER (LEN=80) err_text

      INCLUDE 'netcdf.inc'
      INTEGER ncstat,ncvarid
      !INTEGER :: k

! Read NETCDF data

      IF(p_pe==p_io) THEN
         ALLOCATE(arr_g(ie_g,je_g,klev))
      ELSE
         ALLOCATE(arr_g(0,0,0))
      ENDIF

      IF(p_pe==p_io) THEN
        err_text = 'AUFR: Problem reading '//desc
        ncstat = NF_INQ_VARID(ncid,desc,ncvarid)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all(err_text)
        ncstat = NF_GET_VAR_DOUBLE(ncid,ncvarid,arr_g)
        IF ( ncstat .NE. NF_NOERR ) CALL stop_all(err_text)
      ENDIF

!      DO k=1,klev
        CALL scatter(arr_g(:,:,:),arr(:,:,:),p_io)
!      ENDDO
      END
