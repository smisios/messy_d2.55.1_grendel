
MODULE messy_mpiom_tools_e5

  ! this submodel contains some tools that are used for MPIOM
  ! Author Pozzer Andrea, October 2007
  USE messy_main_constants_mem,  ONLY: DP


  IMPLICIT NONE
  PRIVATE

   PUBLIC :: nullify_borders
   PUBLIC :: read_file_2d
   PUBLIC :: read_file_3d


CONTAINS

!   ! ---------------------------------------------------------------------------

   SUBROUTINE nullify_borders(status, arr)

   USE messy_mpiom_mem_e5,   ONLY : ke, ie, je, ie_g, je_g, icycli,             &
                                    have_g_is, have_g_ie, have_g_js, have_g_je, &
                                    nprocxy,nprocx,nprocy
   USE messy_main_mpi_bi,   ONLY:  p_pe

   IMPLICIT NONE

   real(dp), INTENT(INOUT) :: arr(:,:,:)
   integer, intent(INOUT) :: status
   integer :: i,nm,np

   ! this subroutine nullify the exchanged border in the
   ! MPIOM grid. Sometime useful.

     status = 1

    DO i=1,SIZE(arr,DIM=3)

        IF (nprocx>1 .AND. p_pe<nprocxy ) THEN
          ! x-direction
            IF(icycli/=0 .OR. .NOT. have_g_is) THEN
              IF (MOD(p_pe,nprocx)/=0) arr(1,:,i) = 0._dp
            ENDIF
            IF(icycli/=0 .OR. .NOT. have_g_ie) THEN
              IF ( MOD(p_pe,nprocx)/=nprocx-1) arr(ie,:,i) = 0._dp
            ENDIF
        ENDIF
        
        IF (nprocy>1 .AND. p_pe<nprocxy) THEN
          nm = MOD(p_pe+nprocxy-nprocx,nprocxy)
          np = MOD(p_pe+nprocxy+nprocx,nprocxy)
          ! y-direction
          ! Note that there is no action required if nprocy==1
          IF (nm/=np) THEN
            IF(.NOT. have_g_js) THEN
              arr(:,1,i) = 0._dp
            ENDIF
            IF(.NOT. have_g_je) THEN
              arr(:,je,i) = 0._dp
            ENDIF
          ENDIF
        ENDIF

    ENDDO
   
   END SUBROUTINE nullify_borders

  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_2d(fname, varname,  data_file)

   USE netcdf

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname         ! variable name
    REAL(dp), DIMENSION(:,:), INTENT(OUT):: data_file ! INTENT(OUT)
    
    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_file_2d'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        STOP
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0
    
  END SUBROUTINE read_file_2d
  ! ------------------------------------------------------------------------
  SUBROUTINE read_file_3d(fname, varname,  data_file)

   USE netcdf

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*),          INTENT(IN)  :: fname       ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: varname     ! variable name
    REAL(dp), DIMENSION(:,:,:),INTENT(OUT) :: data_file   ! INTENT(OUT)
    
    ! LOCAL
    INTEGER                     :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_start_netcdf'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL NFERR( status, &
         nf90_inq_varid(ncid, TRIM(varname), varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable ",varname," not found in NetCDF file!"
        STOP
    ENDIF

    CALL NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0
    
  END SUBROUTINE read_file_3d

  ! ------------------------------------------------------------------
  SUBROUTINE NFERR(status,command,pos)

    USE netcdf

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN) :: command
    INTEGER,          INTENT(IN) :: pos

    status=command
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF ERROR at position: ', pos
       WRITE(*,*) 'netCDF ERROR status     : ',status
       WRITE(*,*) 'netCDF ERROR            : ',nf90_strerror(status)
    END IF

  END SUBROUTINE NFERR
  ! ------------------------------------------------------------------


! ***************************************************************************
END MODULE messy_mpiom_tools_e5
! ***************************************************************************
