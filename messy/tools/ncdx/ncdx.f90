  program ncdx

  ! modifies netCDF file for openDX import
  ! written by Patrick Joeckel, MPICH, Dez 2001
  ! Version 3.8 (April 2002)

    USE f2kcli
    USE netcdf
    IMPLICIT NONE

    ! GLOBAL (netCDF-file)
    INTEGER :: status                  ! netCDF status
    INTEGER :: ncid                    ! netCDF file - ID
    INTEGER :: ndims                   ! number of dimensions in file
    INTEGER :: udimid                  ! dimension ID of unlimited dimension
    INTEGER :: nvar                    ! number of variables
    CHARACTER (LEN=100), ALLOCATABLE :: dimname(:) ! name of dimensions
    INTEGER, ALLOCATABLE :: dimlen(:)  ! dimension lenghts
    INTEGER, ALLOCATABLE :: dimvid(:)  ! corresponding variable ID
    INTEGER, ALLOCATABLE :: dimtype(:) ! type of corresponding coordinate var.
    LOGICAL, ALLOCATABLE :: axdim(:,:,:) ! dimensionality required of dim-axis

    ! VARIABLE SPECIFIC
    INTEGER :: varid                   ! variable ID
    CHARACTER (LEN=100) :: varname     ! variable name
    INTEGER :: nvdim                   ! number of variable dimensions
    INTEGER :: nvsdim                  ! number of 'spatial' var. dimensions
    INTEGER, ALLOCATABLE :: dimids(:)  ! dimension IDs of variable
    CHARACTER (LEN=200)  :: fieldatt   ! additional variable attribute
    CHARACTER (LEN=200)  :: posatt     ! additional variable attribute

    ! NEW VARIABLES / DIMENSIONS
    INTEGER, ALLOCATABLE :: newvid(:,:,:)! IDs of new variables
    INTEGER, ALLOCATABLE :: dxdimid(:) ! new additional dimension IDs
    CHARACTER (LEN=100) :: dxatt       ! ncdx global attribute
    CHARACTER(LEN=8)  :: date          ! current date
    CHARACTER(LEN=10) :: time          ! current time
    CHARACTER(LEN=5)  :: zone          ! time zone
    CHARACTER(LEN=3)  :: dimstr        ! axis dimensionality as string
    CHARACTER(LEN=3)  :: posstr        ! axis position as string

    ! INTERNAL
    LOGICAL :: lex                     ! file exists?
    INTEGER :: i, j, k                 ! counter
    INTEGER :: dcount                  ! dimension counter
    INTEGER :: idx1, idx2              ! string indices
    INTEGER :: nvsdimcnt               ! nvsdim counter

    CHARACTER, ALLOCATABLE :: dc(:)    ! character data
    CHARACTER, ALLOCATABLE :: dch(:)   ! character data
    INTEGER,   ALLOCATABLE :: di(:)    ! integer data
    REAL,      ALLOCATABLE :: dr(:)    ! real (float) data
    DOUBLE PRECISION, ALLOCATABLE :: dd(:)    ! double precision data

    ! FOR COMMAND LINE
    CHARACTER(LEN=256)  :: EXE         ! program name
    CHARACTER (LEN=100) :: file        ! netCDF filename
    INTEGER             :: NARG        ! number of arguments

    ! PARSE COMMAND LINE
    NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
    CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name
    IF (NARG /= 1) THEN 
       WRITE(*,*) 'COMMAND LINE ERROR!'
       CALL USAGE(TRIM(EXE)) 
       STOP
    END IF
    ! GET FILENAME
    CALL GET_COMMAND_ARGUMENT(NARG,file)

    ! CHECK IF FILE EXISTS AND OPEN IT
    INQUIRE(file=TRIM(file), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) 'File '//TRIM(file)//' not found !'
       STOP
    END IF
    CALL nf(nf90_open(TRIM(file),NF90_WRITE,ncid))
    CALL nf(nf90_redef(ncid))

    ! CHECK GLOBAL ATTRIBUTE FOR openDX COMPLIANCE
    dxatt = ''
    status = nf90_get_att(ncid, NF90_GLOBAL,'ncdx',dxatt)
    IF (status == NF90_NOERR) THEN
       WRITE(*,*) 'File '//TRIM(file)//' is already openDX compliant:'
       WRITE(*,*) TRIM(dxatt)
       STOP
    END IF

    ! GET NUMBER OF DIMENSIONS AND NUMBER OF VARIABLES
    CALL nf(nf90_Inquire(ncid, nDimensions=ndims, nVariables=nvar, &
                         unlimitedDimId=udimid))

    WRITE(*,*)
    WRITE(*,*) '---------------------------------------------------------'
    WRITE(*,*) 'Scanning file '//TRIM(file)
    WRITE(*,*) '... found ',ndims,' dimensions'
    WRITE(*,*) '... found ',nvar,' variables'
    WRITE(*,*) '... found unlimited dimension with ID ',udimid

    ! ALLOCATE MEMORY FOR DIM-related INFO
    ALLOCATE(dimlen(ndims),dimname(ndims))  ! name and length of dimension
    ALLOCATE(dimvid(ndims),dimtype(ndims))  ! var-ID and -type of dim.var
    ALLOCATE(axdim(ndims,ndims,ndims))      
                   ! dimension ,dimensionality = 1D 2D ..., position 
    axdim(:,:,:) = .false.
    ALLOCATE(newvid(ndims,ndims,ndims))     ! new variable dimensions
    newvid(:,:,:) = -1
    ALLOCATE(dxdimid(ndims))                ! IDs of new dimensions
    dxdimid(:) = -1

    ! GET DIMENSION LENGTHS, NAMES, ID AND TYPE OF CORRESPONDING VARIABLE
    WRITE(*,*)
    WRITE(*,*) 'Scanning dimensions ...'
    DO i=1, ndims
       CALL nf(nf90_Inquire_Dimension(ncid, i, name=dimname(i), &
                                      len=dimlen(i)))
       status = nf90_inq_varid(ncid, TRIM(dimname(i)), varid)
       ! special treatment for time axis (UNLIMITED DIMENSION)
       IF ((status == NF90_NOERR).AND.(i /= udimid)) THEN
          WRITE(*,*) '... dimension '''//TRIM(dimname(i))//''', ID: ',i
          CALL nf(nf90_Inquire_Variable(ncid, varid, xtype=dimtype(i)))
          dimvid(i) = varid
          WRITE(*,*) '    -> variable ID: ',varid
          WRITE(*,*) '    -> variable type: ',dimtype(i)
       ELSE
          WRITE(*,*) '... dimension '''//TRIM(dimname(i))//''', ID: ',i
          dimtype(i) = -1
          dimvid(i) = -1
          IF (status /= NF90_NOERR) THEN
             WRITE(*,*) '    -> no variable '''//TRIM(dimname(i))//''''
          ELSE
             WRITE(*,*) '    -> UNLIMITED DIMENSION !'
          END IF
       END IF
    END DO

    WRITE(*,*)
    WRITE(*,*) 'Scanning variables ...'

    ! LOOP OVER VARIABLES
    DO i=1, nvar
       ! GET variable name and dimensions
       CALL nf(nf90_Inquire_Variable(ncid, i, name=varname, ndims=nvdim))
       IF (nvdim > 0) THEN
          ALLOCATE(dimids(nvdim))
!          CALL nf(nf90_Inquire_Variable(ncid, i, name=varname, dimids=dimids))
          status = (nf90_Inquire_Variable(ncid, i, name=varname, dimids=dimids))
          WRITE(*,*) '... found variable ''',TRIM(varname),''' ',&
                     'with dim-ID(s):',dimids
       ELSE
          WRITE(*,*) '... found variable '''//TRIM(varname),&
                     ''' with 0 dimensions !'
       ENDIF

       ! INITIALIZATION
       fieldatt = TRIM(varname)//', scalar'   ! field - attribute
       posatt = ''                            ! positions - attribute
       nvsdim = 0                             ! number of 'spatial' dimensions

       ! 1'st LOOP OVER DIMENSIONS OF VARIABLE: count 'spatial' dimensions
       DO j=nvdim, 1, -1
          ! CHECK FOR UNLIMITED DIMENSION-ID
          IF (.NOT.((dimids(j)==udimid).AND.(udimid /= -1))) THEN
             nvsdim = nvsdim + 1                    ! 'spatial dimension'
          END IF
       END DO  ! 1'st LOOP OVER DIMENSIONS OF VARIABLE

       ! positions_attribute string component for this variable
       WRITE(dimstr,'(i3.3)') nvsdim

       ! 2'nd LOOP OVER DIMENSIONS OF VARIABLE
       nvsdimcnt = nvsdim+1
       DO j=nvdim, 1, -1
          ! CHECK FOR UNLIMITED DIMENSION-ID
          IF ((dimids(j)==udimid).AND.(udimid /= -1)) THEN
             fieldatt = TRIM(varname)//', scalar, series'
             WRITE(*,*) '     -> time series'
          ELSE
             nvsdimcnt = nvsdimcnt - 1
             WRITE(posstr,'(i3.3)') nvsdimcnt
             ! BUILD position attribute
             idx1 = LEN_TRIM(posatt)+1
             idx2 = idx1 + LEN_TRIM(dimname(dimids(j)))
             posatt(idx1:idx2) = TRIM(dimname(dimids(j)))
             idx1 = idx2
             idx2 = idx1
             posatt(idx1:idx2) = '_'
             idx1 = idx2+1
             idx2 = idx1+3
             posatt(idx1:idx2) = TRIM(dimstr)
             idx1 = idx2
             idx2 = idx1
             posatt(idx1:idx2) = '_'
             idx1 = idx2+1
             idx2 = idx1+3
             posatt(idx1:idx2) = TRIM(posstr)
             idx1 = idx2
             idx2 = idx1+10
             posatt(idx1:idx2) = ', product;'
          END IF

          ! REQUEST FOR 'nvsdim'-dimensional axis of this dimension ...
          ! ... only IF spatial DIMENSION > 1 ...
          ! ... AND valid dim. (e.g., not time dimension) ...
          ! current position of axis is 
          IF ((nvsdim > 1).AND.(dimtype(dimids(j)) /= -1)) THEN
             axdim(dimids(j),nvsdim,nvsdimcnt) = .true.
             ! ELSE
             ! 0D/1D variable OR only time dimension
          END IF

       END DO  ! 2'nd LOOP OVER DIMENSIONS OF VARIABLE
       
       ! CHECK NUMBER OF SPATIAL DIMENSIONS
       ! nvsdim = 0 -> 0D variable or only time dimension
       IF (nvsdim == 0) THEN
          posatt = ''
          fieldatt = ''
          ! HERE: RE-DEFINE VARIABLE WITH AT LEAST ONE SPATIAL DIMENSION !!!
       END IF

       ! nvsdim = 1 -> 1D variables: no positions attribute required
       IF (nvsdim == 1) THEN
          posatt = ''
       END IF

       ! ADD ATTRIBUTES
       IF (TRIM(fieldatt) /= '') THEN
          WRITE(*,*) '    adding attribute -> ',TRIM(varname),&
                     ':field = "',TRIM(fieldatt),'" ;'
          CALL nf(nf90_put_att(ncid, i, "field", TRIM(fieldatt)))
       END IF

       IF (TRIM(posatt) /= '') THEN
          WRITE(*,*) '    adding attribute -> ',TRIM(varname),&
                     ':positions = "',TRIM(posatt),'" ;'
          CALL nf(nf90_put_att(ncid, i, "positions", TRIM(posatt)))
       END IF

       ! CLEAN UP
       IF (nvdim > 0) THEN
          DEALLOCATE(dimids)
       ENDIF
    END DO  ! LOOP OVER VARIABLES

    ! DEFINE NEW DIMENSIONS
    WRITE(*,*)
    WRITE(*,*) 'Defining new dimensions ...'
    DO j=1, ndims     ! LOOP OVER DIMENSIONALITY
       IF (ANY(axdim(:,j,:))) THEN
          WRITE(dimstr,'(i3.3)') j
          WRITE(*,*) '... '//'''dim_'//TRIM(dimstr)//''''
          CALL nf(nf90_def_dim(ncid, 'dim_'//TRIM(dimstr), j, dxdimid(j)))
          WRITE(*,*) '     -> dim-ID: ',dxdimid(j)
       END IF       
    END DO

    ! DEFINE new VARIABLES for AXES
    WRITE(*,*)
    WRITE(*,*) 'Defining new variables ...'
    DO i=1, ndims      ! LOOP OVER DIMENSIONS
       DO j=1,ndims    ! LOOP OVER DIMENSIONALITY
          WRITE(dimstr,'(i3.3)') j
          DO k=1,ndims ! LOOP OVER POSITIONS
             IF (axdim(i,j,k)) THEN
                WRITE(posstr,'(i3.3)') k
                WRITE(*,*) '... '''//TRIM(dimname(i))//'_'//TRIM(dimstr)&
                          &//'_'//TRIM(posstr)//''''
                CALL nf(nf90_def_var(ncid,                                   &
                     TRIM(dimname(i))//'_'//TRIM(dimstr)//'_'//TRIM(posstr), &
                            NF90_FLOAT, (/ dxdimid(j), i /), newvid(i,j,k)))
                WRITE(*,*) '     -> var-ID : ',newvid(i,j,k)
                WRITE(*,*) '     -> dim-IDs: ',dxdimid(j), i
             END IF
          END DO
       END DO
    END DO

    ! WRITE DATA TO NEW VARIABLES
    ! data-mode
    CALL nf(nf90_enddef(ncid))
    WRITE(*,*)
    WRITE(*,*) 'Writing new variables ...'
    
    DO i=1, ndims         ! LOOP OVER DIMENSIONS
       DO j=1, ndims      ! LOOP OVER DIMENSIONALITY
          DO k=1, ndims   ! LOOP OVER POSITION
             IF (axdim(i,j,k)) THEN
                WRITE(dimstr,'(i3.3)') j
                WRITE(posstr,'(i3.3)') k
                WRITE(*,*) '... '''//TRIM(dimname(i))//&
                          &'_'//TRIM(dimstr)//'_'//TRIM(posstr)//''''
                SELECT CASE (dimtype(i))
                   CASE (NF90_FLOAT)
                      ALLOCATE(dr(dimlen(i)))
                      CALL nf(nf90_get_var(ncid, dimvid(i), dr))
                      DO dcount=1, j
                         IF (dcount == k) THEN
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dr                                &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         ELSE
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dr*0.0                            &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         END IF
                      END DO
                      DEALLOCATE(dr)
                   CASE (NF90_BYTE, NF90_CHAR)
                      ALLOCATE(dc(dimlen(i)),dch(dimlen(i)))
                      dch(:) = ''
                      CALL nf(nf90_get_var(ncid, dimvid(i), dc))
                      DO dcount=1, j
                         IF (dcount == k) THEN
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dc                                &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         ELSE
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dch                               &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         END IF
                      END DO
                      DEALLOCATE(dc)                      
                   CASE (NF90_SHORT, NF90_INT)
                      ALLOCATE(di(dimlen(i)))
                      CALL nf(nf90_get_var(ncid, dimvid(i), di))
                      DO dcount=1, j
                         IF (dcount == k) THEN
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,di                                &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         ELSE
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,di*0                              &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         END IF
                      END DO
                      DEALLOCATE(di)
                   CASE (NF90_DOUBLE)
                      ALLOCATE(dd(dimlen(i)))
                      CALL nf(nf90_get_var(ncid, dimvid(i), dd))
                      DO dcount=1, j
                         IF (dcount == k) THEN
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dd                                &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         ELSE
                            CALL nf(nf90_put_var(ncid, newvid(i,j,k)   &
                                    ,dd*0.0                            &
                                    ,start=(/dcount, 1/)               &
                                    ,count=(/1, dimlen(i)/) ))
                         END IF
                      END DO
                      DEALLOCATE(dd)
                   CASE DEFAULT              ! == -1 
                END SELECT
             END IF
          END DO
       END DO
    END DO


    ! redef-mode
    CALL nf(nf90_redef(ncid))

    WRITE(*,*)
    WRITE(*,*) '... adding global attributes'

    ! ADD GLOBAL ATTRIBUTES
    WRITE(*,*) ':ncdx = "modified for openDX import by ncdx version 3.8" ;'
    CALL nf(nf90_put_att(ncid, NF90_GLOBAL, "ncdx", &
            "modified for openDX import by ncdx version 3.8"))
    CALL DATE_AND_TIME(date, time, zone)
    WRITE(*,*) ':ncdx_date = "',date,'" ;'
    CALL nf(nf90_put_att(ncid, NF90_GLOBAL, "ncdx_date", date))
    WRITE(*,*) ':ncdx_time = "',TRIM(time//zone),'"'
    CALL nf(nf90_put_att(ncid, NF90_GLOBAL, "ncdx_time", TRIM(time//zone)))

    ! CLOSE netCDF file
    CALL nf(nf90_close(ncid))

    ! CLEAN UP
    DEALLOCATE(dimlen, dimname, dimtype, dimvid)
    DEALLOCATE(axdim, dxdimid, newvid)

    WRITE(*,*)
    WRITE(*,*) '... done!'
    WRITE(*,*) '---------------------------------------------------------'
    
CONTAINS

! ------------------------------------------------------------------
! SUBROUTINE nf
!   netCDF error message
! INPUT  : STATUS
! OUTPUT : ---
!
! Author: Rolf Sander, MPICH, SEP 2001

SUBROUTINE nf(status)
  INTEGER, INTENT(IN):: status
  IF (status /= NF90_NOERR) THEN
     WRITE(*,*) 'netCDF-ERROR: ',nf90_strerror(status)
     STOP 'STOPPED !'
  ENDIF
END SUBROUTINE nf
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! SUBROUTINE USAGE
!   usage message
! INPUT  : program name
! OUTPUT : ---
!
! Author: Patrick Joeckel, MPICH, Dec 2001

  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <netCDF-filename>'
    WRITE(*,*)
  END SUBROUTINE USAGE

! ------------------------------------------------------------------

END program ncdx
! ------------------------------------------------------------------
! ------------------------------------------------------------------
