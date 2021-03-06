MODULE mo_nudging_sst
!BOP
  ! !MODULE: mo_nudging_sst (layer 2)

  ! !DESCRIPTION: 
  ! update SST field for nudging

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI Hamburg, April-2001
  ! R. Johanni, IPP Garching, May-2002, parallel version
  ! I. Kirchner, MPI Hamburg, Aug-2002, revision

  ! !LAST CHANGE: 
  ! Joachim Buchholz, MPI Mainz, 29. July 2003

  ! !USES:
  USE mo_kind,              ONLY: dp
  USE mo_decomposition,     ONLY: ldc => local_decomposition, &
                                  gdc => global_decomposition
  ! op_pj_20131112+
  USE mo_nudging_constants, ONLY: inudgformat
  USE mo_io,            ONLY: io_read, io_open, IO_close
  USE mo_netCDF,        ONLY: IO_inq_dimid, IO_INQ_DIMLEN,    &
                              IO_INQ_VARID,IO_GET_VAR_DOUBLE, &
                              IO_GET_VARA_DOUBLE, FILE_INFO,  &
                              NF_INQ_VARID, NF_NOERR, IO_info_construct
  USE mo_filename,      ONLY: NETCDF
  ! op_pj_20131112-
  ! op_pj_20180704+
  USE mo_netcdf,        ONLY: NF_MAX_NAME
  ! op_pj_20180704-

!BOX
  IMPLICIT NONE
  SAVE ! op_pj_20131115

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: NudgingReadSST
  PUBLIC :: NudgingSSTnew
  PUBLIC :: NudgingSSTClose

  ! !PUBLIC DATA MEMBERS:

  !* namelist parameters
  CHARACTER(len=256), PUBLIC :: ndg_file_sst        ! template name of sst file

  REAL(kind=dp), PUBLIC      :: ndg_freez  = 271.65_dp ! correction of freezing point
                                                    ! with ERA15  -1.5 C recommended
  INTEGER, PUBLIC            :: nsstinc = 24        !  hours between new SST fields
                                                    !  sst data given at every 24 hour
  INTEGER, PUBLIC            :: nsstoff = 12        ! offset in hours to 00UTC
!EOP

  ! local variables

  REAL(kind=dp), ALLOCATABLE         :: zts(:,:)  ! local buffer
  REAL(kind=dp), ALLOCATABLE, TARGET :: sstn(:,:) ! global buffer
  ! op_pj_20131112+
  TYPE (FILE_INFO)                   :: ndgfile    ! netcdf-file with nudg data 
  INTEGER                            :: ndtsc = 0
  INTEGER                            :: nts = 0
  INTEGER                            :: nvarid_sst = 0
  INTEGER                            :: nvarid_sic = 0
  REAL(dp), ALLOCATABLE              :: timevals(:)
  REAL(kind=dp), ALLOCATABLE         :: zsic(:,:)  ! local buffer
  REAL(kind=dp), ALLOCATABLE, TARGET :: sicn(:,:)  ! global buffer
  LOGICAL                            :: lsic = .FALSE.
  ! op_pj_20131112-

! op_pg_20120713+
!!$  INTEGER     :: sstblock = -1    ! number of sst file
!!$  INTEGER, SAVE  :: sstblock = -1    ! number of sst file
  INTEGER     :: sstblock = -1    ! number of sst file
! op_pg_20120713-
  INTEGER     :: kfiles
  LOGICAL     :: lsstn            ! true if new sst needed
  INTEGER     :: iheads           ! YYYYMMDDHH sst time step
! op_pj_20131115+
!!$  INTEGER, SAVE   :: ipos_sst     ! position in sst file
  INTEGER     :: ipos_sst         ! position in sst file
! op_pj_20131115-

  CHARACTER(len=256) :: mess

  ! op_pj_20180704+
  CHARACTER(LEN=NF_MAX_NAME) :: old_file = ''
  LOGICAL                    :: lreopen
  ! op_pj_20180704-

CONTAINS

  !======================================================================
!BOP
  !
  ! !IROUTINE:  NudgingReadSST
  ! !INTERFACE:

  SUBROUTINE NudgingReadSST

    ! !DESCRIPTION: 

    ! Update global SST field for nudging experiments. For first
    ! timestep, SST is taken from history files, whereas nudging
    ! fields will have to be read in immediately after model start.
    ! Update sst every time step in time step loop.

    ! !USES:
    !* general modules
    USE mo_mpi,               ONLY: p_parallel_io, p_io, p_bcast !!$&
!!$                                  , p_parallel ! op_pj_20180704
    USE mo_time_control,      ONLY: &
         get_date_components, next_date, current_date, &
         out_convert_date, nudg_sst_time, add_date, start_date
    USE mo_exception,         ONLY: message, finish
    USE mo_transpose,         ONLY: scatter_gp
    USE mo_filename,          ONLY: str_filter
    ! op_pg_20120724+
    USE mo_time_conversion, ONLY: time_days, TC_set, TC_get
    ! op_pg_20120724-

    !* nudging modules
    USE mo_nudging_constants, ONLY: lnudgcli, lnudg_run
#ifdef LITTLE_ENDIAN
    USE mo_nudging_utils,     ONLY: swap64, cpbread, WORD_LEN, HEAD_LEN
#else
    USE mo_nudging_utils,     ONLY: cpbread, WORD_LEN, HEAD_LEN
#endif
!EOP

    REAL(kind=dp), POINTER :: gl_sstn(:,:) !
    REAL(kind=dp), POINTER :: gl_sicn(:,:) ! op_pj_20131112
    REAL(kind=dp), ALLOCATABLE :: zhbuf(:)
    INTEGER :: ierr !mz_rs_20040326
    INTEGER       :: &
         krtim, iday, isec, kret, jr, iheads_need, iheads_new, &
         ihead(HEAD_LEN), iymd, ihms, yr, mo, dy, hr, mi, se
    LOGICAL            :: found
    CHARACTER (len=WORD_LEN) :: yhead(HEAD_LEN)
    CHARACTER(len=256) :: cfile
    TYPE(time_days) :: file_date ! op_pg_20120724
    TYPE(time_days) :: help_date
    LOGICAL         :: lstartup = .FALSE.
    ! op_pj_20131112+
    INTEGER                :: ndimid, nvarid
    INTEGER                :: ihead_nc(2)
    INTEGER                :: start(4), count(4)
    INTEGER                :: status
    ! op_pj_20131112-
    INTEGER                :: ncheck

    ! mz_jb_20040610+
    EXTERNAL pbopen, pbread, pbseek, &
             util_cray2ieee, util_i8toi4
    INTRINSIC allocated, mod, trim
    ! mz_jb_20040610-

    IF (.NOT. lnudg_run .OR. nsstinc == 0) RETURN

    NULLIFY(gl_sstn)
    NULLIFY(gl_sicn) ! op_pj_20131112

    IF (.NOT. ALLOCATED(zts)) THEN
       ALLOCATE(zts(ldc%nproma,ldc%ngpblks))
       CALL message('',' Nudging SST local memory initialized')
    END IF

    IF (.NOT.ALLOCATED(sstn)) THEN
       ALLOCATE(sstn(ldc%nlon,ldc%nlat))
       CALL message('',' Nudging SST global memory initialized')
    END IF

    ! op_pj_20131112+
    IF (inudgformat == 2) THEN
       IF (.NOT. ALLOCATED(zsic)) THEN
          ALLOCATE(zsic(ldc%nproma,ldc%ngpblks))
          CALL message('',' Nudging SIC local memory initialized')
       END IF
       
       IF (.NOT.ALLOCATED(sicn)) THEN
          ALLOCATE(sicn(ldc%nlon,ldc%nlat))
          CALL message('',' Nudging SIC global memory initialized')
       END IF
    END IF
    ! op_pj_20131112-

    IF (p_parallel_io) THEN

! op_pg_20120713+
!!$       CALL get_date_components(next_date,yr,mo,dy,hr,mi,se)
       CALL TC_get(current_date, dy, se)
       CALL TC_set(dy,se,file_date)
       ! CALCULATE DIFFERENCE TO START DATE
       CALL TC_set(dy,se,help_date)
       CALL TC_get(start_date, dy, se)
       CALL add_date(-dy, -se, help_date)
       CALL TC_get(help_date, dy, se)
       ncheck = nsstoff
       IF (nsstoff == 0) ncheck=nsstinc
       IF (dy > 0 .or. se > ncheck*3600) THEN
          lstartup = .FALSE.
       ELSE
          lstartup = .TRUE.
       END IF
       ! enable start without sst of previous month
       IF (.NOT.lstartup) THEN
          CALL add_date(0,(-nsstinc+MOD(nsstoff,nsstinc))*3600,file_date)
       ENDIF
       CALL get_date_components(file_date,yr,mo,dy,hr,mi,se)
! op_pg_20120713-

       ! get actually date/time of NSTEP
       nudg_sst_time = current_date
       CALL out_convert_date(nudg_sst_time, iymd, ihms)

       ! find actually sst /date/time
       krtim = nsstoff + 24  ! in hours of the day
       DO
          krtim = krtim - nsstinc
          IF (krtim*10000 <= ihms) EXIT
       END DO
       IF (krtim < 0) THEN
          krtim = krtim + 24
          iday = -1
          isec = 0
          CALL add_date(iday, isec, nudg_sst_time)
          CALL out_convert_date(nudg_sst_time, iymd, ihms)
       END IF

       iheads_need = iymd*100 + krtim              ! YYYYMMDDHH
       IF (lnudgcli) iheads_need = MOD(iheads_need,1000000)    ! remove year

       ! op_pj_20131112+
       SELECT CASE (inudgformat)   

       CASE(0)       ! nudging data is in cray-format
       ! op_pj_20131112-
                    
       lsstn = .FALSE.

       IF (sstblock < 0) THEN
          ! call at first time, no initialization before
          sstblock = 1

          cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
          WRITE(mess,*) 'use SST file : ',TRIM(cfile)
          CALL message('',mess)
          INQUIRE(file=cfile,exist=found)
          IF (.NOT.found) &
               CALL finish('NudgingReadSST','Nudging SST data file not found.')
          CALL pbopen(kfiles,cfile,'r',kret)
          IF (kret/=0) CALL finish('NudgingReadSST','nudging SST files available?')

          ipos_sst   = 0
! op_pg_20120713+
!!$          iheads     = -99
          IF (lstartup) THEN
             iheads = -99
          ELSE
             iheads = -1
          ENDIF
! op_pg_20120713-
          lsstn      = .TRUE.
          
       ELSE IF (iheads < iheads_need) THEN
          ! the last record is older than needed, new sst will be read
          lsstn = .TRUE.

       ENDIF

       IF (lsstn) THEN

          DO
             ! search next sst data record

             CALL cpbread(kfiles,yhead,WORD_LEN*HEAD_LEN,kret)
             IF (kret/=WORD_LEN*HEAD_LEN) THEN
                ! open next SST data file

                CALL NudgingSSTClose(.false.) ! mz_jb_20060316
                sstblock = sstblock + 1
                IF ((sstblock > 3) .AND. lnudgcli) THEN
                   sstblock = 1
! mz_pj_20040322+
!                ELSE IF (sstblock == 5) THEN
!                   CALL finish('NudgingReadSST','Stop reading SST files.')
! mz_pj_20040322-
                END IF
                
                cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
                WRITE(mess,*) 'use SST file : ',TRIM(cfile)
                CALL message('',mess)
                INQUIRE(file=cfile,exist=found)
                IF (.NOT.found) &
                     CALL finish('NudgingReadSST','Nudging SST data file not found.')
                CALL pbopen(kfiles,cfile,'r',kret)

                ipos_sst = 0
                CYCLE
             ENDIF

             CALL util_i8toi4(yhead(1), ihead(1), HEAD_LEN)
             iheads_new = ihead(3)*100+ihead(4)
             IF (lnudgcli) iheads_new = MOD(iheads_new,1000000)  ! remove year

             IF (iheads_new == iheads_need) THEN                 ! next SST record found
                iheads = iheads_new
                EXIT

             ELSE IF ( (iheads_new > iheads_need) .AND. (iheads == -99)) THEN
                ! at initialization time the first SST may be not available
                CALL out_convert_date(current_date, iymd, ihms)
                iheads = iymd*100 + ihms/10000
                IF (lnudgcli) iheads = MOD(iheads,1000000)
                EXIT

             ELSE IF (iheads > iheads_need) THEN
                CALL finish('NudgingReadSST','sst date fault')

             END IF

             ! skip data set
             ipos_sst = ipos_sst + HEAD_LEN + ldc%nlon*ldc%nlat
             CALL pbseek(kfiles,ipos_sst*WORD_LEN,0,kret)
             WRITE (mess,*) 'skip record ',ihead
             CALL message('NudgingReadSST',mess)

          ENDDO

          ipos_sst = ipos_sst + HEAD_LEN + ldc%nlon*ldc%nlat
          WRITE (mess,*) 'USE SST record ',ihead
          CALL message('NudgingReadSST',mess)

          ! read new sst field
          ALLOCATE(zhbuf(ldc%nlat*ldc%nlon)); zhbuf(:) = 0.0_dp
          CALL message('NudgingReadSST',' Attention convert SST DATA to IEEE')

          CALL pbread(kfiles,zhbuf(1),ldc%nlon*ldc%nlat*WORD_LEN,kret)
#ifdef LITTLE_ENDIAN
          CALL swap64(zhbuf,ldc%nlon*ldc%nlat)
#endif
          CALL util_cray2ieee(zhbuf,sstn(1,1),ldc%nlon*ldc%nlat)

          DEALLOCATE (zhbuf, stat=ierr)
          IF (ierr /= 0) CALL finish('NudgingReadSST', &
               'Error while deallocating zhbuf.')

       ENDIF

    ! op_pj_20131112+
    CASE(2)        ! nudging data is in netcdf-format

       lsstn = .FALSE.

       IF (sstblock < 0) THEN
          ! call at first time, no initialization before
          sstblock = 1

          cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
          cfile = TRIM(cfile)//'.nc'
          WRITE(mess,*) 'use SST file : ',TRIM(cfile)
          CALL message('',mess)
          lreopen = (TRIM(cfile) == TRIM(old_file)) ! op_pj_20180704
          INQUIRE(file=cfile,exist=found)
          IF (.NOT.found) THEN
             CALL finish('NudgingReadSST','Nudging SST data file not found.')
          ELSE
             CALL IO_info_construct(ndgfile)
             ndgfile%format = NETCDF
             CALL IO_open (cfile, ndgfile, IO_READ)
          END IF

          ! op_pj_20180704+
          !!$IF (p_parallel) CALL p_bcast(lreopen, p_io)
          IF (lreopen) CALL finish('OpenOneBlock',&
               'file that has just been closed is reopened again')
          ! op_pj_20180704-

          CALL IO_INQ_DIMID(ndgfile%file_id, 'time', ndimid)
          CALL IO_INQ_DIMLEN(ndgfile%file_id, ndimid, nts)
          CALL IO_INQ_VARID(ndgfile%file_id, 'time', nvarid)
          IF (ALLOCATED(timevals)) DEALLOCATE(timevals)
          ALLOCATE (timevals(nts))
          CALL IO_GET_VAR_DOUBLE (ndgfile%file_id, nvarid, timevals)
          ndtsc=0 

          ! sst must be present
          CALL IO_INQ_VARID (ndgfile%file_id, 'stl1', nvarid_sst)
          ! sic optional
          lsic = .FALSE.
          status = NF_INQ_VARID (ndgfile%file_id, 'sic', nvarid_sic)
          IF (status == NF_NOERR) THEN
             lsic = .TRUE.
             WRITE(mess,*) 'use SIC from file : ',TRIM(cfile)
             CALL message('',mess)
          END IF

! op_pg_20120713+
!!$          iheads     = -99
          IF (lstartup) THEN
             iheads = -99
          ELSE
             iheads = -1
          ENDIF
! op_pg_20120713-
          lsstn      = .TRUE.
          
       ELSE IF (iheads < iheads_need) THEN
          ! the last record is older than needed, new sst will be read
          lsstn = .TRUE.

       ENDIF

       IF (lsstn) THEN

          DO

             ndtsc=ndtsc + 1 

             IF (ndtsc > nts) THEN

                ! close file
                CALL NudgingSSTClose(.false.)

                ! open next SST data file
                sstblock = sstblock + 1
                IF ((sstblock > 3) .AND. lnudgcli) THEN
                   sstblock = 1
! mz_pj_20040322+
!                ELSE IF (sstblock == 5) THEN
!                   CALL finish('NudgingReadSST','Stop reading SST files.')
! mz_pj_20040322-
                END IF
                
                cfile = str_filter(ndg_file_sst,yr,mo,dy,hr,mi,se,sstblock)
                cfile = TRIM(cfile)//'.nc'
                WRITE(mess,*) 'use SST file : ',TRIM(cfile)
                CALL message('',mess)
                INQUIRE(file=cfile,exist=found)
                IF (.NOT.found) &
                     CALL finish('NudgingReadSST','Nudging SST data file not found.')
                
                ndgfile%format = NETCDF
                CALL IO_open (cfile, ndgfile, IO_READ)

                CALL IO_INQ_DIMID(ndgfile%file_id, 'time', ndimid)
                CALL IO_INQ_DIMLEN(ndgfile%file_id, ndimid, nts)
                CALL IO_INQ_VARID(ndgfile%file_id, 'time', nvarid)
                IF (ALLOCATED(timevals)) DEALLOCATE(timevals)
                ALLOCATE (timevals(nts))
                CALL IO_GET_VAR_DOUBLE (ndgfile%file_id, nvarid, timevals)
                ndtsc=0 

                ! sst must be present
                CALL IO_INQ_VARID (ndgfile%file_id, 'stl1', nvarid_sst)
                ! sic optional
                lsic = .FALSE.
                status = NF_INQ_VARID (ndgfile%file_id, 'sic', nvarid_sic)
                IF (status == NF_NOERR) THEN 
                   lsic = .TRUE.
                   WRITE(mess,*) 'use SIC from file : ',TRIM(cfile)
                   CALL message('',mess)
                ENDIF

                CYCLE
             ENDIF

             ! search next sst data record
             ihead_nc(1) = FLOOR(timevals(ndtsc))                  ! YYYYMMDD
             ihead_nc(2) = INT((timevals(ndtsc)-ihead_nc(1))*24._dp) ! HH

             iheads_new = ihead_nc(1)*100+ihead_nc(2) ! YYYYMMDDHH
             IF (lnudgcli) iheads_new = MOD(iheads_new,1000000)  ! remove year

             IF (iheads_new == iheads_need) THEN                 ! next SST record found
                iheads = iheads_new
                EXIT

             ELSE IF ( (iheads_new > iheads_need) .AND. (iheads == -99)) THEN
                ! at initialization time the first SST may be not available
                CALL out_convert_date(current_date, iymd, ihms)
                iheads = iymd*100 + ihms/10000
                IF (lnudgcli) iheads = MOD(iheads,1000000)
                EXIT

             ELSE IF (iheads > iheads_need) THEN
                CALL finish('NudgingReadSST','sst date fault')

             END IF

             ! skip data set
             WRITE (mess,*) 'skip record ',ihead_nc,' (',ndtsc,' of ',nts,')'
             CALL message('NudgingReadSST',mess)

          ENDDO

          WRITE (mess,*) 'USE SST record ',ihead_nc,' (',ndtsc,' of ',nts,')'
          CALL message('NudgingReadSST',mess)

          ! read new sst field
          count(:) = (/ ldc%nlon, ldc%nlat, 1, 1 /)
          start(:) = (/ 1, 1, 1, ndtsc /)
          CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid_sst &
               , start, count, sstn(:,:))
         
          ! read new sic field (optional)
          IF (lsic) THEN
             count(:) = (/ ldc%nlon, ldc%nlat, 1, 1 /)
             start(:) = (/ 1, 1, 1, ndtsc /)
             CALL IO_GET_VARA_DOUBLE (ndgfile%file_id, nvarid_sic &
                  , start, count, sicn(:,:))
          ENDIF

       ENDIF

    CASE default
       
       WRITE (mess,*) 'inudgformat =', inudgformat, ' in ndgctl is not supported'
       CALL message('NudgingReadSST',mess)
       CALL finish('NudgingReadSST','Run terminated because of invalid value of inudgformat.')
     
    ! op_pj_20131112-

    END SELECT ! op_pj_20131112

       gl_sstn => sstn

    ENDIF

    CALL scatter_gp(gl_sstn,zts,gdc)

    ! op_pj_20131112+
    CALL p_bcast(lsic, p_io)
    IF (lsic) THEN
       gl_sicn => sicn
       CALL scatter_gp(gl_sicn,zsic,gdc)
    ENDIF
    ! op_pj_20131112-

  END SUBROUTINE NudgingReadSST

  !======================================================================
!BOP
  ! !IROUTINE:  NudgingSSTnew
  ! !INTERFACE:

  SUBROUTINE NudgingSSTnew(krow)

    ! !DESCRIPTION: 

    ! Replace the SST at every timestep by the field provided in the nudging 
    ! forcing files. Since no information about fractional seaice is provided 
    ! in this dataset, seaice=1 is set for SSTs lower or equal to the 
    ! freezing/melting temperature of seaice as it is done in earlier versions 
    ! of the model. Furthermore the new temperature replaces the ice temperature 
    ! tsi, while the water temperature is set to ctfreez (=271.38 K) (the latter 
    ! might not be necessary).
    ! 
    ! For SSTs larger than ctfreez, the water temperature tsw is replaced and 
    ! the seaice fraction is set to zero. The sea ice detection
    ! temperature can be modified using {\it NDG\_FREEZ}.

    ! !USES:
    USE mo_kind,              ONLY: dp
    USE mo_memory_g3b,        ONLY: tsi, tsw, slm, seaice
#ifndef MESSY
    USE mo_physc2,            ONLY: ctfreez
#else
    USE messy_main_constants_mem, ONLY: ctfreez
#endif
    USE mo_nudging_constants, ONLY: lnudg_run &
                                  , inudgformat     ! op_pj_20131112
    USE mo_exception,         ONLY: message
!EOP

    INTEGER                   :: krow
    INTEGER                  :: jrow, jp, kproma
    REAL(kind=dp), PARAMETER :: zdt=0.01_dp  ! correct temperature 
                                             ! near freezing point

    INTRINSIC MAX, MIN
    EXTERNAL clsst

    jrow = krow
    IF (nsstinc==0 .OR. .NOT. lnudg_run) THEN
      CALL clsst      ! use the standard sst for nudging

    ELSE               ! update sea surface temperatures at every time step

     ! op_pj_20131112+
     SELECT CASE (inudgformat)   

     CASE(0)       ! nudging data is in cray-format
     ! op_pj_20131112-

!!$   IF (jrow == 1) &
!!$      CALL message('NudgingSSTnew','Update sst from nudging data set')
         
      IF (jrow==ldc%ngpblks) THEN
        kproma = ldc%npromz
      ELSE
        kproma = ldc%nproma
      END IF
         
      DO jp =1,kproma
        IF (slm(jp,jrow) < 1.0_dp) THEN             ! sea fraction present

          IF (zts(jp,jrow) <= ndg_freez) THEN
            ! below freezing level
            seaice(jp,jrow) = 1.0_dp
            tsi(jp,jrow)  = MIN(zts(jp,jrow),ctfreez-zdt)
            tsw(jp,jrow)  = ctfreez

          ELSE
            ! above freezing level
            seaice(jp,jrow) = 0.0_dp
            tsi(jp,jrow)  = ctfreez
            tsw(jp,jrow)  = MAX(zts(jp,jrow),ctfreez+zdt)

          ENDIF

        ENDIF
      ENDDO

     ! op_pj_20131112+
     CASE(2)       ! nudging data is in netCDF format

      IF (jrow==ldc%ngpblks) THEN
        kproma = ldc%npromz
      ELSE
        kproma = ldc%nproma
      END IF
         
      IF (lsic) THEN
         ! if sea-ice fraction is present, use consistently with sst
!!$      IF (jrow == 1) &
!!$           CALL message('NudgingSSTnew' &
!!$           ,'Update sst and sic from nudging data set')
         seaice(1:kproma,jrow) = zsic(1:kproma,jrow)
         tsi(1:kproma,jrow) = zts(1:kproma,jrow)
         tsw(1:kproma,jrow) = zts(1:kproma,jrow)
      ELSE
         ! old setup as above (with IEEE-data)
         ! (can be switched off with nsstinc=0 in namelist)
         IF (jrow == 1) &
              CALL message('NudgingSSTnew' &
              ,'Update sst from nudging data set')
         DO jp =1,kproma
            IF (slm(jp,jrow) < 1.0_dp) THEN             ! sea fraction present

               IF (zts(jp,jrow) <= ndg_freez) THEN
                  ! below freezing level
                  seaice(jp,jrow) = 1.0_dp
                  tsi(jp,jrow)  = MIN(zts(jp,jrow),ctfreez-zdt)
                  tsw(jp,jrow)  = ctfreez
                  
               ELSE
                  ! above freezing level
                  seaice(jp,jrow) = 0.0_dp
                  tsi(jp,jrow)  = ctfreez
                  tsw(jp,jrow)  = MAX(zts(jp,jrow),ctfreez+zdt)
                  
               ENDIF
               
            ENDIF
         ENDDO
      ENDIF

     END SELECT
     ! op_pj_20131112-

    ENDIF

  END SUBROUTINE NudgingSSTnew

  !======================================================================
!BOP
  ! !IROUTINE:  NudgingSSTClose
  ! !INTERFACE:

  SUBROUTINE NudgingSSTClose(lmem)

    ! !DESCRIPTION: 
    ! close nudging sst file and deallocate memory (optional)
    
    USE mo_exception, ONLY: finish
    
    ! !INPUT PARAMETERS: 
    USE mo_mpi,            ONLY: p_parallel_io
    LOGICAL, INTENT(in) :: lmem
!EOP

    INTEGER :: kret, ierr

    EXTERNAL pbclose

    ! op_pj_20131115+
    SELECT CASE (inudgformat)   
    CASE(0)       ! nudging data is in cray-format
    ! op_pj_20131115-
       IF (p_parallel_io) CALL pbclose(kfiles,kret)
    ! op_pj_20131115+
    CASE(2)       ! nudging data is in netCDF format
       ! op_pj_20180704+
       ! save name of closed file
       old_file = ndgfile%file_name
       ! op_pj_20180704-
       IF (p_parallel_io) CALL IO_close(ndgfile)
       IF (ALLOCATED(timevals)) DEALLOCATE(timevals)
    END SELECT
    ! op_pj_20131115-

    IF (lmem) THEN
      DEALLOCATE(zts,sstn, stat=ierr)
      IF (ierr /= 0) CALL finish('NudgingSSTClose', &
           'Error while deallocating zts, sstn.')
      ! op_pj_20131112+
      IF (ALLOCATED(zsic)) THEN
         DEALLOCATE(zsic, stat=ierr)
         IF (ierr /= 0) CALL finish('NudgingSSTClose', &
              'Error while deallocating zsic.')
      END IF
      IF (ALLOCATED(sicn)) THEN
         DEALLOCATE(sicn, stat=ierr)
         IF (ierr /= 0) CALL finish('NudgingSSTClose', &
              'Error while deallocating sicn.')
      END IF
      ! op_pj_20131112-
    END IF

  END SUBROUTINE NudgingSSTClose

  !======================================================================

END MODULE mo_nudging_sst
