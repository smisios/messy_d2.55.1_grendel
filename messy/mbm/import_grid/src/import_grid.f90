! ******************************************************************
! ------------------------------------------------------------------
PROGRAM IMPORT_GRID
! ------------------------------------------------------------------
! Author: Astrid Kerkweg, Uni Mainz, July 2013 
!  based on NCREGRID by Patrick Joeckel, MPICH, Mainz, June 2002
! ******************************************************************

  USE mo_f2kcli                    ! command line interface
  USE messy_main_import_grid,     ONLY: READ_CONTROL, RG_CTRL,  RG_NML         &
                                      , RG_STATUS, RGSTAT_STOP, RG_PROC        &
                                      , RG_SCAN, NML_NEXT, NML_STAY, RG_STOP
  USE messy_main_grid,            ONLY: t_geohybgrid, INIT_GEOHYBGRID          &
                                      , NEW_GEOHYBGRID, COPY_GEOHYBGRID        &
                                      , EXPORT_GEOHYBGRID, PRINT_GEOHYBGRID
  USE messy_main_grid_netcdf,     ONLY: t_ncvar, QDEF_NCVAR, INIT_NCVAR        &
                                      , COPY_NCVAR, EXPORT_NCVAR               &
                                      , RGMLE, RGMLWC, ERRMSG, GRD_MAXSTRLEN   &
                                      , NULL_XTYPE, RGMLI, RGMLW, RGMSG        &
                                      , DOUBLE_NARRAY, POSITION, ELEMENT       &
                                      , INIT_NARRAY, NF90_FLOAT, NULL_DIMID    &
                                      , NULL_VARID, VTYPE_REAL, VTYPE_DOUBLE   &
                                      , INIT_NCDIM, COPY_NCDIM                 &
                                      , t_ncatt, EXPORT_NCATT, INIT_NCATT      &
                                      , IMPORT_NCATT, NF90_GLOBAL, NF90_CHAR   &
                                      , VTYPE_CHAR, VTYPE_UNDEF, CAT_NARRAY    &
                                      , MAIN_GRID_SET_MESSAGEMODE              &
                                      , MSGMODE_VL, MSGMODE_S, MSGMODE_E       &
                                      , MSGMODE_W, MSGMODE_VM, MSGMODE_I
  USE messy_main_grid_trafo,      ONLY: GTRF_NONE, GTRF_NRGD, GTRF_SCRP        &
                                      , BALANCE_GEOHYBGRID, COMPLETE_GEOHYBGRID&
                                      , GORD_LON, CHECK_NCVAR_ON_GEOHYBGRID    &
                                      , EXPAND_INPUT_GRID_LON
  USE messy_main_grid_trafo_nrgd, ONLY: REGRID_CONTROL, ncregridvers
  USE messy_main_grid_trafo_nrgd_base, ONLY: SNREGRID
  USE messy_main_grid_trafo_scrp, ONLY: t_scrip_data, SCRIP_CONTROL            &
                                      , CALC_SCRIP_WEIGHTS, CALC_SCRIPDATA     &
                                      , CONSTRUCT_VERTICAL_AXIS                &
                                      , INTERPOL_GEOHYBGRID_PS
  USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT                         &
                                      , RGTOOL_CONVERT_DAT2VAR
  USE messy_main_tools,           ONLY: int2str, read_nml_open, read_nml_check &
                                      , read_nml_close
  USE messy_main_constants_mem,   ONLY: dp, STRLEN_ULONG

#ifdef PNCREGRID
  USE messy_main_import_grid_par, ONLY: INIT_PARALLEL
#endif

  IMPLICIT NONE

  INTRINSIC :: TRIM

#ifdef PNCREGRID
#include <mpif.h>
#endif
 
  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

  CHARACTER(LEN=*), PARAMETER :: substr = 'import_grid'

  LOGICAL,  SAVE                        :: lint    ! use input time ?
  INTEGER, DIMENSION(:),  POINTER       :: RGT => NULL() ! regridding type
  TYPE (t_ncvar), DIMENSION(:), POINTER :: rvar  => NULL() ! list of variables
  TYPE (t_ncvar), DIMENSION(:), POINTER :: ovar  => NULL() ! list of variables
  TYPE (t_ncvar), DIMENSION(:), POINTER :: sovar => NULL() ! list of variables
  TYPE (t_ncvar), DIMENSION(:), POINTER :: sovartmp => NULL() ! list of variables
  TYPE (t_geohybgrid)                   :: igrid   ! output grid info
  TYPE (t_geohybgrid)                   :: ogrid   ! output grid info
  TYPE (t_geohybgrid)                   :: ogridtmp  ! output grid info
  TYPE (t_geohybgrid)                   :: intgrid ! intermediate grid info
  TYPE (t_geohybgrid)                   :: int2grid! intermediate grid info
  TYPE (t_geohybgrid)                   :: ipsgrid ! intermediate grid info
  TYPE (t_geohybgrid)                   :: intzogrid  ! output grid info
  INTEGER                               :: tt      ! time step
  INTEGER                               :: count = 0 ! procedure counter
  INTEGER                               :: i       ! counter
  TYPE(t_scrip_data), POINTER           :: PSD
  INTEGER                               :: igridid
  INTEGER                               :: ogridid
  CHARACTER(LEN=4)                      :: cstr = ''
  INTEGER                               :: status
  INTEGER                               :: SCRIP_ID     ! SCRIP DATA ID   
  LOGICAL                               :: lrgz         ! regrid in z direction
  INTEGER                               :: iipol
  LOGICAL                               :: lex
  INTEGER                               :: fstat
  INTEGER                               :: tc, tmin, tmax, tret, tstep
  INTEGER,           PARAMETER          :: t_undef  = -99
  INTEGER                               :: ix
  LOGICAL                               :: lrgxy
  LOGICAL                               :: lpres
  LOGICAL                               :: lwork
  LOGICAL                               :: lpar
  INTEGER                               :: iloc = 0

  REAL(DP), DIMENSION(:,:,:,:), POINTER :: dat => NULL()
  TYPE (t_ncvar), DIMENSION(:), POINTER :: var => NULL() ! list of variables
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: vdat => NULL()
  REAL(dp), DIMENSION(:), ALLOCATABLE   :: help
  REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_in  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER   :: vax_out => NULL()
  INTEGER                               :: xsize, ysize
  LOGICAL                               :: lpresax
  LOGICAL                               :: lrot = .FALSE.
  INTEGER                               :: nvar, nx, ny, nn, nz
  INTEGER                               :: zdim
  TYPE (t_ncatt), SAVE                  :: rggatt   ! RG global attribute
  CHARACTER(LEN=GRD_MAXSTRLEN)          :: timename
  CHARACTER(LEN=STRLEN_ULONG)           :: infostr
  INTEGER                               :: cntrgd
  INTEGER                               :: cntt

#ifdef PNCREGRID
  INTEGER         :: ierr, root
  INTEGER         :: rank, nproc
#endif

  NAMELIST /GTRF/ iipol

#ifdef PNCREGRID
  ! MPI Setup
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank,  ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
 
  CALL INIT_PARALLEL(rank, nproc, MPI_COMM_WORLD)
  lpar = .TRUE.
  write (*,*) 'SPECIES PARALLELSATION ON', nproc
#else
  lpar = .FALSE.
  write (*,*) 'SPECIES PARALLELSATION OFF'
#endif

  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name

  IF (NARG >= 2) THEN 
     WRITE(*,*) 'Too many arguments !'
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  IF (NARG == 0) THEN 
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1,CMD)  

  ! SET DEFAULT
  iipol = 1

  INQUIRE(file=TRIM(CMD), exist=lex)
  IF (.not.lex)  THEN  ! CMD does not exist
     WRITE(*,*) '*** WARNING: FILE '''//TRIM(CMD)//'''  NOT FOUND !'
     STOP    
  END IF

  ! OPEN FILE
  OPEN(100,file=TRIM(CMD))
  WRITE(*,*) 'Reading namelist GTRF'
  REWIND(100)
  READ(100, NML=GTRF, IOSTAT=fstat)
  IF (fstat /= 0) THEN
     WRITE(*,*) '*** ERROR: READ ERROR in NAMELIST GTRF in FILE '''//TRIM(CMD)//''' !', fstat
  ELSE
     WRITE(*,*) ' ... OK !'
  END IF
  CLOSE(100)

  IF (iipol == GTRF_SCRP) THEN
     write (*,*) '*** INFO: SCRIP INTERPOLATION (iipol =',iipol,')'
  ELSE IF (iipol == GTRF_NRGD) THEN
     write (*,*) '*** INFO: NCREGRID INTERPOLATION (iipol =',iipol,')'
  ELSE IF (iipol == GTRF_NONE) THEN
     write (*,*) '*** INFO: NO INTERPOLATION (iipol =',iipol,')'
     STOP
  ELSE
     write (*,*) ' WRONG INTERPOLATION OPTION: ', iipol
     STOP
  ENDIF

  ! REGRIDDING ...
  RG_CTRL = RG_SCAN
  RG_NML  = NML_NEXT
  cntrgd  = 0
  regrid_loop: DO ! endless DO loop (must be terminated with EXIT)

     cntrgd = cntrgd + 1
     
     CALL INIT_GEOHYBGRID(igrid)
     CALL INIT_GEOHYBGRID(ogrid)
     CALL INIT_GEOHYBGRID(ogridtmp)
     igridid = -99
     ogridid = -99

     NULLIFY(rvar)

     tmin  = t_undef
     tmax  = t_undef
     tstep = t_undef
     tret  = t_undef

     CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS      &
            , rvar, igrid, ogrid, RGT, lint, TRIM(CMD) &
            , infostr                                  &
            , tc=1, tmin=tmin, tmax=tmax, tstep=tstep  &
            , iounit=100, lvarparallel=lpar            &
            , lpresaxis=lpres, lwork=lwork )

     IF (RG_STATUS == RGSTAT_STOP) EXIT ! leave endless DO loop
     infostr= ' '

     time_loop: DO tc=tmin,tmax,tstep
        RG_CTRL=RG_STOP
        RG_NML=NML_STAY
        infostr= ' '
        ! CLOSE NAMELIST
        CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS    &
             , rvar, igrid, ogrid, RGT, lint, TRIM(CMD) &
             , infostr                                  &
             , tc=tc ,  iounit=100, lvarparallel=lpar   )

        cntt = 0
        DO 
           RG_CTRL=RG_SCAN 
           RG_NML=NML_NEXT
           infostr= ' '
           cntt = cntt + 1
           CALL READ_CONTROL(RG_CTRL, RG_NML, RG_STATUS    &
                , rvar, igrid, ogrid, RGT, lint, TRIM(CMD) &
                , infostr                                  &
                , tc=tc ,  iounit=100, lvarparallel=lpar   )
        
           IF (RG_STATUS == RGSTAT_STOP) EXIT regrid_loop ! leave endless DO loop
           IF (cntrgd == cntt) EXIT
        END DO
        RG_CTRL=RG_PROC
        RG_NML=NML_STAY

        SELECT CASE(iipol)
        CASE(GTRF_NONE)
           ! NOT IMPLEMENTED FOR BOX MODEL
        CASE(GTRF_NRGD)

           IF (lwork) &
           CALL REGRID_CONTROL(igrid, ogrid, rvar, ovar, RGT, lint &
                , lfirsto=(tc==tmin)) 
           
        CASE(GTRF_SCRP)
           CALL COMPLETE_GEOHYBGRID(igrid)
           CALL COMPLETE_GEOHYBGRID(ogrid)

           CALL NEW_GEOHYBGRID(status, igridid, igrid) 
           IF (status /= 0 .AND. status /= 1) THEN
              CALL RGMSG(substr, RGMLE &
                   , 'NEW_GEOHYBGRID status ', status &
                   ,' for input grid', .TRUE.)
              EXIT
           END IF
           CALL NEW_GEOHYBGRID(status, ogridid, ogrid) 
           IF (status /= 0 .AND. status /= 1)  THEN
              CALL RGMSG(substr, RGMLE &
                   , 'NEW_GEOHYBGRID status ', status &
                   ,' for output grid', .TRUE.)
              EXIT
           END IF
           
           lrgxy = QDEF_NCVAR(ogrid%lonm)  .OR. QDEF_NCVAR(ogrid%loni)  .OR. &
                   QDEF_NCVAR(ogrid%clonm) .OR. QDEF_NCVAR(ogrid%cloni) .OR. &
                   QDEF_NCVAR(ogrid%rlonm) .OR. QDEF_NCVAR(ogrid%rloni) .OR. &
                   QDEF_NCVAR(ogrid%ulonm) .OR. QDEF_NCVAR(ogrid%uloni) .OR. &
                   QDEF_NCVAR(ogrid%latm)  .OR. QDEF_NCVAR(ogrid%lati)  .OR. &
                   QDEF_NCVAR(ogrid%clatm) .OR. QDEF_NCVAR(ogrid%clati) .OR. &
                   QDEF_NCVAR(ogrid%rlatm) .OR. QDEF_NCVAR(ogrid%rlati) .OR. &
                   QDEF_NCVAR(ogrid%ulatm) .OR. QDEF_NCVAR(ogrid%ulati)

           lrgz  = QDEF_NCVAR(ogrid%hybm)  .OR. QDEF_NCVAR(ogrid%hybi)  .OR. &
                   QDEF_NCVAR(ogrid%hyam)  .OR. QDEF_NCVAR(ogrid%hyai)

           write (*,*) ' Horizontal / Vertical Interpolation ?', lrgxy, lrgz
              
           IF (.NOT. lrgxy .AND. .NOT. lrgz) THEN
              write (*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              write (*,*) '   NEITHER horizontal                        '
              write (*,*) '      NOR  vertical interpolation            ' 
              write (*,*) '              requested:                     '
              write (*,*) '                CYCLE                        '
              write (*,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              EXIT
           ENDIF

           CALL COPY_GEOHYBGRID(ogridtmp, ogrid)
           IF (lrgz .AND. .NOT. lrgxy) THEN
              CALL BALANCE_GEOHYBGRID(igrid, ogridtmp)
           ENDIF
           IF ((igrid%lonm%xtype /= NULL_XTYPE).AND. &
                  (igrid%lonm%dim(1)%len==1)) THEN
              CALL EXPAND_INPUT_GRID_LON(igrid,rvar)
           END IF

           IF (lrgxy .OR. lrgz) THEN
              CALL CALC_SCRIPDATA(status, igrid, ogridtmp, RGT, SCRIP_ID &
                   , PSD=PSD)
              IF (status /= 0 ) THEN
                 IF (status /= 01 )  THEN ! status == 01 : PSD exists already
                    CALL RGMSG(substr, RGMLE, 'calc_scrip_data ', status &
                         ,' ', .TRUE.)
                    EXIT
                 ENDIF
              ELSE
                 ! CALCULATE WEIGHTS
                 CALL CALC_SCRIP_WEIGHTS(status, PSD)
                 IF (status /= 0) THEN
                    CALL RGMSG(substr, RGMLE, 'calc_scrip_weights ' &
                         , status ,' ', .TRUE.)
                    EXIT
                 END IF
              END IF
           ENDIF

           IF (lrgxy) THEN
              CALL SCRIP_CONTROL(status, SCRIP_ID, igrid, ogrid, RGT, lint&
                   , rvar, sovar,intgrid, llrgz=lrgz, lfirsto=(tc==tmin))
              IF (status /= 0) THEN
                 CALL RGMSG(substr, RGMLE, 'SCRIP_CONTROL ', status &
                      ,' ', .TRUE.)
                 EXIT
              END IF
           ELSE
              ALLOCATE(sovar(SIZE(rvar)))
              DO ix = 1, SIZE(rvar)
                 CALL COPY_NCVAR(sovar(ix), rvar(ix))
              END DO
              CALL COPY_GEOHYBGRID(intgrid, igrid)
           ENDIF

           ! vertical interpolation
           IF (lrgz) THEN

              CALL RGMSG(substr, RGMLI, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
              CALL RGMSG(substr, RGMLI, 'START VERTICAL REGRIDDING >>>')

              ! CONVERT VARIABLE into 4D SPACE
              CALL RGMSG(substr, RGMLI, 'CONVERT SOVAR for hori. dims')
              CALL RGTOOL_CONVERT(sovar(1), dat, intgrid,order='xyzn') 
              xsize = SIZE(dat,1)
              ysize = SIZE(dat,2)
              DEALLOCATE(dat)

              ! a DEFINITION OF VERTICAL AXIS: IN-GRID
              CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL InGrid Axis')
              ! 
              CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_in &
                   , igrid, lpresax, SCRIP_ID, ogrid, RGT, lint)
              IF (status /= 0) &
                   CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,1)

              ! a DEFINITION OF VERTICAL AXIS: IN-GRID
              CALL RGMSG(substr, RGMLI, 'DEFINE VERTICAL OutGrid Axis')
              ! 
              ! CHECK FOR SPECIAL CASE, vertical axis defined with hybrid 
              ! coefficients without surface pressure => get surface pressure
              ! from input grid (if available)
              IF (.NOT.(QDEF_NCVAR(ogrid%pressi) &
                   .OR. QDEF_NCVAR(ogrid%pressm))) THEN
                 IF (QDEF_NCVAR(ogrid%hybi)  &
                      .AND.( .NOT. QDEF_NCVAR(ogrid%ps))) THEN
                    IF (QDEF_NCVAR(igrid%ps)) THEN
                        CALL INTERPOL_GEOHYBGRID_PS(status, igrid, ogrid &
                             , SCRIP_ID)
                    ELSE
                       CALL ERRMSG('WRONG INPUT FOR VERTICAL AXIS DEFINITION: '&
                         ,42,1)
                    END IF
                 END IF
              END IF
              CALL CONSTRUCT_VERTICAL_AXIS(status, xsize, ysize, vax_out &
                   , ogrid, lpresax)
              IF (status /= 0) &
                   CALL ERRMSG('CONSTRUCT_VERTICAL_AXIS: ' ,status,3)
              ! CHECK if both vertical axis are orientated in the same way
              ! assume axis orientation is equal for all columns, check only
              ! point (1,1)
              IF ( ( (vax_in(1,1,1)-vax_in(1,1,SIZE(vax_in,3))) * &
                   (vax_out(1,1,1)-vax_out(1,1,SIZE(vax_out,3)))) < 0) THEN
                 ! V-AXIS orientation differs
                 ! => AXIS ROTATION in VAX_IN and DAT required
                 lrot = .TRUE.
              ELSE
                 lrot = .FALSE.
              END IF
              ALLOCATE(var(SIZE(sovar)))
              numvar: DO nvar =1, SIZE(sovar)
                 ! CONVERT VARIABLE ON INTGRID TO 3D
                 ! NOTE: vertical coordinate is in intgrid 'n'
                 CALL RGTOOL_CONVERT(sovar(nvar), dat, intgrid,order='xyzn') 
                 ! ALLOCATE MEMORY FOR VERTICAL REMAPPED FIELD
                 ALLOCATE(vdat(SIZE(dat,1),SIZE(dat,2) &
                      ,SIZE(vax_out,3)-1,SIZE(dat,3)))
                 vdat = 0._dp
                 ! ROTATE VAX-IN and DATA IF REQUIRED
                 IF (lrot) THEN
                    zdim = SIZE(VAX_IN,3)
                    ALLOCATE(help(zdim))
                    DO nx = 1, SIZE(dat,1)
                       DO ny = 1, SIZE(dat,2)
                          help(:) = VAX_IN(nx,ny,:)
                          DO nz = 1, SIZE(help)
                             VAX_IN(nx,ny,nz) = help(zdim-nz+1)
                          END DO
                          DO nn = 1, SIZE(dat,3)
                             help(1:zdim-1) = dat(nx,ny,nn,:)
                             DO nz = 1, zdim-1
                                dat(nx,ny,nn,nz) = help(zdim-1-nz+1)
                             END DO
                          END DO
                       END DO
                    END DO
                    DEALLOCATE(help)
                 END IF
                 DO nx = 1, SIZE(dat,1)
                    DO ny = 1, SIZE(dat,2)
                       DO nn = 1, SIZE(dat,3)
                          ! 
                          CALL SNREGRID(status                             &
                               , vax_in(nx,ny,:), vax_out(nx,ny,:)         &
                               , dat(nx,ny,nn,:), vdat(nx,ny,:,nn), .FALSE.)
                          IF (status /= 0) &
                               CALL ERRMSG('SNREGRID ERROR: ' ,status,24)
                       END DO
                    END DO
                 END DO
                 ! convert vdat to var ..
                 CALL RGTOOL_CONVERT_DAT2VAR(var(nvar), vdat &
                      , var(nvar)%name, ogrid, 'xyzn')
                 var(nvar)%name = sovar(nvar)%name 
                 ! free memory
                 DEALLOCATE(dat)
                 DEALLOCATE(vdat)

              END DO numvar
              DEALLOCATE(vax_in)
              NULLIFY(vax_in)
              DEALLOCATE(vax_out)
              NULLIFY(vax_out)

              CALL MAIN_GRID_SET_MESSAGEMODE(MSGMODE_S + MSGMODE_E &
                   + MSGMODE_VL ) !&
!                   + MSGMODE_W + MSGMODE_VM + MSGMODE_I)
              ! EXPORT VARIABLES:
                 ! OUTPUT GRID AND VARIABLES
              IF (TRIM(ogrid%file) /= '') THEN
                 CALL RGMSG(substr, RGMLI,'  WRITE OUTFILE')
                 
                 CALL PRINT_GEOHYBGRID(ogrid, 'OUT ')
                 CALL EXPORT_GEOHYBGRID(ogrid)
                 ! EXPORT GLOBL ATTRIBUTE, ONLY AT FIRST TIME STEP
                 IF (tc==tmin) THEN
             
                    ! GET OLD ATTRIBUTE
                    CALL IMPORT_NCATT(rggatt, varid=NF90_GLOBAL   &
                         ,attname = 'RG_HISTORY'     &
                         ,file = TRIM(ogrid%file), lnostop=.true.)
                    ! CHECK IF EMPTY OR CHAR-ATT
                    IF ((rggatt%dat%type == VTYPE_UNDEF).OR.  &
                         (rggatt%dat%type == VTYPE_CHAR)) THEN
                       ! APPEND NEW DATA
                       CALL CAT_NARRAY(rggatt%dat, ogrid%att%dat)
                       rggatt%xtype = NF90_CHAR
                       rggatt%len   = rggatt%dat%dim(1)
                    ELSE  ! EXISTS WITH NON-CHAR TYPE
                       CALL RGMSG(substr, RGMLW, &
                            'ATTRIBUTE '''//TRIM(rggatt%name)//'''')
                       CALL RGMSG(substr, RGMLWC, &
                            'EXISTS ALREADY, BUT IS NOT OF TYPE CHAR !')
                    END IF
                    CALL EXPORT_NCATT(rggatt, file=TRIM(ogrid%file)  &
                         ,clobber=.true.)
                    ! EXPORT VARIABLES
                 END IF
                 DO i=1, SIZE(var)
                    DO ix = 1, var(i)%ndims
                       IF (var(i)%dim(ix)%fuid) THEN
                          timename = TRIM(var(i)%dim(ix)%name)
                       END IF
                       EXIT
                    END DO
                    var(i)%ustep = ogrid%t
                    CALL EXPORT_NCVAR(var(i), file=TRIM(ogrid%file))
                 END DO
                 CALL INIT_NCATT(rggatt)
              END IF 
           ENDIF
        END SELECT

        CALL INIT_GEOHYBGRID(intgrid)
        
        IF (ASSOCIATED(rvar)) THEN
           DO i=1,SIZE(rvar)
              CALL INIT_NCVAR(rvar(i))
           END DO
           DEALLOCATE(rvar)
           NULLIFY(rvar)
        ENDIF

        IF (ASSOCIATED(ovar)) THEN
           DO i=1, SIZE(ovar)
              CALL INIT_NCVAR(ovar(i))
           END DO
           DEALLOCATE(ovar)
           NULLIFY(ovar) 
        ENDIF
        
        IF (ASSOCIATED(sovar)) THEN
           DO i=1, SIZE(sovar)
              CALL INIT_NCVAR(sovar(i))
           END DO
           DEALLOCATE(sovar)
           NULLIFY(sovar) 
        ENDIF
        
     END DO time_loop
     RG_CTRL = RG_SCAN
     RG_NML  = NML_NEXT
        
  END DO regrid_loop
  ! END OF REGRIDDING

#ifdef PNCREGRID
  CALL MPI_Finalize(ierr)
#endif

CONTAINS

  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'IMPORT_GRID Version 0.9b'
    WRITE(*,*) 'Author: Astrid Kerkweg, Uni Mainz, July 2013'
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'based on NCREGRID Version ',NCREGRIDVERS
    WRITE(*,*) 'Author: Patrick Joeckel, MPICH, June 2002'
#ifdef PNCREGRID
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'MPI parallelisation (number of variables)'
!!$    WRITE(*,*) 'Author: Klaus Ketelsen, MPICH, Dec 2007'
#endif
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '--------------------------------------------'
  END SUBROUTINE USAGE

END PROGRAM IMPORT_GRID
! ------------------------------------------------------------------
