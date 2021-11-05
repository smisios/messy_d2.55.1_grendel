!+ 3DVAR/COSMO source module for feedback file interface to COSMO
!
! $Id: mo_fdbk_cosmo.f90,v 4.29 2013-10-04 08:51:39 for0adm Exp $
!-------------------------------------------------------------------------------

MODULE mo_fdbk_cosmo

!-------------------------------------------------------------------------------
!
! Description:
!   COSMO interface to write NetCDF feedobs (or feedback) file (FOF).
!   (Common format for FOF in 3DVAR and COSMO.)
!   This module is temporarily included in the 3DVAR source tree to check
!   consistent handling with 3D-VAR.
!
!   This module contains the following module procedures:
!     - write_report
!     - write_report_radar_1
!     - write_report_radar_2
!     - add_data : interface for: add_inte_vala, add_real_vala_1D,
!                                 add_text_vala, add_real_vala_2D
!     - fill_bodybuf_int, fill_bodybuf_real
!   It uses from:
!
! Current Code Owners:
!    For DWD 3DVAR:                        For COSMO:
!    DWD, Andreas Rhodin                   DWD, Christoph Schraff
!    phone: +49 69 8062 2722               phone: +49 69 8062 2725
!    fax:   +49 69 8062 3721               fax:   +49 69 8062 3721
!    email: andreas.rhodin@dwd.de          email: christoph.schraff@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_22        2012/01/31 Christoph Schraff
!  Initial release, based on original code from Marek Lazanowicz, IMGW Poland,
!  and introduced in 3DVAR version V1_1.
!  Code restructured, error handling adapted to COSMO needs.
!  Adapted to updates of feedback file routines (e.g. from 3DVAR, sun_zenith
!  etc.) and redefinitions of feedback file format, plus bug fixes.
! V4_26        2012/12/06 Andreas Messer
!  Modification for using also RTTOV10: introduction of 'sat_zenit'.
! V4_28        2013/07/12 Ulrich Blahak, Christoph Schraff
!  - Added routines 'write_report_radar_1', 'write_report_radar_2' to write
!    radar reports to feedback files efficiently.
!  - Added component 'spec_index' to type 't_acc_body'.
!  - 'veri_data' in 't_acc_body' with fixed size 1 instead of pointer, i.e.
!    simulated obs from at most 1 model run can be stored.
!
! CAUTION: This module is used commonly by the 3DVAR and COSMO source trees. !!!
!!!        Therefore, anybody wanting to introduce a modification to this    !!!
!!!        module in the context of either of these programs must consult    !!!
!!!        the 'current code owner' of this module for the other program,    !!!
!!!        in order to allow for checking that the modification will comply  !!!
!!!        with both program packages. This must be done before the          !!!
!!!        modification is put into the Version Control System (VCS).        !!!
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================

USE mo_fdbk,          ONLY: t_fdbk

USE mo_t_netcdf_file, ONLY: t_netcdf_file, nlen

USE mo_netcdf_param   ! provides access to the parameters defined
                      ! in 'netcdf.inc' of the NetCDF package

IMPLICIT  NONE
!================
! public entities
!================
PRIVATE
!-------------------------
! derived type definitions
!-------------------------
PUBLIC :: t_account
PUBLIC :: t_acc_header
PUBLIC :: t_acc_body
!------------
! subroutines
!------------
PUBLIC :: write_report, write_report_radar_1, write_report_radar_2

!Interface block
 INTERFACE add_data
   MODULE PROCEDURE              &
     add_inte_vala,              &
     add_real_vala_1D,           &
     add_real_vala_2D,           &
     add_text_vala
END INTERFACE


TYPE t_acc_header
  INTEGER           ::   i_body            ! body offset in netcdf file
  INTEGER           ::   l_body            ! number of body elements in report
  INTEGER           ::   n_level           ! numbers of level in report
  INTEGER           ::   data_category     !
  INTEGER           ::   sub_category      !
  INTEGER           ::   center            !
  INTEGER           ::   sub_center        !
  INTEGER           ::   obstype           !
  INTEGER           ::   codetype          !
  INTEGER           ::   ident             !
  CHARACTER(LEN=10) ::   statid            !
  REAL              ::   lat               !
  REAL              ::   lon               !
  REAL              ::   sun_zenit         !
  REAL              ::   sat_zenit         !
  INTEGER           ::   time              !
  INTEGER           ::   time_nomi         !
  INTEGER           ::   time_dbase        !
  INTEGER           ::   z_station         !
  INTEGER           ::   z_modsurf         !
  INTEGER           ::   r_state           !
  INTEGER           ::   r_flags           !
  INTEGER           ::   r_check           !
  INTEGER           ::   sta_corr          !
  INTEGER           ::   index_x           !
  INTEGER           ::   index_y           !
  INTEGER           ::   mdlsfc            !
  INTEGER           ::   instype           !
  INTEGER           ::   retrtype          !
  INTEGER           ::   phase             !
  INTEGER           ::   tracking          !
  INTEGER           ::   meas_type         !
  INTEGER           ::   rad_corr          !
  INTEGER           ::   surftype          !
  INTEGER           ::   flg_1dvar         !
  INTEGER           ::   flg_cld           !
  INTEGER           ::   obs_id            !
  INTEGER           ::   source            !
  INTEGER           ::   record            !
  INTEGER           ::   subset            !
  INTEGER           ::   dbkz              !
  INTEGER           ::   index_d           !
END TYPE t_acc_header

TYPE t_acc_body
  INTEGER           ::   varno             !
  REAL              ::   obs               !
  REAL              ::   bcor              !
  REAL              ::   e_o               !
  REAL              ::   level             !
  INTEGER           ::   level_typ         !
  INTEGER           ::   level_sig         !
  INTEGER           ::   state             !
  INTEGER           ::   flags             !
  INTEGER           ::   check             !
  INTEGER           ::   qual              !
  INTEGER           ::   spec_index        !
  REAL              ::   plevel            !
  REAL              ::   accuracy          !
  REAL              ::   w_qc              !
  REAL              ::   veri_data(1) !
  !   assume that length of 'veri_data' is <= 1
  !   i.e. simulated obs from at most 1 model run are stored
! REAL,    POINTER  ::   veri_data(:) => NULL() !
END TYPE t_acc_body

TYPE t_account
  INTEGER                         :: len                ! report length
  INTEGER                         :: offset             ! report offset in FOF
  TYPE ( t_acc_header )           :: header             ! report header
  TYPE ( t_acc_body )  , POINTER  :: body(:) => NULL()  ! report body
END TYPE t_account


CONTAINS

!===============================================================================

SUBROUTINE add_inte_vala ( file, ivala, in, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk),     INTENT(in)    :: file       ! meta file
INTEGER,          INTENT(in)    :: ivala(:)   ! 1D integer buffer to be stored
INTEGER,          INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_int ( file% nc% ncid,                &
                             file% nc% vars(in)% varid,     &
                             start,                         &
                             count,                         &
                             ivala(1:count)             )

END SUBROUTINE add_inte_vala

!===============================================================================

SUBROUTINE add_real_vala_1D (file, rvala1, in, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk),     INTENT(in)    :: file       ! meta file
REAL,             INTENT(in)    :: rvala1(:)  ! 1D real buffer to be stored
INTEGER,          INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_real (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             start,                       &
                             count,                       &
                             rvala1(1:count)              )

END SUBROUTINE add_real_vala_1D

!===============================================================================

SUBROUTINE add_real_vala_2D (file,  rvala2, in, nr, start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk)    , INTENT(in)    :: file         ! meta file
REAL,             INTENT(in)    :: rvala2(:,:)  ! 2D real buffer to be stored
INTEGER,          INTENT(in)    :: in        ,& ! variable number in meta file
                                                !   to be stored
                                   nr           ! position data in rvala2
INTEGER,          INTENT(in)    :: start        ! start position in netcdf file
INTEGER,          INTENT(in)    :: count        ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror       ! netcdf error
!-------------------------------------------------------------------------------

  ierror = nf_put_vara_real (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             (/start, nr/),               &
                             (/count,  1/),               &
                             rvala2(:,:)                )

END SUBROUTINE add_real_vala_2D

!===============================================================================

SUBROUTINE add_text_vala (file,  cvala, in,  start, count, ierror )

!-------------------------------------------------------------------------------
TYPE(t_fdbk)    , INTENT(in)    :: file       ! meta file
CHARACTER(LEN=*), INTENT(in)    :: cvala(:)   ! character buffer to be stored
INTEGER         , INTENT(in)    :: in         ! variable number in meta file
                                              !   to be stored
INTEGER,          INTENT(in)    :: start      ! start position in netcdf file
INTEGER,          INTENT(in)    :: count      ! number of data to be stored
INTEGER,          INTENT(inout) :: ierror     ! netcdf error
!-----------------------------------------------------------
INTEGER   ::  cn      ! string length to be stored
!-------------------------------------------------------------------------------
  cn = file% nc% vars(in)% p(1)% dim% len

  ierror = nf_put_vara_text (file% nc% ncid,              &
                             file% nc% vars(in)% varid,   &
                             (/1,    start/),             &
                             (/cn,   count/),             &
                             cvala(1:count)(1:cn)     )

END SUBROUTINE add_text_vala

!===============================================================================


SUBROUTINE write_report ( file, report, nrep, nbody, ihoff, iboff              &
                        , imdi, rmdich, jerr, yerr )
 
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------
  TYPE(t_fdbk),        INTENT (inout) :: file       ! metafile
  TYPE(t_account),     INTENT (in)    :: report(:)  ! reports to be stored
  INTEGER,             INTENT (in)    :: nrep       ! number of reports
  INTEGER,             INTENT (in)    :: nbody      ! number of observations
  INTEGER,             INTENT (in)    :: ihoff      ! header offset for time box
  INTEGER,             INTENT (in)    :: iboff      ! body offset for time box
  INTEGER,             INTENT (in)    :: imdi       ! missing data indicator for
                                                    !   integer (2^31-1)
  REAL   ,             INTENT (in)    :: rmdich     ! check value for missing
                                                    !   real data (-1.E30)
  INTEGER,             INTENT (inout) :: jerr       ! error status variable
  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= * ), INTENT (inout) :: yerr       ! error message

!---------------------------------------------------------
  INTEGER                       ::  &
    irep    ,& ! loop indices
    jj, iob ,& ! loop indices
    kcase   ,& ! type of variable (int/real, header/body, etc)
!   kd_hdr  ,& ! length of header dimension
!   kd_body ,& ! length of body   dimension
    varid   ,& ! variable id in netcdf file
    start   ,& ! start position in netcdf file
    count   ,& ! number of elements written
    nvdim   ,& ! nr of dimension veri_data field in meta file
    in      ,& ! variable number to be stored in meta file
    nobs    ,& ! number of observations (body length)
    kos     ,& ! index offset obs from current report
    kerr       ! error status variable

  CHARACTER (LEN=nlen) :: yvar    ! NetCDF variable name
  CHARACTER (LEN=nlen) :: ydim    ! NetCDF dimension name

  INTEGER ,           ALLOCATABLE :: ivalh(:) ,& ! buffer for int header elem.
                                     ivalb(:)    ! buffer for int body elements
  REAL    ,           ALLOCATABLE :: rvalh(:) ,& ! buffer for real header elem.
                                     rvalb(:)    ! buffer for real body   elem.
  REAL    ,           ALLOCATABLE :: rvala2(:,:) ! buffer for real body   elem.
  CHARACTER(len=100), ALLOCATABLE :: cvala(:)    ! buffer for char header elem.
!-------------------------------------------------------------------------------

!CS
! DO jj = 1 , file% nc% ndim
!   ydim = file% nc% dims(jj)% name
!   IF (ydim (1:LEN_TRIM(ydim)) == 'd_hdr' )  kd_hdr  = file% nc% dims(jj)% len
!   IF (ydim (1:LEN_TRIM(ydim)) == 'd_body')  kd_body = file% nc% dims(jj)% len
!   IF (file% nc% dims(jj)% name(1:8) == 'd_hdr   ')                           &
!     kd_hdr  = file% nc% dims(jj)% len
!   IF (file% nc% dims(jj)% name(1:8) == 'd_body  ')                           &
!     kd_body = file% nc% dims(jj)% len
! ENDDO
! WHERE (file% nc% dims% name == 'd_hdr' )  kd_hdr  = file% nc% dims(pos)% len
! WHERE (file% nc% dims% name == 'd_body')  kd_body = file% nc% dims(pos)% len
!CS
 
  jerr = 0
  kerr = 0
  IF (jerr == 0)  ALLOCATE ( ivalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1
  IF (jerr == 0)  ALLOCATE ( rvalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2
  IF (jerr == 0)  ALLOCATE ( cvala (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3
  IF (jerr == 0)  ALLOCATE ( ivalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 4
  IF (jerr == 0)  ALLOCATE ( rvalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 5
  IF (jerr /= 0) THEN
    IF (kerr == 1)  WRITE( yerr,'("ERROR in alloc: ivalh (nrep ) ",I8)' ) nrep
    IF (kerr == 2)  WRITE( yerr,'("ERROR in alloc: rvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 3)  WRITE( yerr,'("ERROR in alloc: cvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 4)  WRITE( yerr,'("ERROR in alloc: ivalb (nbody) ",I8)' ) nbody
    IF (kerr == 5)  WRITE( yerr,'("ERROR in alloc: rvalb (nbody) ",I8)' ) nbody
                                                                          RETURN
  ENDIF
  ivalh(:) = imdi
  rvalh(:) = rmdich *1.1
  cvala(:) = ''
  ivalb(:) = imdi
  rvalb(:) = rmdich *1.1

  kerr = nf_enddef(file% nc% ncid )
! IF (kerr /= NF_NOERR) THEN
!   jerr = kerr
!   WRITE( yerr,'("ERROR in nf_enddef in write_reports")' )
!                                                                         RETURN
! ENDIF

  DO in = 1, file% nc% nvar

    IF (.NOT. file% nc% vars(in)% opt_used)                                CYCLE
    kcase = 0

    yvar(:) = ' '
    yvar  = TRIM(ADJUSTL(file% nc% vars(in)% name))

! 1-dimensional arrays: integer or real, header or body length
! ------------------------------------------------------------

    IF (file% nc% vars(in)% nvdims == 1) THEN

      kcase = 1
      IF (     (file% nc% vars(in)% xtype == NF_FLOAT)                         &
          .OR. (file% nc% vars(in)% xtype == NF_FILL_DOUBLE))  kcase = 2
!     IF (file% nc% vars(in)% xtype == NF_CHAR)  kcase = 5
      ydim  =  file% nc% vars(in)% p(1)% dim% name
      IF (ydim (1:LEN_TRIM( ydim )) == 'd_body') THEN
        kcase = kcase + 2
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_veri') THEN
        kcase = 5
      ELSEIF (ydim (1:LEN_TRIM( ydim )) /= 'd_hdr') THEN
        kcase = -6
      ENDIF

      ! integer header elements
      IF (kcase == 1) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('i_body')
          ivalh (1:nrep) = report(1:nrep)% offset
        CASE ('l_body')
          ivalh (1:nrep) = report(1:nrep)% len
        CASE ('n_level')
          ivalh (1:nrep) = report(1:nrep)% header% n_level
        CASE ('data_category')
          ivalh (1:nrep) = report(1:nrep)% header% data_category
        CASE ('sub_category')
          ivalh (1:nrep) = report(1:nrep)% header% sub_category
        CASE ('center')
          ivalh (1:nrep) = report(1:nrep)% header% center
        CASE ('sub_center')
          ivalh (1:nrep) = report(1:nrep)% header% sub_center
        CASE ('obstype')
          ivalh (1:nrep) = report(1:nrep)% header% obstype
        CASE ('codetype')
          ivalh (1:nrep) = report(1:nrep)% header% codetype
        CASE ('ident')
          ivalh (1:nrep) = report(1:nrep)% header% ident
        CASE ('time')
          ivalh (1:nrep) = report(1:nrep)% header% time
        CASE ('time_nomi')
          ivalh (1:nrep) = report(1:nrep)% header% time_nomi
        CASE ('time_dbase')
          ivalh (1:nrep) = report(1:nrep)% header% time_dbase
        CASE ('z_station')
          ivalh (1:nrep) = report(1:nrep)% header% z_station
        CASE ('z_modsurf')
          ivalh (1:nrep) = report(1:nrep)% header% z_modsurf
        CASE ('r_state')
          ivalh (1:nrep) = report(1:nrep)% header% r_state
        CASE ('r_flags')
          ivalh (1:nrep) = report(1:nrep)% header% r_flags
        CASE ('r_check')
          ivalh (1:nrep) = report(1:nrep)% header% r_check
        CASE ('sta_corr')
          ivalh (1:nrep) = report(1:nrep)% header% sta_corr
        CASE ('mdlsfc')
          ivalh (1:nrep) = report(1:nrep)% header% mdlsfc
        CASE ('instype')
          ivalh (1:nrep) = report(1:nrep)% header% instype
        CASE ('retrtype')
          ivalh (1:nrep) = report(1:nrep)% header% retrtype
        CASE ('phase')
          ivalh (1:nrep) = report(1:nrep)% header% phase
        CASE ('tracking')
          ivalh (1:nrep) = report(1:nrep)% header% tracking
        CASE ('meas_type')
          ivalh (1:nrep) = report(1:nrep)% header% meas_type
        CASE ('rad_corr')
          ivalh (1:nrep) = report(1:nrep)% header% rad_corr
        CASE ('surftype')
          ivalh (1:nrep) = report(1:nrep)% header% surftype
        CASE ('flg_1dvar')
          ivalh (1:nrep) = report(1:nrep)% header% flg_1dvar
        CASE ('flg_cld')
          ivalh (1:nrep) = report(1:nrep)% header% flg_cld
        CASE ('index_x')
          ivalh (1:nrep) = report(1:nrep)% header% index_x
        CASE ('index_y')
          ivalh (1:nrep) = report(1:nrep)% header% index_y
        CASE ('obs_id')
          ivalh (1:nrep) = report(1:nrep)% header% obs_id
        CASE ('source')
          ivalh (1:nrep) = report(1:nrep)% header% source
        CASE ('record')
          ivalh (1:nrep) = report(1:nrep)% header% record
        CASE ('subset')
          ivalh (1:nrep) = report(1:nrep)% header% subset
        CASE ('dbkz')
          ivalh (1:nrep) = report(1:nrep)% header% dbkz
        CASE ('index_d')
          ivalh (1:nrep) = report(1:nrep)% header% index_d
        CASE DEFAULT
          ivalh (1:nrep) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (ivalh == imdi)  ivalh = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalh, in, start, count, jerr )
!       =============

      ! real header elements
      ELSEIF (kcase == 2) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('lat')
          rvalh(1:nrep) = report(1:nrep)% header% lat
        CASE ('lon')
          rvalh(1:nrep) = report(1:nrep)% header% lon
        CASE ('sun_zenit')
          rvalh(1:nrep) = report(1:nrep)% header% sun_zenit
        CASE ('sat_zenit')
          rvalh(1:nrep) = report(1:nrep)% header% sat_zenit
        CASE DEFAULT
          rvalh(1:nrep) = file% nc% vars(in)% rinvalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (rvalh < rmdich)  rvalh = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalh, in, start, count, jerr )
!       =============

      ! integer body elements
      ELSEIF (kcase == 3) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = report(irep)% len
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('varno')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% varno
          CASE ('level_typ')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% level_typ
          CASE ('level_sig')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% level_sig
          CASE ('state')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% state
          CASE ('flags')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% flags
          CASE ('check')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% check
          CASE ('qual')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% qual
          CASE ('spec_index')
            ivalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% spec_index
          CASE DEFAULT
            ivalb (kos+1:kos+nobs) = file% nc% vars(in)% invalid
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start = iboff
        count = nbody
        WHERE (ivalb == imdi)  ivalb = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalb, in, start, count, jerr )
!       =============

      ! real body elements
      ELSEIF (kcase == 4) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = report(irep)% len
!         SELECT CASE (file% nc% vars(in)% name)
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('obs')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% obs
          CASE ('bcor')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% bcor
          CASE ('e_o')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% e_o
          CASE ('level')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% level
          CASE ('plevel')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% plevel
          CASE ('accuracy')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% accuracy
          CASE ('w_qc')
            rvalb (kos+1:kos+nobs) = report(irep)% body(1:nobs)% w_qc
          CASE DEFAULT
            rvalb (kos+1:kos+nobs) = file% nc% vars(in)% rinvalid
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start  = iboff
        count  = nbody
        WHERE (rvalb < rmdich)  rvalb = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalb, in, start, count, jerr )
!       =============

      ! integer 1-d verification meta data elements
      ELSEIF (kcase == 5) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('veri_run_type')
        CASE ('veri_run_class')
        CASE ('veri_forecast_time')
        CASE ('veri_ens_member')
        CASE ('veri_exp_id')
        CASE DEFAULT
          kcase = -kcase
        END SELECT
      ENDIF

! multi-dimensional arrays: 'statid', 'veri_data'
! -----------------------------------------------

    ! 'statid': the only character string (internally: 2-dim array)

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'statid') THEN

      kcase = 11
      DO irep = 1, nrep
        cvala(irep) =               report(irep)% header% statid               &
                       (1:LEN_TRIM( report(irep)% header% statid))
      ENDDO
      start = ihoff
      count = nrep
      CALL add_data ( file, cvala, in, start, count, jerr )
!     =============

    ! 'veri_data': the only truly 2-dimensional array

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'veri_data') THEN

      kcase = 12
!     file% nc% vars(in)% p(2)% dim% len = file% n_veri
      ALLOCATE( rvala2 (nbody, file% n_veri), STAT = jerr )
      IF (jerr /= 0) THEN
        WRITE( yerr,'("ERROR in memory alloc: rvala2 (nbody, file% n_veri) "   &
                    &,2I7)' )  nbody, file% n_veri
                                                                          RETURN
      ENDIF
      rvala2(:,:) = rmdich *1.1
      start = iboff
      count = nbody
      !    it is assumed (in t_acc_body) that (file% n_veri <= 1)
      DO jj = 1, file% n_veri
        kos = 0
        DO irep = 1, nrep
          nobs = report(irep)% len
          DO iob = 1, nobs
            rvala2 (kos+iob,jj) = report(irep)% body(iob)% veri_data(jj)
          ENDDO
          kos = kos + nobs
        ENDDO
        WHERE (rvala2(:,jj)< rmdich) rvala2(:,jj) = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvala2, in, jj, start, count, jerr )
!       =============

      ENDDO
      DEALLOCATE(rvala2)

    ! other multi-dimensional arrays (or scalars)
    ELSE
      kcase = 10
      SELECT CASE (yvar (1:LEN_TRIM( yvar )))
      CASE ('veri_model')
      CASE ('veri_initial_date')
      CASE ('veri_resolution')
      CASE ('veri_domain_size')
      CASE ('veri_description')
      CASE DEFAULT
        kcase = -kcase
      END SELECT
    ENDIF

! caution message for unknown variables, dimensions etc.
! ------------------------------------------------------

    IF (kcase <= 0) THEN
      PRINT '("CAUTION writing feedobs file:",I3," unknown variable: ",A )'    &
             , kcase, yvar(1:LEN_TRIM( yvar ))
    ENDIF

!   IF ((kcase == 0) .OR. (kcase == 9)) THEN
!     jerr = 27
!     WRITE( yerr,'("ERROR in write_report:",I2,": unknown var: ",A)' )        &
!            kcase, file% nc% vars(in)% name
!                                                                         RETURN
!   ENDIF
  ENDDO

  DEALLOCATE(ivalh)
  DEALLOCATE(rvalh)
  DEALLOCATE(cvala)
  DEALLOCATE(ivalb)
  DEALLOCATE(rvalb)

  END SUBROUTINE write_report

!===============================================================================


! Version with rep_header(nrep), rep_body(nrep,:) instead of container report(:)
! otherwise similar to write_report above.
!    Optimized on the NEC for the case that number of reports "nrep" is much 
!    smaller than the typical body length of one report.

SUBROUTINE write_report_radar_1 ( file, rep_header, rep_body, rep_offset       &
                                , rep_len, nrep, dim2_body, nbody, ihoff       &
                                , iboff, imdi, rmdich , jerr, yerr )
 
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------
  TYPE(t_fdbk),       INTENT (inout) :: file       ! metafile
!UB On NEC: serious performance problems with allocatable pointer components!
!UB TYPE(t_account),    INTENT (in)    :: report(:)  ! reports to be stored
  INTEGER,            INTENT (in)    :: dim2_body  ! dim2 of rep_body
  INTEGER,            INTENT (in)    :: nrep       ! number of reports
  TYPE(t_acc_header), INTENT (in)    :: rep_header(nrep) ! report headers
                                                         !   to be stored
  TYPE(t_acc_body),   INTENT (in)    :: rep_body  (nrep,dim2_body) ! rep. bodies
  INTEGER,            INTENT (in)    :: rep_offset(nrep) ! offset of each report
  INTEGER,            INTENT (in)    :: rep_len   (nrep) ! len of each rep. body
  INTEGER,            INTENT (in)    :: nbody      ! number of obs for time box
  INTEGER,            INTENT (in)    :: ihoff      ! header offset for time box
  INTEGER,            INTENT (in)    :: iboff      ! body offset for time box
  INTEGER,            INTENT (in)    :: imdi       ! missing data indicator for
                                                   !   integer (2^31-1)
  REAL   ,            INTENT (in)    :: rmdich     ! check value for missing
                                                   !   real data (-1.E30)
  INTEGER,            INTENT (inout) :: jerr       ! error status variable
  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= *), INTENT (inout) :: yerr       ! error message

!---------------------------------------------------------
  INTEGER                       ::  &
    irep    ,& ! loop indices
    jj, iob ,& ! loop indices
    kcase   ,& ! type of variable (int/real, header/body, etc)
!   kd_hdr  ,& ! length of header dimension
!   kd_body ,& ! length of body   dimension
    varid   ,& ! variable id in netcdf file
    start   ,& ! start position in netcdf file
    count   ,& ! number of elements written
    nvdim   ,& ! nr of dimension veri_data field in meta file
    in      ,& ! variable number to be stored in meta file
    nobs    ,& ! number of observations (body length)
    kos     ,& ! index offset obs from current report
    kerr    ,& ! error status variable
    iu         ! helper variable

  CHARACTER (LEN=nlen) :: yvar    ! NetCDF variable name
  CHARACTER (LEN=nlen) :: ydim    ! NetCDF dimension name

  INTEGER ,           ALLOCATABLE :: ivalh(:) ,& ! buffer for int header elem.
                                     ivalb(:)    ! buffer for int body elements
  REAL    ,           ALLOCATABLE :: rvalh(:) ,& ! buffer for real header elem.
                                     rvalb(:)    ! buffer for real body   elem.
  REAL    ,           ALLOCATABLE :: rvala2(:,:) ! buffer for real body   elem.
  CHARACTER(len=100), ALLOCATABLE :: cvala(:)    ! buffer for char header elem.
!-------------------------------------------------------------------------------
 
  jerr = 0
  kerr = 0
  IF (jerr == 0)  ALLOCATE ( ivalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1
  IF (jerr == 0)  ALLOCATE ( rvalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2
  IF (jerr == 0)  ALLOCATE ( cvala (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3
  IF (jerr == 0)  ALLOCATE ( ivalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 4
  IF (jerr == 0)  ALLOCATE ( rvalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 5
  IF (jerr /= 0) THEN
    IF (kerr == 1)  WRITE( yerr,'("ERROR in alloc: ivalh (nrep ) ",I8)' ) nrep
    IF (kerr == 2)  WRITE( yerr,'("ERROR in alloc: rvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 3)  WRITE( yerr,'("ERROR in alloc: cvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 4)  WRITE( yerr,'("ERROR in alloc: ivalb (nbody) ",I8)' ) nbody
    IF (kerr == 5)  WRITE( yerr,'("ERROR in alloc: rvalb (nbody) ",I8)' ) nbody
                                                                          RETURN
  ENDIF
  ivalh(:) = imdi
  rvalh(:) = rmdich *1.1
  cvala(:) = REPEAT(' ',SIZE(cvala,1))
  ivalb(:) = imdi
  rvalb(:) = rmdich *1.1

  kerr = nf_enddef(file% nc% ncid )
! IF (kerr /= NF_NOERR) THEN
!   jerr = kerr
!   WRITE( yerr,'("ERROR in nf_enddef in write_reports")' )
!                                                                         RETURN
! ENDIF


  write_loop: DO in = 1, file% nc% nvar

    IF (.NOT. file% nc% vars(in)% opt_used)                                CYCLE
    kcase = 0

    yvar(:) = ' '
    yvar  = TRIM(ADJUSTL(file% nc% vars(in)% name))

! 1-dimensional arrays: integer or real, header or body length
! ------------------------------------------------------------

    IF (file% nc% vars(in)% nvdims == 1) THEN

      kcase = 1
      IF (     (file% nc% vars(in)% xtype == NF_FLOAT)                         &
          .OR. (file% nc% vars(in)% xtype == NF_FILL_DOUBLE))  kcase = 2
!     IF (file% nc% vars(in)% xtype == NF_CHAR)  kcase = 5
      ydim  =  file% nc% vars(in)% p(1)% dim% name
      IF (ydim (1:LEN_TRIM( ydim )) == 'd_body') THEN
        kcase = kcase + 2
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_veri') THEN
        kcase = 5
      ELSEIF (ydim (1:LEN_TRIM( ydim )) /= 'd_hdr') THEN
        kcase = -6
      ENDIF

      ! integer header elements
      IF (kcase == 1) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('i_body')
          ivalh (1:nrep) = rep_offset(1:nrep)
        CASE ('l_body')
          ivalh (1:nrep) = rep_len(1:nrep)
        CASE ('n_level')
          ivalh (1:nrep) = rep_header(1:nrep)% n_level
        CASE ('data_category')
          ivalh (1:nrep) = rep_header(1:nrep)% data_category
        CASE ('sub_category')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_category
        CASE ('center')
          ivalh (1:nrep) = rep_header(1:nrep)% center
        CASE ('sub_center')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_center
        CASE ('obstype')
          ivalh (1:nrep) = rep_header(1:nrep)% obstype
        CASE ('codetype')
          ivalh (1:nrep) = rep_header(1:nrep)% codetype
        CASE ('ident')
          ivalh (1:nrep) = rep_header(1:nrep)% ident
        CASE ('time')
          ivalh (1:nrep) = rep_header(1:nrep)% time
        CASE ('time_nomi')
          ivalh (1:nrep) = rep_header(1:nrep)% time_nomi
        CASE ('time_dbase')
          ivalh (1:nrep) = rep_header(1:nrep)% time_dbase
        CASE ('z_station')
          ivalh (1:nrep) = rep_header(1:nrep)% z_station
        CASE ('z_modsurf')
          ivalh (1:nrep) = rep_header(1:nrep)% z_modsurf
        CASE ('r_state')
          ivalh (1:nrep) = rep_header(1:nrep)% r_state
        CASE ('r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% r_flags
        CASE ('r_check')
          ivalh (1:nrep) = rep_header(1:nrep)% r_check
        CASE ('sta_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% sta_corr
        CASE ('mdlsfc')
          ivalh (1:nrep) = rep_header(1:nrep)% mdlsfc
        CASE ('instype')
          ivalh (1:nrep) = rep_header(1:nrep)% instype
        CASE ('retrtype')
          ivalh (1:nrep) = rep_header(1:nrep)% retrtype
        CASE ('phase')
          ivalh (1:nrep) = rep_header(1:nrep)% phase
        CASE ('tracking')
          ivalh (1:nrep) = rep_header(1:nrep)% tracking
        CASE ('meas_type')
          ivalh (1:nrep) = rep_header(1:nrep)% meas_type
        CASE ('rad_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% rad_corr
        CASE ('surftype')
          ivalh (1:nrep) = rep_header(1:nrep)% surftype
        CASE ('flg_1dvar')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_1dvar
        CASE ('flg_cld')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_cld
        CASE ('index_x')
          ivalh (1:nrep) = rep_header(1:nrep)% index_x
        CASE ('index_y')
          ivalh (1:nrep) = rep_header(1:nrep)% index_y
        CASE ('obs_id')
          ivalh (1:nrep) = rep_header(1:nrep)% obs_id
        CASE ('source')
          ivalh (1:nrep) = rep_header(1:nrep)% source
        CASE ('record')
          ivalh (1:nrep) = rep_header(1:nrep)% record
        CASE ('subset')
          ivalh (1:nrep) = rep_header(1:nrep)% subset
        CASE ('dbkz')
          ivalh (1:nrep) = rep_header(1:nrep)% dbkz
        CASE ('index_d')
          ivalh (1:nrep) = rep_header(1:nrep)% index_d
        CASE DEFAULT
          ivalh (1:nrep) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (ivalh == imdi)  ivalh = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalh, in, start, count, jerr )
!       =============

      ! real header elements
      ELSEIF (kcase == 2) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('lat')
          rvalh(1:nrep) = rep_header(1:nrep)% lat
        CASE ('lon')
          rvalh(1:nrep) = rep_header(1:nrep)% lon
        CASE ('sun_zenit')
          rvalh(1:nrep) = rep_header(1:nrep)% sun_zenit
        CASE DEFAULT
          rvalh(1:nrep) = file% nc% vars(in)% rinvalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (rvalh < rmdich)  rvalh = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalh, in, start, count, jerr )
!       =============

      ! integer body elements
      ELSEIF (kcase == 3) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = rep_len(irep)
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('varno')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% varno
          CASE ('level_typ')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% level_typ
          CASE ('level_sig')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% level_sig
          CASE ('state')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% state
          CASE ('flags')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% flags
          CASE ('check')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% check
          CASE ('qual')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% qual
          CASE ('spec_index')
            ivalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% spec_index
          CASE DEFAULT
            ivalb (kos+1:kos+nobs) = file% nc% vars(in)% invalid
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start = iboff
        count = nbody
        WHERE (ivalb == imdi)  ivalb = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalb, in, start, count, jerr )
!       =============

      ! real body elements
      ELSEIF (kcase == 4) THEN
        kos = 0
        DO irep = 1, nrep
          nobs = rep_len(irep)
!         SELECT CASE (file% nc% vars(in)% name)
          SELECT CASE (yvar (1:LEN_TRIM( yvar )))
          CASE ('obs')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% obs
          CASE ('bcor')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% bcor
          CASE ('e_o')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% e_o
          CASE ('level')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% level
          CASE ('plevel')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% plevel
          CASE ('accuracy')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% accuracy
          CASE ('w_qc')
            rvalb (kos+1:kos+nobs) = rep_body(irep,1:nobs)% w_qc
          CASE DEFAULT
            rvalb (kos+1:kos+nobs) = file% nc% vars(in)% rinvalid
            kcase = -kcase
          END SELECT
          kos = kos + nobs
        ENDDO
        start  = iboff
        count  = nbody
        WHERE (rvalb < rmdich)  rvalb = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalb, in, start, count, jerr )
!       =============

      ! integer 1-d verification meta data elements
      ELSEIF (kcase == 5) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('veri_run_type')
        CASE ('veri_run_class')
        CASE ('veri_forecast_time')
        CASE ('veri_ens_member')
        CASE ('veri_exp_id')
        CASE DEFAULT
          kcase = -kcase
        END SELECT
      ENDIF

! multi-dimensional arrays: 'statid', 'veri_data'
! -----------------------------------------------

    ! 'statid': the only character string (internally: 2-dim array)

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'statid') THEN

      kcase = 11
      DO irep = 1, nrep
        cvala(irep) =               rep_header(irep)% statid               &
                       (1:LEN_TRIM( rep_header(irep)% statid))
      ENDDO
      start = ihoff
      count = nrep
      CALL add_data ( file, cvala, in, start, count, jerr )
!     =============

    ! 'veri_data': the only truly 2-dimensional array

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'veri_data') THEN

      kcase = 12
!     file% nc% vars(in)% p(2)% dim% len = file% n_veri
      ALLOCATE( rvala2 (nbody, file% n_veri), STAT = jerr )
      IF (jerr /= 0) THEN
        WRITE( yerr,'("ERROR in memory alloc: rvala2 (nbody, file% n_veri) "   &
                    &,2I7)' )  nbody, file% n_veri
                                                                          RETURN
      ENDIF
      rvala2(:,:) = rmdich *1.1
      start = iboff
      count = nbody
      DO jj = 1, file% n_veri
        kos = 0
        DO irep = 1, nrep
          nobs = rep_len(irep)
          DO iob = 1, nobs
            rvala2 (kos+iob,jj) = rep_body(irep,iob)% veri_data(jj)
          ENDDO
          kos = kos + nobs
        ENDDO
        WHERE (rvala2(:,jj)< rmdich) rvala2(:,jj) = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvala2, in, jj, start, count, jerr )
!       =============

      ENDDO
      DEALLOCATE(rvala2)

    ! other multi-dimensional arrays (or scalars)
    ELSE
      kcase = 10
      SELECT CASE (yvar (1:LEN_TRIM( yvar )))
      CASE ('veri_model')
      CASE ('veri_initial_date')
      CASE ('veri_resolution')
      CASE ('veri_domain_size')
      CASE ('veri_description')
      CASE DEFAULT
        kcase = -kcase
      END SELECT
    ENDIF

! caution message for unknown variables, dimensions etc.
! ------------------------------------------------------

    IF (kcase <= 0) THEN
      PRINT '("CAUTION writing feedobs file:",I3," unknown variable: ",A )'    &
             , kcase, yvar(1:LEN_TRIM( yvar ))
    ENDIF

!   IF ((kcase == 0) .OR. (kcase == 9)) THEN
!     jerr = 27
!     WRITE( yerr,'("ERROR in write_report:",I2,": unknown var: ",A)' )        &
!            kcase, file% nc% vars(in)% name
!                                                                         RETURN
!   ENDIF
  ENDDO write_loop

  DEALLOCATE(ivalh)
  DEALLOCATE(rvalh)
  DEALLOCATE(cvala)
  DEALLOCATE(ivalb)
  DEALLOCATE(rvalb)

END SUBROUTINE write_report_radar_1

!===============================================================================

! .. Version similar to write_report_radar_1, but optimized on the NEC for
!    the case that number of reports "nrep" is much larger than the 
!    typical body length of one report.

SUBROUTINE write_report_radar_2 ( file, rep_header, rep_body, rep_offset       &
                                , rep_len, nrep, dim2_body, nbody, ihoff       &
                                , iboff, imdi, rmdich , jerr, yerr )
 
!-------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------
  TYPE(t_fdbk),       INTENT (inout) :: file       ! metafile
!UB On NEC: serious performance problems with allocatable pointer components!
!UB TYPE(t_account),    INTENT (in)    :: report(:)  ! reports to be stored
  INTEGER,            INTENT (in)    :: dim2_body  ! dim2 of rep_body
  INTEGER,            INTENT (in)    :: nrep       ! number of reports
  TYPE(t_acc_header), INTENT (in)    :: rep_header(nrep) ! report headers
                                                         !   to be stored
  TYPE(t_acc_body),   INTENT (in)    :: rep_body  (nrep,dim2_body) ! rep. bodies
  INTEGER,            INTENT (in)    :: rep_offset(nrep) ! offset of each report
  INTEGER,            INTENT (in)    :: rep_len   (nrep) ! len of each rep. body
  INTEGER,            INTENT (in)    :: nbody      ! number of obs for time box
  INTEGER,            INTENT (in)    :: ihoff      ! header offset for time box
  INTEGER,            INTENT (in)    :: iboff      ! body offset for time box
  INTEGER,            INTENT (in)    :: imdi       ! missing data indicator for
                                                   !   integer (2^31-1)
  REAL   ,            INTENT (in)    :: rmdich     ! check value for missing
                                                   !   real data (-1.E30)
  INTEGER,            INTENT (inout) :: jerr       ! error status variable
  ! it is assumed here that LEN >= 72 !
  CHARACTER (LEN= *), INTENT (inout) :: yerr       ! error message

!---------------------------------------------------------
  INTEGER                       ::  &
    irep    ,& ! loop indices
    jj, iob, iobu ,& ! loop indices
    kcase   ,& ! type of variable (int/real, header/body, etc)
!   kd_hdr  ,& ! length of header dimension
!   kd_body ,& ! length of body   dimension
    varid   ,& ! variable id in netcdf file
    start   ,& ! start position in netcdf file
    count   ,& ! number of elements written
    nvdim   ,& ! nr of dimension veri_data field in meta file
    in      ,& ! variable number to be stored in meta file
    nobs    ,& ! number of observations (body length)
    kos     ,& ! index offset obs from current report
    kerr    ,& ! error status variable
    iu         ! helper variable

  CHARACTER (LEN=nlen) :: yvar    ! NetCDF variable name
  CHARACTER (LEN=nlen) :: ydim    ! NetCDF dimension name

  INTEGER ,           ALLOCATABLE :: ivalh(:) ,& ! buffer for int header elem.
                                     ivalb(:)    ! buffer for int body elements
  REAL    ,           ALLOCATABLE :: rvalh(:) ,& ! buffer for real header elem.
                                     rvalb(:)    ! buffer for real body   elem.
  REAL    ,           ALLOCATABLE :: rvala2(:,:) ! buffer for real body   elem.
  CHARACTER(len=100), ALLOCATABLE :: cvala(:)    ! buffer for char header elem.
  INTEGER :: ifillbuf(nrep,dim2_body)
  REAL    :: rfillbuf(nrep,dim2_body)
!-------------------------------------------------------------------------------
 
  jerr = 0
  kerr = 0
  IF (jerr == 0)  ALLOCATE ( ivalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 1
  IF (jerr == 0)  ALLOCATE ( rvalh (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 2
  IF (jerr == 0)  ALLOCATE ( cvala (nrep ) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 3
  IF (jerr == 0)  ALLOCATE ( ivalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 4
  IF (jerr == 0)  ALLOCATE ( rvalb (nbody) , STAT = jerr )
  IF ((jerr /= 0) .AND. (kerr == 0))  kerr = 5
  IF (jerr /= 0) THEN
    IF (kerr == 1)  WRITE( yerr,'("ERROR in alloc: ivalh (nrep ) ",I8)' ) nrep
    IF (kerr == 2)  WRITE( yerr,'("ERROR in alloc: rvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 3)  WRITE( yerr,'("ERROR in alloc: cvalh (nrep ) ",I8)' ) nrep
    IF (kerr == 4)  WRITE( yerr,'("ERROR in alloc: ivalb (nbody) ",I8)' ) nbody
    IF (kerr == 5)  WRITE( yerr,'("ERROR in alloc: rvalb (nbody) ",I8)' ) nbody
                                                                          RETURN
  ENDIF
  ivalh(:) = imdi
  rvalh(:) = rmdich *1.1
  cvala(:) = REPEAT(' ',SIZE(cvala,1))
  ivalb(:) = imdi
  rvalb(:) = rmdich *1.1

  kerr = nf_enddef(file% nc% ncid )
! IF (kerr /= NF_NOERR) THEN
!   jerr = kerr
!   WRITE( yerr,'("ERROR in nf_enddef in write_reports")' )
!                                                                         RETURN
! ENDIF

  
  write_loop: DO in = 1, file% nc% nvar

    IF (.NOT. file% nc% vars(in)% opt_used)                                CYCLE
    kcase = 0

    yvar(:) = ' '
    yvar  = TRIM(ADJUSTL(file% nc% vars(in)% name))

! 1-dimensional arrays: integer or real, header or body length
! ------------------------------------------------------------

    IF (file% nc% vars(in)% nvdims == 1) THEN

      kcase = 1
      IF (     (file% nc% vars(in)% xtype == NF_FLOAT)                         &
          .OR. (file% nc% vars(in)% xtype == NF_FILL_DOUBLE))  kcase = 2
!     IF (file% nc% vars(in)% xtype == NF_CHAR)  kcase = 5
      ydim  =  file% nc% vars(in)% p(1)% dim% name
      IF (ydim (1:LEN_TRIM( ydim )) == 'd_body') THEN
        kcase = kcase + 2
      ELSEIF (ydim (1:LEN_TRIM( ydim )) == 'd_veri') THEN
        kcase = 5
      ELSEIF (ydim (1:LEN_TRIM( ydim )) /= 'd_hdr') THEN
        kcase = -6
      ENDIF

      ! integer header elements
      IF (kcase == 1) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('i_body')
          ivalh (1:nrep) = rep_offset(1:nrep)
        CASE ('l_body')
          ivalh (1:nrep) = rep_len(1:nrep)
        CASE ('n_level')
          ivalh (1:nrep) = rep_header(1:nrep)% n_level
        CASE ('data_category')
          ivalh (1:nrep) = rep_header(1:nrep)% data_category
        CASE ('sub_category')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_category
        CASE ('center')
          ivalh (1:nrep) = rep_header(1:nrep)% center
        CASE ('sub_center')
          ivalh (1:nrep) = rep_header(1:nrep)% sub_center
        CASE ('obstype')
          ivalh (1:nrep) = rep_header(1:nrep)% obstype
        CASE ('codetype')
          ivalh (1:nrep) = rep_header(1:nrep)% codetype
        CASE ('ident')
          ivalh (1:nrep) = rep_header(1:nrep)% ident
        CASE ('time')
          ivalh (1:nrep) = rep_header(1:nrep)% time
        CASE ('time_nomi')
          ivalh (1:nrep) = rep_header(1:nrep)% time_nomi
        CASE ('time_dbase')
          ivalh (1:nrep) = rep_header(1:nrep)% time_dbase
        CASE ('z_station')
          ivalh (1:nrep) = rep_header(1:nrep)% z_station
        CASE ('z_modsurf')
          ivalh (1:nrep) = rep_header(1:nrep)% z_modsurf
        CASE ('r_state')
          ivalh (1:nrep) = rep_header(1:nrep)% r_state
        CASE ('r_flags')
          ivalh (1:nrep) = rep_header(1:nrep)% r_flags
        CASE ('r_check')
          ivalh (1:nrep) = rep_header(1:nrep)% r_check
        CASE ('sta_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% sta_corr
        CASE ('mdlsfc')
          ivalh (1:nrep) = rep_header(1:nrep)% mdlsfc
        CASE ('instype')
          ivalh (1:nrep) = rep_header(1:nrep)% instype
        CASE ('retrtype')
          ivalh (1:nrep) = rep_header(1:nrep)% retrtype
        CASE ('phase')
          ivalh (1:nrep) = rep_header(1:nrep)% phase
        CASE ('tracking')
          ivalh (1:nrep) = rep_header(1:nrep)% tracking
        CASE ('meas_type')
          ivalh (1:nrep) = rep_header(1:nrep)% meas_type
        CASE ('rad_corr')
          ivalh (1:nrep) = rep_header(1:nrep)% rad_corr
        CASE ('surftype')
          ivalh (1:nrep) = rep_header(1:nrep)% surftype
        CASE ('flg_1dvar')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_1dvar
        CASE ('flg_cld')
          ivalh (1:nrep) = rep_header(1:nrep)% flg_cld
        CASE ('index_x')
          ivalh (1:nrep) = rep_header(1:nrep)% index_x
        CASE ('index_y')
          ivalh (1:nrep) = rep_header(1:nrep)% index_y
        CASE ('obs_id')
          ivalh (1:nrep) = rep_header(1:nrep)% obs_id
        CASE ('source')
          ivalh (1:nrep) = rep_header(1:nrep)% source
        CASE ('record')
          ivalh (1:nrep) = rep_header(1:nrep)% record
        CASE ('subset')
          ivalh (1:nrep) = rep_header(1:nrep)% subset
        CASE ('dbkz')
          ivalh (1:nrep) = rep_header(1:nrep)% dbkz
        CASE ('index_d')
          ivalh (1:nrep) = rep_header(1:nrep)% index_d
        CASE DEFAULT
          ivalh (1:nrep) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (ivalh == imdi)  ivalh = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalh, in, start, count, jerr )
!       =============

      ! real header elements
      ELSEIF (kcase == 2) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('lat')
          rvalh(1:nrep) = rep_header(1:nrep)% lat
        CASE ('lon')
          rvalh(1:nrep) = rep_header(1:nrep)% lon
        CASE ('sun_zenit')
          rvalh(1:nrep) = rep_header(1:nrep)% sun_zenit
        CASE DEFAULT
          rvalh(1:nrep) = file% nc% vars(in)% rinvalid
          kcase = -kcase
        END SELECT
        start = ihoff
        count = nrep
        WHERE (rvalh < rmdich)  rvalh = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalh, in, start, count, jerr )
!       =============

      ! integer body elements
      ELSEIF (kcase == 3) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('varno')
          ifillbuf = rep_body(:,:)% varno
          CALL fill_bodybuf_int (ifillbuf, ivalb)
        CASE ('level_typ')
          ifillbuf = rep_body(:,:)% level_typ
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('level_sig')
          ifillbuf = rep_body(:,:)% level_sig
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('state')
          ifillbuf = rep_body(:,:)% state
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('flags')
          ifillbuf = rep_body(:,:)% flags
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('check')
          ifillbuf = rep_body(:,:)% check
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('qual')
          ifillbuf = rep_body(:,:)% qual
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE ('spec_index')
          ifillbuf = rep_body(:,:)% spec_index
          CALL fill_bodybuf_int (ifillbuf , ivalb)
        CASE DEFAULT
          ivalb (:) = file% nc% vars(in)% invalid
          kcase = -kcase
        END SELECT
        start = iboff
        count = nbody
        WHERE (ivalb == imdi)  ivalb = file% nc% vars(in)% invalid
        CALL add_data ( file, ivalb, in, start, count, jerr )
!       =============

      ! real body elements
      ELSEIF (kcase == 4) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('obs')
          rfillbuf = rep_body(:,:)%obs
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('bcor')
          rfillbuf = rep_body(:,:)%bcor
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('e_o')
          rfillbuf = rep_body(:,:)%e_o
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('level')
          rfillbuf = rep_body(:,:)%level
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('plevel')
          rfillbuf = rep_body(:,:)%plevel
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('accuracy')
          rfillbuf = rep_body(:,:)%accuracy
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE ('w_qc')
          rfillbuf = rep_body(:,:)%w_qc
          CALL fill_bodybuf_real (rfillbuf , rvalb)
        CASE DEFAULT
          rvalb (:) = file% nc% vars(in)% rinvalid
          kcase = -kcase
        END SELECT
        start  = iboff
        count  = nbody
        WHERE (rvalb < rmdich)  rvalb = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvalb, in, start, count, jerr )
!       =============

      ! integer 1-d verification meta data elements
      ELSEIF (kcase == 5) THEN
        SELECT CASE (yvar (1:LEN_TRIM( yvar )))
        CASE ('veri_run_type')
        CASE ('veri_run_class')
        CASE ('veri_forecast_time')
        CASE ('veri_ens_member')
        CASE ('veri_exp_id')
        CASE DEFAULT
          kcase = -kcase
        END SELECT
      ENDIF

! multi-dimensional arrays: 'statid', 'veri_data'
! -----------------------------------------------

    ! 'statid': the only character string (internally: 2-dim array)

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'statid') THEN

      kcase = 11
      DO irep = 1, nrep
        cvala(irep) =               rep_header(irep)% statid               &
                       (1:LEN_TRIM( rep_header(irep)% statid))
      ENDDO
      start = ihoff
      count = nrep
      CALL add_data ( file, cvala, in, start, count, jerr )
!     =============

    ! 'veri_data': the only truly 2-dimensional array

    ELSEIF (yvar (1:LEN_TRIM( yvar )) == 'veri_data') THEN

      kcase = 12
!     file% nc% vars(in)% p(2)% dim% len = file% n_veri
      ALLOCATE( rvala2 (nbody, file% n_veri), STAT = jerr )
      IF (jerr /= 0) THEN
        WRITE( yerr,'("ERROR in memory alloc: rvala2 (nbody, file% n_veri) "   &
                    &,2I7)' )  nbody, file% n_veri
                                                                          RETURN
      ENDIF
      rvala2(:,:) = rmdich *1.1
      start = iboff
      count = nbody
      iobu  = rep_offset(1)
      DO jj = 1, file% n_veri
        DO iob = 1, dim2_body
!CDIR NODEP
          DO irep = 1, nrep
            IF (iob <= rep_len(irep)) THEN
              iu = rep_offset(irep) - iobu
              rvala2 (iu+iob,jj) = rep_body(irep,iob)% veri_data(jj)
            END IF
          END DO
        END DO
        WHERE (rvala2(:,jj)< rmdich) rvala2(:,jj) = file% nc% vars(in)% rinvalid
        CALL add_data ( file, rvala2, in, jj, start, count, jerr )
!       =============

      ENDDO
      DEALLOCATE(rvala2)

    ! other multi-dimensional arrays (or scalars)
    ELSE
      kcase = 10
      SELECT CASE (yvar (1:LEN_TRIM( yvar )))
      CASE ('veri_model')
      CASE ('veri_initial_date')
      CASE ('veri_resolution')
      CASE ('veri_domain_size')
      CASE ('veri_description')
      CASE DEFAULT
        kcase = -kcase
      END SELECT
    ENDIF

! caution message for unknown variables, dimensions etc.
! ------------------------------------------------------

    IF (kcase <= 0) THEN
      PRINT '("CAUTION writing feedobs file:",I3," unknown variable: ",A )'    &
             , kcase, yvar(1:LEN_TRIM( yvar ))
    ENDIF

!   IF ((kcase == 0) .OR. (kcase == 9)) THEN
!     jerr = 27
!     WRITE( yerr,'("ERROR in write_report:",I2,": unknown var: ",A)' )        &
!            kcase, file% nc% vars(in)% name
!                                                                         RETURN
!   ENDIF
  ENDDO write_loop

  DEALLOCATE(ivalh)
  DEALLOCATE(rvalh)
  DEALLOCATE(cvala)
  DEALLOCATE(ivalb)
  DEALLOCATE(rvalb)

CONTAINS

  SUBROUTINE fill_bodybuf_int (ibodydata, ibuf)

!!$    INTEGER, INTENT(in)  :: ibodydata (nrep,dim2_body)
!!$    INTEGER, INTENT(inout) :: ibuf (nbody)
    INTEGER, INTENT(in)  :: ibodydata (:,:)
    INTEGER, INTENT(inout) :: ibuf (:)

    INTEGER :: iu, iob, iobu, irep

    ibuf = -9999
    iobu = rep_offset(1)
    DO iob = 1, dim2_body
!CDIR NODEP
      DO irep = 1, nrep
        IF (iob <= rep_len(irep)) THEN
          iu = rep_offset(irep) - iobu
          ibuf (iu+iob) = ibodydata(irep,iob)
        END IF
      END DO
    END DO
  
  END SUBROUTINE fill_bodybuf_int

  SUBROUTINE fill_bodybuf_real (rbodydata, rbuf)

!!$    REAL, INTENT(in)  :: rbodydata (nrep,dim2_body)
!!$    REAL, INTENT(inout) :: rbuf (nbody)
    REAL, INTENT(in)  :: rbodydata (:,:)
    REAL, INTENT(inout) :: rbuf (:)

    INTEGER :: iu, iob, iobu, irep

    rbuf = -9999.99
    iobu = rep_offset(1)
    DO iob = 1, dim2_body
!CDIR NODEP
      DO irep = 1, nrep
        IF (iob <= rep_len(irep)) THEN
          iu = rep_offset(irep) - iobu
          rbuf (iu+iob) = rbodydata(irep,iob)
        END IF
      END DO
    END DO
  
  END SUBROUTINE fill_bodybuf_real

END SUBROUTINE write_report_radar_2

!-------------------------------------------------------------------------------

END MODULE mo_fdbk_cosmo
