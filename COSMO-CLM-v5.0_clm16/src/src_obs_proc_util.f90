!+ Source module for utilities for the observation processing, data assimilation
!-------------------------------------------------------------------------------

MODULE src_obs_proc_util

!-------------------------------------------------------------------------------
! Description:
!   This module provides utilities used in the (pre-)processing of observation
!   reports read from an AOF file for data assimilation purposes.
!   This includes all the tasks which require (possibly time-dependent) model
!   fields within the observation processing based on AOF reading (except
!   possibly for fields in 'src_obs_1dvar_org').
!   Specific tasks are:
!    - assigning horizontally an observational report to a model grid point
!    - gathering surface fields (orography, land/sea mask) from global domain
!    - assigning position in diagnostic array
!    - assigning observation errors
!    - deleting reports too old to be used any more from internal arrays (ODR)
!    - determining orographic height at a specific grid point
!    - determining pressure from height using model fields
!    - determining temperature from virtual temperature using model fields
!    - determining model layer thickness from pressure using model fields
!   This module contains the following module procedures:
!      name                  : used by (described only if model fields are used)
!    - obs_aof_assign_gridpt
!    - get_global_surf_aof   : obs_read_distribute: to assign obs to grid pt.
!                              obs_read_gps       : to assign obs to grid pt.
!    - obs_pointrs
!    - obs_fix_error
!    - obs_error
!    - obs_del_old_report
!    - zsurf             : obs_del_store_retriev: orographic height at retrieval
!    - psurf             : obs_single_level: to reject AIREP /SATOB obs if close
!                                            to orography + height not reported
!    - tupair            : obs_del_store_retriev: model temperature, to add o.i.
!    - rhupair           : obs_del_store_retriev: model rel humidity to add o.i.
!    - z2p               : obs_multi_level  \  to assign pressure to PILOT /
!                          obs_single_level /  SATOB obs when it is not reported
!                          obs_sort_levels : to sort levels (if p not reported)
!    - tv2t              : obs_RASS_tv2t   : to get T from Tv observed by RASS
!    - p2dp              : obs_air_org_mult: for adapting vertical correlation
!                                            scales for aircraft reports
!  Note: This module is not used any more when reading obs from NetCDF files !
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 3.12       2004/09/15 Christoph Schraff
!  Initial release
! 3.13       2004/12/03 Ulrich Schaettler
!  Adapted interface of routine gather_field
! 3.18       2006/03/03 Christoph Schraff
!  Introduction of routine 'zsurf', cancellation of 'obs_del_old_retriv'.
!  Preparation for use of real-data 1DVar satellite retrievals (MSG, NOAA15-18).
! V4_5         2008/09/10 Christoph Schraff
!  obs_assign_gridpt: Use of radiosonde ship and 'lseaobs' assignment modified.
! V4_8         2009/02/16 Ulrich Schaettler
!  Included routine obs_cdf_distrib_reports from src_obs_proc_cdf
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Routine 'obs_cdf_distrib_reports' moved to new module 'src_obs_cdfin_util'.
!  - Variables moved from 'data_nudge_all' to 'data_obs_lib_cosmo'.
!  - Routine names modified (obs_assign_gridpt  --> obs_aof_assign_gridpt,
!                            get_global_surface --> get_global_surf_aof  ).
!  - 'obs_pointrs' adjusted to cope with obs type = 12 (GPS).
!  Note: This module is not used any more when reading obs from NetCDF files !
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  Introduced MESSy interface
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
! Modules used:
!
!-------------------------------------------------------------------------------

USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers, & ! KIND-type parameter for standard integer variables
    irealgrib    ! KIND-type parameter for real variables in the grib library

!-------------------------------------------------------------------------------

USE data_obs_lib_cosmo, ONLY :  &

! 1. General parameters
! ---------------------

    c0         ,& ! standard real constant 0.0
    c1         ,& ! standard real constant 1.0
    c2         ,& ! standard real constant 2.0
    c05        ,& ! standard real constant 0.5
    c3600      ,& ! standard real constant 3600.0
    rmdi       ,& ! =-1.E31_ireals : commonly used missing data indicator
    rmdich     ,& ! =-1.E30_ireals : commonly used check value for miss data
    epsy       ,& ! = 1.E-8_ireals : commonly used very small value > 0
    i0         ,& ! standard integer constant 0
    i1         ,& ! standard integer constant 1

! 2. Variables and parameters obtained by calling 'obs_cdf_interface'
! -------------------------------------------------------------------

    ! report dimensions of ODR (Observation Data Record, for observation storage
    ! on local sub-domains), and other variables related to namelist parameters
    maxmll     ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl     ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl     ,& ! size (report dimension) of the (ground-based) GPS  ODR
    maxtvl     ,& ! size (report dimension) of the satellite retrieval ODR

! 4. I/O device numbers for obs processing / nudging and file names
! -----------------------------------------------------------------

    ! file names
    yucautn    ,& ! caution messages if too many obs for ODR size
    yuquctl    ,& ! data rejected by threshold QC at obs time
    yustats    ,& ! statistics of processed reports
    yurejct    ,& ! direct reporting of rejected obs. reports
    yuobsdr    ,& ! observations stored in the observation data record
    yuverif    ,& ! VOF (output): verification observation file (obs.
                  !      incl. quality control flag for verification)
    yuprint    ,& ! all the remaining information
    ! device numbers
    nucautn    ,& ! caution messages if too many obs for current ODR size
    nuqc       ,& ! data rejected by threshold quality control at obs time
    nustat     ,& ! statistics of processed reports
    nurej      ,& ! direct reporting of rejected obs. reports
    nuodr      ,& ! observations stored in the observation data record
    nuverif    ,& ! VOF (output): verification observation file (observations
                  !   incl. quality control flag for verification)
    nupr          ! all the remaining information

USE data_obs_lib_cosmo, ONLY :  &

! 5. CMA observation type and code type numbers
! ---------------------------------------------

    nsynop     ,& ! SYNOP reports
    nairep     ,& ! AIREP reports (all aircraft reports) 
    nsatob     ,& ! SATOB reports
    ndribu     ,& ! DRIBU reports
    ntemp      ,& ! TEMP  reports
    npilot     ,& ! PILOT reports
    nsatem     ,& ! SATEM reports 
    nsattv     ,& ! SATEM reports 
    ngps       ,& ! GPS   reports
    nsrscd     ,& !   synop surface report 
    natscd     ,& !   automatic synop surface report
    nshscd     ,& !   ship synop report
    nabscd     ,& !   ship synop abbreviated report
    nshred     ,& !   shred report
    natshs     ,& !   automatic ship synop report
    nmetar     ,& !   Metar
    naircd     ,& !   aircraft report
    ncodar     ,& !   codar report
    ncolba     ,& !   colba report
    namdar     ,& !   amdar report
    nacar      ,& !   acar  report
    nstbcd     ,& !   satob report
    nhrvis     ,& !   high-res. VIS wind report
    namv       ,& !   AMV   report
    nsst       ,& !   sst report
    ndrbcd     ,& !   dribu report 
    nbathy     ,& !   bathy report
    ntesac        !   tesac report

USE data_obs_lib_cosmo, ONLY :  &
    nldtcd     ,& !   temp land   report
    nshtcd     ,& !   temp ship   report
    nmotcd     ,& !   temp mobile report
    ntdrop     ,& !   temp drop   report
    nrocob     ,& !   rocob      report
    nrocsh     ,& !   rocob ship report
    nldpcd     ,& !   pilot land   report
    nshpcd     ,& !   pilot ship   report
    nmopcd     ,& !   pilot mobile report
    nwp_eu     ,& !   European wind profiler report
    nra_eu     ,& !   European SODAR/RASS report
    nravad     ,& !   Radar VAD wind report
    npr_us     ,& !   US Wind Profiler/RASS report
    nstmcd     ,& !   satem report
    nstovs     ,& !   high resolution ATOVS satellite data
    nsmsg1     ,& !   MSG_1  satellite retrieval
    nsmsg2     ,& !   MSG_2  satellite retrieval
    nnoa15     ,& !   NOAA15 satellite retrieval
    nnoa16     ,& !   NOAA16 satellite retrieval
    nnoa17     ,& !   NOAA17 satellite retrieval
    nnoa18     ,& !   NOAA18 satellite retrieval
    ngpgfz        !   GPS report processed by GFZ

! end of data_obs_lib_cosmo 

!-------------------------------------------------------------------------------

USE data_obs_process, ONLY :  &

!         1.5   AOF extracted parameters and variables
!               --------------------------------------
    llship      ,& ! ship switch
    nobtyp      ,& ! observation type
    ncdtyp      ,& ! code type
    rppp        ,& ! reported pressure level

!         1.6   AOF time boxes
!               --------------
    naofbox     ,& ! length of AOF time periods to be read
    naoflbx     ,& ! initial time of last  AOF time box to be read currently

!         3.4   Variables' numbering for obs. error extraction
!               ----------------------------------------------
    nvrpoe      ,& !   variables numbering for obs. error calculation:

!         4.3   Diagnostic arrays and pointers
!               ------------------------------
    nobtpp      ,& ! diagnostic arrays pointer
    ncdtpp      ,& !                 "
    nmxcdt      ,& ! max. no. of code type for any obs. type
    ntotob      ,& ! total number of obs. types

!         6     Report counters
!               ---------------
    nmlob       ,& ! current number of  multi-level report
    nsgob       ,& ! current number of single-level report
    ngpob          ! current number of GPS report

! end of data_obs_process

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!         4.1  Observation error levels
!              ------------------------
    nerlev      ,& ! number of standard error levels
    rlevel      ,& ! error levels
    rolnlv      ,& ! ln(rlevel(15))

!         4.2  Observation error constants
!              ---------------------------
    oevsond     ,& ! radiosonde wind errors
    oezsond     ,& ! radiosonde height errors
    oetsond     ,& ! radiosonde temperature errors
    oeairep     ,& ! airep wind errors
    oetairp     ,& ! airep temperature errors
    oevsynp     ,& ! synop wind errors
    oezsynp     ,& ! synop height errors
    oevdrib     ,& ! dribu wind errors
    oezdrib     ,& ! dribu height errors
    oezship     ,& ! ship height error
    oesatob     ,& ! satob wind errors
!   rherr1      ,& ! (root of) fixed    / normal conditions
!   rherr2      ,& ! relative humidity <  if rel. humidity below 20%
!   rherr3      ,& ! error variances    \ if temperature below 233K

!         8    Temporary global model fields
!              -----------------------------
    hsurf_tot   ,& ! total array of model surface height
    fland_tot      ! total array of fraction of land in each grid element

! end of data_obs_process

!-------------------------------------------------------------------------------

USE data_nudge_all, ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl       ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2      ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot      ,& ! .TRUE if (my_cart_id  == 0) (for print on std. output)
    onl_cart_id ,& ! 'my_cart_id'  of node with area containing (ionl,jonl)
    lfirst      ,& ! .TRUE if 'organize_nudging' is called the first time
    acthr       ,& ! actual forecast hour

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

! 2.1  temporal variables
!      ------------------
    lverif      ,& ! on - off switch for verification
!   lverpas     ,& ! .t. : on - off switch for verif. also of passive reports
    hversta     ,& ! start of verification period in 'model integration hours'
    hverend     ,& ! end of verification period in 'model integration hours'
    tconbox     ,& ! 6*dt: timestep [s] for computing analysis increments
                   !       (i.e. time box of constant analysis increments)
    tipolmx     ,& ! max. time span  for linear interpolat. for upper-air data
    tipmxsu     ,& ! max. time span  for linear interpolat.of surf.-level data
    wtukrsa     ,& ! temporal radius of infl. towards the past for TEMP/PILOT
    wtukrse     ,& ! temporal radius of infl. towards the future for TEMP/PILOT
    wtukara     ,& ! temporal radius of infl. towards the past for aircraft data
    wtukare     ,& ! temporal radius of infl. towards the future for aircr.data
    wtuksua     ,& ! temporal radius of infl. towards the past for surface data
    wtuksue     ,& ! temporal radius of influence towards the future

! 2.2  quality control
!      ---------------
    dtqc        ,& ! timestep (in [s]) for the threshold quality control

! 2.3  observation processing
!      ----------------------
    doromx      ,& !  cut-off and gaussian radius of height differences 
                   !  between model orography and station height
    nolbc          !  number of grid rows at lateral boundaries
                   !  where obs are neglected

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &
 
!       1.1    ODR header format
!              -----------------
    mxrhed      ,& ! header length of multi-level reports
    mxshed      ,& ! header length of single-level reports
    mxghed      ,& ! header length of GPS reports
    mxthed      ,& ! header length of satellite retrieval reports
    nhtime      ,& ! time of observat. in forecast hours
    mxrhdf      ,& ! header length of multi-level reports
    mxshdf      ,& ! header length of single-level reports
    mxghdf      ,& ! header length of GPS reports
    mxthdf      ,& ! header length of satellite retrieval reports
    nhio        ,& ! (local) x-coord. of grid pt. assigned to obs
    nhjo        ,& ! (local) y-coord. of grid pt. assigned to obs
    nhobtp      ,& ! observation type
!   nhcode      ,& ! code type
!   nhpass      ,& ! flag for report being set to 'passive'
    nhqcfw      ,& ! threshold quality control flags for pressure, status of
                   ! verification output
!   nhnlev      ,& ! number of obs. levels (for multi-level reports)

!       1.2    ODR body format
!              ---------------
    maxrsl      ,& ! max. number of levels in multi-level ODR
    maxrtv      ,& ! max. number of levels in satellite retrieval  reports
    mxrbdy      ,& ! body length of multi-level reports
    mxrbdf      ,& ! body length of multi-level reports
    nbtp        ,& ! pressure [Pa]
    mxsbdy      ,& ! body length of single-level reports
    mxsbdf      ,& ! body length of single-level reports
    nbsp        ,& ! pressure                                           [Pa]
    mxgbdy      ,& ! body length of GPS reports
    mxgbdf      ,& ! body length of GPS reports
    nbgp        ,& ! pressure [Pa]
    mxtbdy      ,& ! body length of multi-level reports
    mxtbdf      ,& ! body length of sat retrieval reports
    nbvp           ! pressure [Pa]

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------

    imdi        ,& ! missing data indicator for ODR integers (2^31-1)
!   ntotml      ,& ! tot. number of stored multi-level reports
!   ntotsg      ,& ! tot. number of stored single-level reports
!   ntotgp      ,& ! tot. number of stored GPS reports
    ntottv      ,& ! tot. number of stored satellite retrievals
    fdoro       ,& ! scaling factor to vertical distances betw. model


!       2.     Observation data records (ODR)
!       -------------------------------------

    omlbdy      ,& ! body   of multi-level ODR
    omlhed      ,& ! header of multi-level ODR
    osgbdy      ,& ! body   of single-level ODR
    osghed      ,& ! header of single-level ODR
    ogpbdy      ,& ! body   of GPS ODR
    ogphed      ,& ! header of GPS ODR
    otvbdy      ,& ! body   of satellite retrieval ODR
    otvhed      ,& ! header of satellite retrieval ODR
    momlbd      ,& ! body   of multi-level ODR 
    momlhd      ,& ! header of multi-level ODR 
    mosgbd      ,& ! body   of single-level ODR 
    mosghd      ,& ! header of single-level ODR 
    mogpbd      ,& ! body   of GPS ODR 
    mogphd      ,& ! header of GPS ODR 
    motvbd      ,& ! body   of satellite retrieval ODR
    motvhd      ,& ! header of satellite retrieval ODR
    yomlhd      ,& ! header of multi-level ODR
    yosghd      ,& ! header of single-level ODR
    yogphd      ,& ! header of GPS ODR
    yotvhd         ! header of satellite retrieval ODR

! end of data_obs_record

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie           ,& ! number of grid points in zonal direction
    je           ,& ! number of grid points in meridional direction
    ke           ,& ! number of grid points in vertical direction
    ie_tot       ,& ! number of grid points in zonal direction
    je_tot       ,& ! number of grid points in meridional direction

  
! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    dlon         ,& ! grid point distance in zonal direction (in degrees)
    dlat         ,& ! grid point distance in meridional direction (in degrees)
    startlat_tot ,& ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlon_tot ,& ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    degrad       ,& ! factor for transforming degree to rad


! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt           ,& ! long time-step
    dtdeh        ,& ! dt / 3600 seconds

! 8. Organizational variables to handle the COSMO humidity tracers
! ----------------------------------------------------------------
    idt_qv,  idt_qc

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 2. physical constants and related variables
! -------------------------------------------

    r_d          ,& ! gas constant for dry air
    rdocp        ,& ! r_d / cp_d
    rvd_m_o      ,& ! r_v/r_d - 1
    rdv          ,& ! r_d / r_v
    o_m_rdv      ,& ! 1 - r_d/r_v
    g            ,& ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------

    b1           ,& ! variables for computing the saturation vapour pressure
    b2w, b2i     ,& ! over water (w) and ice (i)
    b3           ,& !               -- " --
    b4w, b4i        !               -- " --

! end of data_constants

!-------------------------------------------------------------------------------

USE data_fields     , ONLY :   &

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

    dp0          ,& ! reference pressure thickness of layers        ( Pa)
    p0           ,& ! reference pressure at main levels             ( Pa)
    hhl          ,& ! geometical height of half model levels        ( m )

! 2. external parameter fields                                        (unit)
! ----------------------------

    hsurf        ,& ! height of surface topography                  ( m   )
    fr_land      ,& ! fraction of land in a grid element              --

! 3. prognostic variables
! -----------------------
    pp           ,& ! deviation from the reference pressure         ( pa  )
    t            ,& ! temperature

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------
    qrs             ! precipitation water (water loading)           (kg/kg)

! end of data_fields

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep       ,& ! actual time step
                    ! indices for permutation of three time levels
    nnew            ! corresponds to ntstep + 1

! end of data_runcontrol

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute  ,& ! number of compute PEs
    nboundlines  ,& ! number of boundary lines of the domain for which
                    ! no forecast is computed = overlapping boundary
                    ! lines of the subdomains
    my_cart_id      ! rank of this subdomain in the cartesian communicator

! end of data_parallel

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    gather_field       ! gathers the parts of a total field from all subdomains

!-------------------------------------------------------------------------------

 USE src_tracer,               ONLY :  &
    trcr_get, trcr_errorstr

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "obs_aof_assign_gridpt.incf"
! INCLUDE "get_global_surf_aof.incf"
! INCLUDE "obs_pointrs.incf"
! INCLUDE "obs_fix_error.incf"
! INCLUDE "obs_error.incf"
! INCLUDE "obs_del_old_report.incf"
! INCLUDE "zsurf.incf"
! INCLUDE "psurf.incf"
! INCLUDE "tupair.incf"
! INCLUDE "rhupair.incf"
! INCLUDE "z2p.incf"
! INCLUDE "tv2t.incf"
! INCLUDE "p2dp.incf"


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for assigning reports to grid points
!-------------------------------------------------------------------------------

SUBROUTINE obs_aof_assign_gridpt ( zio_tot, zjo_tot, isthgh, kobtyp , lseaobs  &
                                 ,  io_tot,  jo_tot,  zsurf, lsfcob )

!-------------------------------------------------------------------------------
! Description:
!   Assign observation location horizontally to a model mass grid point.
!   Remark:  nolbc: number of grid rows at the lateral boundary,
!                   where no observations are used.
!
! Method:
!   For sea observations, for reports with missing station height, or if the
!   search radius specified is negative, then the report is assigned to model
!   (sea) grid point, which is closest in the horizontal.
!   Otherwise, the report is assigned to the grid point with the smallest
!   weighted height difference to the observation station within the positive
!   search radius, unless there is a grid point within half the mesh width
!   with a height difference of less than 40 m.
!   The horizontal search radius is  (SQRT(2) * mesh width)  for TEMP, PILOT,
!   SYNOP and DRIBU. This choice will often yield small height differences, 
!   which is desirable for the use of both surface-level data and of the
!   vertical structure of sounding data, while the horizontal offset will not
!   exceed the maximum offset when the search was limited to the 4 surrounding
!   grid points.
!   (Note: Code is extended to T2m and to search at latitude > 60 deg.) 04.11.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Christoph Schraff
!  Initial release
! 1.19       1998/12/11 Christoph Schraff
!  Reports of stations with too large a height difference to the model orography
!  are prepared to be set to passive (rather than being rejected).
!  ANSI violations removed.
! 1.27       1999/03/29 Christoph Schraff
!  32-bit missing data indicator in the ODR.
! 1.31       1999/07/01 Christoph Schraff
!  Bug correction ('je_tot' has been used to limit 'ieseek').
! 1.39       2000/05/03 Ulrich Schaettler
!  Changed names for variables concerned to latitude or longitude.
! 1.40       2000/05/23 Christoph Schraff
!  Bug correction for rejecting data at the boundary of the inner model domain.
! 3.7        2004/02/18 Ulrich Schaettler
!  Replaced zcphi by zcrlat (due to global renaming)
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Declarations:
!
!===============================================================================

IMPLICIT NONE

!===============================================================================

! Subroutine arguments
! --------------------

  REAL    (KIND=ireals)    , INTENT (IN)    ::       & 
    zio_tot, zjo_tot    ! exact obs. location in grid point units

  INTEGER (KIND=iintegers) , INTENT (IN)    ::       &
    kobtyp              ! observation type

  INTEGER (KIND=iintegers) , INTENT (INOUT) ::       &
    isthgh              ! station height in metres (may be modified for Sat obs)

  LOGICAL                  , INTENT (INOUT) ::       &
    lseaobs             ! true, if obs. released from ship or dribu (IN)
                        ! true, if obs. is treated as sea obs. (OUT)
! output
  INTEGER (KIND=iintegers) , INTENT (OUT)   ::       &
    io_tot  , jo_tot    ! grid point assigned to observation.
                        ! - if io_tot = jo_tot = 0 THEN obs. report is outside
                        !     horizontal limits
                        ! - if io_tot , jo_tot < 0 THEN obs. report is outside
                        !     vertical limits or land / sea type (for synops) 

  REAL    (KIND=ireals)    , INTENT (OUT)   ::       & 
    zsurf               ! height of model grid pt. to which obs. is assigned

  LOGICAL                  , INTENT (OUT)   ::       &
    lsfcob              ! false, if no surface data are to be used from profile

! Local parameters:
! ---------------- 

  REAL    (KIND=ireals)    , PARAMETER  :: &
    rtempmx =  1.414_ireals ,& ! horizontal search radius to assign a TEMP /
                               !   PILOT to a grid pt. (if (rtempmx < 0), search
                               !   is limited to the 4 neighbouring grid pts.)
    rsypmx  =  1.414_ireals    ! as 'rtempmx', but for surface-level reports

! Local variables
! ---------------

  INTEGER (KIND=iintegers)               ::       &
    iolb    , jolb   ,& ! observation location in total grid point units
    nlandgp          ,& ! number of nearby land grid points
    iaseek  , ieseek ,& ! grid point search range
    jaseek  , jeseek ,& !
    i,j                 ! loop indices          

  REAL    (KIND=ireals)                  ::       &  
    zfsea            ,& ! maximum sea fraction of candidate grid points
    rseekmx          ,& ! horizontal search radius for grid pts. to assign 
                        ! a observation to a gridpoint
    vfac             ,& ! modification factor for height difference
    zfio   , zfjo    ,& ! factors for bilinear interpolation of sat data
    zdo              ,& ! horizontal distance between obs and grid point
    zdomx            ,& ! maximum distance
    zrlats           ,& ! gridpoint latitude in degrees
    zcrlat           ,& ! cosine of grid point latitude
    fisdmn           ,& ! minimum height difference
    fisd             ,& ! height difference
    fisdps           ,& ! weighted height difference for surface pressure obs.
    fisdrh           ,& ! weighted height difference for rel. humidity observat.
    fisdtt           ,& ! weighted height difference for temperature observat.
    fisduv              ! weighted height difference for wind observations

  LOGICAL                                ::       &
    linarea          ,& ! .true, if obs inside inner model area
    lupair              ! .true, if upper-air report

  CHARACTER (LEN=20) yroutine ! name of subroutine

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE obs_aof_assign_gridpt
!-------------------------------------------------------------------------------
  yroutine = 'obs_aof_assign_gridpt'

!  Number of lateral grid rows where obs are neglegted
   nolbc   = MAX( nolbc , INT( 1+nboundlines ,iintegers) )
 
   iolb    = INT( zio_tot )
   jolb    = INT( zjo_tot )
   io_tot  = 0
   jo_tot  = 0
   zsurf   = c0
   linarea =      (zio_tot > nolbc+epsy) .AND.(zio_tot < ie_tot-nolbc+c1-epsy) &
            .AND. (zjo_tot > nolbc+epsy) .AND.(zjo_tot < je_tot-nolbc+c1-epsy)
   IF (.NOT. linarea)                                                     RETURN

   lsfcob  = .TRUE.
   lupair  = (kobtyp == ntemp) .OR. (kobtyp == npilot) .OR. (kobtyp == nairep) &
                               .OR. (kobtyp == nsatem) .OR. (kobtyp == nsatob)

! Re-determine 'lseaobs' (i.e. whether land or sea observation),
! to account e.g. for ships on rivers or synops/temps on tiny islands
! (input value of 'lseaobs' depends on the obs type (sea obs if ship or buoy);
!  - if a surface-level sea obs has to be assigned to a land grid point
!    then 'lseaobs' is set false to prevent its further use
!  - if an upper-air sea obs (TEMPSHIP) is assigned to a land grid point with
!    < 5 % water fraction, 'lseaobs' is set false to prevent its further use
!  - if an upper-air sea obs (TEMPSHIP) is assigned to a land grid point with
!    > 5 % water fraction, only 'lsfcob' is set false (surface level not used)
!  - land obs types may be assigned to a sea grid point but will always be used
!    (e.g. Ekofisk) )
   nlandgp = NINT( fland_tot(iolb,jolb  ) ) + NINT( fland_tot(iolb+1,jolb  ) ) &
           + NINT( fland_tot(iolb,jolb+1) ) + NINT( fland_tot(iolb+1,jolb+1) )

   zfsea   = c1 - MIN( fland_tot(iolb,jolb  ) , fland_tot(iolb+1,jolb  )       &
                     , fland_tot(iolb,jolb+1) , fland_tot(iolb+1,jolb+1) )
!  PRINT *,'zfsea ', zfsea, zio_tot, zjo_tot, nlandgp, lseaobs, lupair
   IF (      (lseaobs) .AND. (nlandgp == 4)                                    &
       .AND. (.NOT. ((lupair) .AND. (zfsea > 0.01_ireals))))  lseaobs = .FALSE.

! If one adjacent model grid point is a land point, then the observation
! is considered as land observation
   IF ((.NOT. lseaobs) .AND. (nlandgp == 0)) lseaobs = .TRUE.

! A Sat Retrieval over land should be assigned a station height according to
! obs point given by the bilinear interpolation used for the first guess
   IF (kobtyp == nsattv) THEN
     zfio   = zio_tot - INT( zio_tot )
     zfjo   = zjo_tot - INT( zjo_tot )
     isthgh = NINT(   (c1-zfio)*(c1-zfjo) *hsurf_tot(iolb  ,jolb  )            &
                    +     zfio *(c1-zfjo) *hsurf_tot(iolb+1,jolb  )            &
                    + (c1-zfio)*    zfjo  *hsurf_tot(iolb  ,jolb+1)            &
                    +     zfio *    zfjo  *hsurf_tot(iolb+1,jolb+1) )
   ENDIF

   IF ((isthgh == imdi) .OR. (lseaobs) .OR. (kobtyp == nsattv)) THEN
     rseekmx = c0
   ELSEIF ((kobtyp == ntemp ) .OR. (kobtyp == npilot)) THEN
     rseekmx = rtempmx
     vfac    = c1
   ELSEIF ((kobtyp == nsynop) .OR. (kobtyp == ndribu)) THEN
     rseekmx = rsypmx
     vfac    = fdoro(5)
   ENDIF

! determine the index ranges containing all candidate grid points
   IF (rseekmx < epsy) THEN
! if (rseekmx <= 0) then take into account only adjacent grid pts
     jaseek = MAX( jolb - 0 , nolbc + 1 )
     jeseek = MIN( jolb + 1 , je_tot- nolbc )
     iaseek = MAX( iolb - 0 , nolbc + 1 )
     ieseek = MIN( iolb + 1 , ie_tot- nolbc )
     rseekmx = 1.42_ireals
   ELSE
     jaseek = MAX( jolb - INT(rseekmx)     , nolbc + 1 )
     jeseek = MIN( jolb + INT(rseekmx) + 1 , je_tot- nolbc )
     zcrlat = COS( MAX( ABS( startlat_tot + (jaseek-1) *dlat )                 &
                      , ABS( startlat_tot + (jeseek-1) *dlat ) ) * degrad )
     IF (zcrlat > epsy) THEN
       iaseek = MAX( iolb - INT(rseekmx /zcrlat)     , nolbc + 1 )
       ieseek = MIN( iolb + INT(rseekmx /zcrlat) + 1 , ie_tot- nolbc )
     ELSE
       iaseek = nolbc + 1
       ieseek = ie_tot - nolbc
     ENDIF
   ENDIF

   IF ((lseaobs) .AND. (nlandgp == 4)) THEN
! If a radiosonde ship has to be assigned to a land grid pt. with non-zero water
! fraction, it is used actively as a sea obs, but surface level is set passive
     zfsea = -c1
     iseekc:   DO  i = iaseek , ieseek
     jseekc:   DO  j = jaseek , jeseek
       IF (c1-fland_tot(i,j) > zfsea) THEN
         zfsea  = c1 - fland_tot(i,j)
         io_tot = i
         jo_tot = j
       ENDIF
     END DO jseekc
     END DO iseekc
     lsfcob = .FALSE.

   ELSEIF ((isthgh == imdi) .OR. (lseaobs)) THEN
! If there is no station height, assign the observation to the grid
! point which is closest in the horizontal
! the same applies for observations over the sea, if there are model
! sea grid points within the search radius.
     zdomx    = c2
     iseeks:   DO  i = iaseek , ieseek
     jseeks:   DO  j = jaseek , jeseek
       zrlats   = startlat_tot + (j-1) * dlat
       zcrlat   = COS ( zrlats * degrad )
       zdo      = ( (i-zio_tot)*zcrlat *dlon/dlat )**2 + (j-zjo_tot)**2
       zdo      = SQRT( zdo )
       IF (zdo <= zdomx) THEN
         IF ((isthgh == imdi) .OR. (NINT(fland_tot(i,j)) == 0)) THEN
           zdomx  = zdo
           io_tot = i
           jo_tot = j
           lsfcob = .TRUE.
! assign a sea report to a land grid point if necessary
         ELSEIF (io_tot == 0) THEN
           io_tot = i
           jo_tot = j
           lsfcob = .FALSE.
         ENDIF
       ENDIF
     END DO jseeks
     END DO iseeks

     IF (isthgh == imdi) lsfcob = .FALSE.
   ELSE
     fisdmn  = 10000.0_ireals
! ensure that search radius is large enough to find a land grid pt.
     rseekmx = MAX( ABS(rseekmx) , SQRT( c05 ) )
     IF (nlandgp <= 3) rseekmx = MAX( rseekmx , c1 )
     IF (nlandgp <= 2) rseekmx = MAX( rseekmx , SQRT(1.25_ireals) )
     IF (nlandgp <= 1) rseekmx = MAX( rseekmx , SQRT(2.0_ireals ) )
     iseekl:   DO  i = iaseek , ieseek
     jseekl:   DO  j = jaseek , jeseek
       zrlats   = startlat_tot + (j-1) * dlat
       zcrlat   = COS ( zrlats * degrad )
       zdo = ( (i-zio_tot) *zcrlat *dlon/dlat )**2 + (j-zjo_tot)**2
       zdo = SQRT( zdo )
       IF ((zdo <= rseekmx) .AND. (NINT(fland_tot(i,j)) == 1)) THEN
         fisd = isthgh - hsurf_tot(i,j)
         IF ((zdo < 0.25_ireals) .AND. (ABS(fisd) < 40.0_ireals)) THEN
           io_tot = i
           jo_tot = j
           fisdmn = 0._ireals
           lsfcob = .TRUE.
         ENDIF
! height difference 'fisd' is weighted with factor 'fdoro(5)' if the
! station is below model grid point. 'fisd' is always positive.
! for an extrapolation of 100 m and a temperature error of 12 k,
! the resulting pressure error is about 0.5 hpa.
         fisd = ((vfac-c1)/c2 + SIGN( (vfac+c1)/c2 , fisd )) * fisd
         IF (fisd < fisdmn) THEN
           io_tot = i
           jo_tot = j
           fisdmn = fisd
           lsfcob = .TRUE.
         ENDIF
! assign a land report to a sea grid point if necessary
       ELSEIF ((zdo <= rseekmx) .AND. (io_tot == 0)) THEN
         io_tot = i
         jo_tot = j
         lsfcob = .FALSE.
       ENDIF
     END DO jseekl
     END DO iseekl
   ENDIF

   IF (io_tot  > 0)    zsurf = hsurf_tot(io_tot,jo_tot)

   IF (((kobtyp == nsynop) .OR. (kobtyp == ndribu)) .AND. (.NOT. lsfcob)) THEN
     io_tot = - io_tot
     jo_tot = - jo_tot
   ENDIF
   IF ((io_tot <= 0) .OR. (isthgh == imdi) .OR. (kobtyp == nsattv))       RETURN
!  IF (     (((io_tot == 0) .OR. (jo_tot == 0)) .AND. (.NOT. linarea))         &
!      .OR. (isthgh == imdi))                                             RETURN

   fisd   = isthgh - zsurf
   fisdps = ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd) ) * fisd
   fisdrh = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd) ) * fisd
   fisdtt = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd) ) * fisd
   fisduv = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd) ) * fisd
   IF (      (fisdps > doromx(2)) .AND. (fisdrh > doromx(4))                   &
       .AND. (fisdtt > doromx(3)) .AND. (fisduv > doromx(1))) THEN
     IF ((kobtyp == nsynop) .OR. (kobtyp == ndribu)) THEN
       io_tot = - io_tot
       jo_tot = - jo_tot
     ELSEIF ((kobtyp == ntemp) .OR. (kobtyp == npilot)) THEN
       lsfcob = .FALSE.
     ENDIF
   ENDIF

!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE obs_aof_assign_gridpt

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for deletion of old observations
!-------------------------------------------------------------------------------

SUBROUTINE get_global_surf_aof

!-------------------------------------------------------------------------------
! Description:
!   Delete past observations, which are not used for nudging any more,
!   and put forward the current observations in the data records
!
! Method:
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       20.08.04   Original code.
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Local variables:
! ----------------
  INTEGER (KIND=iintegers) :: &
    irm                  ! error status variables
  CHARACTER (LEN=20)       :: &
    yroutine             ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg              ! error message
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin subroutine get_global_surf_aof
!-------------------------------------------------------------------------------

  yroutine = 'get_global_surf_aof'

  IF (num_compute > 1) THEN

    CALL gather_field ( hsurf     , ie     , je                                &
                      , hsurf_tot , ie_tot , je_tot , 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_cart_id, 11141, yerrmsg, yroutine, irm )

    CALL gather_field ( fr_land   , ie     , je                                &
                      , fland_tot , ie_tot , je_tot , 0, irm )
!   =================

    IF (irm /= 0) WRITE( yerrmsg,'(''ERROR at gather_field'',I6)' ) irm
    IF (irm /= 0) CALL model_abort (my_cart_id, 11142, yerrmsg, yroutine, irm )

  ELSE
    hsurf_tot = hsurf
    fland_tot = fr_land
  ENDIF

!-------------------------------------------------------------------------------
! End subroutine get_global_surf_aof
!-------------------------------------------------------------------------------
 
END SUBROUTINE get_global_surf_aof


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for finding the diag. array position
!-------------------------------------------------------------------------------

SUBROUTINE obs_pointrs ( kobtyp , kcdtyp )

!-------------------------------------------------------------------------------
! Description:
!   Find the position within the diagnostic arrays for a given
!   observation and code type.
!
! Method:
!   The 2-d diagnostic arrays pointers are being found as a function
!   of the observation type and the code type.
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       13.11.97   Original code.    Michael Buchhold
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------
! Modules used: These are declared in the module declaration section
! -------------

  IMPLICIT NONE

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT (IN)    ::         &
    kobtyp          ,&  ! observation type
    kcdtyp              ! code type

! Local variables
! ---------------
  INTEGER (KIND=iintegers) ::  &
    i1,i2               ! intermediate storage for obs type and code type

  CHARACTER (LEN=20)                    :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerrmsg             ! error message

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_pointrs
!-------------------------------------------------------------------------------
  yroutine = 'obs_pointrs'


  !----------------------------------------------------------------------
  ! Section 1.   Find out pointers of the statistic arrays
  !----------------------------------------------------------------------

  i1  =  kobtyp
  i2  =  MAX( kcdtyp , nmxcdt+1 )

! IF (i1 > ntotob) THEN
!   WRITE( yerrmsg,'(''wrong observation type:'',I10)' ) i1
! In fact, this type of error should never occur at PE /= 0 because for each
! report, 'obs_pointrs' is called the first time at reading from file by PE == 0
!   IF (my_cart_id /= 0)  PRINT '(A)', yerrmsg
!   CALL model_abort (my_cart_id, 1021, yerrmsg, yroutine)
! ENDIF

  IF (i1 == 1 ) THEN
!     obs. type 1
!     -----------
      IF (kcdtyp == nsrscd)      i2   =   1
      IF (kcdtyp == natscd)      i2   =   2
      IF (kcdtyp == nshscd)      i2   =   3
      IF (kcdtyp == nabscd)      i2   =   4
      IF (kcdtyp == nshred)      i2   =   5
      IF (kcdtyp == natshs)      i2   =   6
 
  ELSEIF (i1 == 2 ) THEN
!     obs. type 2
!     -----------
      IF (kcdtyp == naircd)      i2   =   2
      IF (kcdtyp == namdar)      i2   =   4
      IF (kcdtyp == ncodar)      i2   =   1
      IF (kcdtyp == ncolba)      i2   =   3
      IF (kcdtyp == nacar )      i2   =   5

  ELSEIF (i1 == 3 ) THEN
!     obs. type 3
!     -----------
      IF (kcdtyp == nstbcd)      i2   =   1
      IF (kcdtyp == nsst  )      i2   =   2

  ELSEIF (i1 == 4 ) THEN
!     obs. type 4
!     -----------
      IF (kcdtyp == ndrbcd)      i2   =   3
      IF (kcdtyp == nbathy)      i2   =   1
      IF (kcdtyp == ntesac)      i2   =   2

  ELSEIF (i1 == 5 ) THEN
!     obs. type 5
!     -----------
      IF (kcdtyp == nldtcd)      i2   =   1
      IF (kcdtyp == nshtcd)      i2   =   2
      IF (kcdtyp == ntdrop)      i2   =   3
      IF (kcdtyp == nrocob)      i2   =   4
      IF (kcdtyp == nrocsh)      i2   =   5

  ELSEIF (i1 == 6 ) THEN
!     obs. type 6
!     -----------
      IF (kcdtyp == nldpcd)      i2   =   1
      IF (kcdtyp == nshpcd)      i2   =   2
      IF (kcdtyp == nwp_eu)      i2   =   3
      IF (kcdtyp == nra_eu)      i2   =   4
      IF (kcdtyp == npr_us)      i2   =   5
      IF (kcdtyp == nravad)      i2   =   6
 
  ELSEIF (i1 == 7 ) THEN
!     obs. type 7
!     -----------
      IF (kcdtyp == nstmcd)      i2   =   1
      IF (kcdtyp == nstovs)      i2   =   2
      IF (kcdtyp == nsmsg1)      i2   =   1
      IF (kcdtyp == nsmsg2)      i2   =   2
      IF (kcdtyp == nnoa15)      i2   =   3
      IF (kcdtyp == nnoa16)      i2   =   4
      IF (kcdtyp == nnoa17)      i2   =   5
      IF (kcdtyp == nnoa18)      i2   =   6
 
  ELSEIF ((i1 == 8 ) .OR. (i1 == 12)) THEN
!     obs. type 8 or 12 (ngps)
!     -----------
      i1 = 8
      IF (kcdtyp == ngpgfz)      i2   =   1

  ELSE
    WRITE( yerrmsg,'(''wrong observation type:'',I10)' ) i1
! In fact, this type of error should never occur at PE /= 0 because for each
! report, 'obs_pointrs' is called the first time at reading from file by PE == 0
    IF (my_cart_id /= 0)  PRINT '(A)', yerrmsg
    CALL model_abort (my_cart_id, 1021, yerrmsg, yroutine)

  ENDIF
 
  !----------------------------------------------------------------------
  ! Section 2.   Assign pointers
  !----------------------------------------------------------------------
 
  IF (i2 > nmxcdt) THEN
    WRITE( yerrmsg,'(''wrong observation code type:'',I10)' ) i2
    IF (my_cart_id /= 0)  PRINT '(A)', yerrmsg
    CALL model_abort (my_cart_id, 1022, yerrmsg, yroutine)
  ENDIF

  nobtpp =  i1
  ncdtpp =  i2

!------------------------------------------------------------------------
 
END SUBROUTINE obs_pointrs

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for observation error evaluation
!-------------------------------------------------------------------------------

SUBROUTINE obs_fix_error ( poberr , klevel , perror )

!-------------------------------------------------------------------------------
! Description:
!   Observation error evaluation
!
! Method:
!   For given pressure obs. error is evaluated at that pressure by
!     the linear interpolation in lnp.
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       17.12.97   Original code.    Michael Buchhold
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT (IN)    ::         &
    klevel             ! dimension of poberr

  REAL    (KIND=ireals)    , INTENT (IN)    ::         &
    poberr (klevel)    ! input array of observation errors

  REAL    (KIND=ireals)    , INTENT (OUT)   ::         &
    perror             ! output interpolated observation error

! Local variables
! ---------------
  INTEGER  (KIND=iintegers)            ::         &
    ilevu              ! level

  REAL (KIND=ireals)                   ::         &
    zlnlev             ! ln (rppp)
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin subroutine obs_fix_error
!-------------------------------------------------------------------------------

!    1.       Calculate observation error
!             ---------------------------
 
!       1.1   Keep observation error constant below lowest error level
  IF (rppp >= rlevel(1)) THEN
     perror = poberr(1)

!       1.2   Keep observation error constant above highest error level
  ELSEIF (rppp <= rlevel(nerlev)) THEN
     perror = poberr(nerlev)
 
!       1.3   find nearest analysis level above observation
  ELSE
     zlnlev = LOG(rppp)
     ilevu  = 2
     DO WHILE ((rolnlv(ilevu) > zlnlev) .AND. (ilevu <= nerlev))
        ilevu = ilevu + 1
     ENDDO
     IF (ABS(zlnlev-rolnlv(ilevu)) <  epsy) THEN
        perror = poberr(ilevu)

!       1.4   Interpolate error linearly in the vertical with respect to ln(p)
     ELSE
        perror = poberr(ilevu-1) + (poberr(ilevu) - poberr(ilevu-1))           &
                                  /(rolnlv(ilevu) - rolnlv(ilevu-1))           &
                                  *(zlnlev        - rolnlv(ilevu-1))
     ENDIF
  ENDIF

  END SUBROUTINE obs_fix_error

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for assigning an observation error
!-------------------------------------------------------------------------------

FUNCTION obs_error (kvar)

!-------------------------------------------------------------------------------
! Description:
!   Retrieve observation error
!
! Method:
!   The observation error is assigned depending on the observation
!   type and its vertical position.
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       17.12.97   Original code.    Michael Buchhold
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Function arguments:
! ------------------
!
  INTEGER (KIND=iintegers) , INTENT (IN)    ::    &
    kvar               ! observation type

! Local variables
! ---------------

  LOGICAL                                   ::    &
    lerr               ! .true. if variable type invalid
 
  REAL (KIND=ireals)                        ::    &
    zerror          ,& ! assigned obs error
    obs_error          ! return value of this funktion

  CHARACTER (LEN=20)                        ::    &
    yroutine           ! name of this subroutine
  CHARACTER (LEN=80)                        ::    &
    yerrmsg            ! error message
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function obs_error
!-------------------------------------------------------------------------------
  yroutine = 'obs_error  '

  lerr = .FALSE.

!        1.       Observation type 1 (SYNOP/SHIP)
!                 -------------------------------
 
!           1.1   Winds
  IF (nobtyp == 1) THEN
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
        CALL obs_fix_error ( oevsynp,nerlev,zerror )

!           1.2   Height (ships or land stations)
     ELSE IF (kvar == nvrpoe(3)) THEN
        IF (llship) THEN
           zerror=oezship
        ELSE
           CALL obs_fix_error ( oezsynp,nerlev,zerror )
        ENDIF
     ELSE
        lerr = .TRUE.
     ENDIF
 
!        2.       Observation type 2 (AIREP)
!                 --------------------------
 
  ELSE IF (nobtyp == 2) THEN
 
!           2.2   Winds
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
        CALL obs_fix_error ( oeairep,nerlev,zerror )
 
!           2.3   Temperature
     ELSEif (kvar == nvrpoe(6)) THEN
        CALL obs_fix_error ( oetairp,nerlev,zerror )
     ELSE
        lerr = .TRUE.
     ENDIF
 
!        3.       Observation type 3 (SATOB)
!                 --------------------------
 
!           3.1   Winds
  ELSE IF (nobtyp == 3) THEN
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
        CALL obs_fix_error ( oesatob,nerlev,zerror )
     ELSE
        lerr = .TRUE.
     ENDIF
  
!        4.       Observation type 4 (DRIBU)
!                 --------------------------
 
!           4.1   Winds
  ELSE IF (nobtyp == 4) THEN
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
         zerror=oevdrib
 
!           4.2   Height
     ELSE IF (kvar == nvrpoe(3)) THEN
         zerror=oezdrib
     ELSE
        lerr = .TRUE.
     ENDIF
 
!        5.       Observation type 5 (TEMP)
!                 -------------------------
 
  ELSE IF (nobtyp == 5) THEN
 
!           5.3   Winds
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
        CALL obs_fix_error ( oevsond,nerlev,zerror )
 
!           5.4   Height
     ELSE IF (kvar == nvrpoe(3)) THEN
           CALL obs_fix_error ( oezsond,nerlev,zerror )
 
!           5.5   temperature
     ELSE IF (kvar == nvrpoe(6)) THEN
        CALL obs_fix_error ( oetsond,nerlev,zerror )
     ELSE
!!      WRITE( nupr,'(''obs_error '', 2i4)' ) kvar, nvrpoe(6),nvrpoe(1)
        PRINT       '(''obs_error '', 2i4)' , kvar, nvrpoe(6),nvrpoe(1)
        lerr = .TRUE.
     ENDIF
 
!        6.       Observation type 6 (PILOT/PROFILER/RASS)
!                 ----------------------------------------
 
!           6.1   Winds
  ELSE IF (nobtyp == 6) THEN
     IF (kvar == nvrpoe(1) .OR. kvar == nvrpoe(2)) THEN
        CALL obs_fix_error ( oevsond,nerlev,zerror )

!           6.2   Temperature
     ELSE IF (kvar == nvrpoe(6)) THEN
        CALL obs_fix_error ( oetsond,nerlev,zerror )
     ELSE
        lerr = .TRUE.
     ENDIF
  ENDIF
 
!        9.       Final processing
!                 ----------------
 
!           9.1   Error processing
  IF (lerr)  THEN
!    PRINT      '('' VARIABLE NUMBER  = '',i10)', kvar
     WRITE( yerrmsg,'(''INVALID VARIABLE TYPE: OBS/CODE TYPE'',I10             &
                    &,''/'',I10,'', VAR. NUMBER'',I10)' )                      &
            nobtyp, ncdtyp, kvar
     IF (my_cart_id /= 0)  PRINT '(A,'', by PE'',I5)', yerrmsg, my_cart_id
     CALL model_abort (my_cart_id, 1023, yerrmsg, yroutine)
  ENDIF

!           9.2   Normal return
  obs_error  =   zerror

END FUNCTION obs_error

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_util" for deletion of old observations
!-------------------------------------------------------------------------------

SUBROUTINE obs_del_old_report

!-------------------------------------------------------------------------------
! Description:
!   Delete past observations, which are not used for nudging any more,
!   and put forward the current observations in the data records
!
! Method:
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       01.11.97   Original code.
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
!
! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Local variables:
! ----------------
  INTEGER (KIND=iintegers) ::  &
    nsypo, nraso, nchk,& ! loop indices
    ngpso             ,& ! loop indices
    ndiff, narest        ! indices

  REAL (KIND=ireals)       ::  &
    wtuke             ,& ! temporal radius of influence
    tmaxbox              ! max. timestep for computing analysis increments

  LOGICAL                  ::  &
    ldel                 ! delete past multi-level report
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin subroutine obs_del_old_report
!-------------------------------------------------------------------------------

  tmaxbox = INT( tconbox/dt   + c1-epsy ) * dtdeh
  tmaxbox = 0

!-------------------------------------------------------------------------------
!  Section 1: Multi-level reports
!-------------------------------------------------------------------------------

  ndiff  = 0
  DO nraso = 1 , nmlob
    wtuke = MAX( wtukrse , tipolmx )
    IF (momlhd(nraso,nhobtp) == nairep) wtuke = wtukare
    wtuke = MAX( wtuke , tmaxbox )
    ldel = (omlhed(nraso,nhtime) + wtuke  <=  ntstep*dtdeh)
    IF (momlhd(nraso,nhobtp) /= nairep) wtuke = MAX( wtukrse , tmaxbox )
! If - a past report could only be still needed for temporal interpolation, and
!    - all reports within the time range for potential partners are already read
! then delete the report unless a partner for temporal interpolation is found
! then check for partners for linear temporal interpolation
    IF (      (.NOT. ldel) .AND. (tipolmx > wtuke)                             &
        .AND. (omlhed(nraso,nhtime) + wtuke  <=  ntstep*dtdeh)                 &
        .AND. (omlhed(nraso,nhtime) + tipolmx <  naoflbx+naofbox-epsy)) THEN
      ldel = .TRUE.
! Note: TEMP / PILOT at the same location are sorted according to the obs. time
      DO nchk = nraso+1 , nmlob
        IF (      (momlhd(nchk,nhio)   == momlhd(nraso,nhio))                  &
            .AND. (momlhd(nchk,nhjo)   == momlhd(nraso,nhjo))                  &
            .AND. (omlhed(nchk,nhtime) <= omlhed(nraso,nhtime)+tipolmx)        &
            .AND. (yomlhd(nchk)        == yomlhd(nraso)     ))  ldel = .FALSE.
      ENDDO
    ENDIF
    IF ((lverif) .AND. (momlhd(nraso,nhqcfw) > -99)                            &
                 .AND. (omlhed(nraso,nhtime) >= hversta-epsy)                  &
                 .AND. (omlhed(nraso,nhtime) <= hverend+epsy)) THEN
      IF (omlhed(nraso,nhtime) + (dtqc+tconbox)/c3600  <=  ntstep*dtdeh) THEN
        PRINT '(''CAUTION: m-l report '',A ,'' ,'',F5.2,'' [hrs],'',F7.0       &
              &,'' [Pa] not verified at timestep'',I4)' ,                      &
          yomlhd(nraso), omlhed(nraso,nhtime), omlbdy(nraso,1,nbtp), ntstep
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) THEN
      ndiff = ndiff + 1
    ELSE
      omlbdy (nraso-ndiff,1:maxrsl,1:mxrbdy) = omlbdy (nraso,1:maxrsl,1:mxrbdy)
      momlbd (nraso-ndiff,1:maxrsl,1:mxrbdf) = momlbd (nraso,1:maxrsl,1:mxrbdf)
      omlhed (nraso-ndiff,1:mxrhed)          = omlhed (nraso,1:mxrhed)
      momlhd (nraso-ndiff,1:mxrhdf)          = momlhd (nraso,1:mxrhdf)
      yomlhd (nraso-ndiff)                   = yomlhd (nraso)
    ENDIF
  ENDDO

  narest = nmlob - ndiff + 1
  omlbdy (narest:nmlob,1:maxrsl,1:mxrbdy) = rmdi
  momlbd (narest:nmlob,1:maxrsl,1:mxrbdf) = imdi
  omlhed (narest:nmlob,1:mxrhed)          = rmdi
  momlhd (narest:nmlob,1:mxrhdf)          = imdi
  yomlhd (narest:nmlob)                   = '        '
 
  nmlob = nmlob - ndiff

!-------------------------------------------------------------------------------
!  Section 2: Single-level reports
!-------------------------------------------------------------------------------
 
  ndiff  = 0
  DO  nsypo = 1 , nsgob
    wtuke = MAX( wtuksue , tipmxsu )
    IF ((mosghd(nsypo,nhobtp) == nairep) .OR.(mosghd(nsypo,nhobtp) == nsatob)) &
      wtuke = wtukare
    wtuke = MAX( wtuke , tmaxbox )
    ldel = (osghed(nsypo,nhtime) + wtuke  <=  ntstep*dtdeh)
    IF ((lverif) .AND. (mosghd(nsypo,nhqcfw) > -99)                            &
                 .AND. (osghed(nsypo,nhtime) >= hversta-epsy)                  &
                 .AND. (osghed(nsypo,nhtime) <= hverend+epsy)) THEN
      IF (osghed(nsypo,nhtime) + (dtqc+tconbox)/c3600  <=  ntstep*dtdeh) THEN
        PRINT '(''CAUTION: s-l report '',A ,'' ,'',F5.2,'' [hrs],'',F7.0       &
              &,'' [Pa] not verified at timestep'',I4)' ,                      &
              yosghd(nsypo), osghed(nsypo,nhtime), osgbdy(nsypo,nbsp), ntstep
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) THEN
      ndiff = ndiff + 1
    ELSE
      osgbdy (nsypo-ndiff,1:mxsbdy) = osgbdy (nsypo,1:mxsbdy)
      mosgbd (nsypo-ndiff,1:mxsbdf) = mosgbd (nsypo,1:mxsbdf)
      osghed (nsypo-ndiff,1:mxshed) = osghed (nsypo,1:mxshed)
      mosghd (nsypo-ndiff,1:mxshdf) = mosghd (nsypo,1:mxshdf)
      yosghd (nsypo-ndiff)          = yosghd (nsypo)
    ENDIF
  ENDDO
  narest = nsgob - ndiff + 1
  osgbdy (narest:nsgob,1:mxsbdy) = rmdi
  mosgbd (narest:nsgob,1:mxsbdf) = imdi
  osghed (narest:nsgob,1:mxshed) = rmdi
  mosghd (narest:nsgob,1:mxshdf) = imdi
  yosghd (narest:nsgob)          = '        '
 
  nsgob = nsgob - ndiff

!-------------------------------------------------------------------------------
!  Section 3: GPS reports
!-------------------------------------------------------------------------------
 
  ndiff  = 0
  DO  ngpso = 1 , ngpob
    wtuke = MAX( wtuksue , tipmxsu , tmaxbox )
    ldel = (ogphed(ngpso,nhtime) + wtuke  <=  ntstep*dtdeh)
    IF ((lverif) .AND. (mogphd(ngpso,nhqcfw) > -99)                            &
                 .AND. (ogphed(ngpso,nhtime) >= hversta-epsy)                  &
                 .AND. (ogphed(ngpso,nhtime) <= hverend+epsy)) THEN
      IF (ogphed(ngpso,nhtime) + (dtqc+tconbox)/c3600  <=  ntstep*dtdeh) THEN
        PRINT '(''CAUTION: GPS report '',A ,'' ,'',F5.2,'' [hrs],'',F7.0       &
              &,'' [Pa] not verified at timestep'',I4)' ,                      &
              yogphd(ngpso), ogphed(ngpso,nhtime), ogpbdy(ngpso,nbgp), ntstep
      ELSE
        ldel = .FALSE.
      ENDIF
    ENDIF
    IF (ldel) THEN
      ndiff = ndiff + 1
    ELSE
      ogpbdy (ngpso-ndiff,1:mxgbdy) = ogpbdy (ngpso,1:mxgbdy)
      mogpbd (ngpso-ndiff,1:mxgbdf) = mogpbd (ngpso,1:mxgbdf)
      ogphed (ngpso-ndiff,1:mxghed) = ogphed (ngpso,1:mxghed)
      mogphd (ngpso-ndiff,1:mxghdf) = mogphd (ngpso,1:mxghdf)
      yogphd (ngpso-ndiff)          = yogphd (ngpso)
    ENDIF
  ENDDO
  narest = ngpob - ndiff + 1
  ogpbdy (narest:ngpob,1:mxgbdy) = rmdi
  mogpbd (narest:ngpob,1:mxgbdf) = imdi
  ogphed (narest:ngpob,1:mxghed) = rmdi
  mogphd (narest:ngpob,1:mxghdf) = imdi
  yogphd (narest:ngpob)          = '        '
 
  ngpob = ngpob - ndiff

!-------------------------------------------------------------------------------
 
END SUBROUTINE obs_del_old_report


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for approx. model surface pressure
!-------------------------------------------------------------------------------

FUNCTION zsurf ( iob , job )

!-------------------------------------------------------------------------------
! Description:
!   Get the orographic height at a specified grid point.
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
! 
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)      ::        &
    iob     , job       ! grid point for which pressure is sought

! Local variables:
! ---------------
  REAL (KIND=ireals)      ::      &
    zsurf               ! return value of function
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function zsurf
!-------------------------------------------------------------------------------
 
  zsurf  =  hsurf(iob,job)

!-------------------------------------------------------------------------------
! End of module function zsurf
!-------------------------------------------------------------------------------
 
END FUNCTION zsurf


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for approx. model surface pressure
!-------------------------------------------------------------------------------

FUNCTION psurf ( iob , job )

!-------------------------------------------------------------------------------
! Description:
!   Get the approximate model surface pressure at a specified grid points
!   (in order to reject AIRCRAFT or SATOB observations that are close to the
!    ground - called if height is not reported).
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
! 
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)      ::        &
    iob     , job       ! grid point for which pressure is sought

! Local variables:
! ---------------
  REAL (KIND=ireals)      ::      &
    psurf               ! return value of function
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function psurf
!-------------------------------------------------------------------------------
 
  psurf  =  p0(iob,job,ke) + pp(iob,job,ke,nnew) + c05* dp0(iob,job,ke)

!-------------------------------------------------------------------------------
! End of module function psurf
!-------------------------------------------------------------------------------
 
END FUNCTION psurf


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for upper-air model temperature
!-------------------------------------------------------------------------------

FUNCTION tupair ( iob , job , kob )

!-------------------------------------------------------------------------------
! Description:
!   Get the model temperature at a specified grid point and model level.
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)      ::        &
    iob , job , kob     ! 3-D grid point for which temperature is wanted

! Local variables:
! ---------------
  REAL (KIND=ireals)      ::      &
    tupair              ! return value of function
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function tupair
!-------------------------------------------------------------------------------

  tupair  =  t(iob,job,kob,nnew)

!-------------------------------------------------------------------------------
! End of module function tupair
!-------------------------------------------------------------------------------

END FUNCTION tupair


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for upper-air model relative humidity
!-------------------------------------------------------------------------------

FUNCTION rhupair ( iob , job , kob )

!-------------------------------------------------------------------------------
! Description:
!   Get the model relative humdity at a specified grid point and model level.
!
! Current Code Owner:  Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
!
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)      ::        &
    iob , job , kob     ! 3-D grid point for which relative humidity is wanted

! Local variables:
! ---------------
  REAL (KIND=ireals)      ::      &
    rhupair             ! return value of function

  REAL (KIND=ireals)       ::  &
    fspw  , fq2pv ,& ! statement functions for comput. of saturation humidity
    zt, zp, zqv      ! dummy arguments of statement functions

  INTEGER ( KIND=iintegers)   :: izerror
  CHARACTER (LEN=80)          :: yzerrmsg

  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:) => NULL()      ! QV at nnew

  CHARACTER (LEN=25)          :: yroutine = 'rhupair'
!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function rhupair
!-------------------------------------------------------------------------------

! Statement - functions for computation of saturation humidity
! ------------------------------------------------------------

! magnus formela for water
  fspw(zt)        = b1 * exp( b2w*(zt - b3)/(zt - b4w) )
! vapour pressure (at T = saturation vapour press. at Td) from specific humidity
  fq2pv (zqv,zp)  = MAX( epsy , zqv ) * zp / (rdv + zqv*o_m_rdv)

  ! retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF

  rhupair  =   fq2pv( qv(iob,job,kob)                                          &
                    , pp(iob,job,kob,nnew) + p0(iob,job,kob) )                 &
             / fspw ( t (iob,job,kob,nnew) )  

!-------------------------------------------------------------------------------
! End of module function rhupair
!-------------------------------------------------------------------------------

END FUNCTION rhupair


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for model pressure at a specified height
!-------------------------------------------------------------------------------

FUNCTION z2p ( zzz , iob , job )

!-------------------------------------------------------------------------------
! FUNCTION src_obs_proc_util   (<InputArguments>) &
!   RESULT (<ResultName>) ! The use of RESULT is recommended
!-------------------------------------------------------------------------------
! Description:  <Say what this routine does>
! Method:    <Say how it does it: include references to external documentation>
!            <If this routine is divided into sections, be brief here,
!             and put Method comments at the start of each section>
! Input files:  <Describe these, and say in which routine they are read>
! Output files: <Describe these, and say in which routine they are written>
!-------------------------------------------------------------------------------
! Description:
!   Get the model pressure for a specified height (in order to assign PILOT or 
!   SATOB reports to a pressure level, if pressure itself is not reported).
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       05.12.97   Original code.    Michael Buchhold
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
! 
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)   :: &
    iob     , job       ! grid point for which pressure is sought
 
  REAL    (KIND=ireals)    , INTENT (IN)   :: &
    zzz                 ! height for which pressure is sought
 

! Local variables:
! ---------------

  INTEGER (KIND=iintegers) ::      &
    k                   ! loop index

  REAL (KIND=ireals)      ::      &
    rppp             ,& ! pressure
    hhke             ,& ! height of lowest full model level
    zzo     , zzu    ,& ! height of upper / lower model level
    zfrac            ,& ! weight
    aloppp           ,& ! log of pressure
    z2p                 ! return value of function

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function z2p
!-------------------------------------------------------------------------------
 
  rppp  =  rmdi
  hhke  =  0.5_ireals*( hhl(iob,job,ke) + hhl(iob,job,ke+1) )
  zzo   =  hhke

  IF (zzz >= hsurf(iob,job)) THEN

! Find appropriate model layer for the specified height
    DO  k = ke, 2, -1
      zzu = zzo
      zzo = 0.5_ireals*( hhl(iob,job,k) + hhl(iob,job,k-1) )
      IF ((zzu <= zzz) .AND. (zzo > zzz)) THEN

! interpolate ln(p) linear in z
        zfrac = ( zzz - zzu ) / ( zzo - zzu )
        aloppp= (1._ireals-zfrac)* LOG( p0(iob,job,k  )+ pp(iob,job,k  ,nnew)) &
                         + zfrac * LOG( p0(iob,job,k-1)+ pp(iob,job,k-1,nnew))
        rppp  = EXP( aloppp )
        EXIT
      ENDIF
    ENDDO
    IF ( rppp == rmdi)   THEN
      IF (zzz < hhke .AND. zzz >= hsurf(iob,job) ) THEN

! reduce pressure from lowest model level to the specified height
        rppp = EXP (  LOG( p0(iob,job,ke)+ pp(iob,job,ke,nnew) )               &
                    + g/(r_d*t(iob,job,ke,nnew)) *(hhke-zzz) )
      ENDIF
    ENDIF
  ENDIF

  z2p  =  rppp

!-------------------------------------------------------------------------------
! End of module function z2p
!-------------------------------------------------------------------------------
 
END FUNCTION z2p

!-------------------------------------------------------------------------------
!+ Module function in "src_obs_proc_util" for temperature from virtual temp.
!-------------------------------------------------------------------------------

FUNCTION tv2t ( iob, job, ztv, zzz )

!-------------------------------------------------------------------------------
! Description:
!   Convert virtual temperature (from RASS / SODAR) to air temperature
!   by use of the model humidity
!
! Method:
!   Find nearest model levels for the observed height and interpolate
!   model humidity. Convert temperature.
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       19.08.04   Original code based on a routine of Michael Buchhold.
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------

IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT (IN)   :: &
    iob    , job        ! horiz. coordinates of grid point

  REAL    (KIND=ireals)    , INTENT (IN)   :: &
    ztv              ,& ! virtual temperature
    zzz                 ! height

! Local parameters: none
! -----------------

! Local variables:
! ----------------
  INTEGER (KIND=iintegers)                :: &
    k                   ! loop index

  REAL (KIND=ireals)                      :: &
    ztt              ,& ! temperature
    zfq              ,& ! interpolated model humidity factor
    hhke             ,& ! height of lowest full model level
    zzo     , zzu    ,& ! height of upper / lower model level
    zfrac            ,& ! weight
    aloqc            ,& ! log of specific humidity
    tv2t                ! return value of function

  INTEGER ( KIND=iintegers) :: izerror
  CHARACTER (LEN=80)        :: yzerrmsg

  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:)  => NULL()      ,& ! QV at tlev=nnew
    qc  (:,:,:)  => NULL()         ! QC at tlev=nnew

  CHARACTER (LEN=25)        :: yroutine = 'tv2t'

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine tv2t
!-------------------------------------------------------------------------------

  ! retrieve the required microphysics tracers (at nnew)
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnew, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnew, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yroutine)
  ENDIF

! get interpolated model humidity
! -------------------------------

  zfq   =  rmdi
  hhke  =  0.5_ireals*( hhl(iob,job,ke) + hhl(iob,job,ke+1) )
  zzo   =  hhke

  IF (zzz >= hsurf(iob,job)) THEN

! Find appropriate model layer for the specified height
    DO  k = ke, 2, -1
      zzu = zzo
      zzo = 0.5_ireals*( hhl(iob,job,k) + hhl(iob,job,k-1) )
      IF ((zzu <= zzz) .AND. (zzo > zzz)) THEN

! interpolate ln(qc) linear in z
        zfrac = ( zzz - zzu ) / ( zzo - zzu )
        aloqc = (c1-zfrac)*LOG( MAX( rvd_m_o*qv(iob,job,k)                    &
                                  -qc(iob,job,k) -qrs(iob,job,k),epsy))       &
               +    zfrac *LOG( MAX( rvd_m_o*qv(iob,job,k-1)                  &
                                  -qc(iob,job,k-1) -qrs(iob,job,k-1),epsy))
        zfq   = EXP( aloqc )
        EXIT
      ENDIF
    ENDDO
  ENDIF
  IF ((zfq < rmdich) .AND. (zzz < hhke)) THEN
!   use specific humidity at lowest model level
    zfq  = qv(iob,job,ke)
    PRINT '("use humidity at lowest model level: zzz, hhke, hsurf: ",3F10.0)', &
                                                 zzz, hhke, hsurf(iob,job)
  ENDIF

! calculate temperature
! ---------------------

  IF (zfq > rmdich) THEN
    ztt = ztv / (c1 + zfq)
  ELSE
    PRINT '("tv2t: no model level found: zzz, hhke, hsurf: ",3F10.0)',         &
                                         zzz, hhke, hsurf(iob,job)
    ztt = rmdi
  ENDIF

  tv2t  =  ztt

!-------------------------------------------------------------------------------
! End of module function tv2t
!-------------------------------------------------------------------------------

END FUNCTION tv2t


!-------------------------------------------------------------------------------
!+ Module funct. in "src_obs_proc_util" for model pressure at a specified height
!-------------------------------------------------------------------------------

FUNCTION p2dp ( zpp , iob , job )

!-------------------------------------------------------------------------------
! Description:
!   Get the approximate model layer thickness at a specified pressure level.
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       07.09.04   Original code.    Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!-------------------------------------------------------------------------------
! 
! Modules used:         These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Parameter list:
! --------------
  INTEGER (KIND=iintegers) , INTENT (IN)   :: &
    iob     , job       ! grid point for which pressure is sought
 
  REAL    (KIND=ireals)    , INTENT (IN)   :: &
    zpp                 ! height for which pressure is sought
 

! Local variables:
! ---------------

  INTEGER (KIND=iintegers) ::      &
    kdp              ,& ! index of model layer encompassing pressure zpp
    kk                  ! loop index

  REAL (KIND=ireals)      ::      &
    p2dp                ! return value of function

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function p2dp
!-------------------------------------------------------------------------------

  kdp = ke - 1
  DO kk = ke-2, 2, -1
    IF (p0(iob,job,kk)+pp(iob,job,kk,nnew)+ c05*dp0(iob,job,kk) > zpp) kdp = kk
  ENDDO

  p2dp  =  dp0(iob,job,kdp) + c05 * (  pp(iob,job,kdp+1,nnew)                  &
                                     - pp(iob,job,kdp-1,nnew))

!-------------------------------------------------------------------------------
! End of module function p2dp
!-------------------------------------------------------------------------------
 
END FUNCTION p2dp

!-------------------------------------------------------------------------------

END MODULE src_obs_proc_util
