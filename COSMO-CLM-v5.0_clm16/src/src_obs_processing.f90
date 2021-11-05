!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_processing

!-------------------------------------------------------------------------------
! Description:
!   This module performs the observation processing. Special tasks are:
!    - reading reports from AOF and distributing them to the PEs        
!    - extraction of data and saving observation header information and
!      observation body data in data arrays for single and multi-level reports
!    - removal of redundant information
!    - production of vertical profiles from single-level aircraft reports.
!   This module contains the following module procedures:
!    - organize_obs_proc : master subroutine,
!    - obs_org_print_stat
!    - obs_proc_init
!    - obs_gps_init
!    - obs_read_gps
!    - obs_bias_corr_iwv
!    - obs_print_events
!    - obs_print_statist
!   It uses from:
!    - src_obs_proc_aof: - obs_aof_init
!                        - obs_read_aof_org: obs_read_distribute, obs_save_head,
!                                            obs_single_level, obs_multi_level,
!                                            obs_redundancy, obs_sort_levels,
!                                            obs_RASS_tv2t, obs_option_groups,
!                                            obs_extract, obs_flag6b
!    - src_obs_proc_air: - obs_air_org_mult: obs_air_make_mult, obs_air_correl,
!                                            obs_air_correl
!    - src_obs_cdfin_org:   - obs_cdf_mult_qualicheck
!                           - obs_cdf_raso_rh_bias
!    - src_obs_cdfin_print: - obs_cdf_print_reject
!                           - obs_cdf_print_odr
!                           - obs_print_number_rep
!    - src_obs_proc_util:   - obs_aof_assign_gridpt
!                           - get_global_surf_aof
!                           - obs_pointrs
!                           - obs_del_old_report
!                         and used by a subroutine from another module as above:
!                           obs_error, obs_fix_error, psurf, z2p, tv2t, p2dp
!  Note: This module is not used any more when reading obs from NetCDF files !
!        (It might be used again in the future if GPS data should also be read
!         (from formatted ASCII files).)
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.13       1998/10/22 Michael Buchhold
!  Initial release
! 1.15       1998/11/02 Christoph Schraff
!  Global declaration of allocatable arrays moved to module 'data_obs_process'.
! 1.19       1998/12/11 Christoph Schraff
!  Verification mode: complete observation processing also for data not used for
!  data assimilation. Revised report and data events. ANSI violations removed.
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables
! 1.26       1999/03/25 Guenther Doms
!  Decleration of three local variables (declerated twice) corrected.
! 1.27       1999/03/29 Christoph Schraff
!  Revised VOF format and ODR flag formats. 6 bit hollerith station identity 
!  replaced by characters. Revised (32-bit) missing data indicator in ODR.
!  Extreme temperature group introduced in AOF optional groups.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules, prepared use of MPE_IO and corrected
!  calls to intrinsics for non-Cray machines
! 1.31       1999/07/01 Christoph Schraff
!  MPI calls related to fixed-length arrays replaced by parallel_utility calls;
!  quantities related to MPI communicator 'icomm_world' removed.
! 1.34       1999/12/10 Ulrich Schaettler
!  Use new timing routines
! 1.36       2000/02/24 Christoph Schraff
!  Observation increment format of VOF included.
!  Introduction of ACAR aircraft reports
!  adaptation to 32 bit AOF word length
! 1.38       2000/04/06 Christoph Schraff
!  Index for pressure tendency, redundant report flag, aircraft roll angle code,
!  etc., added to ODR and VOF.
! 1.39       2000/05/03 Ulrich Schaettler
!  Some global variable names for grid specification changed.
! 1.40       2000/05/23 Christoph Schraff
!  Correction of temporal criterion for writing the observation statistics.
! 1.43       2000/07/18 Jean-Marie Bettems.
!  Use of 'irealgrib' (for obs_single_level).
! 2.4        2001/01/29 Christoph Schraff
!  'nmxbln' reduced.
! 2.5        2001/06/01 Christoph Schraff
!  Introduction of lapse rate and wind shear check for multi-level data.
!  Adaptation of aircraft report listing to cope with flight track checks.
!  Additional data events, and savety test at array allocations.
! 2.6        2001/06/12 Christoph Schraff
!  Variable for message on thinning sequences of aircraft reports.
! 2.12       2001/11/07 Christoph Schraff
!  Flight track check done also without production of multi-level aircraft rep.
! 2.13       2002/01/18 Michael Buchhold
!  Introduction of Wind Profiler / RASS reports. Function 'ichcon' removed.
!  Flight track check and thinning flags moved to station characteristics word.
! 2.19       2002/10/24 Michael Buchhold
!  Introduction of Radar VAD wind reports.
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Extension for reading from a dedicated file and processing of GPS-derived
!  IWV data.
! 3.6        2003/12/11 Christoph Schraff
!  Introduction of obs type SATOB. Ensuring error message at model_abort.
! 3.12       2004/09/15 Christoph Schraff
!  Extension of array allocations and diagnostic printing to prepare to include
!  the assimilation of satellite retrievals.
!  Outsourcing of various subroutines to new modules 'src_obs_proc_air',
!  'src_obs_proc_aof', and 'src_obs_proc_util'. Inclusion of former include-
!  files 'obs_read_gps' and 'obs_print***'. Subroutine 'organize_obs_proc'
!  cleaned from most explicit computations by moving initialisation steps
!  to new routines 'obs_proc_init', 'obs_aof_init', and 'obs_gps_init'.
! 3.15       2005/03/03 Christoph Schraff
!  (LME:) Fraction of maxmll/maxmlo and of maxsgl/maxsgo slightly increased.
!  Replaced FLOAT by REAL (Ulrich Schaettler)
! 3.18       2006/03/03 Christoph Schraff
!  GPS processing adapted for IWV spatial consistency check. Flushing of output
!  files introduced. Statistics output extended to 1DVar satellite retrievals.
!  Caution messages on insufficient ODR size to specific file unit 'nucautn'.
!  Use additional variable lyear_360 for changed interface to get_utc_date
! 3.21       2006/12/04 Burkhardt Rockel
!  Introduced variable polgam
! V3_23        2007/03/30 Ulrich Schaettler
!  Introduced ldump_ascii for flushing the ASCII files
! V4_4         2008/07/16 Ulrich Schaettler
!  Changed NL parameter lyear_360 to itype_calendar, to have several options
! V4_5         2008/09/10 Christoph Schraff
!  Adaptions to read obs from NetCDF files (changed report number/ index names,
!  'o??vip' arrays replaced by 'nbsvi?' or 'nhtvi?' elements of ODR).
! V4_7         2008/12/12 Christoph Schraff
!  Bug correction, for control output (nfmt24, nfmt25 introduced).
! V4_8         2009/02/16 Oliver Fuhrer
!  Use global values only if num_compute greater 1 
! V4_10        2009/09/11 Davide Cesari
!  Add characters after a backslash in comments; g95 interprets that as
!  continuation line
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Christoph Schraff, Klaus Stephan
!  New subroutine and call "obs_cdf_raso_rh_bias" 
!  - for bias correction of Vaisala RS92 radiosonde humidity.
!  The type of bias correction can be selected by new namelist mqcorr92
! V4_22        2012/01/31 Christoph Schraff
!  - Routine:            --> re-named or split into:    --> moved to new module:
!    obs_print_verif     -->    obs_print_vof           --> src_obs_print_vof
!    obs_print_out       --> 1. obs_cdf_print_reject    --> src_obs_cdfin_print
!                            2. obs_cdf_print_odr       --> src_obs_cdfin_print
!                            3. obs_print_number_rep    --> src_obs_cdfin_print
!    obs_print_events    -->    obs_cdf_print_events    --> src_obs_cdfin_print
!    obs_print_statist   -->    obs_cdf_print_statist   --> src_obs_cdfin_print
!    obs_mult_qualicheck -->    obs_cdf_mult_qualicheck --> src_obs_cdfin_org
!    obs_assign_gridpt   -->    obs_aof_assign_gridpt
!    get_global_surface  -->    get_global_surf_aof
!    obs_cdf_raso_rh_bias       (modified)              --> src_obs_cdfin_org
!  - Variables moved from 'data_nudge_all' and 'data_obs_process' to
!    'data_obs_cdfin' and 'data_obs_lib_cosmo'. Modified interfaces to the
!    routines called.
!  - Call of external routine 'difmin' replaced by new routine 'diff_minutes' 
!    from utilities (type of routine arguments: iintegers instead of intgribf).
!  - Call of routine for bias correction of Vaisala RS92 radiosonde humidity.
!  - Karolin Eichler: reprocessed GPS data can be read from ASCII files in
!    COST-716 V2.0 (or V1.0) format.
!  Note: This module is not used any more when reading obs from NetCDF files !
! V4_23        2012/05/10 Ulrich Schaettler
!  Adapted interface to SR diff_minutes (added itype_calendar)
! V4_24        2012/06/22 Hendrik Reich
!  Read also minutes and seconds from start date
! V4_25        2012/09/28 Davide Cesari
!  Patch in order to run with AOF files: without that, in src_obs_cdfin_print.f90
!  the model tries to reopen YUREJECT and YUOBSDR files with a runtime error.
! V4_28        2013/07/12 Christoph Schraff
!  Interface to 'obs_solar_zenith_angle' modified.
!  Direct call of mpi routine replaced by call of 'scatterv_values'.
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

USE data_obs_process, ONLY :  &

! Section 1 : Analysis Observation File formats                  
!------------------------------------------------------------------------------
 
!         1.1   observation file parameters
!               ---------------------------
   nmxbln,     & ! max. length of buffer containing all reports for all PE's
!         1.2   AOF data record header format
!               -----------------------------
   nbufhed,    & ! report header length in send buffer

!         1.5     Extracted parameters and Variables
!                 ----------------------------------
   nobtyp,     & ! observation type
   ncdtyp,     & ! code type
   ntime,      & ! actual time of observation
   taof,       & ! observation time in forecast hour units
   roblat,     & ! station latitude
   roblon,     & ! station lonitude

!         1.6     AOF time boxes
!                 --------------
   naofbox,    & ! length of AOF time periods to be read
   naoflbx,    & ! initial time of last  AOF time box to be read currently
   naoffbx,    & ! initial time of first AOF time box to be read currently

!         2       Numerical constants and scaling constants related to AOF
!                 --------------------------------------------------------
   atol  ,     & ! area tolerance (because of rounding off problem)

!         4       Diagnostic arrays and pointers for event counters (AOF input)
!                 ------------------------------
   noctpr,     & ! 2-d array of processed observations
   noctac,     & ! 2-d array of active observations
   noctps,     & ! 2-d array of passive observations
   noctrj,     & ! 2-d array of rejected observations
   nobtpp,     & ! diagnostic arrays pointer
   ncdtpp,     & !                 "
   nmxcdt,     & ! max. no. of code type for any obs. type
   ntotob,     & ! total number of obs. types

!         5.1    Character descriptions of report and data events and flags
!                ----------------------------------------------------------
   mxcrev,     & ! max number of report events
   mxcdev,     & ! max number of data events
   crepfl,     & ! Description of report flags
   cdatfl,     & ! Description of data flags
   clevua,     & ! Description of upper air levels
   clevsy,     & ! Description of synop level

!         5.2   Character decriptions of observation types and code types
!               ---------------------------------------------------------
   chobtp,     & ! Description of observation types
   chobcd,     & ! Description of observation code types

!         6     Temporary information buffers and report counters
!               -------------------------------------------------
   nodbuf,     & ! buffer containing all reports for a single node
   nmlob,      & ! current number of  multi-level report
   nsgob,      & ! current number of single-level report
   ngpob         ! current number of GPS report

! end of data_obs_process

!-------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

! Section 3 : Data event counter arrays and diagnostics arrays
! ------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------
!         3.1.1   Report event counter array format
!                 ---------------------------------
   mxreve,     & ! length of report event counter array
   nedbfl,     & ! data base flag on loc/tim/alt high
   netime,     & ! time out of range (too old) (or time missing)
   nenoal,     & ! no station altitude
   neloca,     & ! station location out of range
   nezdif,     & ! distance 'model orography - station altitude' too large
   neblak,     & ! blacklisted ship
   neobct,     & ! observation or code type excluded on area with sta. location
   nesodr,     & ! report number exceeding size of ODR (==> adjust Namelist)
   nenoda,     & ! no accepted data in report
   nenops,     & ! pressure too small (< 20hPa), or missing with aircraft report
   neredn,     & ! redundancy between 2 multi-level, or 2 single-level reports
   netrac,     & ! (flight) track suspicious
   nethin,     & ! thinning of aircraft (flight) track
   neredx,     & ! redundancy between 1 multi- and 1 single-level report
   neslml,     & ! one multi-level report made from  single-level reports
   neslps,     & ! single-level report put in multi-level report and set passive
   nenoml        ! multi-levl report not made due to ODR array size

USE data_obs_cdfin, ONLY :  &

!         3.1.2   Data event counter array format
!                 -------------------------------
   mxdeve,     & ! length of data event counter array
   nelodr,     & ! level rejected: number of levels exceeding ODR size
   nelmis,     & ! level rejected: pressure (PILOT: pressure and height) missing
   nelflg,     & ! level rejected: pressure (PILOT: height) flagged
   nelsfc,     & ! level rejected: too many surface levels
   nelnop,     & ! level rejected: PILOT height level not in range of model lev.
   nelext,     & ! level rejected: pressure < 9hPa or level below station height
   nelsig,     & ! level rejected: significant level above a specified limit
   nelrdn,     & ! level rejected: redundant level in report  (not active yet)
   nepmis,     & ! pressure (TEMP: height): missing
   nepflg,     & ! pressure (TEMP: height): flagged
   neprac,     & ! pressure: bad reporting practice
   nepalt,     & ! pressure: sta. height, or height dist. to orography too large
   nepsdt,     & ! pressure tendency: flagged, or absolute value too large
   netmis,     & ! temperature missing
   netflg,     & ! temperature flagged
   netext,     & ! temperature too low or too high
   netalt,     & ! height (diff.) too large for 2m-temp.
   netlps,     & ! lapse rate of multi-level temperature too large
   neqmis,     & ! humidity missing
   neqflg,     & ! humidity flagged
   neqlow,     & ! humidity too low
   neq300,     & ! humidity above 300 hpa
   neqbig,     & ! humidity over allowed value (120%)
   neqsap,     & ! humidity forced to be saturated (t>0)
   neqsam,     & ! humidity forced to be saturated (t<0)
   neqclp,     & ! humidity forced to be <=100% (t>0)
   neqclm,     & ! humidity forced to be <=100% (t<0)
   neqalt,     & ! height (diff.) too large for 2m-humid
   nedmis,     & ! wind direction missing
   nefmis,     & ! wind speed missing
   nedflg,     & ! wind direction flagged
   nefflg,     & ! wind speed flagged
   nefneg,     & ! wind speed too small  ( < 0 ; DRIBU: <= 0)
   nevalt,     & ! height (diff.) too large for 10m-wind
   nefshr,     & ! wind speed shear too large
   nedshr,     & ! directional wind shear too large
   nerlim,     & ! prec.amount exceeds threshold limit
   negmis

USE data_obs_cdfin, ONLY :  &

!         3.2     Event counter arrays
!                 --------------------
   neventr,    & ! counter array of report events
   neventd,    & ! counter array of data events

!         3.3     Character descriptions of events and flags
!                 ------------------------------------------
   crepev,     & ! description of report events
   cdatev,     & ! Description of data events

! Section 4 : Observation errors
! ------------------------------

!         4.1  Observation error levels
!              ------------------------
   nerlev,     & ! number of standard error levels
   rlevel,     & ! error levels
   rolnlv,     & ! ln(rlevel(15))

! Section 5 : Different kinds of limits
! -------------------------------------

!         5.1    Redundancy check limits
!                -----------------------
   rtmlim ,    & !  time limit for all reports except airep    [hrs]
   rhzlim ,    & !  horiz.  dist. limit for obs sta. (\ airep)  [km]
   rvtlim ,    & !  vertic. dist. limit for obs sta. (\ airep)   [m]

!         5.4    Limits for the directional wind shear check
!                -------------------------------------------
   nnqcdd,     & ! number of levels in the quality control threshold tables
   nqcdd,      & ! thresholds for directional shear, as funct. of speed
   nqcddff,    & ! limit values for sum of wind speeds (for 'nqcdd')
   qcfddff       ! upper limit of wind speed, for directional shear check

USE data_obs_cdfin, ONLY :  &

! Section 7 :  For reporting rejection of data: Output buffer, size and formats
! -----------------------------------------------------------------------------

   outbuf,     & ! buffer containing output for a single node
   nacout,     & !  actual number of records stored in the output buffer
   nmxoln,     & !  maximum length of output buffer
   istrej,     & !  length of strings (station id) in output buffer
   nfmt1 ,     & ! no pressure
   nfmt2 ,     & ! excess of precipitation 
   nfmt3 ,     & ! no accepted data
   nfmt4 ,     & ! excess of levels
   nfmt5 ,     & ! several surface levels
   nfmt6 ,     & ! excess of pressure tendency 
   nfmt7 ,     & ! excess of lapse rate
   nfmt8 ,     & ! excess of wind speed shear
   nfmt9 ,     & ! excess of directional shear
   nfmt10,     & ! redundancy of surface-level report
   nfmt11,     & ! redundancy of multi-level report
   nfmt12,     & ! redundancy of aircraft report
   nfmt13,     & ! redundancy of wind
   nfmt14,     & ! redundancy of temperature
   nfmt15,     & ! redundancy of humidity
   nfmt16,     & ! redundancy of pressure / height
   nfmt17,     & ! thinning of aircraft reports 
   nfmt18,     & ! exaggerated flight colocation
   nfmt19,     & ! flight track error
   nfmt20,     & ! message only: fog and precipitation
   nfmt21,     & ! message only: fog and invisible sky
   nfmt22,     & ! message only: fog and no cloud base
   nfmt23,     & ! message only: cloud and no cloud base or fog
   nfmt24,     & ! report (partly) blacklisted
   nfmt25,     & ! report not on whitelist

! Section 8 : Temporary global model fields
! -----------------------------------------

   hsurf_tot,  & ! total array of model surface height
   fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_cdfin

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

    maxmll     ,& ! size (report dimension) of the  multi-level (m-l)  ODR
    maxsgl     ,& ! size (report dimension) of the single-level (s-l)  ODR
    maxgpl     ,& ! size (report dimension) of the (ground-based) GPS  ODR
    maxtvl     ,& ! size (report dimension) of the satellite retrieval ODR
    r_degrad   ,& ! factor for transforming degree to rad
    r_g        ,& ! acceleration due to gravity

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
    lopen_odr  ,& ! .true. if yuobsdr is open
    lopen_rej  ,& ! .true. if yurejct is open
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
    ngpgfz     ,& !   GPS report processed by GFZ
    
! 6. Data type with rules for CMA obs and code types
! --------------------------------------------------

    t_cmatyp   ,& ! data type for information on CMA observation and code types
    cma        ,& ! array of meta data on CMA observation and code types

! 7. Functions  
! ------------  

    i_cma         ! function to determine the index of 'cma'
                  ! referring to a given CMA observation and code type

! end of data_obs_lib_cosmo

!-------------------------------------------------------------------------------

USE data_nudge_all, ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id  == 0) (for print on std. output)
    onl_cart_id  ,& ! 'my_cart_id'  of node with area containing (ionl,jonl)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    lvofoi       ,& ! .TRUE if observation increments are also written to VOF
    acthr           ! actual forecast hour

USE data_nudge_all, ONLY :   &

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

! 2.1  temporal variables
!      ------------------

    lverif ,    & ! on - off switch for verification
    lverpas,    & ! .t. : on - off switch for verif. also of passive reports
    mruntyp,    & ! -1  : type of current model run used for increments in VOF
    nudgend,    & ! end of nudging period in 'model integration hours'
    nversta,    & ! start of verification period in timesteps
    nverend,    & ! end of verification period in timesteps
    hversta,    & ! start of verification period in 'model integration hours'
    hverend,    & ! end of verification period in 'model integration hours'
!   tconbox,    & ! 6*dt: timestep [s] for computing analysis increments
                  !       (i.e. time box of constant analysis increments)
    gnudggp,    & ! 0*10^-4: nudging coefficient for GPS-derived IWV [1/s]
    tipolmx,    & ! max. time span  for linear interpolat. for upper-air data
    tipmxsu,    & ! max. time span  for linear interpolat.of surf.-level data
    wtukrsa,    & ! temporal radius of infl. towards the past for TEMP/PILOT
    wtukrse,    & ! temporal radius of infl. towards the future for TEMP/PILOT
    wtukara,    & ! temporal radius of infl. towards the past for aircraft data
    wtukare,    & ! temporal radius of infl. towards the future for aircr.data
    wtuksua,    & ! temporal radius of infl. towards the past for surface data
    wtuksue,    & ! temporal radius of influence towards the future

! 2.2  variables used for adjustment of vertical correlation scales
!      ------------------------------------------------------------

    lsvcorl,    & ! adjustment of vertical correlation scales
    msprpar,    & ! switch specifying the surfaces along which observation
                  ! increments of upper-air data are (primarily) spread
    vcorls ,    & ! square of the vertical correlation scale
    rhinfl ,    & ! constant part of the 'correlation scale of the autore-
                  ! gressive horizontal correlation function' (='COSAC') [km]
    rhvfac ,    & ! multipl. factor to the vertically varying part of 'COSAC'

! 2.3  observation processing
!      ----------------------

    thairh ,    & !  maximum horizontal distance for combining single level
                  !  aircraft reports to a multi-level report
    lgpsbias,   & !  bias correction to GPS IWV applied
    mqcorr92   ,& !  switch for bias correction for Vaisala RS92 sonde humidity
    maxmlo ,    & !  max. number of multi-level reports within the total domain
    maxsgo ,    & !  max. number of (surface-level and upper-air single-level
                  !  reports within the total domain
    maxgpo ,    & !  max. number of GPS reports within total domain
    maxmlv ,    & !  max. number of obs levels in multi-level reports in ODR
    nolbc         !  number of grid rows at lateral boundaries
                  !  where obs are neglected

USE data_nudge_all, ONLY :   &

    lsynop,     & !  if SYNOP data is used
    laircf,     & !  if AIREP data is used (aircraft)
    lsatob,     & !  if SATOB data is used
    ldribu,     & !  if DRIBU data is used (drifting buoy)
    ltemp ,     & !  if TEMP  data is used
    lpilot,     & !  if PILOT data is used
    lsatem,     & !  if SATEM data is used
    lgps  ,     & !  if GPS   data is used
    lprodr,     & !  .t. for diagnostic print of obs data records ODR
    lcd096,     & !  .t. if gps   code  96 data is used

    obnlat,     & !  northern boundary of observation area
    obslat,     & !  southern boundary of observation area
    obwlon,     & !  western boundary of observation area
    obelon,     & !  eastern boundary of observation area
    exnlat,     & !  northern boundary for exclusion area
    exslat,     & !  southern boundary for exclusion area
    exwlon,     & !  western boundary for exclusion area
    exelon,     & !  eastern boundary for exclusion area
    mcdmsg1,    & !  processing / use of MSG1   code  71 data
    mcdmsg2,    & !  processing / use of MSG2   code  72 data
    mcdno15,    & !  processing / use of NOAA15 code 206 data
    mcdno16,    & !  processing / use of NOAA16 code 207 data
    mcdno17,    & !  processing / use of NOAA17 code 208 data
    mcdno18,    & !  processing / use of NOAA18 code 209 data
    ionl  ,     & !  / grid point coordinates
    jonl          !  \ for standard output

USE data_nudge_all, ONLY :   &

! 5. I/O device numbers and file names for nudging
! ------------------------------------------------

    nuaof  ,     & ! analysis observation file
    nuaofex,     & ! print out of expanded analysis observation file
    nugps  ,     & ! GPS observation file unit
    yugps  ,     & ! input gps observation file

! 6. Switch for observation processing and varia
! ----------------------------------------------

    lexigps,     & ! .TRUE if GPS input file exists
    tobnext,     & ! next time when observational data have to be read again
    tmaxbox,     & ! maximum interval for computing analysis increments
    liwvssc        ! .t. : spatial consistency check of IWV performed

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &
 
!       1.1    ODR header format
!       ------------------------

!       1.1.1  Header formats of ODR reports:'omlhed','osghed','ogphed','otvhed'
!              -----------------------------------------------------------------
    mxrhed,     & ! header length of multi-level reports
    mxshed,     & ! header length of single-level reports
    mxghed,     & ! header length of GPS reports
    nhilon,     & ! longitude of observing station
    nhjlat,     & ! latitude  of observing station
    nhalt ,     & ! station altitude [m]
    nhtime,     & ! time of observat. in forecast hours
    nhsurf,     & ! height of model grid pt. to which obs. is assigned
    nhzio ,     & ! longitude of obs. station (or lowest datum) in grid pt. unit
    nhzjo ,     & ! latitude  of obs. station in grid pt. units
    nhsolz,     & ! solar zenith angle [deg]
    nhvcbu,     & ! correction factor to vertical correlation scale for wind
                  ! at base of report
    nhvcbt,     & ! as 'nhvcbu', but for temperature
    nhvcbq,     & ! as 'nhvcbu', but for humidity
    nhvctu,     & ! correction factor to vertical correlation scale for wind
                  ! at top of report
    nhvctt,     & ! as 'nhvctu', but for temperature
    nhvctq,     & ! as 'nhvctu', but for humidity

!       1.1.2  Header formats of ODR reports:'momlhd','mosghd','mopghd','motvhd'
!              -----------------------------------------------------------------
    mxrhdf,     & ! header length of multi-level reports
    mxshdf,     & ! header length of single-level reports
    mxghdf,     & ! header length of GPS reports
    nhio  ,     & ! (local) x-coord. of grid pt. assigned to obs
    nhjo  ,     & ! (local) y-coord. of grid pt. assigned to obs
    nhitot,     & ! global x-coord. of grid pt. assigned to obs
    nhjtot,     & ! global y-coord. of grid pt. assigned to obs
    nhobtp,     & ! observation type
    nhcode,     & ! code type
    nhschr,     & ! station characteristics                      (see 1.1.4)
    nhqofl,     & ! report flags (rds) on lat/long/date/time/alt (see 1.1.4)
    nhpass,     & ! flag for report being set to 'passive'       (see 1.1.4)
    nhqcfw,     & ! threshold quality control flags for pressure, status of
                  ! verification output
    nhnlev,     & ! number of obs. levels (for multi-level reports)
    nhvqcf        ! for satellite retrieval: threshold quality control flags

USE data_obs_record, ONLY :   &
    nhflag   ,& ! report flags (obs type, surf., alt., sta ID) (see 1.2.1)
    nhcorr   ,& ! update sequence number (station correction indicator)
    nhcat    ,& ! data     category (from BUFR Section 1) 
    nhcats   ,& ! data sub-category (from BUFR Section 1)
    nhkz     ,& ! DWD internal classification number (observation type)
    nhcent   ,& ! originating centre  +  (1000* sub-centre)
    nhstid   ,& ! station identity number 
    nhdate   ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn   ,& ! absolute exact observation time [hhmm]
    nhsyhr   ,& ! absolute nominal (synoptic) observation time [yymmddhh]
    mxghdf      ! header length of GPS reports

USE data_obs_record, ONLY :   &
!       1.1.3  Header formats of ODR reports:'yomlhd','yosghd','yopghd','yotvhd'
!              -----------------------------------------------------------------
    ilstid,     & ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in ODR (and VOF) header
!       ------------------------------------------------------------------
    nvpabp,     & ! bit pos. for report set passive since it is     nvschr
                  !              used in a multi-level pseudo report
    nvpsbp,     & ! bit pos. for report set passive since at least    "
                  !              1 of the following 4 flags apply
    nvobbp,     & ! bit pos. for flag: 'station location out of       "
                  !                     user-specified area'
    nvalbp,     & ! bit pos. for flag: 'distance model orography -    "
                  !                     station altitude too large'
    nvbkbp,     & ! bit pos. for flag: 'blacklisted station (ship)'   "
    nvexbp,     & ! bit pos. for flag: 'observation or code type      "
                  !                     excluded at station location'
    nvrdbp,     & ! bit pos. for flag: 'redundant report'             "
    nvsebp,     & ! bit pos. for report located at sea grid pt.       "
    nvscbp,     & ! bit pos. for station correction indicator         "
    nvssbp,     & ! bit pos. for station suspicion indicator          "
    nvsibp,     & ! bit pos. for important station indicator          "
    nvscoc,     & ! no. of bits occ. by each indicator or flag        "
    nvinbp,     & ! bit pos. for instrument specification indicator   "
    nvinoc,     & ! no. of bits occ. by instrument specif. ind.       "
    nvhtbp,     & ! bit pos. for flight track error flag              "
    nvhtoc,     & ! no. of bits occ. by flight track error flag       "
    nvhhbp,     & ! bit pos. for flight thinning flag                 "
    nvhhoc,     & ! no. of bits occ. by flight thinning flag          "
    nvapbp,     & ! bit pos. for phase of flight (aircraft)           "
    nvapoc,     & ! no. of bits occ. by phase of flight               "
    nvaabp,     & ! bit pos. for aircraft roll angle (code)           "
    nvaaoc,     & ! no. of bits occ. by aircraft roll angle           "
    nvhbps,     & ! bit pos. for flags on lat/lon/date/time/alt.    nhqofl
    nvhboc,     & ! no. of bits occ. by each flag                     "
    nvhibp,     & ! inner bit structure for each flag                 "
    nvhioc        ! no. bits occ. within flag strucure                "

USE data_obs_record, ONLY :   &

!       1.3    ODR body format
!              ---------------

!       1.3.0  Number of levels in multi-level ODR 'omlbdy', 'momlbd'
!              ------------------------------------------------------
    maxrsl,     & ! max. number of levels in multi-level ODR


!       1.3.1  Body format of ODR of multi-level reports: 'omlbdy'
!              ---------------------------------------------------
    mxrbdy,     & ! body length of multi-level reports
    nbtu  ,     & ! u wind component [m/s]
    nbtv  ,     & ! v wind component [m/s]
    nbtt  ,     & ! temperature [K]
    nbtrh ,     & ! relative humidity [/]
    nbtp  ,     & ! pressure [Pa]
    nbtz  ,     & ! height [m]
    nbtuer,     & ! error of observed wind component
    nbtter,     & ! error of observed temperature
    nbtqer,     & ! error of observed rel. humidity
    nbtzer,     & ! error of observed height
                  !  Note: nbt?er are set to the negative rms errors, if the
                  !  observations have not passed the threshold quality control
    nbtzio,     & ! longitude in grid pt. units
    nbtzjo,     & ! latitude  in grid pt. units
    nbtlop,     & ! LOG( pressure )


!       1.3.2  Body format of ODR of multi-level report flags: 'momlbd'
!              --------------------------------------------------------
    mxrbdf,     & ! body length of multi-level reports
    nbtflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbtqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
    nbterr,     & ! status flag word        (bit pattern, see below: 'nb?err')
!   nbtlsg,     & ! level id (bit pattern, as in NetCDF statistics file)
    nbtlid        ! level identity          (bit pattern, see below: 'nb?lid')

USE data_obs_record, ONLY :   &

!       1.3.3  Body format of ODR of surface reports: 'osgbdy'
!              -----------------------------------------------
    mxsbdy,     & ! body length of single-level reports
    nbsu  ,     & ! u wind component                                   [m/s]
    nbsv  ,     & ! v wind component                                   [m/s]
    nbst  ,     & ! temperature                                        [K]
    nbsrh ,     & ! relative humidity                                  [/]
    nbsp  ,     & ! pressure                                           [Pa]
    nbsz  ,     & ! height                                             [m]
    nbsuer,     & ! error of observed wind component
    nbster,     & ! error of observed temperature
    nbsqer,     & ! error of observed relative humidity
    nbszer,     & ! error of observed height
    nbspst,     & ! (3-hourly) pressure tendency                       [Pa/3h]
    nbspr ,     & ! precipitation amount                               [mm]
    nbscl ,     & ! low cloud cover                                    [octas]
    nbsper,     & ! error of observed precipitation amount
    nbscer,     & ! error of low cloud cover
    nbsvis,     & ! (horizontal) visibility                            [m]
    nbstn ,     & ! minimum temperature (at 2m during past 12 hrs)     [K]
    nbstx ,     & ! maximum temperature (at 2m during past 12 hrs)     [K]
    nbstg ,     & ! ground temperature (min. T-5cm during past 12 hrs) [K]
    nbsfgu,     & ! max. wind speed of gusts                           [m/s]
    nbsfme,     & ! max. wind speed of 10 minute mean wind             [m/s]
    nbsrad,     & ! global radiation, sum over 1 hour                  [J/m2]

!       1.3.4  Body format of ODR of surface report flags: 'mosgbd'
!              ----------------------------------------------------
    mxsbdf,     & ! body length of single-level reports
    nbsflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbsqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
    nbserr,     & ! status flag word        (bit pattern, see below: 'nb?err')
    nbslid,     & ! SYNOP: pressure code (SYNOP)   (code, see below: 'nbslid')
                  ! else : level identity   (bit pattern, see below: 'nb?lid')
    nbscwg,     & ! combined cloud and weather group (set of classes, s below)
    nbscfw,     & ! AOF read: flags for cloud, weather, precip, extreme temp.
                  !     (never set except for 'nvqfbp', 'nvtrbp')  (see below)
    nbswwe,     & ! NetCDF read, SYNOP: weather and ground group word  (below)
    nbstur        ! NetCDF read, Aircraft: degree of turbulence WMO Tab 011031
                  !   (not contained in merged multi-level aircraft reports !)

USE data_obs_record, ONLY :   &

!       1.3.5  Body format of ODR of GPS reports: 'ogpbdy'
!              -------------------------------------------
    mxgbdy,     & ! body length of GPS reports
    nbgtze,     & ! error in total zenith delay [mm]
    nbgzpd,     & ! zenith path delay (total zenith delay)             [mm]
    nbgzwd,     & ! zenith wet delay [mm]
    nbgiwv,     & ! integrated water vapour [mm]
    nbgp  ,     & ! pressure [Pa]
    nbgt  ,     & ! temperature [K]
    nbgrh ,     & ! relative humidity [/]
    nbgbia,     & ! bias correction to integrated water vapour [mm]
    nbgiwa,     & ! adjusted (bias corrected) integrated water vapour [mm]

!       1.3.6  Body format of ODR of GPS report flags: 'mogpbd'
!              -------------------------------------------------
    mxgbdf,     & ! body length of GPS reports
    nbgflg,     & ! main flag word          (bit pattern, see below: 'nb?flg')
    nbgerr,     & ! status flag word        (bit pattern, see below: 'nb?err')
    nbgqcf,     & ! threshold quality control flags      (see below: 'nb?qcf')
    nbglid,     & ! level identity          (bit pattern, see below: 'nb?lid')

!       1.3.7  Body format of ODR of sat retrieval reports: 'otvbdy'
!              -----------------------------------------------------
    mxtbdy,     & ! body length of multi-level reports
    nbvt  ,     & ! temperature [K]
    nbvrh ,     & ! relative humidity [/]
    nbvp  ,     & ! pressure [Pa]  (mandatory)

!       1.3.8  Body format of ODR of sat retrieval report flags: 'motvbd'
!              ----------------------------------------------------------
    mxtbdf,     & ! body length of sat retrieval reports
    nbvflg        ! main flag word          (bit pattern, see below: 'nb?flg')

USE data_obs_record, ONLY :   &

!       1.4.2  Other bit patt. for packed info in ODR (VOF) body, general words
!              ----------------------------------------------------------------
    nvru  ,     & ! bit pos. for status 'active' for horiz. wind    nb?err
    nvrt  ,     & ! bit pos. for status 'active' for temperature      "
    nvrq  ,     & ! bit pos. for status 'active' for humidity         "
    nvrz  ,     & ! bit pos. for status 'active' for pressure/height  "
    nvrw  ,     & ! bit pos. for status 'active' for vertical wind    "
    nvriwv,     & ! bit pos. for status 'active' for IWV            nb?err
    nvfubp,     & ! bit pos. for main flag on wind                  nb?flg
    nvftbp,     & ! bit pos. for main flag on temperature             "
    nvfqbp,     & ! bit pos. for main flag on humidity                "
    nvfzbp,     & ! bit pos. for main flag on pressure / geopot.      "
    nvfaoc,     & ! no. of bits occ. by each main flag                "
    nvfbps,     & ! bit pattern for main flags:                       "
    nvfboc,     & ! no. of bits occ. for main flags                   "
    nvflbp,     & ! bit pos. for level flag: level below surface      "
    nvfloc,     & ! no. of bits occ. by level flag                    "
    nvlidp,     & ! level id. bit pattern                           nb?lid
    nvlido        ! no. bits occ. by each indicator in level id.      "

USE data_obs_record, ONLY :   &

!       1.3    Further quantities related to ODR
!              ---------------------------------
    imdi  ,     & ! missing data indicator for ODR integers (2^31-1)
    ntotml,     & ! tot. number of stored multi-level reports
    ntotsg,     & ! tot. number of stored single-level reports
    ntotgp,     & ! tot. number of stored GPS reports
    ntottv,     & ! tot. number of stored satellite retrievals

!       2.     Observation data records (ODR)
!       -------------------------------------
    omlbdy,     & ! body of multi-level ODR
    momlbd,     & ! body of multi-level ODR
    omlhed,     & ! header of multi-level ODR
    momlhd,     & ! header of multi-level ODR
    osgbdy,     & ! body of single-level ODR
    mosgbd,     & ! body of single-level ODR
    osghed,     & ! header of single-level ODR
    mosghd,     & ! header of single-level ODR
    ogpbdy,     & ! body of GPS ODR
    mogpbd,     & ! body of GPS ODR
    ogphed,     & ! header of GPS ODR
    mogphd,     & ! header of GPS ODR
    otvbdy,     & ! body of satellite retrieval ODR
    motvbd,     & ! body of satellite retrieval ODR
    otvhed,     & ! header of satellite retrieval ODR
    motvhd,     & ! header of satellite retrieval ODR
    yomlhd,     & ! header of multi-level ODR
    yosghd,     & ! header of single-level ODR
    yogphd,     & ! header of GPS ODR
    yotvhd,     & ! header of satellite retrieval ODR

!       3.     Masking constants
!       ------------------------
    nibits        ! masking constants

! end of data_obs_record

!-------------------------------------------------------------------------------

  USE mo_fdbk_tables, ONLY :   &

    FL_OBSTYPE ,& ! passive report type (at obs. location)
    FL_AREA    ,& ! location not in valid area (i.e. outside user-defined area)
    FL_HEIGHT     ! location not in valid height range

! end of mo_fdbk_tables

!-------------------------------------------------------------------------------

USE data_modelconfig, ONLY :   &

! 1. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
  
! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    polgam,       & ! angle between the north poles of the systems
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    degrad,       & ! factor for transforming degree to rad

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

    dt,           & ! long time-step
    dtdeh           ! dt / 3600 seconds

! end of data_modelconfig

!-------------------------------------------------------------------------------

USE data_constants  , ONLY :   &

! 1. mathematical constants
! -------------------------

!   pi,           & ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

!   t0_melt,      & ! melting temperature of ice   
!   r_d,          & ! gas constant for dry air
    g,            & ! acceleration due to gravity
    rdocp           ! r_d / cp_d

! end of data_constants

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels

! 3. controlling the physics
! --------------------------
!   itype_gscp,   & ! type of grid-scale precipitaiton physics

! 7. additional control variables
! -------------------------------

    ldump_ascii,  & ! for flushing (close and re-open) the ASCII files

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------

    itype_calendar   ! for specifying the calendar used
                     !  = 0: gregorian calendar (default)
                     !    (but this needs a bug fix in get_utc_date,
                     !    because up to now we only have julian calendar)
                     !  = 1: every year has 360 days
                     !  = 2: every year has 365 days

! end of data_runcontrol

!-------------------------------------------------------------------------------

USE data_parallel,      ONLY :  &
    num_compute,     & ! number of compute PEs
    nboundlines,     & ! number of boundary lines of the domain for which
                       ! no forecast is computed = overlapping boundary
                       ! lines of the subdomains
    my_cart_id,      & ! rank of this subdomain in the cartesian communicator
    isubpos,         & ! positions of the subdomains in the total domain. Given
                       ! are the i- and the j-indices of the lower left and the
                       ! upper right grid point in the order
                       !                  i_ll, j_ll, i_ur, j_ur.
                       ! Only the interior of the domains are considered, not
                       ! the boundary lines.
    icomm_cart,      & ! communicator for the virtual cartesian topology
    imp_reals,       & ! determines the correct REAL type used in the model
                       ! for MPI
    imp_integers       ! determines the correct INTEGER type used in the model
                       ! for MPI

! end of data_parallel

!-------------------------------------------------------------------------------

USE data_io,            ONLY :  &

    ydate_ini          ! start of the forecast
                       ! yyyymmddhhmmss (year, month, day, hour, min., sec.)

! end of data_io

!-------------------------------------------------------------------------------

 USE environment,              ONLY :  &
    model_abort        ! aborts the program in case of errors

!-------------------------------------------------------------------------------

 USE parallel_utilities,       ONLY :  &
    global_values,   & ! computes global values by operating on local arrays
    gather_values,   & ! gathers a set of values from all nod. to 1 or all nodes
    scatterv_values, & ! scatters batches of variable length to other nodes
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    convert_month,   & ! number of month converted from a 3-character string
    get_utc_date,    & ! actual date of the forecast in different forms
    phi2phirot,      & ! converts phi from the real to the rotated system
    rla2rlarot,      & ! converts lambda from the real to the rotated system
    to_upper,        & ! converts alphabetic characters from lower to upper case
    diff_minutes       ! compute difference in minutes between 2 dates/times

!-------------------------------------------------------------------------------

 USE src_obs_proc_util,        ONLY :  &
    obs_aof_assign_gridpt ,& ! horizontal assignment of obs. report to a grid pt
    get_global_surf_aof   ,& ! gather orography/land fraction from global domain
    obs_pointrs       ,& ! finding the diagnostic array position
    obs_del_old_report,& ! deletion of obs. reports too old to be used any more
    p2dp                 ! get approx. model layer thickness at given pressure

!-------------------------------------------------------------------------------

 USE src_obs_proc_air,         ONLY :  &
    obs_air_org_mult  ,& ! production of multi-level aircraft reports from
                         !   single-level reports
    obs_air_interface    ! provides values for model-dependent input variables

!-------------------------------------------------------------------------------

 USE src_obs_proc_aof,         ONLY :  &
    obs_read_aof_org  ,& ! reading and processing of obs data from AOF
    obs_aof_init         ! initialisation of AOF obs processing

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    obs_solar_zenith_angle  ! get solar zenith angle for an array of reports

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_print,      ONLY :  & 
    obs_cdf_print_reject ,& ! prints messages on report / data rejection
    obs_cdf_print_odr    ,& ! prints a part of the ODR
    obs_print_number_rep    ! prints statistics of processed reports per node

!-------------------------------------------------------------------------------
 
 USE src_obs_cdfin_org,        ONLY :  & 
    obs_cdf_raso_rh_bias ,& ! bias correction of Vaisala RS92 r-sonde humidity
    obs_cdf_mult_qualicheck ! model-independent multi-level gross error checking

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_processing" for organizing the obs processing
!-------------------------------------------------------------------------------

SUBROUTINE organize_obs_proc

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_processing" is the master routine 
!   for the observation reading and processing.
!
! Method:
!   Reading data records from analysis observation file 'AOF' or other files
!   and distributing to the appropriate nodes.
!   Control of data extraction, processing and insertion in the 
!   observation data record 'ODR'.
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       01.11.97   Original code.    Michael Buchhold
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

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Local parameters:
! ----------------

  REAL    (KIND=ireals   ) , PARAMETER  ::  &
    thdlnpt = 28.6_ireals    ! mean vertical potential temperature gradient in
                             ! the troposphere in  log(p) units

! Local variables:
! ----------------
  LOGICAL, SAVE            ::  &
    lobpfrs = .TRUE.    ! variable for 'first time of obs. processing'

  INTEGER (KIND=iintegers) ::  &
    io      ,jo      ,& ! local  horizontal coordinates of observation
    irep    ,ivar    ,& ! loop indices
    nodrtot  (3)     ,& ! total number of multi-level / single-level / GPS rep.
    nodrold  (3)     ,& ! number of multi-level / single-level / GPS reports
                        ! before having read new data at the current timestep
                        !   (can be modified in 'obs_cdf_del_old_rep')
    imaxl    (4)     ,& ! size (report dimension) of the (4 types of) ODR
    nmloba           ,& ! number of multi-level reports  \  before having read
    nsgoba           ,& ! number of single-level reports  > new data at the
    ngpoba           ,& ! number of GPS reports          /  current timestep
    maxmla, nairll   ,& ! 
    nexceair     (2) ,& ! number of multi-level reports in excess of array size
    nstat   ,ierr       ! error status variable

  REAL (KIND=ireals)       ::  &
    zpob             ,& ! pressure at observation
    wtmaxa           ,& ! upper limit to the temporal radius of influence
                        ! towards the past of obs. time for any observation
    tnulast          ,& ! time [hrs] of last call of this subroutine
    tairobox         ,& ! time box for reading aircraft observations [h]
    zfsize           ,& ! 
    vscale       (4) ,& ! scale (in [ln(p)] units) of the static Gaussian
                        ! vertical correlation function (see obs_air_org_mult)
                        ! --> usual and reasonable values: 0.577,0.577, 0.2, 0.2
    tscale              ! scale [hrs] of the temporal Gaussian weight
                        ! --> usual and reasonable value : 0.25

  LOGICAL                  ::  &
    lhflag              ! 'nhflag' entry in ODR header has been set

  CHARACTER (LEN=20)       ::  &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=60)       ::  &
    yerrmsg             ! error message
  CHARACTER (LEN=40)       ::  &
    yerr                ! error message

! Local arrays:
! ------------

  REAL (KIND=ireals)       , ALLOCATABLE  ::       &
    odp (:)             ! approx. model layer thickness at the obs level

  INTEGER (KIND=iintegers) , ALLOCATABLE  ::       &
    i1evnt (:)       ,& ! first  index \ for data event
    i2evnt (:)          ! second index / counters
!
!------------ End of header ----------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Begin Subroutine organize_obs_proc
  !-----------------------------------------------------------------------------

  yroutine = 'organize_obs_proc'

  !-----------------------------------------------------------------------------
  !  Section  1.   Initialise observation processing, including opening of files
  !                and allocation of long term storage arrays and of short term
  !                storage rejection messaging arrays ('outbuf').
  !-----------------------------------------------------------------------------

  CALL obs_proc_init ( lobpfrs , nmxcdt , ntotob )
! ==================

  IF (lobpfrs)  THEN

    CALL obs_aof_init
!   =================

    CALL obs_gps_init
!   =================

! Set initial number of processed reports to zero
    ntotml = 0
    ntotsg = 0
    ntotgp = 0

    nmloba = 0
    nsgoba = 0
    ngpoba = 0

  ENDIF ! lobpfrs == .true.
 
! Set number of reports in the ODR
  nmlob  = ntotml
  nsgob  = ntotsg
  ngpob  = ntotgp

  !-----------------------------------------------------------------------------
  !  Section  2.   Get rid of the old observations in the data records,
  !-----------------------------------------------------------------------------

  IF (.NOT.lobpfrs) THEN

    CALL obs_del_old_report
!   =======================
    naoflbx = naoflbx + naofbox
    naoffbx = naoflbx

! Store number of reports before reading new reports
    nmloba = nmlob
    nsgoba = nsgob
    ngpoba = ngpob
  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  3.   Read GPS reports from dedicated file and fill ODR arrays
  !-----------------------------------------------------------------------------

  IF ((lobpfrs) .AND. (lexigps)) THEN

    CALL obs_read_gps ( .FALSE. )
!   =================

    IF (my_cart_id == 0)  CLOSE (nugps)
  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  4.   Read, extract, select, distribute to appropriate nodes,
  !                process, fill into ODR, and check for redundancy of
  !                observational reports from the AOF
  !-----------------------------------------------------------------------------

  CALL obs_read_aof_org ( nsgoba )
! =====================

  lobpfrs = .FALSE.

! Flush YUPRINT file
  IF (lwonl .AND. ldump_ascii) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1008, yerrmsg, yroutine)
  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  4b.  bias correction dependent on solar zenith angle
  !                for Vaisala RS92 radiosonde humidity
  !-----------------------------------------------------------------------------

  r_degrad = degrad
  r_g      = g

  IF (nmlob > nmloba) THEN

    CALL obs_solar_zenith_angle ( nmlob - nmloba , rmdich                      &
                                , momlhd(nmloba+1:nmlob,nhdate)                &
                                , momlhd(nmloba+1:nmlob,nhhrmn)                &
                                , omlhed(nmloba+1:nmlob,nhjlat)                &
                                , omlhed(nmloba+1:nmlob,nhilon)                &
                                , omlhed(nmloba+1:nmlob,nhsolz) )
!   ===========================

    CALL obs_cdf_raso_rh_bias ( mqcorr92, nmloba, nmlob )
!   =========================

  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  5.   Produce vertical profiles (multi-level reports)
  !                from single-level aircraft reports,
  !                and save total number of stored ODR reports
  !-----------------------------------------------------------------------------

  tairobox = REAL( naofbox , ireals )
  ntotml = nmlob
  ntotsg = nsgob
  ntotgp = ngpob
! (do not use 'nmlob', 'nsgob', 'ngpob' hereafter as total number of reports)
  nexceair = 0

! 'obs_air_interface' has to be called prior to a call to 'obs_air_org_mult'
! unless 'obs_cdf_interface' has already been called within this routine

  CALL obs_air_interface ( acthr , maxmll, maxsgl, maxmlv, lwonl, degrad       &
                         , nboundlines, num_compute, my_cart_id, icomm_cart    &
                         , imp_reals, imp_integers )
! ======================

  nodrold (1)  =  nmloba
  nodrold (2)  =  nsgoba
  nodrold (3)  =  ngpoba

  CALL obs_copy_err_status ( 2, nodrold, ntotml, ntotsg )
! ========================

  IF (lsvcorl) THEN
!   lhflag = (itype_obfile /= 1)
    lhflag = .FALSE.

    ALLOCATE (odp (MAX(ntotsg,1)) , STAT=nstat )
    DO irep = 1 , ntotsg
      io   = mosghd (irep,nhio)
      jo   = mosghd (irep,nhjo)
      zpob = osgbdy (irep,nbsp)
      odp (irep)  =  p2dp ( zpob , io , jo )
!                    ====
    ENDDO
    tscale   = MAX( MIN( wtukara , wtukare )                                   &
                  , (wtukara + wtukare)/ 10.0_ireals ) / c2

! get vertical correlation scale as by the namelist,
! convert it to ln(p) units

    DO ivar = 1 , 4 
      vscale (ivar) = SQRT( vcorls(ivar) )
      IF (msprpar == 2) vscale (ivar) = vscale(ivar) / thdlnpt
    ENDDO

!   PRINT *,'callair1 ', lsvcorl, nsgoba, nmloba, thairh, lhflag, ntotml, ntotsg
!   PRINT *,'callair2 ', MAXVAL( odp ), MAXVAL( -odp )
!   PRINT *,'callair3 ', vscale
!   PRINT *,'callair4 ', tscale

    zfsize = MIN( 4._ireals, (  MAX( wtukrsa , wtukara , tipolmx )             &
                              + MAX( wtukrse , wtukrse , tipolmx ) + tairobox) &
                            /(wtukara + wtukare + epsy) )
    IF (num_compute == 1) zfsize = 1
    maxmla = MAX( NINT( zfsize * maxmlo ) *2 , i1 )
    nairll = MAX( NINT( zfsize * maxsgo )    , i1 )
!   PRINT *,'callair7 ', maxmla, maxmll, nairll, maxsgl

    CALL obs_air_org_mult ( lsvcorl, nsgoba, nmloba, thairh, lhflag, lfirst    &
                          , isubpos , nexceair                                 &
                          , odp , vscale , tscale , rhinfl , rhvfac  )
!   =====================

!   PRINT *,'callair5 ', nsgoba, nmloba, ntotml, ntotsg

    DEALLOCATE (odp , STAT=nstat )

  ELSE
!   lhflag = (itype_obfile /= 1)
    lhflag = .FALSE.

!   PRINT *,'callair0 ', lsvcorl, nsgoba, nmloba, thairh, lhflag, ntotml, ntotsg

    CALL obs_air_org_mult ( lsvcorl, nsgoba, nmloba, thairh, lhflag, lfirst    &
                          , isubpos , nexceair )
!   =====================

!   PRINT *,'callair6 ', nsgoba, nmloba, ntotml, ntotsg

  ENDIF

! write alerting messages if array size is insufficient to accommodate all
! multi-level reports that can be derived from single-level (aircraft) reports

  IF (num_compute > 1) THEN
    CALL global_values ( nexceair, 2,'MAX',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11013, yerr, yroutine)
!                   ----------------
  ENDIF

  IF ((my_cart_id == 0) .AND. (MAXVAL(nexceair) > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=nstat)
    IF (nstat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1413, yerr, yroutine)
    IF (nexceair(1) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," LOC MULTI-LEV. AIRCR"  &
                     &,"AFT REPORTS BEYOND ARRAY SIZE ")' ) acthr, nexceair(1)
      nexceair(1)  =  nexceair(1) *maxmlo /maxmll + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY (AT LEAS" &
                     &,"T) ABOUT",I6)' )  nexceair(1)
      WRITE( nucautn,'("     OR THESE AIRCRAFT REPORTS ARE ASSIMILATED AS SIN" &
                     &,"GLE-LEVEL REPORTS")' )
    ENDIF 
    IF (nexceair(2) > 0) THEN
      WRITE( nucautn,'("CAUTION !!!!! t=",F6.3,":",I5," NEW MULTI-LEV. AIRCR"  &
                     &,"AFT REPS BEYOND ARRAY SIZE ")' ) acthr, nexceair(2)
      nexceair(2)  =  nexceair(2) *maxmlo /maxmll + 1
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceair(2)
    ENDIF 
    CLOSE (nucautn)
  ENDIF


  !-----------------------------------------------------------------------------
  !  Section  6.   Model-independent quality control check for multi-level data
  !-----------------------------------------------------------------------------

  IF (ntotml > nmloba) THEN
    ALLOCATE ( i1evnt (ntotml+1) , STAT=nstat )
    ALLOCATE ( i2evnt (ntotml+1) , STAT=nstat )
    DO irep = nmloba+1 , ntotml
      CALL obs_pointrs ( momlhd(irep,nhobtp) , momlhd(irep,nhcode) )
      i1evnt (irep) = nobtpp
      i2evnt (irep) = ncdtpp
    ENDDO

    CALL obs_cdf_mult_qualicheck ( nmloba , ntotml , maxmlv , i1evnt , i2evnt  &
                                 , rdocp )
!   ============================

    DEALLOCATE ( i1evnt , i2evnt , STAT=nstat )
  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  7.   Print out ODR and diagnostic output (+ de-allocate 'outbuf')
  !-----------------------------------------------------------------------------

  imaxl (1)    =  maxmll
  imaxl (2)    =  maxsgl
  imaxl (3)    =  maxgpl
  imaxl (4)    =  maxtvl
  nodrold (1)  =  nmloba
  nodrold (2)  =  nsgoba
  nodrold (3)  =  ngpoba
  nodrtot (1)  =  ntotml
  nodrtot (2)  =  ntotsg
  nodrtot (3)  =  ntotgp

  CALL obs_cdf_print_reject ( lwonl , num_compute , my_cart_id , icomm_cart    &
                            , imp_integers )
! =========================

  CALL obs_cdf_print_odr ( nodrtot , nodrold , num_compute , my_cart_id        &
                         , icomm_cart , imp_integers , 1 , chobtp )
! ======================

  CALL obs_print_number_rep ( nodrtot , nodrold , imaxl , lwonl , num_compute  &
                            , my_cart_id , icomm_cart , imp_integers )
! =========================

  !-----------------------------------------------------------------------------
  !  Section  8.   Determine next time, when obs. data have to be processed
  !-----------------------------------------------------------------------------

  wtmaxa   =  MAX( wtukrsa , tipolmx , wtukara , wtuksua , tipmxsu )           &
            + c05* tmaxbox /c3600 - epsy
  tobnext  =  naoflbx + naofbox - wtmaxa + epsy

! Flush or close files
! --------------------

! time [hrs] of last timestep, when proc. 'organize_nudging' is called
  tnulast  =  MIN( INT( MAX( nverend , nudgend ) ,iintegers) , nstop ) * dtdeh 

! if proc. 'obs_processing' is called the last time
  IF ((tobnext > tnulast) .AND. (my_cart_id == 0)) THEN
    CLOSE (nuaof) 
    CLOSE (nuaofex)
!   CLOSE (nuodr)

  ELSEIF (my_cart_id == 0) THEN
    CLOSE (nurej)
    CLOSE (nustat)
    CLOSE (nuodr)

    OPEN   (nurej ,FILE=yurejct,FORM='FORMATTED',STATUS='OLD'                  &
                               ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat == 0)                                                            &
      OPEN (nuodr ,FILE=yuobsdr,FORM='FORMATTED',STATUS='OLD'                  &
                               ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat == 0)                                                            &
      OPEN (nustat,FILE=yustats,FORM='FORMATTED',STATUS='OLD'                  &
                               ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) THEN
      yerrmsg = 'OPENING OF FILE yurejct / yuobsdr / yustats FAILED'
      CALL model_abort (my_cart_id, 1007, yerrmsg, yroutine)
    ENDIF
    lopen_rej = .TRUE.
    lopen_odr = .TRUE.
  ENDIF

! Flush YUPRINT file
  IF (lwonl .AND. ldump_ascii) THEN
    CLOSE (nupr)
    OPEN  (nupr,FILE=yuprint,FORM='FORMATTED',STATUS='OLD'                     &
                            ,POSITION='APPEND',IOSTAT=nstat)
    IF (nstat /= 0) yerrmsg = 'OPENING OF FILE yuprint FAILED'
    IF (nstat /= 0) CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
  ENDIF

! CALL obs_copy_err_status ( 1, nodrold, ntotml, ntotsg )
! ========================

!-------------------------------------------------------------------------------
! End Subroutine organize_obs_proc
!-------------------------------------------------------------------------------
END SUBROUTINE organize_obs_proc


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_processing": organizing printing the statistics
!-------------------------------------------------------------------------------

SUBROUTINE obs_org_print_stat

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_processing" organizes the printing
!   of statistics related to the observation processing.
!
! Method:
!   Print statistics on processed reports and summary of report and data events
!   by calling subroutines (at the last timestep nudging or verification is on).
!
! Current Code Owner:  Christoph Schraff
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

IMPLICIT NONE

!-------------------------------------------------------------------------------

! Local parameters: None
! ----------------

! Local variables:
! ----------------
  LOGICAL                  ::  &
    lrepevt                ! .true. if report events are to be printed

  INTEGER (KIND=iintegers) ::  &
    istat               ,& ! error status variables
    ilasev                 ! last obs. event to be printed
!
!------------ End of header ----------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Begin Subroutine obs_org_print_stat
  !-----------------------------------------------------------------------------

! Print statistics on the processed observational reports
  CALL obs_print_statist
! ======================

! Print summary of report and data events
  lrepevt = .TRUE.
  ilasev  = mxreve
  CALL obs_print_events ( lrepevt , neventr , mxreve , 1 , ilasev )
! =====================

  lrepevt = .FALSE.
  ilasev  = 18
  CALL obs_print_events ( lrepevt , neventd , mxdeve , 1 , ilasev )
! =====================
  ilasev  = mxdeve
  CALL obs_print_events ( lrepevt , neventd , mxdeve , 19, ilasev )
! =====================

  DEALLOCATE ( neventr   , STAT = istat )
  DEALLOCATE ( neventd   , STAT = istat )
 
! Close file
! ----------
  IF (my_cart_id == 0) THEN
    CLOSE (nustat)
    CLOSE (nurej)
    CLOSE (nuodr)
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_org_print_stat
!-------------------------------------------------------------------------------
END SUBROUTINE obs_org_print_stat


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_processing" for initialising the obs processing
!-------------------------------------------------------------------------------

SUBROUTINE obs_proc_init ( lobpfrs , ndimev2 , ndimev3 )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_processing" initialises the
!   observation processing.
!
! Method:
!   Opening of files for control output and statistics.
!   Allocation and initialization of long term storage arrays.
!
! Current Code Owner:  Christoph Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       23.08.04   Original code, part of previous 'organize_obs_proc'.
!                      Christoph Schraff
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

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
  LOGICAL                  , INTENT (IN)     :: &
    lobpfrs                ! variable for 'first time of obs. processing'

  INTEGER (KIND=iintegers) , INTENT (IN)     :: &
    ndimev2             ,& ! length of 2nd dimenstion of event counter arrays
    ndimev3                ! length of 3rd dimenstion of event counter arrays

! Local parameters: None
! ----------------

! Local variables:
! ----------------
  INTEGER (KIND=iintegers) ::  &
    istat, nstat           ! error status variables
!   nmloba              ,& ! number of multi-level reports before having
                           ! read new data at current timestep
!   nsgoba              ,& ! number of single-level reports before having
                           ! read new data at current timestep
!   ngpoba                 ! number of GPS reports before having
                           ! read new data at current timestep

  REAL (KIND=ireals)       ::  &
    zfsize              ,& ! factor used to set array sizes
    zfproc                 ! additional factor dependent on number of processors

  CHARACTER (LEN=20)       ::  &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerrmsg                ! error message
!
!------------ End of header ----------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Begin Subroutine obs_proc_init
  !-----------------------------------------------------------------------------

  yroutine = 'obs_proc_init'

  !-----------------------------------------------------------------------------
  !  Section  1.   Opening of control output files for obervation processing
  !-----------------------------------------------------------------------------

  IF (lobpfrs)  THEN

    IF (my_cart_id == 0) THEN

      OPEN (nurej,FILE=yurejct,FORM='FORMATTED',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yurejct FAILED'
        CALL model_abort (my_cart_id, 1003, yerrmsg, yroutine)
      ENDIF
      REWIND nurej
      lopen_rej = .TRUE.

      OPEN (nuodr,FILE=yuobsdr,FORM='FORMATTED',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yuobsdr FAILED'
        CALL model_abort (my_cart_id, 1004, yerrmsg, yroutine)
      ENDIF
      REWIND nuodr
      lopen_odr = .TRUE.

      OPEN (nustat,FILE=yustats,FORM='FORMATTED',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yustats FAILED'
        CALL model_abort (my_cart_id, 1005, yerrmsg, yroutine)
      ENDIF
      REWIND nustat

! dedicated GPS input observation file is opened further below

    ENDIF

  !-----------------------------------------------------------------------------
  !  Section  2.   Allocation and initialization of long term storage arrays
  !-----------------------------------------------------------------------------

!   report and data event counters
!   ------------------------------
    IF (ALLOCATED( neventr )) THEN
      PRINT '("CAUTION in src_obs_processing: neventr is already allocated "   &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( neventr , neventd , STAT=istat )
    ENDIF
    ALLOCATE ( neventr( mxreve, ndimev2, ndimev3), STAT=istat )
    ALLOCATE ( neventd( mxdeve, ndimev2, ndimev3), STAT=istat )

    neventr (:,:,:) = 0
    neventd (:,:,:) = 0
    
!   Observation data record (ODR)
!   -----------------------------

!   Determine the maximum number of reports in a subdomain as a function of the
!   nodes, the temporal weighting, and the total number of reports as specified
!   in the namelist input.

    zfproc = MIN( c1 , c1/10 + 9/(10 *SQRT(SQRT( c1*num_compute ))) )
!   zfproc = MIN( c1 , c1/5  + 4/(5  *     SQRT( c1*num_compute ) ) )

    zfsize = MIN( (  MAX( wtukrsa, wtukara, wtuksua, tipolmx, tipmxsu )        &
                   + MAX( wtuksue, tipmxsu ) + naofbox)                        &
                 /   MAX( wtuksua+wtuksue, c2*tipmxsu, c1 )  ,  4.0_ireals )
    maxsgl = MAX( i1 , NINT( maxsgo * zfproc * zfsize ) )

    zfsize = MAX( wtukrsa+wtukrse, wtukara+wtukare, c2*tipolmx, c1 )
    zfsize = (zfsize + naofbox) / zfsize
    maxmll = MAX( i1 , NINT( maxmlo * zfproc * zfsize ) )

    zfsize = MAX( c1 , MIN( nstop, nudgend ) *dtdeh )
    maxgpl = MAX( i1 , NINT( maxgpo * zfproc * zfsize ) )

!   zfsize = SQRT( ie_tot*dlon *je_tot*dlat /num_compute )
!   maxsgl = MIN( MAX( INT( zfsize *maxsgo/12 ) , 400 ) , maxsgo )
!   maxmll = MIN( MAX( INT( zfsize *maxmlo/10 ) , 4   ) , maxmlo )

    IF (lwonl) WRITE(nupr,'(" MAXIMUM NUMBER OF REPORTS IN THE TOTAL DOMAIN:", &
                           &"   MAXMLO=",I4,"  MAXSGO=",I5,/,                  &
                           &" MAXIMUM NUMBER OF REPORTS IN SUBDOMAINS:      ", &
                           &"   MAXMLL=",I4,"  MAXSGL=",I5)')                  &
                     maxmlo, maxsgo, maxmll, maxsgl

    IF (lwonl) WRITE(nupr,'(" MAXIMUM NUMBER OF GPS REP IN THE TOTAL DOMAIN:", &
                           &"   MAXGPO=",I4,/,                                 &
                           &" MAXIMUM NUMBER OF GPS REPORTS IN SUBDOMAINS:  ", &
                           &"   MAXGPL=",I4)')                                 &
                     maxgpo, maxgpl

    maxrsl = maxmlv

    IF (ALLOCATED( omlbdy )) THEN
      PRINT '("CAUTION in src_obs_processing: omlbdy is already allocated "    &
            &,"at time / PE",I6,I5)', ntstep, my_cart_id
      DEALLOCATE ( omlbdy , omlhed , momlbd , momlhd , yomlhd , STAT=istat )
      DEALLOCATE ( osgbdy , osghed , mosgbd , mosghd , yosghd , STAT=istat )
      DEALLOCATE ( ogpbdy , ogphed , mogpbd , mogphd , yogphd , STAT=istat )
    ENDIF
    ALLOCATE ( omlbdy( maxmll, maxrsl, mxrbdy), STAT=istat)
    ALLOCATE ( omlhed( maxmll,         mxrhed), STAT=istat)
    ALLOCATE ( momlbd( maxmll, maxrsl, mxrbdf), STAT=istat)
    ALLOCATE ( momlhd( maxmll,         mxrhdf), STAT=istat)
    ALLOCATE ( yomlhd( maxmll                ), STAT=istat)
    ALLOCATE ( osgbdy( maxsgl,         mxsbdy), STAT=istat)
    ALLOCATE ( osghed( maxsgl,         mxshed), STAT=istat)
    ALLOCATE ( mosgbd( maxsgl,         mxsbdf), STAT=istat)
    ALLOCATE ( mosghd( maxsgl,         mxshdf), STAT=istat)
    ALLOCATE ( yosghd( maxsgl                ), STAT=istat)
    ALLOCATE ( ogpbdy( maxgpl,         mxgbdy), STAT=istat)
    ALLOCATE ( ogphed( maxgpl,         mxghed), STAT=istat)
    ALLOCATE ( mogpbd( maxgpl,         mxgbdf), STAT=istat)
    ALLOCATE ( mogphd( maxgpl,         mxghdf), STAT=istat)
    ALLOCATE ( yogphd( maxgpl                ), STAT=istat)

    omlbdy (:,:,:) = rmdi
    omlhed (:,:)   = rmdi
    momlbd (:,:,:) = imdi
    momlhd (:,:)   = imdi
    yomlhd (:)     = '        '
    osgbdy (:,:)   = rmdi
    osghed (:,:)   = rmdi
    mosgbd (:,:)   = imdi
    mosghd (:,:)   = imdi
    yosghd (:)     = '        '
    ogpbdy (:,:)   = rmdi
    ogphed (:,:)   = rmdi
    mogpbd (:,:)   = imdi
    mogphd (:,:)   = imdi
    yogphd (:)     = '        '

! Calculate log of pressure of the observation error levels
    rolnlv(:) = LOG (rlevel(:))

! Set initial number of processed reports to zero
!   ntotml = 0
!   ntotsg = 0
!   ntotgp = 0

!   nmloba = 0
!   nsgoba = 0
!   ngpoba = 0

  ENDIF

  !-----------------------------------------------------------------------------
  !  Section  3.   initialise buffer for report rejection messages
  !-----------------------------------------------------------------------------

! Find maximum length of output buffer
  istrej  = ilstidp
  nmxoln  = maxmll*(100+istrej) + (maxsgl+maxgpl)*(5+(istrej+1)/2)

  IF (ALLOCATED( outbuf )) THEN
    PRINT '("CAUTION in src_obs_processing: outbuf is already allocated "      &
          &,"at time / PE",I6,I5)', ntstep, my_cart_id
    DEALLOCATE ( outbuf , STAT=istat )
  ENDIF
  ALLOCATE ( outbuf(nmxoln)             , STAT=istat )
  outbuf (:)   = 0

! Set initial number of output records to 0 
  nacout       = 0


!-------------------------------------------------------------------------------
! End Subroutine obs_proc_init
!-------------------------------------------------------------------------------
END SUBROUTINE obs_proc_init

!===============================================================================
!+ Module procedure in "src_obs_processing" for copying status 'active/passive'
!-------------------------------------------------------------------------------

SUBROUTINE obs_copy_err_status ( iactio, nodrold, ntotml, ntotsg )

!-------------------------------------------------------------------------------
!
! Description:
!   Procedure copies the status 'active / passive' from entries 'nb?err'
!   into 'nbt?er', 'nbs?er' or vice versa, because later routines currently
!   still ask for the latter entries.
!   When library 2 on observation operators is finished, this routine should
!   become obsolete.
!
! Method:
!   29.01.09
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 272in_org5
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
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

! Subroutine arguments:
! --------------------

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    iactio           ,& ! = 1 : copy status from 'nb?err' to 'nb??er'
                        ! = 2 : vice versa
    nodrold  (3)     ,& ! number of multi-level / single-level / GPS reports
                        ! before having read new data at the current timestep
                        !   (can be modified in 'obs_cdf_del_old_rep')
    ntotml , ntotsg     ! new number of multi-level / single-level reports

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    klev , irep  ,& ! loop indices
    nlev         ,& ! number of vertical levels
    ibits        ,& ! statement function to unpack any bit pattern
    invar        ,& ! word to be unpacked
    ipos , icovr ,& ! position of first bit / no of bits to be unpacked
    ireplace     ,& ! statement function to replace any bit pattern
    iboc         ,& ! no. of bit occ. by the bit structure to be replaced
    irepl        ,& ! word containing the replacing bit structure
    iposr           ! bit position of the replacing bit structure

! Local arrays:
! ------------
!
!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ipos,icovr) = IAND (ISHFT(invar,-ipos),icovr)

! Statement function to replace any bit pattern by another bit pattern
! --------------------------------------------------------------------
  ireplace ( invar, ipos, iboc, irepl, iposr ) =                               &
               IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )          &
                  , ISHFT( IAND( ISHFT( irepl,-iposr ), nibits(iboc) ), ipos ) )

!-------------------------------------------------------------------------------
! Begin Subroutine obs_copy_err_status
!-------------------------------------------------------------------------------

IF (iactio == 1) THEN

! loops have to start at 1, because reports can be written
! to irep <= nodrold in the redundancy checking !

! DO irep = nodrold(1)+1 , ntotml
  DO irep = 1 , ntotml
    nlev = momlhd(irep,nhnlev)
    DO klev = 1 , nlev
      IF (ibits( momlbd(irep,klev,nbterr),nvru,1 ) == 0)                       &
        omlbdy (irep,klev,nbtuer) = rmdi
      IF (ibits( momlbd(irep,klev,nbterr),nvrt,1 ) == 0)                       &
        omlbdy (irep,klev,nbtter) = rmdi
      IF (ibits( momlbd(irep,klev,nbterr),nvrq,1 ) == 0)                       &
        omlbdy (irep,klev,nbtqer) = rmdi
      IF (ibits( momlbd(irep,klev,nbterr),nvrz,1 ) == 0)                       &
        omlbdy (irep,klev,nbtzer) = rmdi
    ENDDO
  ENDDO

  DO irep = 1 , ntotsg
    IF (ibits( mosgbd(irep,nbserr),nvru,1 ) == 0)  osgbdy (irep,nbsuer) = rmdi
    IF (ibits( mosgbd(irep,nbserr),nvrt,1 ) == 0)  osgbdy (irep,nbster) = rmdi
    IF (ibits( mosgbd(irep,nbserr),nvrq,1 ) == 0)  osgbdy (irep,nbsqer) = rmdi
    IF (ibits( mosgbd(irep,nbserr),nvrz,1 ) == 0)  osgbdy (irep,nbszer) = rmdi
  ENDDO

ELSE

  DO irep = 1 , ntotml
    nlev = momlhd(irep,nhnlev)
    DO klev = 1 , nlev
      IF (omlbdy(irep,klev,nbtuer) < rmdich)                                   &
        momlbd(irep,klev,nbterr) = ireplace(momlbd(irep,klev,nbterr),nvru,1,0,0)
      IF (omlbdy(irep,klev,nbtter) < rmdich)                                   &
        momlbd(irep,klev,nbterr) = ireplace(momlbd(irep,klev,nbterr),nvrt,1,0,0)
      IF (omlbdy(irep,klev,nbtqer) < rmdich)                                   &
        momlbd(irep,klev,nbterr) = ireplace(momlbd(irep,klev,nbterr),nvrq,1,0,0)
      IF (omlbdy(irep,klev,nbtzer) < rmdich)                                   &
        momlbd(irep,klev,nbterr) = ireplace(momlbd(irep,klev,nbterr),nvrz,1,0,0)
    ENDDO
  ENDDO

  DO irep = 1 , ntotsg
    IF (osgbdy(irep,nbsuer) < rmdich)                                          &
      mosgbd (irep,nbserr) = ireplace( mosgbd(irep,nbserr), nvru, 1, 0, 0 )
    IF (osgbdy(irep,nbster) < rmdich)                                          &
      mosgbd (irep,nbserr) = ireplace( mosgbd(irep,nbserr), nvrt, 1, 0, 0 )
    IF (osgbdy(irep,nbsqer) < rmdich)                                          &
      mosgbd (irep,nbserr) = ireplace( mosgbd(irep,nbserr), nvrq, 1, 0, 0 )
    IF (osgbdy(irep,nbszer) < rmdich)                                          &
      mosgbd (irep,nbserr) = ireplace( mosgbd(irep,nbserr), nvrz, 1, 0, 0 )
  ENDDO

ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_copy_err_status
!-------------------------------------------------------------------------------

END SUBROUTINE obs_copy_err_status


!-------------------------------------------------------------------------------
!+ Module procedure in OBS_PROCESSING to read GPS reports
!-------------------------------------------------------------------------------

SUBROUTINE  obs_read_gps ( lcdf )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "OBS_PROCESSING" performs the
!   GPS report processing
!
! Method:
!   GPS data extraction from ascii files in COST-716 format.
!   Reference: "Format specification for COST-716 processed data."
!   COST-716 Document. Version 1.0-3, 28 September 2000.    
!   More data files can be joined in the input file. 
!   Values less than 0.0 are considered missing values.  
!   The grid point assignement is perfomed as for synop observations.   
!   A redundancy check is also included.
!
! Current Code Owner: DWD, Maria Tomassini
!  phone:  +49  69  8062 2711
!  fax:    +49  69  8062 3721
!  email:  maria.tomassini@dwd.de
!
! S-story:
! Version    Date       Name
! ---------- ---------- ----
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Initial release
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

! Subroutine arguments:
! --------------------

  LOGICAL                 , INTENT (IN)         ::  &
    lcdf                ! obs have been read from NetCDF files, which implies:
                        ! - 'cma' is used (instead of 'noctps') for statistics

! Local parameters:
! -----------------  

  INTEGER (KIND=iintegers) ::  &
    nlengps = 29        ! report header length in send buffer

  REAL    (KIND = ireals)    , PARAMETER  :: &
    c1000   =  1000.0_ireals     ,& !
    c1000r  =   0.001_ireals        !

! Local scalars:
! -------------

  LOGICAL                  ::  &
    lprcs               ,& ! obs to be processed (if lverif or liwvscc)
    lsfcob              ,& ! false, if no surface data are to be used
                           ! from a temp or pilot
    lseaobs             ,& ! true, if observation taken on sea
    lcoloc                 ! .TRUE. for spatial or temporal colocation

  INTEGER (KIND=iintegers) ::  &
    ianyy     ,imoyy    ,& ! year     / of the date and time
    ianmm     ,imomm    ,& ! month   /  of the current
    iandd     ,imodd    ,& ! day    <   observation
    ianhh     ,imohh    ,& ! hour    \  resp. of the initial
    ianmin    ,imomin   ,& ! minute   \ model date and time
    iansec    ,imosec   ,& ! second    \ model date and time
    imindif             ,& ! time difference [min] between these two times
    ierrf                  ! error status for time difference routine

  INTEGER (KIND=iintegers) ::  &
    icma     , iev2,     & ! indices for statistics ('cma', events)
    istat, irm, nstat,   & ! error status variables
    izmplcode, ierr,     & ! for MPI error code
    izlocstat,           & ! for space allocation error
    npe,                 & ! rank of the receiving PE
    nobbuf,              & ! number of reports in report buffers
    ibuflen,iallen,      & ! total length of report buffers
    inodcnt,             & ! length of local report buffer
    inodlen(num_compute),& ! total length of report buffer at each node
    inodpos(num_compute),& ! displacement of local report buffers in 'iallbuf'
    ipos,                & ! start address - 1 of an observation in 
                           ! report buffer
    irplen,              & ! report length
    ilen,                & ! length of report
    iobdat,              & ! observation date as yymmdd or yyyymodd
    nelevt,              & ! station altitude
    nstachr,             & ! station characteristics
    nkflag,              & ! processing flag
    mpassiv,             & ! passive flag assigned to report
    mrepflg,             & ! report flag word as part of station charcteristics
    nbufbit,             & ! buffer of bit pattern
    ncentre,             & ! data centre       (WMO Common Code Table C11 ???)
    nexceed,             & ! number of reports in excess of array size
    iob,     job,        & ! grid point asigned to observation (sub-domain)
    iob_tot, job_tot,    & ! grid points assigned to obs (total model area)
    nstcor,  nstcor2,    & ! station correction bit
    istcor,              & ! == 1, if active report is a station correction
    ngplo ,              & ! ODR index of GPS report
    nodract, nodrrej,    & ! index of active / rejected ODR report
    nobtypr, ncdtypr,    & ! observation/code type of rejected report
    irnode,              & ! process an observation belongs to
    inode,   inobs,      & ! loop indices
    icl,     i,          & ! loop indices
    insert,              & ! statement function to insert any bit pattern
    invar,               & ! word to be unpacked / to be partly replaced
    inval,               & ! bit pattern to be inserted
    ibit,                & ! bit position
    ibits,               & ! statement function to unpack any bit pattern
    ibp,                 & ! bit position of bit struct. to be unpacked/replaced
    icovr                  ! no. of bit occ. by the bit structure

  REAL (KIND=ireals)       ::  &
    timlim,              & ! limit for time difference in forecast hour units
    rscale,              & ! scaling factor
    rdegkm2,             & ! factor for conversion of degree to km (**2)
    cosolat,             & ! cos (latitude)
    rhzob,               & ! horizontal distance in km
    zsurf,               & ! height of model grid pt. to which obs. is assigned
    zio_tot, zjo_tot,    & ! observation position in g.p. units (total area)
    rlon,    rlat,       & ! latitude and longitude in the rotated system
    tzderr,              & !      / variables
    tzd,                 & !     /  for
    zwd,                 & !    /   temporary storage
    iwv,                 & !   <    of GPS measurements,
    surfpres,            & !    \   as
    surftemp,            & !     \  in
    surfrh,              & !      \ GPSbufr
    hhday,               & ! GPS observation hour
    miday,               & ! GPS observation minutes
    ttday,               & ! GPS observation time of the day
    x2                     ! variable  for computation of ttday
 
 
  CHARACTER (LEN=ilstidp)  ::  &
    ystid   , ystidr       ! station identity
  CHARACTER (LEN=12)       :: &
    yobdat                 ! date and time of current observation
  CHARACTER (LEN=20)       :: &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=30)       :: &
    yerrmsg                ! error message
  CHARACTER (LEN=40)       ::  &
    yerr                   ! error message

                     
! GPS variables: 

! Parameter
  CHARACTER (LEN=4),        PARAMETER :: &
      COSTfmtver = "V1.0" ! COST format version supported here
  CHARACTER (LEN=4),        PARAMETER :: &
      COSTfmtve2 = "V2.0" ! COST format version supported here
  CHARACTER (LEN=10),       PARAMETER :: &
      SOVFmark   = "COST-716 V"    ! Start of vfile marker
  INTEGER (KIND=iintegers), PARAMETER :: &
      EOVFstat    = -999            ! End of vfile status code
  CHARACTER (LEN=100),      PARAMETER :: &
      EOVFmark    = REPEAT("-",100) ! End   of vfile marker
  INTEGER (KIND=iintegers), PARAMETER :: &
      MaxSamples  = 288             ! Max. no. of data samples per vfile
  INTEGER (KIND=iintegers), PARAMETER :: & 
      MaxSlants   = 24              ! Max. no. of slant obs per data sample

! COST File Header info

  TYPE COSTheader
    CHARACTER (LEN=4)  :: FmtVer        ! File format ident code (Vm.n)
    CHARACTER (LEN=4)  :: StnID         ! Station ident code
    CHARACTER (LEN=20) :: RcvrType      ! Receiver type 
    CHARACTER (LEN=20) :: AntType       ! Antenna/Radome types
    REAL (KIND=ireals) :: StnLat        ! Station latitude  wrt WGS-84 (deg)
    REAL (KIND=ireals) :: StnLon        ! Station longitude wrt WGS=84 (deg)
    REAL (KIND=ireals) :: StnHt         ! Station height    wrt WGS-84 (m)
    REAL (KIND=ireals) :: StnHtAMSL     ! Station height    wrt EGM96  (m)
    CHARACTER (LEN=20) :: DateTimeFirst ! Date/time of 1st sample(dd-MMM-yyyy)
    CHARACTER (LEN=20) :: DateTimeProc  ! Date/time of process.  (dd-MMM-yyyy)
    CHARACTER (LEN=20) :: ProcCentre    ! Processing centre name
    CHARACTER (LEN=20) :: ProcMethod    ! Processing software type
    CHARACTER (LEN=20) :: OrbitType     ! GNSS orbit source and type
    CHARACTER (LEN=20) :: MetSource     ! Source of surface meteorological data
    INTEGER (KIND=iintegers) :: TimeInc    ! Sample time incr./resolut. (min.)
    INTEGER (KIND=iintegers) :: UpdateInt  ! Batch update interval (minutes)
    INTEGER (KIND=iintegers) :: BatchLen   ! Tot. length of batch data (minutes)
    INTEGER (KIND=iintegers) :: PCDH       ! Header Product Confidence Data
    INTEGER (KIND=iintegers) :: NumSamples ! Number of data samples in file
  END TYPE
  TYPE (COSTheader) Header 


! Slant delay & associated parameters 

  TYPE COSTslant
    CHARACTER (LEN=4)  :: SatID         ! GNSS satellite idents
    REAL (KIND=ireals) :: TSD           ! Total Slant delays (mm)
    REAL (KIND=ireals) :: TSDErr        ! Estimates errors in TSD (mm)
    REAL (KIND=ireals) :: Azim          ! Slant azimuth angles wrt North(deg)
    REAL (KIND=ireals) :: Elev          ! Slant elevations wrt horizon(deg)
  END TYPE
  TYPE (COSTslant) Slant

! File Data parameters

  TYPE COSTsample
    INTEGER (KIND=iintegers), DIMENSION(3) :: & 
                          TimeStamp     ! Data sample time stamp(hr,min,sec)
    INTEGER (KIND=iintegers)               :: &
                          PCDD          ! Data Product Confidence data
    REAL (KIND=ireals) :: TZD           ! Total Zenith Delay (mm)
    REAL (KIND=ireals) :: TZDErr        ! Estimated error in TZD (mm)
    REAL (KIND=ireals) :: ZWD           ! Zenith Wet Delay (mm)
    REAL (KIND=ireals) :: IWV           ! Integrated Water Vapour (mm)
    REAL (KIND=ireals) :: SurfPres      ! Surface pressure (hPa)
    REAL (KIND=ireals) :: SurfTemp      ! Surface temperature (K)
    REAL (KIND=ireals) :: SurfRH        ! Surface relative humidity (%)
    REAL (KIND=ireals) :: GradientN     ! N/S gradient (mm)
    REAL (KIND=ireals) :: GradientE     ! E/W gradient (mm)
    REAL (KIND=ireals) :: GrdErrN       ! Estimated error in N/S gradient(mm)
    REAL (KIND=ireals) :: GrdErrE       ! Estimated error in E/W gradient(mm)
    REAL (KIND=ireals) :: LogTEC        ! (Log10 of) total electron content
    INTEGER (KIND=iintegers)              :: &
                          NumSlants     ! Number of slants this data sample
    TYPE (COSTslant)   Slant(MaxSlants) ! Define up to MaxSlants slants
  END TYPE
  TYPE (COSTsample)       Sample 
  
! Input format definitions
  CHARACTER (LEN=5),  PARAMETER :: hfmt1="(A20)"
  CHARACTER (LEN=65), PARAMETER :: hfmt2=           &
            "(A4/A20,5X,A20/2F12.5,2F12.3/2(A20,5X)/4(A20,5X)/3I5/Z8.8/I4)"
  CHARACTER (LEN=3),  PARAMETER :: dfmt1="(A)"
  CHARACTER (LEN=27), PARAMETER :: dfmt2="(3I3,Z9.8,7F7.1,4F7.2,F7.3)"
  CHARACTER (LEN=4),  PARAMETER :: dfmt3="(I4)"
  CHARACTER (LEN=10), PARAMETER :: dfmt4="(A4,4F7.1)"

  CHARACTER (LEN=20)            :: hline      ! header line    
  CHARACTER (LEN=110)           :: dline      ! data line 
  INTEGER  (KIND=iintegers)     :: nsamp      ! sample loop counter
  INTEGER  (KIND=iintegers)     :: nslan      ! slant  loop counter

  INTEGER  (KIND=iintegers)     :: nconv      ! used for date conversion

  LOGICAL, SAVE                 :: leof  = .false.   ! end-of-file flag

! Local arrays:
! ------------
  
  INTEGER (KIND=iintegers), ALLOCATABLE :: &
    iallbuf(:)          ,& ! buffer containing observations from all PEs
    iposbuf(:)          ,& ! position of the reports in the report buffer
    inodbuf(:)             ! node the reports belong to
 

!
!------------- End of header ---------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ibp ,icovr) = IAND (ISHFT(invar,-ibp),nibits(icovr))

! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar,ISHFT(inval,ibit) )

!-------------------------------------------------------------------------------
! Begin Subroutine obs_read_gps
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Section 1.0: Initialization of organizational variables
  !-----------------------------------------------------------------------------

  yroutine = 'obs_read_gps'

  lprcs    = (lverif) .OR. (liwvssc)
  ngpob = ntotgp
! Find maximum length of the buffer containing all reports for all PE's
  nmxbln  = ( maxgpo * nlengps )
  nexceed = 0

  ALLOCATE( nodbuf (nmxbln)             , STAT=izlocstat )
  nodbuf (:)   = 0

  ibuflen      = 0
  npe          = 0
  IF (num_compute > 1) THEN
    ALLOCATE ( iallbuf (nmxbln)           , STAT=izlocstat )
    ALLOCATE ( inodbuf (nmxbln/nlengps+1) , STAT=izlocstat )
    ALLOCATE ( iposbuf (nmxbln/nlengps+1) , STAT=izlocstat )
! (nlengps is the minimum length of a report in the buffer 'iallbuf') 
    iallbuf (:)  = 0
    inodbuf (:)  = 0
    iposbuf (:)  = 0
    nobbuf       = 0
  ENDIF 
  
! Initial date and time of the forecast run
  READ (ydate_ini,'(I4,5I2)') imoyy, imomm, imodd, imohh, imomin, imosec

  !-------------------------------------------------------------------------
  ! Section 1.1  Compose full fields of model surface height / land fraction
  !-------------------------------------------------------------------------

  ALLOCATE ( hsurf_tot (ie_tot,je_tot)      , STAT=izlocstat )
  ALLOCATE ( fland_tot (ie_tot,je_tot)      , STAT=izlocstat )

  CALL get_global_surf_aof
! ========================

  !-----------------------------------------------------------------------------
  ! Section 2.    Read GPS data
  !-----------------------------------------------------------------------------

  IF (my_cart_id == 0) THEN
! =========================

    IF (lcdf) THEN
!CS: open files if reading other obs from NetCDF files
      OPEN (nurej  , FILE=yurejct, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yurejct FAILED'
        CALL model_abort (my_cart_id, 1010, yerrmsg, yroutine)
      ENDIF
      IF (nstat == 0)                                                          &
        OPEN (nustat , FILE=yustats, FORM='FORMATTED', STATUS='UNKNOWN'        &
                                   , POSITION='APPEND', IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yustats FAILED'
        CALL model_abort (my_cart_id, 1011, yerrmsg, yroutine)
      ENDIF
    ENDIF

    read_gps : DO WHILE (.NOT.leof)
!   -------------------------------
  
! Read (next) header

    hline   = ' ' 
    istat = 0
    DO WHILE ( hline(1:10) /= SOVFmark .AND. istat == 0 ) 
      READ ( UNIT=nugps, FMT=hfmt1, IOSTAT=istat ) hline
    ENDDO

! If OK, read rest of header info
    IF ( istat == 0 ) THEN
      Header%FmtVer = hline(10:13)
      IF ((Header%FmtVer /= COSTfmtver).AND.(Header%FmtVer /= COSTfmtve2)) THEN
         PRINT *, "COST format version ", Header%FmtVer, " not supported"
         PRINT *, "**** NO GPS DATA READ FROM INPUT FILE ****"
         EXIT read_gps
      ENDIF
    
      READ ( UNIT=nugps, FMT=hfmt2, IOSTAT=istat ) &
                        Header%StnID,         &
                        Header%RcvrType,      &
                        Header%AntType,       &
                        Header%StnLat,        &
                        Header%StnLon,        &
                        Header%StnHt,         &
                        Header%StnHtAMSL,     &
                        Header%DateTimeFirst, &
                        Header%DateTimeProc,  &
                        Header%ProcCentre,    &
                        Header%ProcMethod,    &
                        Header%OrbitType,     &
                        Header%MetSource,     &
                        Header%TimeInc,       &
                        Header%UpdateInt,     &
                        Header%BatchLen,      &
                        Header%PCDH,          &
                        Header%NumSamples
! Ensure strings are upper case

      CALL to_upper ( Header%StnID )
      CALL to_upper ( Header%RcvrType )
      CALL to_upper ( Header%AntType )
      CALL to_upper ( Header%DateTimeFirst )
      CALL to_upper ( Header%DateTimeProc )
      CALL to_upper ( Header%ProcCentre )
      CALL to_upper ( Header%ProcMethod )
      CALL to_upper ( Header%OrbitType )
      CALL to_upper ( Header%MetSource )
          
! If header ok, read samples
      nsamp = 0
      IF ( Header%NumSamples < 0 ) Header%NumSamples = MaxSamples

      DO WHILE ( nsamp < Header%NumSamples .AND. istat==0 )
      READ ( UNIT=nugps, FMT=dfmt1, IOSTAT=istat ) dline

      IF ( istat == 0 ) THEN
        IF ( dline(1:100) /= EOVFmark ) THEN
          READ (dline, FMT=dfmt2, IOSTAT=istat ) &
                     Sample%TimeStamp, &
                     Sample%PCDD,      &
                     Sample%TZD,       &
                     Sample%TZDErr,    & 
                     Sample%ZWD,       & 
                     Sample%IWV,       & 
                     Sample%SurfPres,  &
                     Sample%SurfTemp,  &
                     Sample%SurfRH,    &
                     Sample%GradientN, &
                     Sample%GradientE, &
                     Sample%GrdErrN,   &
                     Sample%GrdErrE,   &
                     Sample%LogTEC
          nsamp = nsamp + 1        ! (Another) sample read successfully
  
! Next read no. of slant samples and input that number of slants
          READ ( UNIT=nugps, FMT=dfmt3, IOSTAT=istat ) &
                      Sample%NumSlants

          DO nslan = 1, MIN(Sample%NumSlants,MaxSlants)
          READ ( UNIT=nugps, FMT=dfmt4, IOSTAT=istat ) &
                     Sample%Slant(nslan)%SatID,  &
                     Sample%Slant(nslan)%TSD,    &
                     Sample%Slant(nslan)%TSDErr, &
                     Sample%Slant(nslan)%Azim,   &
                     Sample%Slant(nslan)%Elev
          END DO ! on slants

! %%%% no correct allignement below
  !-----------------------------------------------------------------------------
  ! Section 2.1   Pre-selection of GPS reports
  !-----------------------------------------------------------------------------

 
  ! Skip empty report, i.e. if negative or zero (?) TZD 
  !-----------------------------------------------------------------------
    IF (Sample%TZD <= c0)  CYCLE read_gps 

    nkflag  = 0

  ! Station id 
  !-----------------------------------------------------------------------
    WRITE( ystid(1:4),'(A4)' ) Header%StnID
    DO icl = 5, ilstidp
      ystid (icl:icl) = ' '
    ENDDO

  ! Observation and code type
  !-----------------------------------------------------------------------
    nobtyp  = ngps
    ncdtyp  = ngpgfz

!   get diagnostic array position (--> nobtpp, ncdtpp)
!   --------------------------------------------------
    IF (lcdf) THEN
      icma  = i_cma ( nobtyp , ncdtyp )
      cma(icma)%cnt_pr = cma(icma)%cnt_pr + 1
      ncdtpp = icma
      nobtpp = 1
    ELSE
      CALL obs_pointrs ( nobtyp , ncdtyp )
      noctpr (nobtpp,ncdtpp) = noctpr (nobtpp,ncdtpp) + 1
    ENDIF
!   pre-set passive flag to zero
    mpassiv = 0
    mrepflg = 0

  ! Black lists
  !-----------------------------------------------------------------------
! to be done
!   IF (blacklisted) THEN
!     IF (mpassiv == 0)                                                        &
!       neventr (neblak,ncdtpp,nobtpp) = neventr(neblak,ncdtpp,nobtpp) + 1
!     WRITE( nurej,'(" STATION ",A ," : BLACKLISTED GPS")') ystid
!     IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
!     IF (.NOT. lprcs ) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
!     IF (.NOT. lprcs ) CYCLE read_gps
!     noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
!     mpassiv = 2
!     mrepflg = insert( mrepflg , 1 , nvbkbp )
!   ENDIF


  ! Observation time
  !-----------------------------------------------------------------------

!   observation time
!   ----------------
! Extract date from header  
! Year 
    READ(Header%DateTimeFirst(8:11),'(I4)') nconv
    iobdat = nconv * 10000
! Month (trasform month name in month number) 
    CALL convert_month ( Header%DateTimeFirst(4:6), nconv , 1) 
    iobdat = iobdat + nconv*100
! Day 
    READ(Header%DateTimeFirst(1:2),'(I2)') nconv 
    iobdat = iobdat + nconv 
! Extract time from sample 
    ntime  = Sample%TimeStamp(1)*100 + Sample%TimeStamp(2)          

    WRITE( yobdat,'(I8, I4)' )  iobdat, ntime
    READ ( yobdat,'(I4,4I2)' )  ianyy, ianmm, iandd, ianhh, ianmin
       
!   Time difference in minutes between observation time and initial model time 
    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                    &
                        ianyy, ianmm, iandd, ianhh, ianmin,                    &
                        itype_calendar, imindif, ierrf )
!   =================

!   observation time in forecast hour units  
    taof   = imindif/60.0_ireals
    
!   reject report if too old  
!   ------------------------ 
    timlim  =  taof + MAX( wtuksue , tipmxsu )   
    IF(ntstep*dt/3600._ireals > timlim ) THEN
      IF (.NOT. lcdf) neventr (netime,ncdtpp,nobtpp) = neventr(netime,ncdtpp,nobtpp) + 1
      IF (      lcdf) neventr (netime,icma,1) = neventr(netime,icma,1) + 1
!CS
      IF (.NOT. lcdf) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (      lcdf) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
!CS
      WRITE(nurej,'(" GPS STA ",A ," : OBSERVATION TOO OLD: ",F6.1,            &
                &" [FORECAST HRS]")' ) ystid ,taof
      CYCLE read_gps
    ENDIF
 
      
  ! Observation position
  !-----------------------------------------------------------------------

!   observation height (above Mean Sea Level) 
!   -----------------------------------------
    nelevt     =   Header%StnHtAMSL
    IF ((nelevt < -400) .OR. (nelevt > 9000)) THEN
      IF (.NOT. lcdf) neventr (nenoal,ncdtpp,nobtpp) = neventr(nenoal,ncdtpp,nobtpp) + 1
      IF (      lcdf) neventr (nenoal,icma,1) = neventr(nenoal,icma,1) +1
!CS
      IF (.NOT. lcdf) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (      lcdf) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
!CS
      WRITE( nurej,'(" STATION ",A ," : GPS STATION HEIGHT MISSING"            &
                   &,F6.1," HRS")' ) ystid, taof
      CYCLE read_gps
    ENDIF

!   observation location
!   --------------------
    roblat     =   Header%StnLat 
    IF ( Header%StnLon > 180._ireals ) THEN
        roblon =  Header%StnLon - 360._ireals
    ELSE 
        roblon =  Header%StnLon
    ENDIF

    rlon = rla2rlarot (roblat, roblon, pollat, pollon, polgam)
    rlat = phi2phirot (roblat, roblon, pollat, pollon)
!   observation position in grid point units (total area)
    zio_tot    = 1._ireals + (rlon - startlon_tot) /dlon
    zjo_tot    = 1._ireals + (rlat - startlat_tot) /dlat


!   Treat GPS as synop observation to  
!   assign observation to grid point
!   --------------------------------
    nobtyp  = nsynop
    lseaobs =    .FALSE.
    CALL obs_aof_assign_gridpt ( zio_tot, zjo_tot, nelevt, nobtyp , lseaobs    &
                               , iob_tot, job_tot, zsurf , lsfcob )
!   ==========================
!   reassign correct observation type 
    nobtyp  = ngps

!   reject report if observation location outside model domain
!   ----------------------------------------------------------
    IF   (iob_tot == 0) THEN
      
      IF (.NOT. lcdf .AND. mpassiv == 0)                                       &
         neventr (neloca,ncdtpp,nobtpp) = neventr (neloca,ncdtpp,nobtpp) + 1
      IF (      lcdf .AND. mpassiv == 0)                                       &
         neventr (neloca,icma,1) = neventr (neloca, icma,1) + 1
      WRITE(nurej,'(" GPS STA ",A ," : OBS. LOCATION OUT OF DOMAIN ",2F7.1     &
                  &," , ",F6.1," HRS")' ) ystid, roblon, roblat, taof
!CS
      IF (lcdf) THEN
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      ENDIF
!CS, etc.
      CYCLE read_gps
    ENDIF

!   PRINT '("ZG ",A,I10,2x,L,2I5,2x,L,F5.0,2F7.2)', ystid, nelevt, lseaobs     &
!                                           , iob_tot,job_tot, lsfcob, zsurf   &
!                                           , roblat, roblon

!   observation location outside user-defined area
!   ----------------------------------------------
    IF (     (roblat > obnlat) .OR. (roblat < obslat)                          &
        .OR. ((obwlon >= obelon) .AND. (roblon < obwlon)                       &
                                 .AND. (roblon > obelon))                      &
        .OR. ((obwlon <  obelon) .AND. (     (roblon < obwlon)                 &
                                        .OR. (roblon > obelon)))) THEN
      IF (.NOT. lcdf .AND. mpassiv == 0)                                      &
          neventr (neloca,ncdtpp,nobtpp) = neventr (neloca,ncdtpp,nobtpp) + 1
      IF (      lcdf .AND. mpassiv == 0)                                      &
          neventr (neloca,icma,1) = neventr (neloca,icma,1) + 1
      WRITE( nurej,'(" GPS STA ",A ," : OBS. LOCATION OUT OF USER-SPECIFIED"   &
                   &," AREA ",2F7.1," , ",F6.1," HRS")' )                      &
             ystid, roblon, roblat, taof
      IF (lcdf) THEN
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        IF (.NOT. lprcs ) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      ENDIF
      IF (.NOT. lprcs ) CYCLE read_gps
 
      IF (      lcdf) cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      IF (.NOT. lcdf) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvobbp )
      nkflag  = insert( nkflag  , 1 , FL_AREA )
    ENDIF

!   distance "model orography - station altitude" too large
!   -------------------------------------------------------
!   as for SYNOP 
    IF (iob_tot <  0) THEN
      
      IF (.NOT. lcdf .AND. mpassiv == 0)                                      &
         neventr (nezdif,ncdtpp,nobtpp) = neventr (nezdif,ncdtpp,nobtpp) + 1
      IF (      lcdf .AND. mpassiv == 0)                                      &
         neventr (nezdif,icma,1) = neventr (nezdif,icma,1) + 1
      WRITE( nurej,'(" GPS STA ",A ," : HEIGHT ",I5," DIFF. TO MODEL ",F5.0    &
                   &," TOO LARGE, ",F7.1," HRS")') ystid, nelevt, zsurf, taof
      IF (lcdf) THEN
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        IF (.NOT. lprcs ) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      ENDIF
      IF (.NOT. lprcs ) CYCLE read_gps
      
      IF (      lcdf) cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
      IF (.NOT. lcdf) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvalbp )
      nkflag  = insert( nkflag  , 1 , FL_HEIGHT )
      iob_tot = ABS( iob_tot )
      job_tot = ABS( job_tot )
    ENDIF

  ! Observation type, code type
  !-----------------------------------------------------------------------

    IF (nobtyp == ngps   .AND. .NOT.lgps      .OR.                             &
        ncdtyp == ngpgfz .AND. .NOT.lcd096        )              THEN
      IF (      (roblat <= exnlat+atol) .AND. (roblat >= exslat-atol)          &
          .AND. (roblon >= exwlon-atol) .AND. (roblon <= exelon+atol))  THEN
   
        IF (.NOT. lcdf .AND. mpassiv == 0)                                    &
            neventr (neobct,ncdtpp,nobtpp) = neventr(neobct,ncdtpp,nobtpp) + 1
        IF (      lcdf .AND. mpassiv == 0)                                    &
            neventr (neobct,icma,1) = neventr(neobct,icma,1)
        IF (lcdf) THEN
          IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
          IF (.NOT. lprcs ) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
        ELSE
          IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
          IF (.NOT. lprcs ) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        ENDIF
        IF (.NOT. lprcs ) CYCLE read_gps
        IF (      lcdf ) cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
        IF (.NOT. lcdf ) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
        mpassiv = 2
        mrepflg = insert( mrepflg , 1 , nvexbp )
        nkflag  = insert( nkflag  , 1 , FL_OBSTYPE )
      ENDIF
    ENDIF

!csc
    IF (TRIM(Header%ProcCentre) == 'GFZ') THEN
      ncentre = 23023
    ELSE
      ncentre = imdi
    ENDIF
    nstachr = 0

  !-----------------------------------------------------------------------
  ! Section  2.2      Find the appropriate node for this report
  !-----------------------------------------------------------------------

    irnode = imdi
    Nodes: DO  inode = 0, num_compute-1
      IF ( isubpos(inode,1) <= iob_tot .AND.                                   &
           isubpos(inode,3) >= iob_tot .AND.                                   &
           isubpos(inode,2) <= job_tot .AND.                                   &
           isubpos(inode,4) >= job_tot )                   THEN
        irnode = inode
        EXIT Nodes
      ENDIF
    ENDDO Nodes
    IF (irnode == imdi) THEN
 
      IF (.NOT. lcdf .AND. mpassiv == 0)                                     &
         neventr (neloca,ncdtpp,nobtpp) = neventr(neloca,ncdtpp,nobtpp) + 1
      IF (      lcdf .AND. mpassiv == 0)                                     &
         neventr (neloca,icma,1) = neventr(neloca,icma,1) + 1
      IF (lcdf) THEN
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      ENDIF
      WRITE( nurej,'(" GPS STA ",a ," : NO APPROPRIATE NODE FOUND ",           &
                   &2F7.1," , ",F6.1," HRS")' ) ystid, roblon, roblat, taof
      CYCLE read_gps
    ENDIF

!   count accepted reports
   IF (lcdf) THEN
     IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac + 1
   ELSE
     IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) + 1
   ENDIF
  !----------------------------------------------------------------------
  ! Section  2.3  Store this report in buffer 'nodbuf' ( --> 'iallbuf')
  !               together with the already extracted header information
  !----------------------------------------------------------------------

! passive reports are rejected, if lverpas = .FALSE.

   IF (lcdf) THEN
     IF ((.NOT. lverpas) .AND. (mpassiv > 0)) THEN
       cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
       cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1  
       CYCLE read_gps
     ENDIF
   ELSE
     IF ((.NOT. lverpas) .AND. (mpassiv > 0)) THEN
       noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
       noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
       CYCLE read_gps
     ENDIF
   ENDIF

!  GPS observations to be saved 
    tzderr    = Sample%TZDErr
    tzd       = Sample%TZD
    zwd       = Sample%ZWD
    iwv       = Sample%IWV
    surfpres  = Sample%SurfPres
    surftemp  = Sample%SurfTemp
    surfrh    = Sample%SurfRH
!CS
!   PRINT '(" OBS_READ_GPS: ",A,2F6.1,4F8.1,3I5,F5.0)',  &
!   ystid,taof,iwv,tzd,surfpres,surftemp,surfrh,iob_tot,job_tot,nelevt,zsurf

!   actual buffer length of array 'nodbuf'
    ipos   = ibuflen
!   report length of current observation: 
!          header information + values in GPS report buffer gpsbuf
    irplen = nlengps 
    IF (ipos+irplen > nmxbln)       THEN
      nexceed  =  nexceed + 1
      WRITE( nustat ,'(A20,": CAUTION: maxgpo TOO SMALL ==> REPORT REJEC"      &
                     &,"TED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:"   &
                     &,I7," / ",I7)' )                  yroutine, ipos, nmxbln
      PRINT          '(A20,": CAUTION: maxgpo TOO SMALL ==> REPORT REJEC"      &
                     &,"TED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:"   &
                     &,I7," / ",I7)' ,                  yroutine, ipos, nmxbln
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!     IF (mpassiv == 0)                                                        &
      IF (.NOT. lcdf) neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
      IF (      lcdf) neventr (nesodr,icma,1) = neventr(nesodr,icma,1) + 1
      IF (lcdf) THEN
        IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
      ELSE
        IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1 
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      ENDIF
      IF (      lcdf) cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      IF (.NOT. lcdf) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      CYCLE read_gps
    ENDIF
    IF (mpassiv == 2) mrepflg = insert( mrepflg , 1 , nvpsbp )
    nodbuf (ipos+ 1) = irplen
    nodbuf (ipos+ 2) = nobtyp
    nodbuf (ipos+ 3) = ncdtyp
    nodbuf (ipos+ 4) = nstachr
    nodbuf (ipos+ 5) = nkflag
    nodbuf (ipos+ 6) = ncentre
    nodbuf (ipos+ 7) = iobdat
    nodbuf (ipos+ 8) = ntime
    nodbuf (ipos+ 9) = imindif
    nodbuf (ipos+10) = nelevt
    nodbuf (ipos+11) = NINT( zsurf   *c1000 )
    nodbuf (ipos+12) = NINT( zio_tot *c1000 )
    nodbuf (ipos+13) = NINT( zjo_tot *c1000 )
    nodbuf (ipos+14) = NINT( roblon  *c1000 )
    nodbuf (ipos+15) = NINT( roblat  *c1000 )
    nodbuf (ipos+16) = iob_tot
    nodbuf (ipos+17) = job_tot
    nodbuf (ipos+18) = mrepflg
    nodbuf (ipos+19) = NINT( tzderr  *c1000 )
    nodbuf (ipos+20) = NINT( tzd     *c1000 )
    nodbuf (ipos+21) = NINT( zwd     *c1000 )
    nodbuf (ipos+22) = NINT( iwv     *c1000 )
    nodbuf (ipos+23) = NINT( surfpres*c1000 )
    nodbuf (ipos+24) = NINT( surftemp*c1000 )
    nodbuf (ipos+25) = NINT( surfrh  *c1000 )
    nodbuf (ipos+26) = ICHAR(ystid(1:1))
    nodbuf (ipos+27) = ICHAR(ystid(2:2))
    nodbuf (ipos+28) = ICHAR(ystid(3:3))
    nodbuf (ipos+29) = ICHAR(ystid(4:4))

    IF (num_compute > 1) THEN
!     Update number of reports
      nobbuf           = nobbuf + 1
!     Node the report belongs to
      inodbuf (nobbuf) = irnode
!     Position of the report in 'nodbuf'
      iposbuf (nobbuf) = ipos
    ENDIF
!   Update buffer length
    ibuflen          = ibuflen + irplen
          
! %%%% no correct alligment above
        
        ELSE
          istat = EOVFstat

        END IF  ! on dline 

      END IF    ! istat on Sample 
      END DO    ! on Samples 


    ELSEIF ( istat < 0 ) THEN 
      leof = .true.         ! end-of-file
  
    ELSEIF ( istat > 0 ) THEN 
!   I/O error
      WRITE( nurej,'(" I/O ERROR DETECTED ON GPS FILE")')
      yerrmsg = 'I/O ERROR DETECTED ON GPS FILE'
      CALL model_abort (my_cart_id, 1001, yerrmsg, yroutine)  

    ENDIF ! istat of header 

  
    ENDDO read_gps
!   --------------

  !----------------------------------------------------------------------
  ! Section  2.4  Sort reports according to the nodes they belong to
  !----------------------------------------------------------------------

    IF (num_compute > 1) THEN
      iallen      = 0
      inodlen (:) = 0
      Sort: DO inode = 0 , num_compute-1
        DO inobs = 1 , nobbuf
          IF (inodbuf(inobs) == inode) THEN
            ipos = iposbuf (inobs)
            ilen = nodbuf (ipos+1)
            DO i = 1 , ilen
              iallbuf (iallen+i) = nodbuf (ipos+i)
            ENDDO
            inodlen (inode+1) = inodlen(inode+1) + ilen
            iallen            = iallen           + ilen
          ENDIF
        ENDDO
      ENDDO Sort
    ENDIF  ! num_compute > 1

    IF (lcdf) THEN
      CLOSE (nurej)
      CLOSE (nustat)
    ENDIF

    IF (nexceed > 0) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=istat)
      IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER LENGTH FOR"   &
                     &," READING GPS:",I9)' ) ntstep, nmxbln
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST" &
                     &,I6)' ) nexceed
      CLOSE (nucautn)
    ENDIF
  ENDIF ! my_cart_id == 0 
 
  !----------------------------------------------------------------------
  ! Section  2.5  Scatter the reports to the nodes they belong to
  !----------------------------------------------------------------------

  IF (num_compute > 1) THEN

! Broadcast the number of elements for each node

    CALL distribute_values ( inodlen, num_compute, 0, imp_integers, icomm_cart &
                           , izmplcode )
!   ======================

! Determine the displacement of the first of the elements for each node
! within array 'iallbuf', the number of elements for the current node

    inodpos (1) = 0
    DO inode = 2 , num_compute
      inodpos (inode) = inodpos(inode-1) + inodlen(inode-1)
    ENDDO
    inodcnt = inodlen(my_cart_id+1)

! Determine 'nodbuf' by scattering the elements to the nodes they belong to

    nodbuf (:) = 0
    IF (inodpos(num_compute)+inodlen(num_compute) > 0) THEN

      CALL scatterv_values ( iallbuf, nmxbln, 0, inodlen                       &
                           , nodbuf , nmxbln, my_cart_id, num_compute          &
                           , imp_integers, icomm_cart, yerrmsg, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11177, yerrmsg, yroutine, izmplcode)
    ENDIF

  ENDIF  ! num_compute > 1


  IF (lwonl)                                                                   &
    WRITE(nupr,'(" OBS_READ_GPS:  length of first report received",     &
               &" at node number ",i3," :  ",i4," words")') my_cart_id,nodbuf(1)

  !--------------------------------------------------------------------------
  ! Section 3.  Deallocate local arrays
  !--------------------------------------------------------------------------

  IF (num_compute > 1) THEN
    DEALLOCATE ( iallbuf   , STAT = istat )
    DEALLOCATE ( iposbuf   , STAT = istat )
    DEALLOCATE ( inodbuf   , STAT = istat )
  ENDIF
  DEALLOCATE ( hsurf_tot   , STAT = istat )
  DEALLOCATE ( fland_tot   , STAT = istat )

  
  !--------------------------------------------------------------------------
  ! Section 4.  Fill GPS observation data arrays
  !--------------------------------------------------------------------------

  nexceed = 0

! start index -1 of first observation in report buffer

  ipos = 0  

  get_gps_report :   DO
! ~~~~~~~~~~~~~~~~~~~~~

! get length of next report
    IF  (ipos +1 > nmxbln)                                   EXIT get_gps_report
    irplen   = nodbuf(ipos+1) 
    IF  (irplen == 0)                                        EXIT get_gps_report

! read GPS buffer
    irplen   = nodbuf(ipos + 1)
    nobtyp   = nodbuf(ipos + 2)
    ncdtyp   = nodbuf(ipos + 3)
    nstachr  = nodbuf(ipos + 4)
    nkflag   = nodbuf(ipos + 5)
    ncentre  = nodbuf(ipos + 6)
    iobdat   = nodbuf(ipos + 7)
    ntime    = nodbuf(ipos + 8)
    taof     = REAL ( nodbuf(ipos + 9) ) /60._ireals
    nelevt   = nodbuf(ipos +10)
    zsurf    = REAL ( nodbuf(ipos +11), ireals ) * c1000r
    zio_tot  = REAL ( nodbuf(ipos +12), ireals ) * c1000r
    zjo_tot  = REAL ( nodbuf(ipos +13), ireals ) * c1000r
    roblon   = REAL ( nodbuf(ipos +14), ireals ) * c1000r
    roblat   = REAL ( nodbuf(ipos +15), ireals ) * c1000r
    iob_tot  = nodbuf(ipos +16)
    job_tot  = nodbuf(ipos +17)
    mrepflg  = nodbuf(ipos +18)
    mpassiv  = ibits( mrepflg , nvpsbp , nibits(nvscoc) ) * 2
    tzderr   = REAL ( nodbuf(ipos +19), ireals ) * c1000r
    tzd      = REAL ( nodbuf(ipos +20), ireals ) * c1000r
    zwd      = REAL ( nodbuf(ipos +21), ireals ) * c1000r
    iwv      = REAL ( nodbuf(ipos +22), ireals ) * c1000r
    surfpres = REAL ( nodbuf(ipos +23), ireals ) * c1000r
    surftemp = REAL ( nodbuf(ipos +24), ireals ) * c1000r
    surfrh   = REAL ( nodbuf(ipos +25), ireals ) * c1000r
    ystid(1:1)=ACHAR( nodbuf(ipos +26) )
    ystid(2:2)=ACHAR( nodbuf(ipos +27) )
    ystid(3:3)=ACHAR( nodbuf(ipos +28) )
    ystid(4:4)=ACHAR( nodbuf(ipos +29) )
    DO icl = 5, ilstidp
      ystid (icl:icl) = ' '
    ENDDO


! Assigned grid point in the processor domain
    iob     = iob_tot - isubpos (my_cart_id,1) + 1 + nboundlines
    job     = job_tot - isubpos (my_cart_id,2) + 1 + nboundlines

! check whether there is space in the ODR
! ---------------------------------------
    IF (ngpob >= maxgpl) THEN 
      nexceed  =  nexceed + 1
      IF (lwonl)                                                               &
        WRITE( nupr,'(" CAUTION !!!!! STATION ",A ,": ",I5,                    &
                     &"th GPS REPORT EXCEEDS ODR SIZE MAXGPL ",I5)')           &
               ystid, ngpob, maxgpl
      PRINT         '(" CAUTION !!!!! STATION ",A ,": ",I5,                    &
                     &"th GPS REPORT EXCEEDS ODR SIZE MAXGPL")' ,              &
               ystid, ngpob
      IF (lcdf) THEN
        icma  = i_cma ( nobtyp , ncdtyp )
!       IF (mpassiv == 0)                                                      &
           neventr (nesodr,icma,1) = neventr(nesodr,icma,1) + 1
        IF (mpassiv == 0) cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
        IF (mpassiv >  0) cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
      ELSE
        CALL obs_pointrs ( nobtyp , ncdtyp )
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!       IF (mpassiv == 0)                                                      &
          neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
        IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      ENDIF
      ipos = ipos + irplen
      CYCLE get_gps_report
    ENDIF

! create ODR station characteristic word, incl. flags and instr. specification
! ----------------------------------------------------------------------------
    nstachr = mrepflg
    IF ((iwv >  epsy) .AND. (.NOT. lgpsbias))  nbufbit = 36
    IF ((iwv >  epsy) .AND. (      lgpsbias))  nbufbit = 37
    IF ((iwv <= epsy) .AND. (.NOT. lgpsbias))  nbufbit = 38
    IF ((iwv <= epsy) .AND. (      lgpsbias))  nbufbit = 39
    nstachr = insert( nstachr , nbufbit , nvinbp )

! fill GPS observation data arrays
!------------------------------------------------------------------------------

    ngpob = ngpob + 1       
    ogphed (ngpob,nhilon) = roblon
    ogphed (ngpob,nhjlat) = roblat
    ogphed (ngpob,nhalt ) = REAL ( nelevt, ireals )
    ogphed (ngpob,nhtime) = taof
    ogphed (ngpob,nhsurf) = zsurf
    ogphed (ngpob,nhzio ) = zio_tot
    ogphed (ngpob,nhzjo ) = zjo_tot
 
    mogphd (ngpob,nhio  ) = iob
    mogphd (ngpob,nhjo  ) = job
    mogphd (ngpob,nhitot) = iob_tot
    mogphd (ngpob,nhjtot) = job_tot
    mogphd (ngpob,nhobtp) = nobtyp
    mogphd (ngpob,nhcode) = ncdtyp
    mogphd (ngpob,nhschr) = nstachr
!   mogphd (ngpob,nhqofl) = 0
    mogphd (ngpob,nhpass) = mpassiv
    mogphd (ngpob,nhqcfw) = -1

    mogphd (ngpob,nhflag) = nkflag
    mogphd (ngpob,nhcorr) = imdi
    mogphd (ngpob,nhcat ) = 0
    mogphd (ngpob,nhcats) = 99
    mogphd (ngpob,nhkz  ) = imdi
    mogphd (ngpob,nhcent) = ncentre
    mogphd (ngpob,nhstid) = imdi
    mogphd (ngpob,nhdate) = iobdat
    mogphd (ngpob,nhhrmn) = ntime
    mogphd (ngpob,nhsyhr) = MOD( iobdat,1000000 ) *100 + ntime/100
 
    yogphd (ngpob) = ystid
    DO icl = ilstidp, ilstid
      yogphd(ngpob) (icl:icl) = ' '
    ENDDO
 
! Write values if  >= 0 on the ODR (tzd has been already checked)
    ogpbdy (ngpob,nbgzpd) = tzd 
    IF ( tzderr   >= c0 )  ogpbdy (ngpob,nbgtze) = tzderr
    IF ( zwd      >= c0 )  ogpbdy (ngpob,nbgzwd) = zwd
    IF ( iwv      >= c0 )  ogpbdy (ngpob,nbgiwv) = iwv
    IF ( surfpres >= c0 )  ogpbdy (ngpob,nbgp)   = surfpres *100._ireals
    IF ( surftemp >= c0 )  ogpbdy (ngpob,nbgt)   = surftemp
! Temperature before March 2001 was in C
!   IF ( surftemp >= -99._ireals ) ogpbdy (ngpob,nbgt)   = surftemp
    IF ( surfrh   > epsy)  ogpbdy (ngpob,nbgrh)  = surfrh /100._ireals
 
    mogpbd (ngpob,nbgflg) = 0
    mogpbd (ngpob,nbgqcf) = 0
    mogpbd (ngpob,nbglid) = insert( 0, 1, nvlidp(7) )
!CS
    mogpbd (ngpob,nbgerr) = 0
    IF ( tzderr   >= c0 )  mogpbd (ngpob,nbgerr) = insert( 0, 1, nvriwv )
!CS

! Update ipos for next report  
    ipos = ipos + irplen 


  !--------------------------------------------------------------------------
  ! Section 5.  Bias correction to GPS Integrated water vapour.
  !             ---------------    
  !             This old version is cancelled !
  !             It is based on statistics GPS minus LM1 (102 sta., Aug. 2002)
  !             depends on daytime, and is approximated by function
  !             obs_bias_corr_iwv: from 18 to 8  UTC :  0.2 mm, 
  !                                from  8 to 14 UTC :  0.12*time - 0.72 mm
  !                                from 14 to 18 UTC : -0.12*time + 2.40 mm .
  !             The corrected 'observed' IWV value is
  !             IWVcorr = IWVorig - obs_bias_corr_iwv
  !--------------------------------------------------------------------------

!   IF (lgpsbias) THEN
!!    READ (ydate_ini,'(i4,i2,i2,i2)') imoyy,imomm,imodd,imohh
!     hhday = imohh + INT (taof + 24._ireals + epsy )
!     x2 = 24._ireals
!     hhday = MOD( hhday, x2 )
!csc --> include term 24. also for fraction of hour (for obs time < 0)
!     miday = (24._ireals + taof) - INT(24._ireals+taof)
!     ttday = hhday + miday
!     ogpbdy (ngpob,nbgbia) = - obs_bias_corr_iwv ( ttday , imomm )
!!                              =============
!   ELSE
      ogpbdy (ngpob,nbgbia) = c0
!   ENDIF


  !--------------------------------------------------------------------------
  ! Section 6.  Redundancy check
  !--------------------------------------------------------------------------

    rscale   =  100._ireals*(1._ireals+epsy)
    rdegkm2  =  111._ireals **2
    nstcor   =  ibits(mogphd(ngpob,nhschr),nvscbp,nibits(nvscoc))
    cosolat  =  COS( (startlat_tot+(job_tot-1)*dlat) * degrad )

!  Preset ODR index for check against all prior reports
    ngplo  = 0

!  Get next report as 'second' report
!  ----------------------------------
    next_gps_report: DO
!   ~~~~~~~~~~~~~~~~~~~
      ngplo  = ngplo + 1
!  Leave 'loop' if no redundancy with present report
      IF (ngplo == ngpob)                                   EXIT next_gps_report
!  Check if current / second report is (partly) redundant
!  ------------------------------------------------------

!  Check time difference
      lcoloc = (ABS( ogphed(ngplo,nhtime) - ogphed(ngpob,nhtime) )  <=  rtmlim)
      IF (.NOT. lcoloc)                                    CYCLE next_gps_report
!  Check horizontal (and vertical) (quasi-) colocation and station id
      rhzob  =   (cosolat *(mogphd(ngplo,nhitot) -iob_tot) * dlon)**2          &
               +          ((mogphd(ngplo,nhjtot) -job_tot) * dlat)**2
      lcoloc =      (rhzob *rdegkm2 <= rhzlim *rhzlim)                         &
              .AND. (ABS( ogphed(ngplo,nhalt)-ogphed(ngpob,nhalt) ) <=  rvtlim)
      IF (.NOT. lcoloc)                                    CYCLE next_gps_report
      IF (.NOT. (yogphd(ngplo) == yogphd(ngpob)))          CYCLE next_gps_report
!  Decide which report is redundant
!  --------------------------------

      nstcor2  =  ibits(mogphd(ngplo,nhschr),nvscbp,nibits(nvscoc))
      istcor   =  0
      IF ((mogphd(ngpob,nhpass) < 2) .AND. (mogphd(ngplo,nhpass) == 2)) THEN
        nodract  =  ngpob
      ELSEIF ((mogphd(ngpob,nhpass) == 2) .AND. (mogphd(ngplo,nhpass) < 2)) THEN
        nodract  =  ngplo
      ELSEIF (nstcor > nstcor2) THEN
        nodract  =  ngpob
        istcor   =  1
      ELSEIF (nstcor < nstcor2) THEN
        nodract  =  ngplo
        istcor   =  1
      ELSEIF ((tzderr > epsy) .AND. (     (ogpbdy(ngplo,nbgtze) > tzderr)      &
                                     .OR. (ogpbdy(ngplo,nbgtze) < rmdich))) THEN
        nodract  =  ngpob
      ELSE
        nodract  =  ngplo
      ENDIF
      nodrrej  =  ngpob + ngplo - nodract

      nobtypr  =  mogphd(nodrrej,nhobtp)
      ncdtypr  =  mogphd(nodrrej,nhcode)
      mpassiv  =  mogphd(nodrrej,nhpass)
      IF (lprodr)   ystidr   =  yogphd(nodrrej)

!  Reject redundant data
!  ---------------------

      IF ((lprcs ) .AND. (lverpas)) THEN
        mogphd (nodrrej,nhpass)  = 2
        mogphd (nodrrej,nhschr)  = insert( mogphd(nodrrej,nhschr) , 1 , nvrdbp )
      ELSE
        IF (nodract == ngpob) THEN
          ogpbdy (ngplo ,1:mxgbdy) = ogpbdy (nodract,1:mxgbdy)
          mogpbd (ngplo ,1:mxgbdf) = mogpbd (nodract,1:mxgbdf)
          ogphed (ngplo ,1:mxghed) = ogphed (nodract,1:mxghed)
          mogphd (ngplo ,1:mxghdf) = mogphd (nodract,1:mxghdf)
          yogphd (ngplo )          = yogphd (nodract)
        ENDIF
        ogpbdy (ngpob ,1:mxgbdy) = rmdi
        mogpbd (ngpob ,1:mxgbdf) = imdi
        ogphed (ngpob ,1:mxghed) = rmdi
        mogphd (ngpob ,1:mxghdf) = imdi
        yogphd (ngpob )          = '        '
      ENDIF

!  Update statistics and report events and print message
!  -----------------------------------------------------

      IF (lcdf) THEN
        icma  = i_cma ( nobtypr , ncdtypr )
        IF (mpassiv < 2) THEN
          cma(icma)%cnt_ac = cma(icma)%cnt_ac - 1
          neventr (neredn,icma,1) = neventr(neredn,icma,1) + 1
        ELSE
          cma(icma)%cnt_ps = cma(icma)%cnt_ps - 1
        ENDIF
        IF ((lprcs ) .AND. (lverpas)) THEN
          cma(icma)%cnt_ps = cma(icma)%cnt_ps + 1
        ELSE
          cma(icma)%cnt_rj = cma(icma)%cnt_rj + 1
        ENDIF
      ELSE
        CALL obs_pointrs ( nobtypr , ncdtypr )
        IF (mpassiv < 2) THEN
          noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
          neventr (neredn,ncdtpp,nobtpp) = neventr(neredn,ncdtpp,nobtpp) + 1
        ELSE
          noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        ENDIF
        IF ((lprcs ) .AND. (lverpas)) THEN
          noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
        ELSE
          noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        ENDIF
      ENDIF

      IF (lprodr) THEN
        ilen = 15 + 2*istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+ 1) = ilen
          outbuf(nacout+ 2) = nfmt10
          outbuf(nacout+ 3) = ncdtypr
          DO icl = 1, 5
            outbuf(nacout+ 3+icl) = imdi
          ENDDO
          IF (ogpbdy(nodract,nbgiwv) > rmdich)                                 &
            outbuf(nacout+ 4) = INT(ogpbdy(nodract,nbgiwv)*rscale)
          IF (ogpbdy(nodract,nbgzwd) > rmdich)                                 &
            outbuf(nacout+ 5) = INT(ogpbdy(nodract,nbgzwd)*rscale)
          IF (ogpbdy(nodract,nbgzpd) > rmdich)                                 &
            outbuf(nacout+ 6) = INT(ogpbdy(nodract,nbgzpd)*rscale)
          IF (ogpbdy(nodract,nbgtze) > rmdich)                                 &
            outbuf(nacout+ 7) = INT(ogpbdy(nodract,nbgtze)*rscale)
          IF (ogpbdy(nodract,nbgp  ) > rmdich)                                 &
            outbuf(nacout+ 8) = INT(ogpbdy(nodract,nbgp  )*rscale)
          outbuf(nacout+ 9) = mogphd(nodract,nhcode)
          outbuf(nacout+10) = INT(ogphed(nodract,nhilon)*rscale)
          outbuf(nacout+11) = INT(ogphed(nodract,nhjlat)*rscale)
          outbuf(nacout+12) = 1 + MAX( i0 , MIN( i1, nodract - nodrrej ) )
          outbuf(nacout+13) = istcor
          outbuf(nacout+14) = 0
          outbuf(nacout+15) = INT(ogphed(nodract,nhtime)*rscale)
 
          DO icl = 1 , istrej
            outbuf(nacout+15       +icl) = ICHAR( ystidr          (icl:icl) )
            outbuf(nacout+15+istrej+icl) = ICHAR( yogphd(nodract) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
! Adjust total number of GPS reports

      IF (.NOT. ((lprcs ) .AND. (lverpas)))  ngpob = ngpob - 1
      IF (.NOT. ((lprcs ) .AND. (lverpas)))                 EXIT next_gps_report

    ENDDO next_gps_report
!   ~~~~~~~~~~~~~~~~~~~~~

  ENDDO get_gps_report
! ~~~~~~~~~~~~~~~~~~~~
ntotgp = ngpob
  DEALLOCATE ( nodbuf    , STAT = izlocstat )


! number of GPS obs not used due to insufficient ODR size
! -------------------------------------------------------

  IF (num_compute > 1) THEN
    CALL global_values ( nexceed, 1,'MAX', imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11014, yerr, yroutine)
!                   ----------------
  ENDIF

  IF ((my_cart_id == 0) .AND. (nexceed > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
    WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL UP-AIR SING-L."      &
                   &," OBS. INCR. BEYOND maxgpl ",I5)' ) ntstep, nexceed, maxgpl
    WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxgpo BY AT LEAST"   &
                   &,I6)' ) nexceed
    CLOSE (nucautn)
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_read_gps
!-------------------------------------------------------------------------------
END SUBROUTINE obs_read_gps


!-------------------------------------------------------------------------------
!+ Extra module procedure in "OBS_PROCESSING" to check for and open GPS file
!-------------------------------------------------------------------------------

SUBROUTINE obs_gps_init

!-------------------------------------------------------------------------------
!
! Description :
!   This module procedure of module "src_obs_processing" initialises the reading
!   from a dedicated GPS data input file:
!   1. Check existence of GPS input file ('lexigps')
!   2: Open GPS input file, if it exists.
!
! Current Code Owner:  Ch. Schraff
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       23.08.04   Original code.    Christoph Schraff
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments: None
! --------------------

! Local variables:
! ---------------

  INTEGER (KIND=iintegers) ::  &
    nstat, ierr, jerr      ! error status variables

  CHARACTER (LEN=70)       ::  &
    yerrmsl                ! error message
  CHARACTER (LEN=20)       ::  &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerr                   ! error message

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_gps_init
!-------------------------------------------------------------------------------

  yroutine = 'obs_gps_init'

  ierr = 0
  jerr = 0
  IF (my_cart_id == 0) THEN
    INQUIRE (FILE=yugps, EXIST=lexigps)

!   IF ((.NOT. lexigps) .AND. ((gnudggp > epsy) .OR. (lgps))) THEN
!     WRITE( yerrmsl,'("ERROR: FILE ",A ," MISSING (or set gnudggp=0 and "     &
!                    &,"lgps=.FALSE.)")' )  yugps
!     PRINT '("!! either produce file ''",A,"'' , or set "                     &
!           &,"''gnudggp = 0'' and ''lgps =.FALSE.'' !!")' , yugps
!     CALL model_abort (my_cart_id, 1013, yerrmsl, yroutine)
!   ENDIF
    IF (lexigps) jerr = 1
  ENDIF
  IF (num_compute > 1)                                                         &
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
!   ==================

  lexigps = (jerr >= 1) 

  IF (lexigps) THEN
! GPS input file exists
    IF (my_cart_id == 0) THEN
      OPEN (nugps,FILE=yugps,FORM='FORMATTED',STATUS='OLD',IOSTAT=nstat)
      IF ((nstat /= 0) .OR. (ierr /= 0)) THEN
        WRITE( yerrmsl,'("OPENING OF FILE yugps FAILED: 0 /= ",I5,",",I5)' )   &
               nstat, ierr
        CALL model_abort (my_cart_id, 1014, yerrmsl, yroutine)
      ENDIF
      REWIND nugps
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! End of module procedure obs_gps_init
!-------------------------------------------------------------------------------
END SUBROUTINE obs_gps_init


!-------------------------------------------------------------------------------
!+ Extra module procedure in "OBS_PROCESSING" to compute IWV bias correction
!-------------------------------------------------------------------------------

FUNCTION obs_bias_corr_iwv ( ttday , month )

!-------------------------------------------------------------------------------
!
! Description :
!   Compute seasonal daytime-dependent bias correction for GPS IWV data
!   IWV corrected = IWV original - obs_bias_corr_iwv
!   The bias correction is an approximation of the bias GPS - LM1MO.
!
! Current Code Owner:  M. Tomassini
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       09.11.02   Original code.    Maria Tomassini
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Declarations:
!
!-------------------------------------------------------------------------------

! Function arguments:
! ------------------

  REAL    (KIND=ireals   ), INTENT (IN)         ::       &
    ttday                                      ! hour of the day

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    month                                      ! month

! Local variables:
! ---------------

  REAL (KIND = ireals)     , PARAMETER  :: &
    cseason (12)  =  (/ 0.00_ireals , 0.07_ireals , 0.25_ireals                &
                      , 0.50_ireals , 0.75_ireals , 0.93_ireals                &
                      , 1.00_ireals , 0.93_ireals , 0.75_ireals                &
                      , 0.50_ireals , 0.25_ireals , 0.07_ireals /)
                !  = 1 for July,  = 0 for January, sin-function in between, i.e.
                !  (1 + SIN( (month - 4) *PI /6 ) /2   .

  REAL (KIND=ireals)                  ::      &
    obs_bias_corr_iwv                          ! return value of function


!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Function obs_bias_corr_iwv
!-------------------------------------------------------------------------------

  IF     ((ttday <   8.0_ireals) .OR.  (ttday >  18.0_ireals)) THEN
    obs_bias_corr_iwv =  0.2_ireals
  ELSEIF ((ttday >=  8.0_ireals) .AND. (ttday <  14.0_ireals)) THEN
    obs_bias_corr_iwv =  0.12_ireals * ttday -0.72_ireals
  ELSEIF ((ttday >= 14.0_ireals) .AND. (ttday <= 18.0_ireals)) THEN
    obs_bias_corr_iwv = -0.12_ireals * ttday +2.4_ireals
  ENDIF

! seasonal scaling;
! 'cseason' has to be scaled again, because unscaled bias correction is valid
! for months MAY - AUGUST rather than only for JULY

  obs_bias_corr_iwv = obs_bias_corr_iwv * cseason(month)                       &
                                *4/(cseason(5)+cseason(6)+cseason(7)+cseason(8))

!-------------------------------------------------------------------------------
! End of module function obs_bias_corr_iwv
!-------------------------------------------------------------------------------
END FUNCTION obs_bias_corr_iwv


!-------------------------------------------------------------------------------
!+ Module procedure in "OBS_PROCESSING" for printing the observation statistics
!-------------------------------------------------------------------------------

SUBROUTINE obs_print_statist

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "obs_processing" prints a summary of the
!   statistics on the processed observational reports.
!
! Method:
!   The diagnostic arrays are inspected to produce a summary of the reports
!   processed and rejected.
!   Formatted printing. (External 'pointrs': find the diagnostic array position)
!   If the program is run in parallel mode, the diagnostic arrays summed up
!   over all nodes.
!   23.10.97
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
!  Description and writing PE changed.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO.
! 1.31       1999/07/01 Christoph Schraff
!  Input to call of 'global_values' adjusted.
! 1.36       2000/02/24 Michael Buchhold
!  Introduction of ACAR aircraft reports.
! 2.13       2002/01/18 Michael Buchhold
!  Introduction of Wind Profiler / RASS reports.
! 2.19       2002/10/24 Michael Buchhold
!  Introduction of Radar VAD wind reports.
! 3.3        2003/04/22 Christoph Schraff
!  New report events and warning messages for of GPS reports.
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

! Subroutine arguments: None
! --------------------

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    iobtpr           ,& ! number of processed reports of current observat. type
    iobtac           ,& ! number of  active   reports of current observat. type
    iobtps           ,& ! number of  passive  reports of current observat. type
    iobtrj           ,& ! number of rejected  reports of current observat. type
    jconto           ,& ! loop index over code types
    jcdt             ,& ! loop index over code types
    jobt             ,& ! loop index over observation types
    jznoct           ,& ! 1-d loop index over observation and code types
    izerror             ! for error status

  CHARACTER (LEN= 5)       ::  yform10, yform11                      ! formats
  CHARACTER (LEN=38)       ::  yform72                               ! format
  CHARACTER (LEN=45)       ::  yform71                               ! format
  CHARACTER (LEN=47)       ::  yform70                               ! format
  CHARACTER (LEN=59)       ::  yform18                               ! format
  CHARACTER (LEN=82)       ::  yform17                               ! format

  CHARACTER (LEN=80)       ::  yzerrmsg

! Local (automatic) arrays:
! -------------------------

  INTEGER (KIND=iintegers) ::  &
    inoctpr (ntotob*nmxcdt) ,& ! 1-d array of processed reports
    inoctac (ntotob*nmxcdt) ,& ! 1-d array of  active   reports
    inoctps (ntotob*nmxcdt) ,& ! 1-d array of  passive  reports
    inoctrj (ntotob*nmxcdt)    ! 1-d array of rejected  reports
!
!------------ End of header ----------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '

!-------------------------------------------------------------------------------
! Begin Subroutine obs_print_statist
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Summing of report counters over all nodes
!-------------------------------------------------------------------------------


    IF (num_compute > 1) THEN

! proc. 'global_values' expects a REAL vector as input
      DO     jcdt   =  1   ,   nmxcdt
        DO   jobt   =  1   ,   ntotob
          jznoct  =  (jobt-1) *nmxcdt + jcdt
          inoctpr (jznoct)  =  noctpr(jobt,jcdt)
          inoctac (jznoct)  =  noctac(jobt,jcdt)
          inoctps (jznoct)  =  noctps(jobt,jcdt)
          inoctrj (jznoct)  =  noctrj(jobt,jcdt)
        ENDDO
      ENDDO

      CALL global_values ( inoctpr, ntotob*nmxcdt, 'SUM', imp_integers         &
                         , icomm_cart, 0, yzerrmsg, izerror )
!     ===================

      CALL global_values ( inoctac, ntotob*nmxcdt, 'SUM', imp_integers         &
                         , icomm_cart, 0, yzerrmsg, izerror )
!     ===================

      CALL global_values ( inoctps, ntotob*nmxcdt, 'SUM', imp_integers         &
                         , icomm_cart, 0, yzerrmsg, izerror )
!     ===================

      CALL global_values ( inoctrj, ntotob*nmxcdt, 'SUM', imp_integers         &
                         , icomm_cart, 0, yzerrmsg, izerror )
!     ===================

      IF (my_cart_id == 0) THEN
        DO     jcdt   =  1   ,   nmxcdt
          DO   jobt   =  1   ,   ntotob
            jznoct  =  (jobt-1) *nmxcdt + jcdt
            noctpr (jobt,jcdt)  =  inoctpr(jznoct)
            noctac (jobt,jcdt)  =  inoctac(jznoct)
            noctps (jobt,jcdt)  =  inoctps(jznoct)
            noctrj (jobt,jcdt)  =  inoctrj(jznoct)
          ENDDO
        ENDDO
      ENDIF ! my_cart_id == 0

    ENDIF ! num_compute > 1


!-------------------------------------------------------------------------------
!  Section 2: Distribution counting
!-------------------------------------------------------------------------------


  IF (my_cart_id == 0) THEN



!          2.0   preliminaries

!          2.0.1 set formats

      yform10 = '("0")'
      yform11 = '("1")'
      yform17 (1 :54) = '("0 ---DISTRIBUTION OF PROCESSED/ACTIVE/PASSIVE/REJECT'
      yform17 (55:78) = 'ED REPORTS FOR NUDGING")'
      yform18 = '("0",42X," processed   active  passive rejected")'
      yform71 = '("0 --- observation type ",I1,A12,8X,4I9)'
      yform72 = '("      code type ",I3,6X,A19,4I9)'
      yform70 = '("0 --- total number of reports  ",13X,4I9)'

!          2.0.2 print headers

      WRITE(nustat,yform11)
      WRITE(nustat,yform10)
      WRITE(nustat,yform17)
      WRITE(nustat,yform18)
      WRITE(nustat,yform10)

!          2.0.3 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.1   observation type 1 (six code types)

      DO   jconto   =   1,   6
        iobtpr           =  noctpr(nsynop,jconto) + iobtpr
        iobtac           =  noctac(nsynop,jconto) + iobtac
        iobtps           =  noctps(nsynop,jconto) + iobtps
        iobtrj           =  noctrj(nsynop,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(nsynop,nsrscd)
      WRITE(nustat,yform71) nsynop,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.1.1 SYNOP code type 11

      WRITE(nustat,yform72) nsrscd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,1),noctac(nsynop,1)                  &
                           ,noctps(nsynop,1),noctrj(nsynop,1)

!          2.1.2 SYNOP code type 14

      CALL obs_pointrs(nsynop,natscd)
      WRITE(nustat,yform72) natscd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,2),noctac(nsynop,2)                  &
                           ,noctps(nsynop,2),noctrj(nsynop,2)

!          2.1.3 SYNOP code type 21

      CALL obs_pointrs(nsynop,nshscd)
      WRITE(nustat,yform72) nshscd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,3),noctac(nsynop,3)                  &
                           ,noctps(nsynop,3),noctrj(nsynop,3)

!          2.1.4 SYNOP code type 22

      CALL obs_pointrs(nsynop,nabscd)
      WRITE(nustat,yform72) nabscd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,4),noctac(nsynop,4)                  &
                           ,noctps(nsynop,4),noctrj(nsynop,4)

!          2.1.5 SYNOP code type 23

      CALL obs_pointrs(nsynop,nshred)
      WRITE(nustat,yform72) nshred,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,5),noctac(nsynop,5)                  &
                           ,noctps(nsynop,5),noctrj(nsynop,5)

!          2.1.6 SYNOP code type 24

      CALL obs_pointrs(nsynop,natshs)
      WRITE(nustat,yform72) natshs,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nsynop,6),noctac(nsynop,6)                  &
                           ,noctps(nsynop,6),noctrj(nsynop,6)

!          2.1.7 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.2   observation type 2 (five code types)

      DO   jconto   =   1,   5
        iobtpr           =  noctpr(nairep,jconto) + iobtpr
        iobtac           =  noctac(nairep,jconto) + iobtac
        iobtps           =  noctps(nairep,jconto) + iobtps
        iobtrj           =  noctrj(nairep,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(nairep,ncodar)
      WRITE(nustat,yform71) nairep,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.2.1 AIREP code type 141

      WRITE(nustat,yform72) ncodar,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nairep,1),noctac(nairep,1)                  &
                           ,noctps(nairep,1),noctrj(nairep,1)

!          2.2.2 AIREP code type 41

      CALL obs_pointrs(nairep,naircd)
      WRITE(nustat,yform72) naircd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nairep,2),noctac(nairep,2)                  &
                           ,noctps(nairep,2),noctrj(nairep,2)

!          2.2.3 AIREP code type 241

      CALL obs_pointrs(nairep,ncolba)
      WRITE(nustat,yform72) ncolba,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nairep,3),noctac(nairep,3)                  &
                           ,noctps(nairep,3),noctrj(nairep,3)

!          2.2.4 AMDAR code type 144

      CALL obs_pointrs(nairep,namdar)
      WRITE(nustat,yform72) namdar,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nairep,4),noctac(nairep,4)                  &
                           ,noctps(nairep,4),noctrj(nairep,4)

!          2.2.5 ACAR  code type 244

      CALL obs_pointrs(nairep,nacar )
      WRITE(nustat,yform72) nacar ,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(nairep,5),noctac(nairep,5)                  &
                           ,noctps(nairep,5),noctrj(nairep,5)

!          2.2.6 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.3   observation type 3 (two code types)

      DO   jconto   =   1,   2
        iobtpr          =  noctpr(nsatob,jconto) + iobtpr
        iobtac          =  noctac(nsatob,jconto) + iobtac
        iobtps          =  noctps(nsatob,jconto) + iobtps
        iobtrj          =  noctrj(nsatob,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(nsatob,nstbcd)
      WRITE(nustat,yform71) nsatob,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj
      IF (lsatob) THEN

!          2.3.1 SATOB code type 88

        WRITE(nustat,yform72) nstbcd,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatob,1),noctac(nsatob,1)                &
                             ,noctps(nsatob,1),noctrj(nsatob,1)

!          2.3.2 SATOB code type 188

        CALL obs_pointrs(nsatob,nsst)
        WRITE(nustat,yform72) nsst,chobcd(nobtpp,ncdtpp)                       &
                             ,noctpr(nsatob,2),noctac(nsatob,2)                &
                             ,noctps(nsatob,2),noctrj(nsatob,2)

!          2.3.4 reset counters

      ENDIF
      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.4   observation type 4 (three code types)

      DO   jconto   =   1,   3
        iobtpr         =  noctpr(ndribu,jconto) + iobtpr
        iobtac         =  noctac(ndribu,jconto) + iobtac
        iobtps         =  noctps(ndribu,jconto) + iobtps
        iobtrj         =  noctrj(ndribu,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(ndribu,nbathy)
      WRITE(nustat,yform71) ndribu,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.4.1 DRIBU code type 165

      WRITE(nustat,yform72) nbathy,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ndribu,1),noctac(ndribu,1)                  &
                           ,noctps(ndribu,1),noctrj(ndribu,1)

!          2.4.2 DRIBU code type 63

      CALL obs_pointrs(ndribu,ntesac)
      WRITE(nustat,yform72) ntesac,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ndribu,2),noctac(ndribu,2)                  &
                           ,noctps(ndribu,2),noctrj(ndribu,2)

!          2.4.3 DRIBU code type 64

      CALL obs_pointrs(ndribu,ndrbcd)
      WRITE(nustat,yform72) ndrbcd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ndribu,3),noctac(ndribu,3)                  &
                           ,noctps(ndribu,3),noctrj(ndribu,3)

!          2.4.4 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.5   observation type 5 (five code types)

      DO   jconto   =   1,   5
        iobtpr         =  noctpr(ntemp ,jconto) + iobtpr
        iobtac         =  noctac(ntemp ,jconto) + iobtac
        iobtps         =  noctps(ntemp ,jconto) + iobtps
        iobtrj         =  noctrj(ntemp ,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(ntemp,nldtcd)
      WRITE(nustat,yform71) ntemp,chobtp(nobtpp)                               &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.5.1 TEMP code type 35

      WRITE(nustat,yform72) nldtcd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ntemp ,1),noctac(ntemp ,1)                  &
                           ,noctps(ntemp ,1),noctrj(ntemp ,1)

!          2.5.2 TEMP code type 36

      CALL obs_pointrs(ntemp,nshtcd)
      WRITE(nustat,yform72) nshtcd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ntemp ,2),noctac(ntemp ,2)                  &
                           ,noctps(ntemp ,2),noctrj(ntemp ,2)

!          2.5.3 TEMP code type 135

      CALL obs_pointrs(ntemp,ntdrop)
      WRITE(nustat,yform72) ntdrop,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ntemp ,3),noctac(ntemp ,3)                  &
                           ,noctps(ntemp ,3),noctrj(ntemp ,3)

!          2.5.4 TEMP code type 39

      CALL obs_pointrs(ntemp,nrocob)
      WRITE(nustat,yform72) nrocob,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ntemp ,4),noctac(ntemp ,4)                  &
                           ,noctps(ntemp ,4),noctrj(ntemp ,4)

!          2.5.5 TEMP code type 40

      CALL obs_pointrs(ntemp,nrocsh)
      WRITE(nustat,yform72) nrocsh,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(ntemp ,5),noctac(ntemp ,5)                  &
                           ,noctps(ntemp ,5),noctrj(ntemp ,5)

!          2.5.6 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.6   observation type 6 (two code types)

      DO   jconto   =   1,   2
        iobtpr          =  noctpr(npilot,jconto) + iobtpr
        iobtac          =  noctac(npilot,jconto) + iobtac
        iobtps          =  noctps(npilot,jconto) + iobtps
        iobtrj          =  noctrj(npilot,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(npilot,nldpcd)
      WRITE(nustat,yform71) npilot,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.6.1 PILOT code type 32

      WRITE(nustat,yform72) nldpcd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,1),noctac(npilot,1)                  &
                           ,noctps(npilot,1),noctrj(npilot,1)

!          2.6.2 PILOT code type 33

      CALL obs_pointrs(npilot,nshpcd)
      WRITE(nustat,yform72) nshpcd,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,2),noctac(npilot,2)                  &
                           ,noctps(npilot,2),noctrj(npilot,2)

!          2.6.3 PILOT code type 132 European wind profiler

      CALL obs_pointrs(npilot,nwp_eu)
      WRITE(nustat,yform72) nwp_eu,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,3),noctac(npilot,3)                  &
                           ,noctps(npilot,3),noctrj(npilot,3)

!          2.6.4 PILOT code type 133 European sodar/rass report

      CALL obs_pointrs(npilot,nra_eu)
      WRITE(nustat,yform72) nra_eu,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,4),noctac(npilot,4)                  &
                           ,noctps(npilot,4),noctrj(npilot,4)

!          2.6.5 PILOT code type 136 US wind profiler / rass report

      CALL obs_pointrs(npilot,npr_us)
      WRITE(nustat,yform72) npr_us,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,5),noctac(npilot,5)                  &
                           ,noctps(npilot,5),noctrj(npilot,5)

!          2.6.6 PILOT code type 137 Radar VAD wind report

      CALL obs_pointrs(npilot,nravad)
      WRITE(nustat,yform72) nravad,chobcd(nobtpp,ncdtpp)                       &
                           ,noctpr(npilot,6),noctac(npilot,6)                  &
                           ,noctps(npilot,6),noctrj(npilot,6)


!          2.6.7 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.7   observation type 7 (two code types)

      DO   jconto   =   1,   2
        iobtpr          =  noctpr(nsatem,jconto) + iobtpr
        iobtac          =  noctac(nsatem,jconto) + iobtac
        iobtps          =  noctps(nsatem,jconto) + iobtps
        iobtrj          =  noctrj(nsatem,jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(nsatem,nsmsg1)
      WRITE(nustat,yform71) nsatem,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj

!          2.7.1 SATEM code type  71 MSG_1

      IF (mcdmsg1 >= 1) THEN
        CALL obs_pointrs(nsatem,nsmsg1)
        WRITE(nustat,yform72) nsmsg1,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,1),noctac(nsatem,1)                &
                             ,noctps(nsatem,1),noctrj(nsatem,1)
      ENDIF

!          2.7.2 SATEM code type  72 MSG_2
  
      IF (mcdmsg1 >= 1) THEN
        CALL obs_pointrs(nsatem,nsmsg2)
        WRITE(nustat,yform72) nsmsg2,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,2),noctac(nsatem,2)                &
                             ,noctps(nsatem,2),noctrj(nsatem,2)
      ENDIF

!          2.7.3 SATEM code type 206 ATOVS_NOAA15
  
      IF (mcdno15 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa15)
        WRITE(nustat,yform72) nnoa15,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,3),noctac(nsatem,3)                &
                             ,noctps(nsatem,3),noctrj(nsatem,3)
      ENDIF

!          2.7.4 SATEM code type 207 ATOVS_NOAA16

      IF (mcdno16 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa16)
        WRITE(nustat,yform72) nnoa16,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,4),noctac(nsatem,4)                &
                             ,noctps(nsatem,4),noctrj(nsatem,4)
      ENDIF

!          2.7.5 SATEM code type 208 ATOVS_NOAA17

      IF (mcdno17 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa17)
        WRITE(nustat,yform72) nnoa17,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,5),noctac(nsatem,5)                &
                             ,noctps(nsatem,5),noctrj(nsatem,5)
      ENDIF

!          2.7.6 SATEM code type 209 ATOVS_NOAA18

      IF (mcdno18 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa18)
        WRITE(nustat,yform72) nnoa18,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(nsatem,6),noctac(nsatem,6)                &
                             ,noctps(nsatem,6),noctrj(nsatem,6)
      ENDIF

!          2.7.7 reset counters

      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0

!          2.8   observation type 8 (one code types)

      DO   jconto   =   1,   1
        iobtpr          =  noctpr(MIN( ngps,8 ),jconto) + iobtpr
        iobtac          =  noctac(MIN( ngps,8 ),jconto) + iobtac
        iobtps          =  noctps(MIN( ngps,8 ),jconto) + iobtps
        iobtrj          =  noctrj(MIN( ngps,8 ),jconto) + iobtrj
      ENDDO

      CALL obs_pointrs(ngps  ,ngpgfz)
      WRITE(nustat,yform71) ngps  ,chobtp(nobtpp)                              &
                           ,iobtpr,iobtac,iobtps,iobtrj
      IF (lexigps) THEN

!          2.8.1 GPS   code type 96

        WRITE(nustat,yform72) ngpgfz,chobcd(nobtpp,ncdtpp)                     &
                             ,noctpr(MIN( ngps,8 ),1), noctac(MIN( ngps,8 ),1) &
                             ,noctps(MIN( ngps,8 ),1), noctrj(MIN( ngps,8 ),1)

!          2.8.2 reset counters

      ENDIF
      iobtpr=0
      iobtac=0
      iobtps=0
      iobtrj=0


!-------------------------------------------------------------------------------
!  Section 3: Total number of reports
!-------------------------------------------------------------------------------


!          3.1   sum them up

      DO     jcdt   =  1   ,   nmxcdt
        DO   jobt   =  1   ,   ntotob
          iobtpr          =  noctpr(jobt,jcdt) + iobtpr
          iobtac          =  noctac(jobt,jcdt) + iobtac
          iobtps          =  noctps(jobt,jcdt) + iobtps
          iobtrj          =  noctrj(jobt,jcdt) + iobtrj
        ENDDO
      ENDDO

      WRITE(nustat,yform70) iobtpr,iobtac,iobtps,iobtrj

!          3.2   print notes.

      WRITE(nustat,'("0")' )
      WRITE(nustat,'("  --- Notes on the table above:")' )
      WRITE(nustat,'(6X,''"Rejected"/"passive" means that the whole ''      &
                      &,''report is rejected /set passive.'')' )
      WRITE(nustat,'(6X,''Partly rejected and partly passive reports are '' &
                      &,''labeled "active".'')' )
      WRITE(nustat,'(6X,''A report can be labeled "active" even if part ''  &
                      &,''of its data is black listed.'')' )
      WRITE(nustat,'("0")' )
      WRITE(nustat,'(6X,''Reports may only be rejected or set passive ''    &
                      &,''for reasons given by report'')' )
      WRITE(nustat,'(6X,''  events 1 - 11 , except for events 3 and 5 ''    &
                      &,''(on station altitude) for'')' )    
      WRITE(nustat,'(6X,''  TEMPs and PILOTs. Hence, the number of the''    &
                      &,''se events must equal the'')' )
      WRITE(nustat,'(6X,''  number of rejected and passive reports.'')' )
      WRITE(nustat,'(6X,''In the verification mode (i.e if data are wr''    &
                      &,''itten to the VOF), reports are'')' )
      WRITE(nustat,'(6X,''  rejected in case of report events 1 - 3, 8''    &
                      &,'' - 11, and event 4 if the re-'')' )    
      WRITE(nustat,'(6X,''  port is outside the model domain. Otherwis''    &
                      &,''e the reports are set passive,'')' )    
      WRITE(nustat,'(6X,''  except that already passive reports are al''    &
                      &,''so rejected if they do not'')' )    
      WRITE(nustat,'(6X,''  contain any data or if they are redundant ''    &
                      &,'' and a subset of an active'')' )    
      WRITE(nustat,'(6X,''  report.'')' )    
      WRITE(nustat,'(6X,''  Without verification, reports are always r''    &
                      &,''ejected for events 1 to 11.'')' )    

  ENDIF ! (my_cart_id == 0)


!-------------------------------------------------------------------------------
! End of module procedure obs_print_statist
!-------------------------------------------------------------------------------

END SUBROUTINE obs_print_statist


!-------------------------------------------------------------------------------
!+ Module procedure in "OBS_PROCESSING" for printing report or data events
!-------------------------------------------------------------------------------

SUBROUTINE obs_print_events ( lrepevt , nevent , mxeve , nfrsev , nlasev )

!-------------------------------------------------------------------------------
!
! Description:
!   This module procedure of module "obs_processing" prints a summary of
!   report or data events. The printed report events are ordered according
!   to the order of the checks.
!   Additionally, WARNING messages are issued, if the size of the observation
!   Data Record (ODR) as specified by namelist parameters is too small to
!   accommodate all observations.
!
! Method:
!   The events table is inspected to produce a summary of the events.
!   Formatted printing. (External 'pointrs': find the diagnostic array position)
!   If the program is run in parallel mode, the diagnostic arrays summed up
!   over all nodes.
!   23.10.97
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
!  Revised report and data event. Formats, description, and writing PE changed.
! 1.29       1999/05/11 Ulrich Schaettler
!  Adapted interfaces to utility-modules and prepared use of MPE_IO.
! 1.31       1999/07/01 Christoph Schraff
!  Input to call of 'global_values' adjusted.
! 1.36       2000/02/24 Michael Buchhold
!  Introduction of ACAR aircraft reports.
! 2.5        2001/06/01 Christoph Schraff
!  New data events for pressure tendency, lapse rate check, wind shear checks.
!  New report event for flight track check and flight report thinning.
! 2.13       2002/01/18 Michael Buchhold
!  New report events for Wind Profiler/RASS reports.
! 2.19       2002/10/24 Michael Buchhold
!  New report events for Radar VAD wind reports.
! 3.3        2003/04/22 Christoph Schraff
!  New report events and warning messages for of GPS reports.
! 3.6        2003/12/11 Christoph Schraff
!  SATOB taken into account for warning message related to array size.
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

! Subroutine arguments:
! --------------------

  LOGICAL                 , INTENT (IN)         ::       &
    lrepevt             ! if .TRUE: print report events, else print data events

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    mxeve            ,& ! length of first dimension of record 'nevent'
    nfrsev              ! number of first event to be printed

  INTEGER (KIND=iintegers), INTENT (INOUT)      ::       &
    nlasev              ! number of last  event to be printed

  INTEGER (KIND=iintegers), INTENT (IN)         ::       &
    nevent (mxeve,nmxcdt,ntotob)  ! record containing the events

! Local parameters: None
! ----------------

! Local scalars:
! -------------

  INTEGER (KIND=iintegers) ::  &
    jrev,jevt        ,& ! loop index over report events
    jdev             ,& ! loop index over  data  events
    nrevc            ,& ! loop index over        events
    jcdt             ,& ! loop index over code types
    jobt             ,& ! loop index over observation types
    jevc             ,& ! adjusted index over events
    jzevent             ! 1-d loop index over events, observation and code types

  CHARACTER (LEN=14)       ::  yevents
  CHARACTER (LEN=62)       ::  y1900 , yform90 , yform91 , yform92   ! formats
  CHARACTER (LEN=12)       ::  yform12                               ! format 
  CHARACTER (LEN= 5)       ::  yform10, yform11                      ! formats


! Local (automatic) arrays:
! -------------------------

  INTEGER (KIND=iintegers) ::  &
    nevtno  (mxeve)               ,& ! event indices
    inevent (mxeve*nmxcdt*ntotob) ,& ! 1-d real array containing the events
    ievent  (mxeve,nmxcdt,ntotob) ,& ! record containing all events globally
    nwarn   (4)                   ,& ! local max of events due to ODR size limit
    nwarn0  (4)                      ! events due to ODR size limit at node 0

  INTEGER (KIND=iintegers) ::  izerror 
  CHARACTER (LEN=80)       ::  yzerrmsg
!
!------------ End of header ----------------------------------------------------

  izerror  = 0
  yzerrmsg = '   '

!-------------------------------------------------------------------------------
! Begin Subroutine obs_print_events
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Section 1: Summing of event counters over all nodes
!-------------------------------------------------------------------------------


  IF (lrepevt) THEN
    nlasev = MIN( nlasev , mxcrev )
  ELSE
    nlasev = MIN( nlasev , mxcdev )
  ENDIF

  ievent = nevent

! prepare printing of warning messages

  IF (lrepevt) THEN
    nwarn (:) = 0
    DO     jcdt   =  1   ,   nmxcdt
      nwarn (1) =   nwarn(1)                   + nevent(nesodr,jcdt,nsynop)    &
                  + nevent(nesodr,jcdt,nairep) + nevent(nesodr,jcdt,nsatob)    &
                  + nevent(nesodr,jcdt,ndribu)
      nwarn (2) =   nwarn(2)                   + nevent(nesodr,jcdt,ntemp )    &
                  + nevent(nesodr,jcdt,npilot)
      nwarn (3) =   nwarn(3)                   + nevent(nesodr,jcdt,MIN(ngps,8))
      nwarn (4) =   nwarn(4)                   + nevent(nesodr,jcdt,nsatem)
    ENDDO
    IF (my_cart_id == 0) THEN
      nwarn0 (:) = nwarn(:)
      nwarn  (:) = 0
    ENDIF
    DO     jcdt   =  1   ,   nmxcdt
      nwarn (2) =   nwarn(2)                   + nevent(nenoml,jcdt,nairep)
    ENDDO
  ENDIF

  IF (num_compute > 1) THEN

    DO       nrevc  =  nfrsev  ,   nlasev
      jevc          =  nrevc - nfrsev + 1
      DO     jcdt   =  1   ,   nmxcdt
        DO   jobt   =  1   ,   ntotob
          jzevent   =  (jevc-1) *ntotob *nmxcdt  +  (jcdt-1) *ntotob  +  jobt
          inevent (jzevent)  =  nevent (nrevc,jcdt,jobt)
        ENDDO
      ENDDO
    ENDDO

    CALL global_values ( inevent, ntotob*nmxcdt*(nlasev-nfrsev+1), 'SUM'       &
                       , imp_integers, icomm_cart, 0, yzerrmsg, izerror)
!   ===================

    IF (lrepevt) CALL global_values ( nwarn, 4, 'MAX', imp_integers            &
                                    , icomm_cart, 0, yzerrmsg, izerror )
!                ==================

    IF (my_cart_id == 0) THEN
      DO       nrevc  =  nfrsev  ,   nlasev
        jevc          =  nrevc - nfrsev + 1
        DO     jcdt   =  1   ,   nmxcdt
          DO   jobt   =  1   ,   ntotob
            jzevent   =  (jevc-1) *ntotob *nmxcdt  +  (jcdt-1) *ntotob  +  jobt
            ievent  (nrevc,jcdt,jobt)  =  inevent (jzevent)
          ENDDO
        ENDDO
      ENDDO
    ENDIF ! my_cart_id == 0

  ENDIF ! num_compute > 1

!-------------------------------------------------------------------------------
!  Section 1b: Printing of warning messages
!-------------------------------------------------------------------------------

  IF ((lrepevt) .AND. (my_cart_id == 0)) THEN
    nwarn (1) = MAX( nwarn(1) , NINT(  nwarn0(1) /                          &
                         SQRT( REAL (num_compute, ireals) ) + 0.49_ireals ) )
    nwarn (2) = MAX( nwarn(2) , NINT(  nwarn0(2) /                          &
                         SQRT( REAL (num_compute, ireals) ) + 0.49_ireals ) )
    nwarn (3) = MAX( nwarn(3) , NINT(  nwarn0(3) /                          &
                         SQRT( REAL (num_compute, ireals) ) + 0.49_ireals ) )
    nwarn (4) = MAX( nwarn(4) , NINT(  nwarn0(4) /                          &
                         SQRT( REAL (num_compute, ireals) ) + 0.49_ireals ) )
    IF (MAX( nwarn(1),nwarn(2),nwarn(3),nwarn(4) ) > 0) THEN
      WRITE( nustat,'(''1'')' )
      WRITE( nustat,'(''0'')' )
      WRITE( nustat,'(''+     !!! CAUTION !!!!! CAUTION !!!!! ''               &
                    &,''CAUTION !!!!! CAUTION !!!!! CAUTION !!!!!'')' )
      WRITE( nustat,'(''+     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''               &
                    &,''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'')' )
    ENDIF
    IF (nwarn(1) > 0)                                                          &
      WRITE( nustat,'(''0     WARNING: array size for single-level ''          &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXSGO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(1)
    IF (nwarn(2) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for multi-level ''         &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXMLO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(2)
    IF (nwarn(3) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for GPS ''                 &
                    &,''observations is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXGPO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(3)
    IF (nwarn(4) > 0)                                                          &
        WRITE( nustat,'(''0     WARNING: array size for SAT ''                 &
                    &,''retrievals   is too small'', /                         &
                    &,''      =======  to accommodate all observations'', /    &
                    &,''       --> Increase "MAXTVO" (namelist) by'',I5        &
                    &,'': usually ok for local obs. array'', /                 &
                    &,''           =================                   ''      &
                    &,'' (possibly still insufficient for'', /                 &
                    &,''                                               ''      &
                    &,'' the global obs increment array!)'')' )        nwarn(4)
  ENDIF


!-------------------------------------------------------------------------------
!  Section 2: Events distribution
!-------------------------------------------------------------------------------


  IF (my_cart_id == 0) THEN



!          2.0   preliminaries: formats and headers

!          2.0.1 set formats

      yform10 = '("0")'
      yform11 = '("0")'
      yform12 = '("0*",A8)'
      yform90 = '("+",A10,I3,I5,I4,I5,I4,I3,I5,I3,I5,I3,I5,2I3,I3,I4,I5,I3) '
      yform91 = '("+",A10,I4,I2,I3,I2,I3,I4,I5,I3,I5,I4,I3,I4,I3,I5,2I3,I5,I3) '
      yform92 = '("+",A10,I4,I2,I2,I3,I3,I4,I4,I2,I4,I4,I4,I4,3I3,I5,I2,2I3,I2)'

!          2.0.1 print definition of events

      WRITE(nustat,yform11)
      WRITE(nustat,yform10)
 
      IF (lrepevt) THEN
        nlasev = MIN( nlasev , mxcrev )
        WRITE( nustat,'(''1 *** REPORT EVENTS DEFINITIONS (THEIR ''            &
                       &,''ORDER MATCHES THE ORDER OF THE CHECKS):'')')
        DO  jrev  = nfrsev , nlasev
          WRITE( nustat,'(8X,A)') crepev(jrev)
        ENDDO
        y1900 = yform90
      ELSE
        nlasev = MIN( nlasev , mxcdev )
        IF (nlasev == 1) THEN
          WRITE( nustat,'(''1 *** DATA EVENTS DEFINITIONS (LEVEL ''            &
                        &,'' EVENTS APPLY TO MULTI-LEVEL DATA ONLY,'')')
          WRITE( nustat,'(''                               THE OR''            &
                        &,''DER OF ALL EVENTS MATCHES THE ORDER OF'')')
          WRITE( nustat,'(''                               THE CH''            &
                        &,''ECKS EXCEPT FOR EVENT 8)'')')
        ELSE
          WRITE( nustat,'(''1 *** DATA EVENTS DEFINITIONS (CONTINUED):'')')
        ENDIF
        DO  jdev  = nfrsev , nlasev
          WRITE( nustat,'(8X,A)') cdatev(jdev)
        ENDDO
        y1900 = yform92
        IF (nfrsev == 1) y1900 = yform91
      ENDIF

!          2.0.2 print header for events

      WRITE(nustat,yform11)
      WRITE(nustat,yform10)
      DO nrevc = nfrsev , nlasev
        nevtno (nrevc) = nrevc
      ENDDO
      yevents = '  events      '
      WRITE(nustat,y1900) yevents, (nevtno(nrevc), nrevc=nfrsev,nlasev)
      WRITE(nustat,'('' '',''  ------ '')')

!          2.1   observation type 1 (six code types)

      CALL obs_pointrs(nsynop,nsrscd)
      WRITE(nustat,yform12) chobtp(nobtpp)

!          2.1.1 SYNOP code type 11

      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,1,nsynop), jevt=nfrsev,nlasev)

!          2.1.2 SYNOP code type 14

      CALL obs_pointrs(nsynop,natscd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,2,nsynop), jevt=nfrsev,nlasev)

!          2.1.3 SYNOP code type 21

      CALL obs_pointrs(nsynop,nshscd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,3,nsynop), jevt=nfrsev,nlasev)

!          2.1.4 SYNOP code type 22

      CALL obs_pointrs(nsynop,nabscd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,4,nsynop), jevt=nfrsev,nlasev)

!          2.1.5 SYNOP code type 23

      CALL obs_pointrs(nsynop,nshred)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,5,nsynop), jevt=nfrsev,nlasev)

!          2.1.6 SYNOP code type 24

      CALL obs_pointrs(nsynop,natshs)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,6,nsynop), jevt=nfrsev,nlasev)

!          2.2   observation type 2 (four code types)

      CALL obs_pointrs(nairep,ncodar)
      WRITE(nustat,yform12) chobtp(nobtpp)

!          2.2.1 AIREP code type 141

      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,1,nairep), jevt=nfrsev,nlasev)

!          2.2.2 AIREP code type 41

      CALL obs_pointrs(nairep,naircd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,2,nairep), jevt=nfrsev,nlasev)

!          2.2.3 AIREP code type 241

      CALL obs_pointrs(nairep,ncolba)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,3,nairep), jevt=nfrsev,nlasev)

!          2.2.4 AIREP code type 144 (amdar)

      CALL obs_pointrs(nairep,namdar)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,4,nairep), jevt=nfrsev,nlasev)

!          2.2.5 AIREP code type 244 (acar)

      CALL obs_pointrs(nairep,nacar)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,5,nairep), jevt=nfrsev,nlasev)


!          2.3   observation type 3 (two code types)

      IF (lsatob) THEN
        CALL obs_pointrs(nsatob,nstbcd)
        WRITE(nustat,yform12) chobtp(nobtpp)

!          2.3.1 SATOB code type 88

        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,1,nsatob), jevt=nfrsev,nlasev)

!          2.3.2 SATOB code type 188

        CALL obs_pointrs(nsatob,nsst)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,2,nsatob), jevt=nfrsev,nlasev)
      ENDIF

!          2.4   observation type 4 (three code types)

      CALL obs_pointrs(ndribu,nbathy)
      WRITE(nustat,yform12) chobtp(nobtpp)

!          2.4.1 DRIBU code type 165

      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,1,ndribu), jevt=nfrsev,nlasev)

!          2.4.2 DRIBU code type 63

      CALL obs_pointrs(ndribu,ntesac)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,2,ndribu), jevt=nfrsev,nlasev)

!          2.4.3 DRIBU code type 64

      CALL obs_pointrs(ndribu,ndrbcd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,3,ndribu), jevt=nfrsev,nlasev)

!          2.5   observation type 5 (five code types)

      CALL obs_pointrs(ntemp,nldtcd)
      WRITE(nustat,yform12) chobtp(nobtpp)

!          2.5.1 TEMP code type 35

      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,1,ntemp ), jevt=nfrsev,nlasev)

!          2.5.2 TEMP code type 36

      CALL obs_pointrs(ntemp,nshtcd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,2,ntemp ), jevt=nfrsev,nlasev)

!          2.5.3 TEMP code type 135

      CALL obs_pointrs(ntemp,ntdrop)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,3,ntemp ), jevt=nfrsev,nlasev)

!          2.5.4 TEMP code type 39

      CALL obs_pointrs(ntemp,nrocob)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,4,ntemp ), jevt=nfrsev,nlasev)

!          2.5.5 TEMP code type 40

      CALL obs_pointrs(ntemp,nrocsh)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,5,ntemp ), jevt=nfrsev,nlasev)

!          2.6   observation type 6 (six code types)

      CALL obs_pointrs(npilot,nldpcd)
      WRITE(nustat,yform12) chobtp(nobtpp)

!          2.6.1 PILOT code type 32

      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,1,npilot), jevt=nfrsev,nlasev)

!          2.6.2 PILOT code type 33

      CALL obs_pointrs(npilot,nshpcd)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,2,npilot), jevt=nfrsev,nlasev)

!          2.6.3 PILOT code type 132  European wind profiler report

      CALL obs_pointrs(npilot,nwp_eu)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,3,npilot), jevt=nfrsev,nlasev)

!          2.6.4 PILOT code type 133  European sodar/rass report

      CALL obs_pointrs(npilot,nra_eu)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,4,npilot), jevt=nfrsev,nlasev)

!          2.6.5 PILOT code type 136  US profiler/rass report

      CALL obs_pointrs(npilot,npr_us)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,5,npilot), jevt=nfrsev,nlasev)

!          2.6.6 PILOT code type 137  Radar VAD wind report

      CALL obs_pointrs(npilot,nravad)
      WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                                &
                         , (ievent(jevt,6,npilot), jevt=nfrsev,nlasev)

!          2.7   observation type 7 (two code types)

      IF (MAX( mcdmsg1,mcdmsg2,mcdno15,mcdno16,mcdno17,mcdno18 ) >= 1) THEN
        CALL obs_pointrs(nsatem,nsmsg1)
        WRITE(nustat,yform12) chobtp(nobtpp)
      ENDIF

!          2.7.1 SATEM code type  71 MSG_1

      IF (mcdmsg1 >= 1) THEN
        CALL obs_pointrs(nsatem,nsmsg1)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,1,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.7.2 SATEM code type  72 MSG_2

      IF (mcdmsg2 >= 1) THEN
        CALL obs_pointrs(nsatem,nsmsg2)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,2,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.7.3 SATEM code type 206 ATOVS NOAA15

      IF (mcdno15 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa15)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,3,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.7.4 SATEM code type 207 ATOVS NOAA16

      IF (mcdno16 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa16)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,4,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.7.5 SATEM code type 208 ATOVS NOAA17

      IF (mcdno17 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa17)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,5,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.7.6 SATEM code type 209 ATOVS NOAA18

      IF (mcdno18 >= 1) THEN
        CALL obs_pointrs(nsatem,nnoa18)
        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,6,nsatem), jevt=nfrsev,nlasev)
      ENDIF

!          2.8   observation type 8 (one code type )

      IF (lexigps) THEN
        CALL obs_pointrs(ngps  ,ngpgfz)
        WRITE(nustat,yform12) chobtp(nobtpp)

!          2.8.1 GPS   code type 96

        WRITE(nustat,y1900) chobcd(nobtpp,ncdtpp)                              &
                           , (ievent(jevt,1,MIN( ngps,8 ) ), jevt=nfrsev,nlasev)
      ENDIF


  ENDIF ! (my_cart_id == 0)


!-------------------------------------------------------------------------------
! End of module procedure obs_print_events
!-------------------------------------------------------------------------------

END SUBROUTINE obs_print_events


!-------------------------------------------------------------------------------

END MODULE src_obs_processing
