!+ Source module for the observation processing in the data assimilation mode
!-------------------------------------------------------------------------------

MODULE src_obs_proc_aof

!-------------------------------------------------------------------------------
! Description:
!   This module performs the basic parts of the observation processing of
!   reports read from the AOF. 
!   Special tasks are:
!    - reading reports from AOF and distributing them to the PEs        
!    - extraction of data and saving observation header information and
!      observation body data in data arrays for single and multi-level reports
!    - removal of redundant information
!
! Method:
!   This module contains the following module procedures:
!    - obs_aof_init           : called by organize_obs_proc (src_obs_processing)
!    - obs_read_aof_org       : called by organize_obs_proc (src_obs_processing)
!    - obs_extract
!    - obs_sort_levels        : called by     obs_multi_level
!    - obs_flag6b             : called by       obs_option_groups
!    - obs_option_groups      : called by     obs_single_level
!    - obs_RASS_tv2t          : called by     obs_multi_level
!    - obs_read_distribute    : called by   obs_read_aof_org
!    - obs_save_head          : called by   obs_read_aof_org
!    - obs_multi_level        : called by   obs_read_aof_org
!    - obs_single_level       : called by   obs_read_aof_org
!    - obs_redundancy         : called by   obs_read_aof_org
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
!  Initial release, based on some previous modules of src_obs_processing.
!  No redundancy checking between PILOTs and profilers.
!  Avoid verification of reports that are completely contained in other reports.
!  Very small correction to avoid sub-saturation due to rounding errors.
! 3.15       2005/03/03 Ulrich Schaettler
!  Replaced FLOAT by REAL
! 3.18       2006/03/03 Christoph Schraff
!  Caution messages on insufficient ODR size to specific file unit 'nucautn'.
!  REAL replaced by REAL of kind ireals.
! 3.21       2006/12/04 Burkhardt Rockel
!  Introduced variable polgam
! V3_24        2007/04/26 Christoph Schraff
!  Special grid point assignment for Lindenberg reports.
! V4_8         2009/02/16 Oliver Fuhrer
!  Use global values only if num_compute greater 1
! V4_9         2009/07/16 Ulrich Schaettler 
!  Replaced some quotes with double quotes (problems with IBM preprocessor)
! V4_12        2010/05/11 Ulrich Schaettler
!  Renamed t0 to t0_melt because of conflicting names
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_22        2012/01/31 Christoph Schraff
!  - Observed 'temperature variable' (e.g. virtual temperature) not converted
!    into temperature any more. Some variables moved from 'data_obs_process' and
!    'data_nudge_all' to new modules 'data_obs_cdfin' and 'data_obs_lib_cosmo'.
!  - Bug correction at selection of surface observations (usage of 'fdoro').
!  - Call of external routine 'atmhp' (libmisc) replaced by 'std_atmosphere'.
!  - Call of external routine 'difmin' replaced by new routine 'diff_minutes
!    from utilities (type of routine arguments: iintegers instead of intgribf).
!  - Convert aircraft roll angle from 4-bit flag into 2-bit flag.
!  - Some (but not all) new ODR entries filled by correct values.
!  - Temporary: Set all radiosondes to Vaisala RS92 except a few specified ones.
! V4_23        2012/05/10 Ulrich Schaettler
!  Adapted interface to SR diff_minutes (added itype_calendar)
! V4_24        2012/06/22 Hendrik Reich
!  Adapted length of strings for date variables
!  Read also minutes and seconds from start date
! V4_28        2013/07/12 Christoph Schraff
!  Direct call of mpi routine replaced by call of 'scatterv_values'.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!       Version  Name      : Old Story (of src_obs_processing)
!       -------  ----        ---------
!       V 1.13   Buchhold (1998/10/22) : Initial release.
!       V 1.15   Schraff   : Global declaration of allocatable arrays moved to
!                            module 'data_obs_process'.
!       V 1.19   Schraff   : Verification mode: complete obs processing also for
!                            data not used for data assimilation. Revised report
!                            and data events. ANSI violations removed.
!       V 1.20   Doms, 1999: Renaming of some global variables
!       V 1.26   Doms      : Bug corr of double declaration of 3 variables.
!       V 1.27   Schraff   : Revised VOF format and ODR flag formats. Extreme
!                            temperatur group introduced in AOF optional groups.
!                            6-bit hollerith station ID replaced by characters.
!                            Revised (32-bit) missing data indicator in ODR.
!       V 1.29   Schaettler: Adapted interfaces to utility-modules, prepared use
!                            of MPE_IO and corr calls to intrinsics on non-Cray.
!       V 1.31   Schraff   : MPI calls replaced by parallel_utility calls;
!                            quantities related to 'icomm_world' removed.
!       V 1.34   Schaettler: Use new timing routines
!       V 1.36   Schraff   : VOF obs increment format included. Adapt to 32-bit
!                 (2000)     AOF. Introduction of ACAR aircraft reports.
!       V 1.38   Schraff   : Index for pressure tendency, redundant report flag,
!                            aircraft roll angle code,... added to ODR and VOF.
!       V 1.39   Schaettler: Changed global var names for grid specification.
!       V 1.40   Schraff   : Corr of time criterion for writing obs statistics.
!       V 1.43   Bettems   : Use of 'irealgrib' (for obs_single_level).
!       V 2.4    Schraff   : 'nmxbln' reduced.
!                 (2001)
!       V 2.5    Schraff   : Introduction of lapse rate and wind shear check for
!                            multi-level data. Adapt aircraft report listing to
!                            cope with flight track checks. Additional data
!                            events, and savety test at array allocations.
!       V 2.6    Schraff   : Messaging on thinning sequences of aircraft reps.
!       V 2.12   Schraff   : Flight track check done also without production of
!                            multi-level aircraft reports.
!       V 2.13   Buchhold  : Introduction of Wind Profiler / RASS reports.
!                 (2002)     Funct. 'ichcon' removed. Thinning and flight track
!                            check flags moved to station characteristics word
!       V 2.19   Buchhold  : Introduction of Radar VAD wind reports.
!       V 3.3    Tomassini + Schraff : Extension for reading from a dedicated
!                 (2003)     file and processing of GPS-derived IWV data.
!       V 3.6    Schraff   : Introduction of obs type SATOB.
!                            Ensuring error message at model_abort.
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

!------------------------------------------------------------------------------
! Section 1 : Analysis Observation File formats                  
!------------------------------------------------------------------------------
!
!               data description record format
!               ------------------------------
   nlmaxd,     & ! maximum data record length
   ntiop ,     & ! initial time of observation period
   nteop ,     & ! end of observation period

!         1.1   observation file parameters
!               ---------------------------
   nmxlob,     & ! maximum length of AOF data record
   nmxbln,     & ! max. length of buffer containing all reports for all PE's
   lastob,     & ! .TRUE. if all observations have been processed

!         1.2   AOF data record header format
!               -----------------------------
   nhdob ,     & ! header len.
   nrecln,     & ! record len.
   nobtp ,     & ! observation type
   ncdtp ,     & ! code type
   nlatit,     & ! observation latitude
   nlongt,     & ! observation longitude
   ndate ,     & ! date of observation
   nsyntm,     & ! synoptic time of observation
   nextim,     & ! exact time of observation
   nstid1,     & ! part 1 of station id
   nstid2,     & ! part 2 of station id
   naltit,     & ! station altitude
   nstcar,     & ! station characteristics
   ninstr,     & ! instrument specifications
   nqllta,     & ! flags on lat./lon./dat./tim./alt.
   nbufhed,    & ! report header length in send buffer

!         1.3.1  Bit patterns for packed information
!                -----------------------------------
!                station characteristics (nstcar)
   nstcbp,     & ! bit pos. for st. correct. indic.  
   nstcoc,     & ! no. of bits occ. by st. corr. ind.
   nstibp,     & ! bit pos. for important st. indic.
   nstioc,     & ! no. of bits occ. by imp. st. ind.
   nscfbp,     & ! bit pos. for st. confidence indic.
   nscfoc,     & ! no. of bits occ. by st. conf. ind.
   narabp,     & ! bit pos. for aircraft roll angle      "
   naraoc,     & ! no. of bits occ. by airc. roll angle  
   npafbp,     & ! bit position for phase of airc. flight
   npafoc        ! no. of bits occ. by phase of flight


!         1.3.2   instrument specifications (ninstr)
!                 ----------------------------------
USE data_obs_process, ONLY :  &

   ninsbp,     & ! bit pos. for instr. spec. indic.
   ninsoc,     & ! no. of bits occ. by ins. spec. in.
   ni4bp ,     & ! bit pos. for I4 parameter             "
                 ! type of measurement equipment used    "
                 ! or data processing tehnique
   ni4oc ,     & ! of bits occ. by I4 param.             "

!         1.3.3   Bit patterns for flag information
!                 ---------------------------------
!                 Bit position for each flag
   nf3bps,     & !  bit pos. for flags on:          
   nf3boc,     & !  no. of bits occ. by each flag  
   nf3ibp,     & !  inner bit struct. for each flag
   nf3ioc,     & !  no. bits occ. within flag struc.

!         1.3.4   Constants added on certain parameters (lat/lon/alt)
!                 ---------------------------------------------------
   nltmcn,     & !  cons. added on lat.
   nlgmcn,     & !  cons. added on lon.
   nalmcn        !  cons. added on altit.

USE data_obs_process, ONLY :  &

!         1.4      AOF data record body format
!                  ---------------------------
   nwrvr ,     & ! body format matrix of var.,ob.tp.

!         1.4.1   Level len. arrays and multi/single level report switch
   nlev0,      & ! level 0 len. array of (obs. type)
   nlevn,      & ! level n len. array of (obs. type)
   lmulti,     & ! multi/single level report switch of (obs. t.)

!         1.4.2   Bit patterns for packed information (weather/cloud group)
!                 --------------------------------------------------------
   nwwbp,      & ! bit position  for ww
   nvibp,      & !         "         v
   nwwoc,      & ! no. of bits occupied by ww
   nvioc,      & !           "             v
   nchbp,      & ! bit position for ch
   ncmbp,      & !         "        cm
   nclbp,      & !         "        cl
   nnbp ,      & !         "        nn
   nhbp ,      & !         "        h
   nbp  ,      & !         "        n
   nchoc,      & ! no. of bits occupied by ch
   ncmoc,      & !           "             cm
   ncloc,      & !           "             cl
   nnoc ,      & !           "             nn
   nhoc ,      & !           "             h
   noc           !           "             n

USE data_obs_process, ONLY :  &

!         1.4.3   Additional group indic. and bit pattern (only SYNOP)
!                 ---------------------------------------------------
   naddgr,     & ! additional groups indicator       
   naclbp,     & ! bit pos. for additional clouds
   nacloc,     & ! no. of bits occ. by above indic.
   ngrnbp,     & ! bit pos. for ground group
   ngrnoc,     & ! no. of bits occ. by above indic.
   ntgbp ,     & ! bit position for tgtg                 "
   ntgoc ,     & ! no. of bits occupied tgtg             "
   nsppbp,     & ! bit pos. for special phenomena
   nsppoc,     & ! no. of bits occ. by above indic.
   nlspbp,     & ! bit position for spsp              calculated
   nuspbp,     & !          "       spsp                 "
   nlspoc,     & ! no. of bits occupied spsp             "
   nuspoc,     & !              "       spsp             "
   nicebp,     & ! bit pos. for ice group 
   niceoc,     & ! no. of bits occ. by above indic.
   nrinbp,     & ! bit pos. for rain group
   nrinoc,     & ! no. of bits occ. by above indic.
   ntrbp ,     & ! bit position for trtr 
   nrrbp ,     & !          "       rr
   ntroc ,     & ! no. of bits occupied trtr
   nrroc ,     & !              "       rr
   nshpbp,     & ! bit pos. for ship group
   nshpoc,     & ! no. of bits occ. by above indic.
   nwavbp,     & ! bit pos. for waves group 
   nwavoc,     & ! no. of bits occ. by above indic.
   nextbp,     & ! bit pos. for extreme temp. group   naddgr
   nextoc,     & ! no. of bits occ. by above indic.      "
   ntnbp ,     & ! bit position for tntn            calculated
   ntxbp ,     & !          "       txtx                 "
   ntnoc ,     & ! no. of bits occupied tntn             "
   ntxoc ,     & !              "       txtx             "
   nradbp,     & ! bit pos. for radiation group       naddgr
   nradoc,     & ! no. of bits occ. by above indic.      "
   nffbp ,     & ! bit position for ffff            calculated
   nffoc         ! no. of bits occupied ffff             "

USE data_obs_process, ONLY :  &

!         1.4.4   Flag formats (set as 2-d arr. of varb.,obs. typ.)
!                 -------------------------------------------------
   nwrvrf,     & ! flag format matrix of var.,obs.types

!         1.4.5    6 bit struct. within flag words (2-d arr.)
!                  ------------------------------------------
   n6vfb,      & ! 6bit flag patterns matrix of varib.,obs. type
!                          def. by nwrvrf matrix

!         1.4.6    Bit pattern within 6 bit structures
!                  ------------------------------------
   nf1aoc,     & ! no. bits in active part of 6bit pattern
   nf1bps,     & ! bit positions within 6bit pattern
   nf1boc,     & ! no. of bits occ. by each infor.

!         1.4.7    Other bit patterns for flags
!                  ----------------------------
   nwwfbp,     & ! bit position for flag on ww          nwrvrf matrix
   nvfbp ,     & !             "            v                "
   nchfbp,     & !             "            ch               "
   ncmfbp,     & !             "            cm               "
   nnhfbp,     & !             "            nh               "
   nclfbp,     & !             "            cl               "
   nnnfbp,     & !             "            nn               "
   nnfbp ,     & !             "            n                "
   ntgfbp,     & !             "            tgtg             "
   nlspfb,     & !             "            spsp             "
   nuspfb,     & !             "            spsp             "
   ntrfbp,     & !             "            trtr             "
   nrrfbp,     & !             "            rr               "
   nttnbp,     & !             "            tntn             "
   nttxbp,     & !             "            txtx             "
   nfffbp,     & !             "            ffff             "
   nf2boc,     & ! no. of bits occupied by each above flag

!         1.4.8   Bit pattern for lev. ind. and press. code
!                 -----------------------------------------
!                  other bit patterns within flag words
   nlinbp,     & ! bit position for level id
   nlinoc,     & ! no. of bits occupied by level id
   npcdbp,     & ! bit position for pressure code
   npcdoc,     & ! no. of bits occupied by press. code
   nlidbp,     & ! bit pattern for level id
   nlidoc        ! no. of bits occupied by each above information

USE data_obs_process, ONLY :  &

!         1.5     Extracted parameters and Variables
!                 ----------------------------------
   llship,     & ! ship switch
   lsurfob,    & ! if true then derive surface report from temp
   nobtyp,     & ! observation type
   ncdtyp,     & ! code type
   ntime,      & ! actual time of observation
   nvarib,     & ! current variable
   nwwgr ,     & ! weather group
   ngclgr,     & ! general cloud group
   naclgr,     & ! additonal cloud groups
   ngrngr,     & ! ground group
   nsppgr,     & ! special phenomena group
   nicegr,     & ! ice group
   nraing,     & ! rain group
   nshpgr,     & ! ship group
   nwavgr,     & ! waves group
   nextgr,     & ! extreme temperature group
   nradgr        ! radiation group flag word   "

USE data_obs_process, ONLY :  &

! flags to be extracted
   nvrfwr,     & ! current variable flag (bit pattern)
   nflag ,     & ! actual flag for current variable (0-3)
   nlevin,     & ! level id. (bit pattern)
   nwwfw ,     & ! weather flag word 
   ngclfw,     & ! cloud flag word
   naclfw,     & ! additonal clouds flag word 
   ngrnfw,     & ! ground group flag word
   nsppfw,     & ! special phenomena flag word
   nicefw,     & ! ice group flag word
   nranfw,     & ! rain group flag word 
   nshpfw,     & ! ship movement flag word
   nwavfw,     & ! wavws group flag word
   nextfw,     & ! extr. temp. group flag word
   nradfw,     & ! radiation group flag word   "
   taof,       & ! observation time in forecast hour units
   rppp,       & ! reported pressure level
   roblat,     & ! station latitude
   roblon,     & ! station lonitude

!         1.6     AOF time boxes
!                 --------------
   naofbox,    & ! length of AOF time periods to be read
   naoflbx,    & ! initial time of last  AOF time box to be read currently
   naoffbx,    & ! initial time of first AOF time box to be read currently

!------------------------------------------------------------------------------
! Section 2 : Numerical constants and scaling constants related to AOF
!------------------------------------------------------------------------------

   nmdi  ,     & ! missing data indicator for AOF data  (2^30-1 --> 2^59-1)
   atol  ,     & ! area tolerance (because of rounding off problem)

   npfact,     & ! pressure multiple (depending on press. code)
   nzfact,     & ! geopot. multiple (depending on press. code)
   npflfc,     & ! press. flag multiple
   nzflfc,     & ! geop. flag multiple
   rfact ,     & ! scaling factors
   rpset ,     & ! set values of press. for diff. press. codes
   rzset         ! set values of geop. for diff. press. codes

USE data_obs_process, ONLY :  &

!------------------------------------------------------------------------------
! Section 3 : Variable extraction and output numbers related to AOF
!------------------------------------------------------------------------------

!         3.1   Variables extraction numbering
!               ------------------------------
   nvar  ,     & !   variables extraction numbering:

!         3.2    Variables' extraction inventory table
!                -------------------------------------
   nzex  ,     & !   geopotential extraction table
   nuex  ,     & !   u-comp. extraction table
   nvex  ,     & !   v-comp. extraction table
   ntex  ,     & !   temperature extraction table
   ntdex ,     & !   dew-point extraction table
   nqex  ,     & !   water content extraction table
   nt2ex ,     & !   2m temperature extraction table
   ntd2ex,     & !   2m dew-point extraction table

!         3.2.1  Variables' extraction inventory table (SYNOP additional groups)
!                ---------------------------------------------------------------
   nwwex ,     & !   weather group extraction table
   ngclex,     & !   general cloud group extraction table
   naclex,     & !   additional cloud groups extraction table
   ngrnex,     & !   ground group extraction table
   nsppex,     & !   special phenomena group table
   niceex,     & !   ice group extraction table
   nraine,     & !   rain group extraction table
   nshipe,     & !   ship group extraction table
   nwavee,     & !   wave groups extraction table
   nextex,     & !   extreme temperature group extraction table
   nradex,     & !   radiation group extraction table

!         3.3    Variables' output inventory table
!                ---------------------------------
   nzot  ,     & !   geopotential output table
   nuot  ,     & !   u-comp. output table
   nvot  ,     & !   v-comp. output table
   ntot  ,     & !   temperature output table
   ntdot ,     & !   dew-point output table
   nt2ot ,     & !   2m temperature output table
   ntd2ot,     & !   dew-point output table

!         3.4    Variables' numbering for obs. error extraction
!                ----------------------------------------------
   nvrpoe        !   variables numbering for obs. error calculation:

USE data_obs_process, ONLY :  &

!-------------------------------------------------------------------------------
! Section 4 : Diagnostic arrays and pointers for event counters (AOF input only)
!-------------------------------------------------------------------------------

   noctpr,     & ! 2-d array of processed observations
   noctac,     & ! 2-d array of active observations
   noctps,     & ! 2-d array of passive observations
   noctrj,     & ! 2-d array of rejected observations
   nobtpp,     & ! diagnostic arrays pointer
   ncdtpp,     & !                 "

!------------------------------------------------------------------------------
! Section 6 : Temporary information buffers and report counters
!------------------------------------------------------------------------------
   aofbuf,     & ! observation report in aof format
   nodbuf,     & ! buffer containing all reports for a single node
   nmlob,      & ! current number of  multi-level report
   nsgob         ! current number of single-level report

! end of data_obs_process

!------------------------------------------------------------------------------

USE data_obs_cdfin, ONLY :  &

!------------------------------------------------------------------------------
! Section 3 : Data event counter arrays and diagnostics arrays
!------------------------------------------------------------------------------

!         3.1     Format of event counters
!                 ------------------------

!         3.1.1   Report event counter array format
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
   neredx,     & ! redundancy between 1 multi- and 1 single-level report
!  neslml,     & ! one multi-level report made from  single-level reports
!  neslps,     & ! single-level report put in multi-level report and set passive
!  nenoml,     & ! multi-levl report not made due to ODR array size
!  netrac,     & ! (flight) track suspicious
   nethin        ! thinning of aircraft (flight) track

USE data_obs_cdfin, ONLY :  &

!         3.1.2   Data event counter array format
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
!  nefshr,     & ! wind speed shear too large
!  nedshr,     & ! directional wind shear too large
   nerlim,     & ! prec.amount exceeds threshold limit

!         3.2     Event counter arrays
!                 --------------------
   neventr,    & ! counter array of report events
   neventd       ! counter array of data events

USE data_obs_cdfin, ONLY :  &

!------------------------------------------------------------------------------
! Section 4 : Observation errors
!------------------------------------------------------------------------------

!         4.1  Observation error levels
!              ------------------------
!  nerlev,     & ! number of standard error levels
!  rlevel,     & ! error levels
!  rolnlv,     & ! ln(rlevel(15))

!
!         4.2  Observation error constants
!              ---------------------------
   rherr1 ,     & ! (root of) fixed    / normal conditions
   rherr2 ,     & ! relative humidity <  if temperature below 233K
   rherr3 ,     & ! error variances    \ if rel. humidity below 20%

!------------------------------------------------------------------------------
! Section 5 : Different kinds of limits
!------------------------------------------------------------------------------

!      Redundancy check limits
!      -----------------------
   rtmlim ,    & !  time limit for all reports except airep    [hrs]
   rhzlim ,    & !  horiz.  dist. limit for obs sta. (\ airep)  [km]
   rvtlim ,    & !  vertic. dist. limit for obs sta. (\ airep)   [m]
   rtmlair,    & !  time limit for reports of obs type 'airep' [hrs]
                 !  (time of lowest level of multi-level ODR)
   rhzlair,    & !  horizont. distance limit for airep reports  [km]
   rvtlair,    & !  vertical  distance limit for airep reports [hpa]
   rdplim ,    & !  vertic. dist. limit within multi-level rep.[hpa]
   rprlim ,    & !  vertic. dist. limit for replacing 'missing
                 !    data' within multi-level reports           [hpa]

!      temperature/humidity/pressure/height limits
!      -------------------------------------------
   rttlim,     & !  temperature limit below which rel. hum. error is increased
   rrhlim,     & !  rel. hum. limit below which rel. hum. error is increased
   rtshlm,     & !  rel. hum. limit for extra level of temp/
                 !  synop humidity data control
   rerrpp,     & !  msl pressure above which observed pressure is
                 !  reckoned to be erroneous
   rpplim,     & !  pressure level obove which observations are not used
   rhtsat,     & !  rel. humidity threshold for saturation with real obs
   pminsig,    & !  significant-level TEMP / PILOT data at pressure
                 !       less than 10000 [pa] are neglected
   pqmin,      & !  pressure [pa] of level above which moisture obs are not used
   vfoglim,    & !  visibility threshold [m] below which the existence of low
                 !       cloud (fog) is assumed in the presence of precipitation

!      height/pressure limits for reporting practice check
!      ---------------------------------------------------
   rh300 ,     & !  300m station height limit
   rh800 ,     & !  800m station height limit
   rh1700,     & ! 1700m station height limit
   rh2300,     & ! 2300m station height limit
   rh3700        ! 3700m station height limit

USE data_obs_cdfin, ONLY :  &

!------------------------------------------------------------------------------
! Section 10: Output buffer size and formats: for reporting rejection of data
!------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
! Section 8 : Temporary global model fields
!------------------------------------------------------------------------------

   hsurf_tot,  & ! total array of model surface height
   fland_tot     ! total array of fraction of land in each grid element

! end of data_obs_process

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
    ngpgfz        !   GPS report processed by GFZ

! end of data_obs_lib_cosmo 

!-------------------------------------------------------------------------------

USE data_nudge_all, ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    mxda         ,& ! max. number of data pts 'always' / never to be nudged
    lwonl        ,& ! .TRUE if grid pt. (ionl ,jonl ) lies in the local domain
    lwonl2       ,& ! .TRUE if grid pt. (ionl2,jonl2) lies in the local domain
    lcroot       ,& ! .TRUE if (my_cart_id  == 0) (for print on std. output)
    onl_cart_id  ,& ! 'my_cart_id'  of node with area containing (ionl,jonl)
    lfirst       ,& ! .TRUE if 'organize_nudging' is called the first time
    acthr        ,& ! actual forecast hour

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

! 2.1  temporal variables
!      ------------------

    lverif ,    & ! on - off switch for verification
    lverpas,    & ! .t. : on - off switch for verif. also of passive reports
    tconbox,    & ! 6*dt: timestep [s] for computing analysis increments
                  !       (i.e. time box of constant analysis increments)
    gnudggp,    & ! 0*10^-4: nudging coefficient for GPS-derived IWV [1/s]
    tipolmx,    & ! max. time span  for linear interpolat. for upper-air data
    tipmxsu,    & ! max. time span  for linear interpolat.of surf.-level data
    wtukrsa,    & ! temporal radius of infl. towards the past for TEMP/PILOT
    wtukrse,    & ! temporal radius of infl. towards the future for TEMP/PILOT
    wtukara,    & ! temporal radius of infl. towards the past for aircraft obs
    wtukare,    & ! temporal radius of infl. towards the future for aircr.data
    wtuksua,    & ! temporal radius of infl. towards the past for surface data
    wtuksue,    & ! temporal radius of influence towards the future

! 2.3  observation processing
!      ----------------------

    yaofpath,   & !  path(name) of input observation file 'AOF'
    doromx ,    & !  cut-off and gaussian radius of height differences 
                  !  between model orography and station height
    altopsu,    & !  SYNOP obs. above height 'altopsu' are not assimilated
    maxmlo ,    & !  max. number of multi-level reports within the total domain
    maxsgo ,    & !  max. number of (surface-level and upper-air single-level
                  !  reports within the total domain
    mrhyes ,    & !  station ID's of SURFACE humidity data that are
                  !  always used (irrespective of 'doromx', 'altopsu)
    mrhno  ,    & !  station ID's of SURFACE humidity data never used
    muvyes ,    & !  station ID's of SURFACE wind data always used
    muvno         !  station ID's of SURFACE wind data never used

USE data_nudge_all, ONLY :   &

    lsynop,     & !  if SYNOP data is used
    laircf,     & !  if AIREP data is used (aircraft)
    lsatob,     & !  if SATOB data is used
    ldribu,     & !  if DRIBU data is used (drifting buoy)
    ltemp ,     & !  if TEMP  data is used
    lpilot,     & !  if PILOT data is used
    lsatem,     & !  if SATEM data is used
    lgps  ,     & !  if GPS   data is used
    lsurfa,     & !  if surface fields are analysed
    lt2m  ,     & !  if 2m temperat. field is analysed
    lrh2m ,     & !  if 2m rel. hum. field is analysed
    lprecp,     & !  if precipitation is analysed
    lpraof,     & !  .t. for diagnostic print of analysis obs file AOF
    lprodr,     & !  .t. for diagnostic print of obs data records ODR
    ldiasa,     & !  .t. for diagnostics of surface analysis
    lcd011,     & !  .t. if synop code 11 data is used
    lcd014,     & !  .t. if synop code 14 data is used
    lcd021,     & !  .t. if synop code 21 data is used
    lcd022,     & !  .t. if synop code 22 data is used
    lcd023,     & !  .t. if synop code 23 data is used
    lcd024,     & !  .t. if synop code 24 data is used
    lcd041,     & !  .t. if airrep code 41 data is used
    lcd141,     & !  .t. if airrep code 141 data is used
    lcd241,     & !  .t. if airrep code 241 data is used
    lcd144,     & !  .t. if airrep code 144 data is used
    lcd244,     & !  .t. if airrep code 244 data is used
    lcd088,     & !  .t. if satob code 88 data is used
    lcd188,     & !  .t. if satob code 188 data is used
    lcd063,     & !  .t. if dribu code 63 data is used
    lcd064,     & !  .t. if dribu code 64 data is used
    lcd165        !  .t. if dribu code 165 data is used

USE data_nudge_all, ONLY :   &

    lcd035,     & !  .t. if temp code  35 data is used
    lcd036,     & !  .t. if temp code  36 data is used
    lcd135,     & !  .t. if temp code 135 data is used
    lcd039,     & !  .t. if temp code  39 data is used
    lcd040,     & !  .t. if temp code  40 data is used
    lcd032,     & !  .t. if pilot code  32 data is used
    lcd033,     & !  .t. if pilot code  33 data is used
    lcd132,     & !  .t. if pilot code 132 data is used
    lcd133,     & !  .t. if pilot code 133 data is used
    lcd136,     & !  .t. if pilot code 136 data is used
    lcd137,     & !  .t. if pilot code 137 data is used
    lcd086,     & !  .t. if satem code  86 data is used
    lcd186,     & !  .t. if tovs  code 186 data is used

    noctrq,     & !  observation/code type of diagnstic print
    hprc  ,     & !  time of prec-ana in hours since model start
    raintp,     & !  time period of precipitation analysis
    obnlat,     & !  northern boundary of observation area
    obslat,     & !  southern boundary of observation area
    obwlon,     & !  western boundary of observation area
    obelon,     & !  eastern boundary of observation area
    exnlat,     & !  northern boundary for exclusion area
    exslat,     & !  southern boundary for exclusion area
    exwlon,     & !  western boundary for exclusion area
    exelon,     & !  eastern boundary for exclusion area
    dinlat,     & !  northern boundary of diagnost. print area
    dislat,     & !  southern boundary of diagnost. print area
    diwlon,     & !  western boundary of diagnost. print area
    dielon,     & !  eastern boundary of diagnost. print area

    ionl  ,     & !  / grid point coordinates
    jonl  ,     & !  \ for standard output
    ionl2 ,     & !  / 2nd grid pt coordinates
    jonl2         !  \ for other standard output

USE data_nudge_all, ONLY :   &

! 3. I/O device numbers and file names for nudging
! ------------------------------------------------

    nuaof  ,    & ! analysis observation file
    nuaofex,    & ! print out of expanded analysis observation file
!   yuaof  ,    & ! input observation file in AOF format
    yuaofex,    & ! expanded analysis observation file

! 4. Surface analysis limits
! --------------------------

    rraint,     & ! rain period limits
    rrainl        ! amount of rain limit

! end of data_nudge_all

!-------------------------------------------------------------------------------

USE data_obs_record, ONLY :   &
 
!       1.1    ODR header format
!       ------------------------

!       1.1.1  Header formats of ODR reports: 'omlhed' and 'osghed'
!              ----------------------------------------------------
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
    nhvcbu,     & ! correction factor to vertical correlation scale for wind
                  ! at base of report
    nhvcbt,     & ! as 'nhvcbu', but for temperature
    nhvcbq,     & ! as 'nhvcbu', but for humidity
    nhvctu,     & ! correction factor to vertical correlation scale for wind
                  ! at top of report
    nhvctt,     & ! as 'nhvctu', but for temperature
    nhvctq        ! as 'nhvctu', but for humidity

USE data_obs_record, ONLY :   &

!       1.1.2  Header formats of ODR reports: 'momlhd' and 'mosghd'
!              ----------------------------------------------------
    mxrhdf,     & ! header length of multi-level reports
    mxshdf,     & ! header length of single-level reports
    mxghdf,     & ! header length of GPS reports
    nhio  ,     & ! (local) x-coord. of grid pt. assigned to obs
    nhjo  ,     & ! (local) y-coord. of grid pt. assigned to obs
    nhitot,     & ! global x-coord. of grid pt. assigned to obs
    nhjtot,     & ! global y-coord. of grid pt. assigned to obs
    nhobtp,     & ! observation type
    nhcode,     & ! code type
    nhschr,     & ! station characteristics (packed as in VOF, see below)
    nhqofl,     & ! report flags (rds) on lat/long/date/time/alt (as in VOF)
    nhpass,     & ! flag for report being set to 'passive' (as in VOF)
    nhqcfw,     & ! threshold quality control flags for pressure, status of
                  ! verification output
    nhdate     ,& ! absolute exact observation date [yyyymmdd]
    nhhrmn     ,& ! absolute exact observation time [hhmm]
    nhrtyp     ,& ! radiosonde type    (NRARA, see WMO common code Table C2)
    nhnlev,     & ! number of obs. levels (for multi-level reports)
    nhvqcf,     & ! for satellite retrieval: threshold quality control flags
    nhaexi,     & ! flag for exist. of wind or temperature in multi-level report
    nhuexi,     & ! flag for existence of wind data        in multi-level report
    nhtexi,     & ! flag for existence of temperature data in multi-level report
    nhqexi,     & ! flag for existence of humidity data    in multi-level report


!       1.1.3  Header formats of ODR reports: 'yomlhd' and 'yosghd'
!              ----------------------------------------------------

    ilstid,     & ! character length of the station identity
    ilstidp       ! char. length used for printing the station ID
                  ! Note: (ilstid >= ilstidg >= ilstidp), cf. data_nudge_gather

USE data_obs_record, ONLY :   &

!       1.2    Bit patterns for packed information in VOF header
!              -------------------------------------------------
    nvpabp,     & ! bit pos. for report set passive since it is     nvhsch
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
    nvhbps,     & ! bit pos. for flags on lat/lon/date/time/alt.    nvhqfl
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
    nbtflg,     & ! main flag word (bit pattern as in VOF body, word 'nvbmfw')
    nbtqcf,     & ! threshold quality control flags
    nbtlid        ! level id (bit pattern, as in VOF)

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
    nbsflg,     & ! main flag word (bit pattern as in VOF body, word 'nvbmfw')
    nbsqcf,     & ! threshold quality control flags (as in VOF)
    nbslid,     & ! pressure code (SYNOP) or (else) level id. (bit pattern, VOF)
    nbscwg,     & ! combined cloud and weather group (set of classes, as in VOF)
    nbscfw        ! flag word for cloud, weather, precip., extreme temp., cf VOF

USE data_obs_record, ONLY :   &

!       1.4.2  Bit patterns for packed information in VOF body, general words
!              --------------------------------------------------------------
    nvfubp,     & ! bit pos. for main flag on wind                  nvbmfw
    nvftbp,     & ! bit pos. for main flag on temperature             "
    nvfqbp,     & ! bit pos. for main flag on humidity                "
    nvfzbp,     & ! bit pos. for main flag on pressure / geopot.      "
    nvfaoc,     & ! no. of bits occ. by each main flag                "
    nvfbps,     & ! bit pattern for main flags:                       "
    nvfboc,     & ! no. of bits occ. for main flags                   "
    nvflbp,     & ! bit pos. for level flag: level below surface      "
    nvfloc,     & ! no. of bits occ. by level flag                    "
    nvlidp,     & ! level id. bit pattern                           nvblid
    nvlido        ! no. bits occ. by each indicator in level id.      "

USE data_obs_record, ONLY :   &

!       1.4.3  Bit patterns for 'optional groups' in VOF body
!              ----------------------------------------------
                  ! combined cloud and weather group                nvbcwg
    nvchbp,     & ! bit position for ch  (type of high cloud)
    nvcmbp,     & !         "        cm  (type of middle cloud)
    nvclbp,     & !         "        cl  (type of low cloud)
    nvnhbp,     & !         "        nh  (cover of low, else of middle cloud)
    nvhbp ,     & !         "        h   (cloud base height)
    nvnbp ,     & !         "        n   (total cloud cover)
    nvwwbp,     & !         "        ww  (present weather)
                  !                      (see VUB WMO Code tables:)
    nvchoc,     & ! no. of bits occupied by ch    [Code table 0509]
    nvcmoc,     & !           "             cm    [Code table 0515]
    nvcloc,     & !           "             cl    [Code table 0513]
    nvnhoc,     & !           "             nh    [Code table 2700]
    nvhoc ,     & !           "             h     [Code table 1600]
    nvnoc ,     & !           "             n     [Code table 2700]
    nvwwoc,     & !           "             ww    [Code table 4677]
                  ! combined optional group flag word               nvbcfw
    nvufbp,     & ! bit position for flag on ch
    nvmfbp,     & !             "            cm
    nvhfbp,     & !             "            h
    nvlfbp,     & !             "            cl
    nvbfbp,     & !             "            nh
    nvnfbp,     & !             "            n
    nvwfbp,     & !             "            ww
    nvvfbp,     & !             "            v     (visibility)
    nvrfbp,     & !             "            rr    (total precipitation)
    nvcfbp,     & !             "            ccl   (low cloud cover)
    nvqfbp,     & ! bit position for refined quality flag on ccl
    nv2foc,     & ! no. of bits occupied by each of above flags
    nvgfbp,     & ! bit position for flag on tgtg  (ground temperature)
    nvtfbp,     & !             "            tntn  (minimum temperature)
    nvxfbp,     & !             "            txtx  (maximum temperature)
    nvffbp,     & !             "            fxgu  (maximum gusts)
    nvefbp,     & !             "            fxme  (maximum 10' mean wind)
    nvsfbp,     & !             "            ffff  (global radiation)
    nv1foc,     & ! no. of bits occupied by each of above flags
    nvtrbp,     & ! bit position for code of precipitation
                  ! measurement duration          [Code table 4019, keys 0-7]
    nvtroc        ! no. of bits occupied by precipitation obs. duration code

USE data_obs_record, ONLY :   &

!       1.5    Further quantities related to ODR
!              ---------------------------------

!       1.5.0  Missing data indicators in ODR's
!              --------------------------------
    imdi  ,     & ! missing data indicator for ODR integers (2^31-1)

!       1.5.1  Total number of stored reports in ODR's
!              ---------------------------------------
    ntotml,     & ! tot. number of stored multi-level reports
    ntotsg,     & ! tot. number of stored single-level reports

!       1.5.2  Use of surface-level data with orography differences
!              ----------------------------------------------------
    fdoro ,     & ! scaling factor to vertical distances betw. model

!       2.     Observation data records (ODR) and associated arrays
!       ----------------------------------------------------------------------

!       2.1    Formats of observation data records
!              -----------------------------------
    omlbdy,     & ! body of multi-level ODR
    omlhed,     & ! header of multi-level ODR
    osgbdy,     & ! body of single-level ODR
    osghed,     & ! header of single-level ODR
    momlhd,     & ! header of multi-level ODR
    momlbd,     & ! body of multi-level ODR
    mosgbd,     & ! body of single-level ODR
    mosghd,     & ! header of single-level ODR
    yomlhd,     & ! header of multi-level ODR
    yosghd,     & ! header of single-level ODR

!       3.     Masking constants
!       ------------------------

    nibits        ! masking constants

! end of data_obs_record

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

! 2. physical constants and related variables
! -------------------------------------------

    t0_melt,      & ! melting temperature of ice   
    r_d,          & ! gas constant for dry air
    g,            & ! acceleration due to gravity

! 3. constants for parametrizations
! ---------------------------------

    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i)
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i             !               -- " --

! end of data_constants

!-------------------------------------------------------------------------------

USE data_runcontrol , ONLY :   &

! 1. start and end of the forecast
! --------------------------------
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels

! 3. controlling the physics
! --------------------------
    itype_gscp,   & ! type of grid-scale precipitaiton physics

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
    scatterv_values, & ! scatters batches of variable length to other nodes
    distribute_values  ! distributes a set of values from one node to all others

!-------------------------------------------------------------------------------

 USE utilities,                ONLY :  &
    phi2phirot,      & ! converts phi from the real to the rotated system
    rla2rlarot,      & ! converts lambda from the real to the rotated system
    uv2uvrot,        & ! converts wind components from normal to rotated system
    diff_minutes       ! compute difference in minutes between 2 dates/times

!-------------------------------------------------------------------------------

 USE src_obs_proc_util,        ONLY :  &
    obs_aof_assign_gridpt,& ! horizontal assignment of obs. report to a grid pt
    get_global_surf_aof  ,& ! gather orography/land fraction from global domain
    obs_pointrs     ,& ! finding the diagnostic array position
    obs_error       ,& ! assignment of an observation error
    psurf           ,& ! computes (approx.) surface pressure at a specific point
    z2p             ,& ! converts height to pressure at a specific grid point
    tv2t               ! converts virtual temperature to temperature at a point

!-------------------------------------------------------------------------------

 USE src_obs_cdfin_util,       ONLY :  &
    std_atmosphere     ! convert variables according to the standard atmosphere

!-------------------------------------------------------------------------------
! Local declarations
!-------------------------------------------------------------------------------

IMPLICIT NONE

! 1. Parameters
! -------------

  REAL (KIND=ireals)       , PARAMETER  :: &
    fbogus  =   0.0_ireals   ! if (fbogus /= 0) then bogus data are produced by
                             !   applying factor 'fbogus' to obs 'corrections
                             !   as def. in *obs_multi_level*

!===============================================================================

!-------------------------------------------------------------------------------
! Public and Private Subroutines
!-------------------------------------------------------------------------------

CONTAINS

! INCLUDE "obs_aof_init.incf"
! INCLUDE "obs_read_aof_org.incf"
! INCLUDE "obs_extract.incf"
! INCLUDE "obs_sort_levels.incf"
! INCLUDE "obs_flag6b.incf"
! INCLUDE "obs_option_groups.incf"
! INCLUDE "obs_RASS_tv2t.incf"
! INCLUDE "obs_read_distribute.incf"
! INCLUDE "obs_save_head.incf"
! INCLUDE "obs_multi_level.incf"
! INCLUDE "obs_single_level.incf"
! INCLUDE "obs_redundancy.incf"


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for initialising the reading from AOF
!-------------------------------------------------------------------------------

SUBROUTINE obs_aof_init

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" initialises the reading
!   from the AOF ('Analysis Observation File') file.
!
! Method:
!   Opening of files. Checking of period (i.e. whether data in AOF are too old
!   for current model time or not).
!   Determination of AOF time boxes to be read ('naoffbx', 'naoflbx').
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

! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE

! Local parameters: None
! ----------------
  INTEGER (KIND=iintegers) , PARAMETER :: &
    nfddrl = 14,         & ! AOF FDR (file description) record length
    nddrrl = 17            ! AOF DDR (data description) record length

! Local variables:
! ----------------
! LOGICAL, SAVE            ::  &
!   lobpfrs = .TRUE.       ! variable for 'first time of obs. processing'

  INTEGER (KIND=iintegers) ::  &
    nandat              ,& ! initial date of AOF
    nantim              ,& ! initial time of AOF
    nendat              ,& ! end date of AOF
    nentim              ,& ! end time of AOF
    iyy                 ,& ! year
    ienyy     ,imoyy    ,& ! year     / of the initial or final
    ienmm     ,imomm    ,& ! month   /  date and time
    iendd     ,imodd    ,& ! day    <   of the observations (AOF)
    ienhh     ,imohh    ,& ! hour    \  resp. of the initial
    ienmin    ,imomin   ,& ! minute   \ model date and time
    iensec    ,imosec   ,& ! seconds   \ model date
    imindif             ,& ! time difference [min] between these two dates
    ierrf                  ! error status for time difference routine

  INTEGER (KIND=iintegers) ::  &
    nfdrao (nfddrl)     ,& ! analysis observation file fdr
    nddrao (nddrrl)     ,& ! analysis observation file ddr
    icd                 ,& ! loop index
    nstat, ierr, jerr   ,& ! error status variables
    implcode               ! MPI status variable

  REAL (KIND=ireals)       ::  &
    t_last              ,& ! time of latest obs to be read at start time
    t_first                ! time of first  obs to be read at start time

  CHARACTER (LEN=20)       ::  &
    yroutine               ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerrmsg  , yerr        ! error message
!
!------------ End of header ----------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Begin Subroutine obs_aof_init
  !-----------------------------------------------------------------------------

  yroutine = 'obs_aof_init'

  !-----------------------------------------------------------------------------
  !  Section  1.   Opening of files used within obervation processing
  !-----------------------------------------------------------------------------

! IF (lobpfrs)     THEN
    IF (my_cart_id == 0) THEN

      OPEN (nuaof,FILE=yaofpath,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yaofpath FAILED'
        CALL model_abort (my_cart_id, 1001, yerrmsg, yroutine)
      ENDIF
      REWIND nuaof

      OPEN (nuaofex,FILE=yuaofex,FORM='FORMATTED',IOSTAT=nstat)
      IF (nstat /= 0) THEN
        yerrmsg = 'OPENING OF FILE yuaofex FAILED'
        CALL model_abort (my_cart_id, 1002, yerrmsg, yroutine)
      ENDIF
      REWIND nuaofex

! dedicated GPS input observation file is opened further below

    ENDIF

  !-----------------------------------------------------------------------------
  !  Section  2.   Initialize analysis observation file processing
  !-----------------------------------------------------------------------------

!   Initial date and time of the forecast run
    READ (ydate_ini,'(I4,5I2)') imoyy, imomm, imodd, imohh, imomin, imosec

    IF (my_cart_id == 0) THEN
!     Read AOF file description record
!     --------------------------------
      READ (nuaof, IOSTAT=nstat) (nfdrao(icd), icd=1,nfddrl)
      IF (nstat /= 0) THEN
        yerrmsg = 'READING FDR OF FILE aof FAILED'
        CALL model_abort (my_cart_id, 1006, yerrmsg, yroutine)
      ENDIF

!     Read AOF data description record
!     --------------------------------
      READ (nuaof, IOSTAT=nstat) (nddrao(icd), icd=1,nddrrl)
      IF (nstat /= 0) THEN
        yerrmsg = 'READING DDR OF FILE aof FAILED'
        CALL model_abort (my_cart_id, 1007, yerrmsg, yroutine)
      ENDIF
    ENDIF ! my_cart_id == 0

    IF (num_compute > 1) THEN
!     Broadcast AOF data description record to all PEs.
      CALL distribute_values ( nddrao, nddrrl, 0, imp_integers, icomm_cart     &
                             , implcode )
!     ======================
      jerr = ABS( implcode )
      CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
      IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
        WRITE( yerrmsg,'(" * Error in distribute_values *",2I5)' ) jerr, ierr
        CALL model_abort (my_cart_id, 1111, yerrmsg, yroutine)
      ENDIF
    ENDIF

!!   Initialize 'nmdi' for 64-bit machine: convert (2^30-1) to (2^59-1)
!!   ------------------------------------
!!    nmdi  =  nmdi+1
!!    nmdi  =  nmdi/2 *nmdi - 1
!   IF (my_cart_id == 0)  PRINT '(''    nmdi : '',o22,i22       )' , nmdi, nmdi


!   Extract initial and last AOF date and time and
!   maximum observation length from AOF ddr
!   -----------------------------------------
    nandat = nddrao ( ntiop ) / 100
    nantim = MOD( nddrao ( ntiop ) , 100 )
    iyy = nandat/10000
    IF (iyy < 1000 )    THEN
       IF (iyy < 70 ) THEN
         nandat = 20000000 + nandat
       ELSE
         nandat = 19000000 + nandat
       ENDIF
    ENDIF
    nendat = nddrao ( nteop ) / 100
    nentim = MOD( nddrao ( nteop ) , 100 )
    iyy = nendat/10000
    IF (iyy < 1000 )    THEN
       IF (iyy < 70 ) THEN
         nendat = 20000000 + nendat
       ELSE
         nendat = 19000000 + nendat
       ENDIF
    ENDIF

    nmxlob = nddrao ( nlmaxd )

    IF (lwonl) THEN
      WRITE(nupr,'(" INITIAL TERM OF FORECAST RUN:     ",A14)') ydate_ini
      WRITE(nupr,'(" INITIAL DATE OF OBSERVATION FILE: ",I10)') nandat
      WRITE(nupr,'(" INITIAL TIME:                     ",I10)') nantim
      WRITE(nupr,'(" LAST DATE OF OBSERVATION FILE: .  ",I10)') nendat
      WRITE(nupr,'(" LAST TIME:                        ",I10)') nentim
      WRITE(nupr,'(" MAXIMUM REPORT LENGTH:            ",I10)') nmxlob
    ENDIF

!   Check time difference between end of AOF observation period
!   and initial model time
!   -----------------------------------------------------------
    ienyy  = nendat/10000
    ienmm  = INT (nendat/100) - ienyy*100
    iendd  = MOD ( nendat, 100 )
    ienhh  = nentim
    ienmin = 0
    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                    &
                        ienyy, ienmm, iendd, ienhh, ienmin,                    &
                        itype_calendar, imindif, ierrf )
    IF (imindif < 0)  THEN
      IF (lwonl) THEN
        WRITE( nupr,'(" Time difference between end of AOF observation ",      &
                     &"period and initial model time: ",i10, "minutes")')      &
                       imindif
        WRITE( nupr,'(" OBSERVATION FILE TOO OLD!       A B O R T ")')
      ENDIF
      yerrmsg ='OBSERVATION FILE TOO OLD'
      CALL model_abort (my_cart_id, 1009, yerrmsg, yroutine)
    ENDIF


!   Determine initial time of AOF time box to be read now
!   -----------------------------------------------------
!   (time in forecast hour units)

!   time of first obs. allowed in the AOF, i.e initial time of first AOF timebox
    ienyy  = nandat/10000
    ienmm  = INT (nandat/100) - ienyy*100
    iendd  = MOD ( nandat, 100 )
    ienhh  = nantim
    ienmin = 0
    CALL diff_minutes ( imoyy, imomm, imodd, imohh, imomin,                    &
                        ienyy, ienmm, iendd, ienhh, ienmin,                    &
                        itype_calendar, imindif, ierrf )
    naoffbx = imindif / 60.0_ireals

!   initial time of last  AOF timebox influencing the forecast at the start time
    t_last  =   ntstep*dtdeh + tconbox/3600.0_ireals                           &
              + MAX( wtukrsa , tipolmx , wtukara , wtuksua , tipmxsu )
    naoflbx = naoffbx + naofbox *INT( (t_last -naoffbx) /naofbox )
!   naoflbx = INT( t_last-epsy ) /naofbox
!   naoflbx = (t_last+naofbox-epsy)/naofbox

!   initial time of first AOF timebox influencing the forecast at the start time
    t_first =   ntstep*dtdeh                                                   &
              - MAX( wtukrse , tipolmx , wtukare , wtuksue , tipmxsu )
    naoffbx = naoffbx + naofbox *INT( (t_first-naoffbx) /naofbox )

    IF (lwonl)                                                                 &
      WRITE (nupr,'(" LAST TIME BOX TO BE READ AT START TIME: ",F4.0)') naoflbx

    lastob = .FALSE.

! ENDIF ! lobpfrs == .true.

!-------------------------------------------------------------------------------
! End Subroutine obs_aof_init
!-------------------------------------------------------------------------------
END SUBROUTINE obs_aof_init


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for organizing the reading from AOF
!-------------------------------------------------------------------------------

SUBROUTINE obs_read_aof_org ( nsgoba )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" is the master routine 
!   for reading, extraction and processing of data from the AOF ('Analysis
!   Observation File').
!
! Method:
!   Reading data records from analysis observation file 'AOF' and 
!   distributing to the appropriate nodes.
!   Control of data extraction, processing and insertions in the 
!   observation data record 'ODR'.
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

! Modules used: These are declared in the module declaration section
! -------------

IMPLICIT NONE


! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT (INOUT)  :: &
    nsgoba             ! number of single-level reports before having
                       ! read new data at current timestep

! Local parameters: None
! ----------------

! Local variables:
! ----------------
  LOGICAL                  ::  &
    laddrep         ,& ! if true then add present report to ODR
    lvprof             ! .true. if number of levels > 1

  INTEGER (KIND=iintegers) ::  &
    nexceed(4)      ,& ! number of reports in excess of array size
    ipos            ,& ! start index - 1 of current observation in report buffer
    irplen          ,& ! report length
    istat, ierr, jerr  ! error status variables

  CHARACTER (LEN=70)       ::  &
    yerrmsl            ! error message
  CHARACTER (LEN=20)       ::  &
    yroutine           ! name of this subroutine
  CHARACTER (LEN=40)       ::  &
    yerr               ! error message

!
!------------ End of header ----------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Begin Subroutine obs_read_aof_org
  !-----------------------------------------------------------------------------

  yroutine = 'obs_read_aof_org'

  !-----------------------------------------------------------------------------
  !  Section  1.   Allocate and initialize buffers to read from AOF
  !-----------------------------------------------------------------------------

! Allocate space for report and output buffers and initialize arrays
! ------------------------------------------------------------------

! Find the maximum length of the buffer containing all reports for all PE's.
  nmxbln  = maxmlo
  IF ((laircf) .OR. (lverpas)) nmxbln  =  MIN( maxmlo , 150 )
! nmxbln  = NINT(   maxmlo * (nmxlob+nbufhed)                                  &
  nmxbln  = NINT(   nmxbln * (nmxlob+nbufhed)                                  &
                  * MAX( c1 , naofbox /MAX( wtukrsa+wtukrse, 2*tipolmx, c1 ) ) &
                 +  maxsgo * (52+nbufhed)                                      &
                  * MAX( c1 , naofbox /MAX( wtuksua+wtuksue, 2*tipmxsu, c1 ) ) )

  IF (ALLOCATED( aofbuf )) THEN
    PRINT '("CAUTION in src_obs_proc_aof: aofbuf is already allocated "        &
          &,"at time / PE",I6,I5)', ntstep, my_cart_id
    DEALLOCATE ( aofbuf , nodbuf , STAT=istat )
  ENDIF

  ALLOCATE ( aofbuf(nmxlob)             , STAT=istat )
  ALLOCATE ( nodbuf(nmxbln)             , STAT=istat )

  ierr = 0
  jerr = ABS( istat )
  IF (num_compute > 1)                                                         &
    CALL global_values( jerr, 1,'MAX',imp_integers,icomm_cart, -1, yerr,ierr )
  IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
    WRITE( yerrmsl,'(''could not allocate nodbuf:'',2I6                        &
                  &,'' ==> more PE, smaller "maxmlo"'')' ) jerr, ierr
    CALL model_abort (my_cart_id, 1017, yerrmsl, yroutine)
  ENDIF

  aofbuf (:)   = 0
  nodbuf (:)   = 0
  nexceed(:)   = 0

  !-----------------------------------------------------------------------------
  !  Section  2.   Read and select reports of the next timebox,
  !                distribute reports to appropriate nodes
  !-----------------------------------------------------------------------------

  Get_reports_of_next_aof_timebox:   DO
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Write info on AOF time boxes for control
    IF (lwonl)                                                                 &
      WRITE(nupr,'("OBS_READ_AOF_ORG: CURRENT AOF TIME BOX TO BE READ:",F4.0,/,&
                   &15X,"LENGTH OF AOF TIME BOX: ",I3)' )       naoffbx, naofbox

    CALL obs_read_distribute
!   ========================
 
  !-----------------------------------------------------------------------------
  !  Section  3.   Observation Processing
  !-----------------------------------------------------------------------------
 
! start index - 1 of first observation in report buffer
    ipos = 0

    Get_next_report:   DO
!   ~~~~~~~~~~~~~~~~~~~~~

! get length of next report
      IF (ipos +1 > nmxbln)          EXIT Get_next_report
      irplen = nodbuf (ipos+1)
      IF (irplen == 0)               EXIT Get_next_report

  !-----------------------------------------------------------------------------
  !  Section  3.1  Fill observation header arrays
  !-----------------------------------------------------------------------------

      CALL obs_save_head ( ipos , laddrep , nexceed )
!     ==================

      ipos   = ipos + irplen

      IF (.NOT.laddrep)  CYCLE Get_next_report
 
  !-----------------------------------------------------------------------------
  !  Section  3.2  Fill observation body arrays
  !-----------------------------------------------------------------------------
 
! Multi-level obs. report (TEMP/PILOT)

      lvprof = .FALSE.
 
      IF (      lmulti(nobtyp)) CALL obs_multi_level ( lvprof )
!                               ====================
 
! Single-level obs. report (SYNOP/AIREP/DRIBU)
 
      IF (.NOT. lmulti(nobtyp)) CALL obs_single_level
!                               =====================
 
  !-----------------------------------------------------------------------------
  !  Section  3.3  Redundancy check (removal of redundant information)
  !-----------------------------------------------------------------------------

      CALL obs_redundancy ( lvprof , nsgoba )
!     ===================
 
! Close loop : fetch next report if necessary
! -------------------------------------------
    ENDDO Get_next_report 
!   ~~~~~~~~~~~~~~~~~~~~~

    naoffbx = naoffbx + naofbox

    IF (naoffbx > naoflbx)   EXIT Get_reports_of_next_aof_timebox

  ENDDO Get_reports_of_next_aof_timebox
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !-----------------------------------------------------------------------------
  !  Section  4.   Deallocate observation buffer, write caution message
  !-----------------------------------------------------------------------------

  DEALLOCATE ( aofbuf    , STAT = istat )
  DEALLOCATE ( nodbuf    , STAT = istat )


! number of obs not used due to insufficient ODR size
! ---------------------------------------------------

  IF (num_compute > 1) THEN
    CALL global_values ( nexceed, 4, 'MAX',imp_integers,icomm_cart, 0,yerr,ierr)
!   ==================
    IF (ierr /= 0)  CALL model_abort (my_cart_id, 11015, yerr, yroutine)
!                   ----------------
  ENDIF

  IF ((my_cart_id == 0) .AND. (MAXVAL(nexceed(1:3)) > 0)) THEN
    OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'            &
                               , POSITION='APPEND', IOSTAT=istat)
    IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
    IF (istat /= 0) CALL model_abort (my_cart_id, 1409, yerr, yroutine)
  ENDIF
  IF ((my_cart_id == 0) .AND. (nexceed(1) > 0)) THEN
    WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL SINGLE-LEVEL OBS."   &
                   &," BEYOND maxsgl ",I5)' ) ntstep, nexceed(1), maxsgl
    nexceed(1)  =  nexceed(1) * maxsgo / maxsgl + 1
    WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST"   &
                   &,I6)' ) nexceed(1)
  ENDIF
  IF ((my_cart_id == 0) .AND. (nexceed(2) > 0)) THEN
    WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," LOCAL MULTI-LEVEL OBS."    &
                   &," BEYOND maxmll ",I5)' ) ntstep, nexceed(2), maxmll
    nexceed(2)  =  nexceed(2) * maxmlo / maxmll + 1
    WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST"   &
                   &,I6)' ) nexceed(2)
  ENDIF
  IF ((my_cart_id == 0) .AND. (nexceed(3) > 0)) THEN
    WRITE( nucautn,'("CAUTION !!!!! t=",I5,":",I5," OBS WITH UNKNOWN OBS TY"   &
                  &,"PE, E.Q. TYPE ",I5)' ) ntstep, nexceed(3), nexceed(4)
    WRITE( nucautn,'("   ==>  CHECK   AOF - FILE ")' )
  ENDIF
  IF ((my_cart_id == 0) .AND. (MAXVAL(nexceed(1:3)) > 0))  CLOSE (nucautn)

!-------------------------------------------------------------------------------
! End Subroutine obs_read_aof_org
!-------------------------------------------------------------------------------
END SUBROUTINE obs_read_aof_org


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for variable and flag extraction
!-------------------------------------------------------------------------------

SUBROUTINE obs_extract ( k, k1, kdata, klen )

!-------------------------------------------------------------------------------
! Description:
!   Extract variable and flags from analysis observation file buffer
!
! Method:
!   Simple extraction and unpacking
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date        Comment
! -------   ----        -------
!   1.0     09.12.97    Original code. 
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
!-------------------------------------------------------------------------------
! Modules used:       These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------
! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT (IN)  ::  &
    k,           & ! variable level pointer
    k1,          & ! variable number
    klen,        & ! report length
    kdata(klen)    ! report

! Local parameters:
! ----------------
  INTEGER (KIND=iintegers) , PARAMETER    ::  &
    i_6  =  63     ! bit 0-5 set to 1

! Local scalars:
! -------------
  INTEGER (KIND=iintegers) ::  &
    ibits ,      & ! statement function to unpack any bit pattern
    inval ,      & ! word to be unpacked
    ipos  ,      & ! bit position of bit struct. to be unpacked
    icovr ,      & ! no. of bit occ. by the bit structure
    iovrid,      & ! override bit
    ibit6 ,      & ! 6bit flag pattern
    invar          ! flag word
!
!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (inval,ipos,icovr)=IAND (ISHFT(inval,-ipos),icovr)

!-------------------------------------------------------------------------------
! Begin Subroutine obs_extract
!-------------------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Section 1.0: Variable and flag word extraction
  !----------------------------------------------------------------------

  IF (nwrvr(k1,nobtyp) > 0) THEN 
    nvarib   =    kdata(k+nwrvr (k1,nobtyp))
  ELSE
    nvarib   =    nmdi
  ENDIF
  ibit6    =    n6vfb(k1,nobtyp)
  invar    =    kdata(k+nwrvrf(k1,nobtyp))
  nvrfwr   =    IAND (ishft (invar,-(ibit6-1)*6), i_6)
  nvrfwr   =    IAND (nvrfwr, nibits(nf1aoc))

  !----------------------------------------------------------------------
  ! Section 2.0:  Flag extraction
  !----------------------------------------------------------------------

  iovrid   =    ibits (nvrfwr,nf1bps( 3),nibits (nf1boc( 3)))
  IF (iovrid == 1)     THEN
     nflag =    0
  ELSE
     nflag =    ibits (nvrfwr,nf1bps( 4),nibits (nf1boc( 4)))
  ENDIF
  
!-----------------------------------------------------------------------
 
END SUBROUTINE obs_extract

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for vertical sorting of reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_sort_levels (kvar, io, jo, zsurf)

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module src_obs_proc_aof sorts multi-level
!   reports by pressure
!
! Method:
!   TEMP/PILOT reports are sorted in vertical in decreasing/
!   increasing order of pressures.
!
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
!
!-------------------------------------------------------------------------------
! Modules used:       These are declared in the module declaration section
!
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
   kvar,        & ! variable which is being used as key for sorting
   io, jo         ! horiz. coordinate of grid point used for computation
                  ! of pressure.

  REAL    (KIND=ireals)    , INTENT (IN)   :: &
   zsurf          ! model surface topography height at (io,jo)

! Local parameters: none
! -----------------

! Local variables:
! ----------------
  INTEGER (KIND=iintegers)                :: &
   iskip , jbuf,        & ! loop indices
   jloop1, jloop2,      &
   istart, iend,        &
   istr  ,              &
   jfill ,              &
   nppp2 ,              & ! integer value of pressure
   ibuff (8 )             ! buffer comtaining data of one level

  REAL (KIND=ireals)                      :: &
   rzzz,                & ! observed height
   rppp                   ! model pressure

  LOGICAL                     :: lcomp

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_sort_levels
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Sort levels by pressure (decreasing order)
  !-----------------------------------------------------------------------------

! Loop over report
!
  istart       =   nhdob + nlev0 (nobtyp)
  iskip        =   nlevn (nobtyp)
  iend         =   aofbuf(nrecln) - nlevn (nobtyp)

  Loop1  : DO    jloop1   =   istart,iend ,iskip

! Check sorting key
 
  IF ((kvar /= nvar(1)) .AND. (kvar /= nvar(7)))                     CYCLE Loop1
  IF (nwrvr(nvar(1),nobtyp) == 0)                                    CYCLE Loop1
  IF ((kvar == nvar(7)) .AND. (nwrvr(nvar(7),nobtyp) == 0))          CYCLE Loop1
 
! Store the level into the auxiliary buffer

  DO      jbuf     =   1,iskip
    ibuff (jbuf) =   aofbuf(jbuf+jloop1)
  ENDDO          
 
! If obs. pressure is missing, get model pressure
 
  IF (kvar == nvar( 7)) THEN
     nvarib =  ibuff (nwrvr (nvar( 1),nobtyp))
     IF (nvarib == nmdi) THEN
        nvarib =  ibuff (nwrvr (nvar( 7),nobtyp))
        IF (nvarib == nmdi)                                          CYCLE Loop1
        rzzz   =  ( nvarib * rfact (nvar(7)) - nalmcn*1._ireals)
        IF (rzzz <  zsurf )                                          CYCLE Loop1
        rppp   =  z2p ( rzzz, io, jo )
        IF (rppp <= rmdich)                                          CYCLE Loop1
        ibuff (nwrvr(nvar(1),nobtyp)) = nint( rppp /rfact(nvar(1)) )
     ENDIF
  ENDIF

! Look for lower/higher pressure
 
  istr         =   jloop1 + iskip
  Loop2 : DO   jloop2   =   istr  ,iend ,iskip

  IF (kvar == nvar  ( 1))                    THEN
    lcomp  =  ibuff (nwrvr(kvar,nobtyp))  <  aofbuf (jloop2+nwrvr(kvar,nobtyp))
  ELSE
    IF (aofbuf (jloop2+nwrvr (nvar( 1) ,nobtyp)) == nmdi) THEN
      nvarib =   aofbuf (jloop2+nwrvr (nvar ( 7) ,nobtyp))
      IF (nvarib == nmdi)                                            CYCLE Loop2
      rzzz   =  ( nvarib * rfact (nvar(7)) - nalmcn*1._ireals )
      IF (rzzz <  zsurf ) THEN
        lcomp  = .TRUE.
      ELSE
        rppp   =  z2p ( rzzz, io, jo )
        IF (rppp <= rmdich) THEN
          lcomp = .FALSE.
        ELSE
          nppp2  =   nint( rppp /rfact(nvar( 1)) )
          lcomp  =   ibuff (nwrvr (nvar(1),nobtyp)) <  nppp2
        ENDIF
      ENDIF
    ELSE
     lcomp = ibuff(nwrvr(nvar(1),nobtyp)) < aofbuf(jloop2+nwrvr(nvar(1),nobtyp))
    ENDIF
  ENDIF

  IF (lcomp) THEN
! Lower pressure found (swap levels)
     DO    jfill    =   1,iskip
        ibuff (       jfill)  =  aofbuf(jloop2+jfill)
     ENDDO         
     DO    jfill    =   1,iskip
        aofbuf(jloop2+jfill)  =  aofbuf(jloop1+jfill)
     ENDDO
     DO    jfill    =   1,iskip
        aofbuf(jloop1+jfill)  =  ibuff (       jfill)
     ENDDO
  ENDIF

  ENDDO Loop2  
  ENDDO Loop1 

!-------------------------------------------------------------------------------
 
END SUBROUTINE obs_sort_levels

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for flag word pointers for add. groups
!-------------------------------------------------------------------------------

SUBROUTINE obs_flag6b ( ktotal, kordgr, kflgwr, koc6b )

!-------------------------------------------------------------------------------
! Description:
!   Determine flag word pointers for additional groups
!
! Method:
!   Simple calculation according to AOF format.
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       15.12.97   Original code.    Michael Buchhold
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
INTEGER   (KIND=iintegers), INTENT(IN)   ::                           &
  ktotal ,    & ! total number of add. groups
  kordgr        ! pointer for actual add. group

INTEGER   (KIND=iintegers), INTENT(OUT)  ::                           &
  kflgwr ,    & ! pointer for flag word
  koc6b         ! pointer for 6 bit pattern within flag word

!
!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_flag6b
!-------------------------------------------------------------------------------

!                Find out pointers (flag word and 6 bit pattern)
!                ------------------------------------------------
 
  kflgwr  =  ktotal + INT(kordgr/4._ireals *(1._ireals+epsy) +0.75_ireals )
  koc6b   =  kordgr - INT(kordgr/4._ireals *(1._ireals+epsy) -0.25_ireals )*4

!-------------------------------------------------------------------------------
 
END SUBROUTINE obs_flag6b

!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for extraction of optional group data
!-------------------------------------------------------------------------------

SUBROUTINE obs_option_groups

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" performs data extraction
!   from SYNOP optional group.
!
! Method:
!   Additional groups which are being extracted are:
!   1. weather group                  (SYNOP)
!   2. general cloud group            (SYNOP)
!   3. additional cloud groups        (SYNOP)
!   4. ground group                   (SYNOP)
!   5. special phenomena group        (SYNOP)
!   6. ice group                      (SYNOP)
!   7. rain group                     (SYNOP)
!   8. ship group                     (SYNOP SHIP)
!   9. waves group                    (SYNOP SHIP)
!  10. extreme temperature group      (SYNOP)
!  11. radiation group                (SYNOP)
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       15.12.97   Original code.    Michael Buchhold
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

! Subroutine arguments: none
! --------------------

! Local parameters:
! ----------------
  INTEGER (KIND=iintegers) , PARAMETER  :: &
    i_6  =  63     ! bit 0-5 set to 1

! Local scalars
! -------------
  INTEGER (KIND=iintegers) ::  &
    i6bit ,      & ! statement function
    invar ,      & ! word to be unpacked
    ibits ,      & ! statement function to unpack any bit pattern
    ipos  ,      & ! bit position of bit struct. to be unpacked
    icovr ,      & ! no. of bit occ. by the bit structure
    ibit6 ,      & ! bit position in 6 bit structure
    iflgwr,      & ! pointer for flag word
    ioc6b ,      & ! pointer for 6 bit pattern within flag word
    jadd, jacl,  & ! loop indices
    jspp, jwav,  & !     "
    jmiss,       & !     "
    itotal,      & ! total number of additional groups
    icrtot,      & ! total number of additional cloud groups
    iaclin,      & ! additional cloud groups indicator
    igrgrp,      & ! ground group indicator
    isppin,      & ! special phenomena indicator
    iicegr,      & ! ice group indicator
    iraing,      & ! rain group indicator
    ishipg,      & ! ship group indicator
    iwavgr,      & ! wave group indicator
    iextgr,      & ! extreme temperature group indicator
    iradgr         ! radiation group indicator
!
!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ipos,icovr)=IAND (ISHFT(invar,-ipos),icovr)

! Statement function to unpack 6 bits from flag word
  i6bit(invar,ibit6)  = IAND (ishft(invar,-(ibit6-1)*6), i_6)

!-------------------------------------------------------------------------------
! Begin Subroutine obs_option_groups
!-------------------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Section 1.0: Synop weather,cloud and additional groups
  !----------------------------------------------------------------------

!     1.1   Weather group
 
  nwwgr   =   nmdi
  IF (      (nobtyp == nsynop)                                                 &
      .AND. ((nwwex (1) /= 0) .OR. (nwwex (2) /= 0) .OR. (nwwex (3) /= 0))) THEN
     CALL obs_extract(nhdob,nvar  ( 9),aofbuf,nmxlob)
     nwwgr   =   nvarib
     nwwfw   =   nvrfwr
  ENDIF

!     1.2   General cloud group
 
  ngclgr  =   nmdi
  IF (      (nobtyp == nsynop)                                                 &
      .AND. (     ngclex(1) /= 0 .OR. ngclex(2) /= 0 .OR. ngclex(3) /= 0       &
             .OR. ngclex(4) /= 0 .OR. ngclex(5) /= 0 .OR. ngclex(6) /= 0)) THEN
!!   CALL obs_extract(nhdob,nvar  (10),aofbuf,nmxlob)
!!   ngclgr  =   nvarib
!!   ngclfw  =   nvrfwr

    IF (nwrvr(nvar(10),nobtyp) > 0) THEN
      ngclgr   =    aofbuf(nhdob+nwrvr (nvar(10),nobtyp))
    ELSE
      ngclgr   =    nmdi
    ENDIF
    ibit6    =    n6vfb(nvar(10),nobtyp)
    invar    =    aofbuf(nhdob+nwrvrf(nvar(10),nobtyp))
    ngclfw   =    IAND (ishft (invar,-(ibit6-1)*6), 4095)
  ENDIF
 
!     1.3   Additonal groups (if any)

  jadd   =   1
  IF  (nobtyp == nsynop)  jadd  =  naddgr
  IF ((nobtyp == nsynop) .AND. (aofbuf (nhdob+jadd) /= 0)) THEN
 
!      1.3.1 Preset counter(s)
 
    itotal   =   0
    icrtot   =   0
 
!      1.3.2 Find out total number of additional groups
 
    DO       jadd   =    1    ,   15
      itotal   =   itotal + ibits (aofbuf (nhdob +naddgr),jadd-1, nibits(1))
    ENDDO
 
!      1.3.3 Additonal cloud groups
 
    iaclin     =  ibits (aofbuf (nhdob+naddgr),naclbp,nibits(nacloc))
    DO      jacl     =    1     ,    4
      naclgr(jacl)  =  nmdi
      IF (ibits (iaclin,jacl-1,nibits( 1)) == 0)                           CYCLE
      icrtot        =  icrtot + 1
      IF (naclex(1) == 0 .AND. naclex(2) == 0 .AND. naclex(3) == 0)        CYCLE
      naclgr(jacl)  =  aofbuf (nhdob +naddgr+icrtot)
      CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
      naclfw(jacl)  =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
    ENDDO


!      1.3.4 Ground group
 
    ngrngr     =  nmdi
    igrgrp     =  ibits (aofbuf (nhdob+naddgr),ngrnbp,nibits(ngrnoc))
    IF (igrgrp /= 0) THEN      
      icrtot     =  icrtot + 1
      IF (ngrnex(1) /= 0 .OR. ngrnex(2) /= 0 .OR. ngrnex(3) /= 0) THEN
        ngrngr     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        ngrnfw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF

!      1.3.5 Special phenomena group
 
    isppin     =  ibits (aofbuf (nhdob+naddgr),nsppbp,nibits(nsppoc))
    DO      jspp     =    1     ,    2
      nsppgr(jspp)  =  nmdi
      IF (ibits (isppin,jspp-1,nibits( 1)) == 0)                           CYCLE
      icrtot        =  icrtot + 1
      IF (nsppex(1) == 0 .AND. nsppex(2) == 0)                             CYCLE
      nsppgr(jspp)  =  aofbuf (nhdob +naddgr+icrtot)
      CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
      nsppfw(jspp)  =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
    ENDDO
 
!      1.3.6 Ice group
 
    nicegr     =  nmdi
    iicegr     =  ibits (aofbuf (nhdob+naddgr),nicebp,nibits(niceoc))
    IF (iicegr /= 0) THEN
      icrtot     =  icrtot + 1
      IF (niceex(1) /= 0 .OR. niceex(2) /= 0 .OR. niceex(3) /= 0) THEN
        nicegr     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        nicefw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF
  
!      1.3.7 Rain group
 
    nraing     =  nmdi
    iraing     =  ibits (aofbuf (nhdob+naddgr),nrinbp,nibits(nrinoc))
    IF (iraing /= 0) THEN
      icrtot     =  icrtot + 1
      IF (nraine(1) /= 0 .OR. nraine(2) /= 0 )         THEN
        nraing     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        nranfw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF

!      1.3.8 Ship group
 
    nshpgr     =  nmdi
    ishipg     =  ibits (aofbuf (nhdob+naddgr),nshpbp,nibits(nshpoc))
    IF (ishipg /= 0) THEN
      icrtot     =  icrtot + 1
      IF (nshipe(1) /= 0 .OR. nshipe(2) /= 0) THEN
        nshpgr     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        nshpfw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF
 
!      1.3.9 Waves group
 
    iwavgr      =  ibits (aofbuf (nhdob+naddgr),nwavbp,nibits(nwavoc))
    DO       jwav    =     1     ,     3
      nwavgr(jwav)=  nmdi
      IF (ibits (iwavgr,jwav-1,nibits( 1)) == 0)                           CYCLE
      icrtot      =  icrtot + 1
      IF (nwavee(1) == 0 .AND. nwavee(2) == 0 .AND. nwavee(3) == 0)        CYCLE
      nwavgr(jwav)=  aofbuf (nhdob +naddgr+icrtot)
      CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
      nwavfw(jwav)=  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
    ENDDO
  
!      1.3.10 Extreme temperature group
 
    nextgr     =  nmdi
    iextgr     =  ibits (aofbuf (nhdob+naddgr),nextbp,nibits(nextoc))
    IF (iextgr /= 0) THEN
      icrtot     =  icrtot + 1
      IF (nextex(1) /= 0 .OR. nextex(2) /= 0) THEN
        nextgr     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        nextfw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF

!      1.3.11 Radiation group
 
    nradgr     =  nmdi
    iradgr     =  ibits (aofbuf (nhdob+naddgr),nradbp,nibits(nradoc))
    IF (iradgr /= 0) THEN
      icrtot     =  icrtot + 1
      IF (nradex(1) /= 0) THEN
        nradgr     =   aofbuf (nhdob +naddgr+icrtot)
        CALL obs_flag6b (itotal,icrtot,iflgwr,ioc6b)
        nradfw     =  i6bit(aofbuf (nhdob +naddgr+iflgwr),ioc6b)
      ENDIF
    ENDIF

  ELSE
 
!       1.4   Set optional groups to missing values

!       1.4.1 Additional cloud groups
    DO       jmiss   =    1    ,   4
      naclgr( jmiss)   =   nmdi
    ENDDO
 
!       1.4.2 ground group
    ngrngr           =   nmdi
 
!       1.4.3 special phenomena groups
    nsppgr(     1)   =   nmdi
    nsppgr(     2)   =   nmdi
 
!       1.4.4 ice group
    nicegr           =   nmdi
 
!       1.4.5 rain group
    nraing           =   nmdi
 
!       1.4.6 ship group
    nshpgr           =   nmdi

!       1.4.7 waves groups
    DO      jmiss   =    1    ,   3
      nwavgr( jmiss)   =   nmdi
    ENDDO

!       1.4.8 extreme temperature group
    nextgr           =   nmdi

!       1.4.9 radiation group
    nradgr           =   nmdi
  ENDIF
   
!-------------------------------------------------------------------------------
! End of module procedure obs_option_groups
!-------------------------------------------------------------------------------
END SUBROUTINE obs_option_groups


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for data assimilation
!-------------------------------------------------------------------------------

SUBROUTINE obs_RASS_tv2t ( io, jo )

!-------------------------------------------------------------------------------
! Description:
!   Convert SODAR/RASS virtual temperature to air temperature
!   by the aid of model humidity
!
! Method:
!   Find nearest model levels for the observed height and interpolate
!   model humidity. Convert temperatue.
!
!
! Current Code Owner:  Michael Buchhold
!
! S-story:
! Version   Date       Comment
! -------   ----       -------
! 1.0       17.05.01   Original code.    Michael Buchhold
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: European Standards for Writing and Documenting
!                       Exchangeable Fortran 90 Code + DWD Extensions.
!
!-------------------------------------------------------------------------------
! Modules used:       These are declared in the module declaration section
!
!-------------------------------------------------------------------------------
IMPLICIT NONE

!-------------------------------------------------------------------------------
!
! Declarations:
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------
  INTEGER (KIND=iintegers) , INTENT(IN)   :: &
   io, jo         ! horiz. coordinate of grid point

! Local parameters: none
! -----------------

! Local variables:
! ----------------
  INTEGER (KIND=iintegers)                :: &
   iskip , jbuf,        & ! loop indices
   jloop1,              &
   istart, iend,        &
   ibuff (8 )             ! buffer containing data of one level

  REAL (KIND=ireals)                      :: &
   zzz,                 & ! height
   ztv,                 & ! virtual temperature
   ztt                    ! temperature

!------------ End of header ----------------------------------------------------

!-------------------------------------------------------------------------------
! Begin Subroutine obs_RASS_tv2t
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Convert observed virtual temperature of SODAR/RASS profile to air
  ! temperature by the aid of model humidity
  !-----------------------------------------------------------------------------

! Loop over report
!
  istart       =   nhdob + nlev0 (nobtyp)
  iskip        =   nlevn (nobtyp)
  iend         =   aofbuf(nrecln) - nlevn (nobtyp)

  Loop1  : DO    jloop1   =   istart,iend ,iskip

! Store the level into the auxiliary buffer

    DO      jbuf     =   1,iskip
      ibuff (jbuf) =   aofbuf(jbuf+jloop1)
    ENDDO

! get observed height and temperature
! -----------------------------------

! get observed height
    nvarib =  ibuff (nwrvr (nvar( 7),nobtyp))
    IF (nvarib == nmdi)  THEN
      PRINT '(" jloop1=",i3," SODAR/RASS height missing!")', jloop1
                                                                     CYCLE Loop1
    ENDIF
    zzz   =  (nvarib *rfact(nvar(7)) - nalmcn*c1)

! get observed virtual temperature
    nvarib = ibuff (nwrvr (nvar( 4),nobtyp))
    IF (nvarib == nmdi)  THEN
      PRINT '(" jloop1=",i3," SODAR/RASS Temp. missing ")', jloop1
                                                                     CYCLE Loop1
    ENDIF
    ztv   =   nvarib *rfact(nvar(4))

! calculate temperature
! ---------------------
    ztt = tv2t ( io , jo , ztv , zzz )
!         ====

    ibuff (nwrvr(nvar(4),nobtyp)) = NINT( ztt /rfact(nvar(4)) )

    DO      jbuf     =   1,iskip
      aofbuf(jbuf+jloop1) = ibuff (jbuf)
    ENDDO

  ENDDO Loop1

!-------------------------------------------------------------------------------
! End of module procedure obs_RASS_tv2t
!-------------------------------------------------------------------------------
END SUBROUTINE obs_RASS_tv2t


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for reading and distributing reports
!-------------------------------------------------------------------------------

SUBROUTINE obs_read_distribute 

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" organizes the reading,
!   pre-selection and distribution of observation reports.
!
! Method:
!   Reading data records from analysis observation file 'AOF',    
!   extraction of header information, selection and sending of
!   observations to the appropriate nodes.
!   28.07.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version  Date       Name
! -------- ---------- ----
! 1.13     1998/10/22 Michael Buchhold  : Initial release
! 1.19     1998/12/11 Christoph Schraff : Introduction of a verification mode
!                                         (reports set passive instead of being
!                                          rejected). Updated report events.
! 1.20     1999/01/07 Guenther Doms     : Renaming of some global variables.
! 1.27     1999/03/29 Christoph Schraff : Report flag word introduced. Revised
!                                         (32-bit) ODR missing data indicator.
! 1.28     1999/04/19 Christoph Schraff : Checking century of observations 
!                                         (AOF format revised due to yy 2000).
! 1.29     1999/05/11 Ulrich Schaettler : Adapted interfaces to utility-modules
!                                         and prepared use of MPE_IO
! 1.31     1999/07/01 Christoph Schraff : MPI calls replaced by parallel_utility
!                                         calls. Passive reports rejected in
!                                         verification mode, if requested.
! 1.36     2000/02/24 Michael Buchhold  : Introduction of ACAR aircraft reports. 
!                                         Adaption to 32 bit AOF word length.
! 1.39     2000/05/03 Ulrich Schaettler : Changed variable names for lat/lon.
! 2.4      2001/01/29 Christoph Schraff : Error message, if 'maxbln' too small.
! 2.11     2001/09/28 Christoph Schraff : Modification for porting to IBM.
! 2.13     2002/01/18 Michael Buchhold  : Introduction of Wind Profiler/RASS.
! 2.14     2002/02/15 Christoph Schraff : AIREP with bad station ID set passive.
! 2.19     2002/10/24 Michael Buchhold  : Introduction of Radar VAD wind report.
! 3.6      2003/12/11 Christoph Schraff : Introduction of obs type SATOB.
!                                         Ensuring error message at model_abort.
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

  LOGICAL                  ::  &
   lsfcob ,              & ! false, if no surface data are to be used
                           ! from a temp or pilot
   lseaobs,              & ! true, if observation taken on sea
   lindenb                 ! true, if sta Lindenberg: special grid pt assignment

  INTEGER (KIND=iintegers) ::  &
   ianyy     ,imoyy     ,& ! year     / of the date and time
   ianmm     ,imomm     ,& ! month   /  of the current
   iandd     ,imodd     ,& ! day    <   observation
   ianhh     ,imohh     ,& ! hour    \  resp. of the initial
   ianmin    ,imomin    ,& ! minute   \ model date and time
   iansec    ,imosec    ,& ! seconds   \ model date
   imindif              ,& ! time difference [min] between these two times
   ierrf                   ! error status for time difference routine

  INTEGER (KIND=iintegers) ::  &
   istat, irm, ierr,jerr,& ! error status variables
   izmplcode,            & ! for MPI error code
   izlocstat,            & ! for space allocation error
   npe,                  & ! rank of the receiving PE
   nobbuf,               & ! number of reports in report buffers
   ibuflen,iallen,       & ! total length of report buffers
   inodcnt,              & ! length of local report buffer
   inodlen(num_compute), & ! total length of report buffer at each node
   inodpos(num_compute), & ! displacement of local report buffers in 'iallbuf'
   ipos,                 & ! start address - 1 of an observation in 
                           ! report buffer
   irplen,               & ! report length
   nexceed(2),           & ! number of reports in excess of array size
   iovrid,               & ! override report flag
   ilen,                 & ! length of AOF report
   ibits ,               & ! statement function to unpack any bit pattern
   invar ,               & ! word to be unpacked
   ibpos ,               & ! bit position of bit struct. to be unpacked
   icovr ,               & ! no. of bit occ. by the bit structure
   insert,               & ! statement function to insert any bit pattern
   inval ,               & ! bit pattern to be inserted
   ibit  ,               & ! bit position
   ist1  , ist2 ,        & ! integer equivalent of station ID (part 1 and 2)
   is    ,               & ! position of a character in ASCII collating seq.
   iobdat,               & ! observation date as yymmdd or yyyymodd
   nelevt, nelevtl,      & ! station altitude
   istcon,               & ! station confidence flag
   iob_tot, job_tot,     & ! grid points assigned to obs (total model area)
   isurfob,              & ! = 1 if lsfcob=.t., else 0
   mrepflg,              & ! report flag word as part of station charcteristics
   mpassiv,              & ! passive flag assigned to report
   irnode,               & ! process an observation belongs to
   i,j,jflg,nc,          & ! loop indices
   inode,inobs,          & ! loop indices
   istr,iend,iend2,      & ! loop indices
   iskip                   !
 
  REAL (KIND=ireals)       ::  &
   timlim,               & ! limit for time difference in forecast hour units
   zsurf,                & ! height of model grid pt. to which obs. is assigned
   zio_tot,zjo_tot,      & !  observation position in g.p. units (total area)
   rlon, rlat              ! latitude and longitude in the rotated system
 
  CHARACTER (LEN=ilstidp)  ::  &
   ystid                   ! station identity
  CHARACTER (LEN=12)       :: &
   yobdat                  ! date and time of current observation
  CHARACTER (LEN=20)       :: &
   yroutine                ! name of this subroutine
  CHARACTER (LEN=30)       :: &
   yerrmsg , yerr          ! error message
  CHARACTER (LEN=70)       :: &
   yerrmsl                 ! error message
 
! Local arrays:
! ------------

  INTEGER (KIND=iintegers), ALLOCATABLE :: &
   iallbuf(:)           ,& ! buffer containing observations from all PEs
   iposbuf(:)           ,& ! position of the reports in the report buffer
   inodbuf(:)              ! node the reports belong to
!
!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ibpos,icovr)=IAND (ISHFT(invar,-ibpos),icovr)

! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar , ishft(inval,ibit))

!-------------------------------------------------------------------------------
! Begin Subroutine obs_read_distribute
!-------------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! Section 1.0: Initialization of organizational variables
  !-----------------------------------------------------------------------------

  yroutine = 'obs_read_distribute'
  ibuflen      = 0
  nodbuf  (:)  = 0
  inodlen (:)  = 0
  npe          = 0
  nexceed (:)  = 0
  IF (num_compute > 1) THEN
    ALLOCATE ( inodbuf (nmxbln/(20+nbufhed)+1) , STAT=izlocstat )
    ALLOCATE ( iposbuf (nmxbln/(20+nbufhed)+1) , STAT=izlocstat )
    ALLOCATE ( iallbuf (nmxbln)                , STAT=izlocstat )
! (20+nbufhed is the minimum length of a report in the buffer 'iallbuf')
    jerr = ABS( izlocstat )

    CALL global_values ( jerr, 1,'MAX', imp_integers,icomm_cart, -1, yerr,ierr )
    IF ((jerr /= 0) .OR. (ierr /= 0)) THEN
      WRITE( yerrmsl,'(''could not allocate nodbuf:'',2I6                      &
                     &,'' ==> more PE, smaller "maxmlo"'')' ) jerr, ierr
      CALL model_abort (my_cart_id, 11007, yerrmsl, yroutine)
    ENDIF
    iallbuf (:)  = 0
    inodbuf (:)  = 0
    iposbuf (:)  = 0
    nobbuf       = 0
  ENDIF

! Initial date and time of the forecast run
  READ (ydate_ini,'(I4,5I2)') imoyy, imomm, imodd, imohh, imomin, imosec

   !---------------------------------------------------------------------------
   ! Section 1.1  Compose full fields of model surface height and land fraction
   !---------------------------------------------------------------------------

  ALLOCATE ( hsurf_tot (ie_tot,je_tot)      , STAT=izlocstat )
  ALLOCATE ( fland_tot (ie_tot,je_tot)      , STAT=izlocstat )

  CALL get_global_surf_aof
! ========================

  !-----------------------------------------------------------------------------
  ! Section 2.    Fetch reports, pre-select and distribute them to the
  !               appropriate nodes
  !-----------------------------------------------------------------------------

  IF ((my_cart_id == 0) .AND. (.NOT. lastob)) THEN
! ================================================

    taof  =  naoffbx - c1

    read_obs:    DO
!   ~~~~~~~~~~~~~~~

!   Check if observation belongs to next timebox
!   --------------------------------------------
    IF (taof >= naoffbx+naofbox-epsy)                              EXIT read_obs

!   read report from AOF
!   --------------------

    READ (UNIT=nuaof, IOSTAT=istat) ilen, (aofbuf(i), i=2,ilen)
    aofbuf(1) =  ilen
    IF (istat < 0)      THEN
!     End of file
      lastob=.TRUE.
      PRINT       '(" LASTOB: LAST OBSERVATION IS READ FROM AOF")'
      EXIT read_obs
    ELSEIF (istat > 0)  THEN
!     I/O error
      WRITE( nurej,'(" I/O ERROR DETECTED ON ANALYSIS OBSERVATION FILE")')
      yerrmsg = 'I/O ERROR DETECTED ON AOF'
      CALL model_abort (my_cart_id, 11001, yerrmsg, yroutine)
    ENDIF

!   Print AOF report
!   ----------------
    IF (lpraof) THEN
!     check if obs./code type to be printed and location inside print region
      IF (     (noctrq == aofbuf(nobtp))                                       &
          .OR. (noctrq == aofbuf(ncdtp)) .OR. noctrq == 9)  THEN
        roblat   =  (aofbuf(nlatit) - nltmcn) * 0.01_ireals
        roblon   =  (aofbuf(nlongt) - nlgmcn) * 0.01_ireals
        IF (      (roblat <= dinlat) .AND. (roblat >= dislat)                  &
            .AND. (roblon >= diwlon) .AND. (roblon <= dielon))  THEN
 
!         observation header
!         ------------------
          WRITE(nuaofex,'("**************************************************" &
                        &,"**** A O F    R E P O R T ************************" &
                        &,"***************")')
          WRITE(nuaofex,'("0",4(5(i11,1x),/," "))') (aofbuf(i),i=1,nhdob)
          WRITE(nuaofex,'(" ",4(5(o11,1X),/," "))') (aofbuf(i),i=1,nhdob)
 
!         observation body
!         ----------------
          iend  =  aofbuf(1)-1
          istr  =  nhdob+1
          IF (    (aofbuf(nobtp) == nsynop) .OR. (aofbuf(nobtp) == nairep)     &
              .OR.(aofbuf(nobtp) == ndribu) .OR. (aofbuf(nobtp) == nsatob)) THEN
!           single-level reports : SYNOP/AIREP/DRIBU/SATOB
            WRITE(nuaofex,'("0",4(5(i11,1x),/," "))') (aofbuf(i),i=istr,iend)
            WRITE(nuaofex,'(" ",4(5(o11,1x),/," "))') (aofbuf(i),i=istr,iend)

          ELSEIF ((aofbuf(nobtp) == ntemp) .OR. (aofbuf(nobtp) == npilot)) THEN
!           multi-level reports  : TEMP/PILOT
            iskip = nlevn (aofbuf ( nobtp))
            DO     j = istr,iend,iskip
              iend2 = iskip + j - 1
              WRITE(nuaofex,'("0",4(5(i11,1x),/," "))') (aofbuf(i),i=j,iend2)
              WRITE(nuaofex,'(" ",4(5(o11,1x),/," "))') (aofbuf(i),i=j,iend2)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  !-----------------------------------------------------------------------
  ! Section  2.1      Pre-selection of reports
  !-----------------------------------------------------------------------

!   observation and code type
    nobtyp   =  aofbuf ( nobtp)
    ncdtyp   =  aofbuf ( ncdtp)
 
!   Get WMO station identity
!   ------------------------
    ist1=aofbuf(nstid1)
    ist2=aofbuf(nstid2)

    DO nc =3,0,-1
     is  =  IAND (ISHFT(ist1,-nc*7) , 127)
     ystid(4-nc:4-nc) =  achar(is)
     is  =  IAND (ISHFT(ist2,-nc*7) , 127)
     ystid(8-nc:8-nc) =  achar(is)
    ENDDO
    DO nc = 9, ilstidp
      ystid(nc:nc) = ' '
    ENDDO

!   get diagnostic array position (--> nobtpp, ncdtpp)
!   --------------------------------------------------
    CALL obs_pointrs ( nobtyp , ncdtyp )
    noctpr (nobtpp,ncdtpp) = noctpr (nobtpp,ncdtpp) + 1

!   pre-set passive flag to zero
    mpassiv = 0
    mrepflg = 0
 
!   reject report if data base flag on lat / long / date / time is high
!   -------------------------------------------------------------------
    DO   jflg =  1,4
      nvrfwr    =  ibits (aofbuf (nqllta),nf3bps(jflg),nibits (nf3boc))
      nflag     =  ibits (nvrfwr,nf3ibp( 4),nibits (nf3ioc( 4)))
      iovrid    =  ibits (nvrfwr,nf3ibp( 3),nibits (nf3ioc( 3)))
      IF (iovrid == 1)           nflag   =  0
      IF (nflag  >= 2)                                          THEN
        neventr (nedbfl,ncdtpp,nobtpp) = neventr(nedbfl,ncdtpp,nobtpp) + 1
        noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        WRITE( nurej,'(" STATION ",A ," : DATA BASE FLAG ",I1,                 &
                      &" TOO HIGH")' ) ystid, jflg
        CYCLE read_obs
      ENDIF
    ENDDO       
 
  ! Observation time
  !-----------------------------------------------------------------------

!   observation time
!   ----------------
!   actual time of observation as hhmm
    iobdat = aofbuf (ndate)
    ntime  = aofbuf (nextim)
!   observation synoptic time
    IF (ntime == nmdi) ntime = aofbuf (nsyntm)
    IF ((iobdat == nmdi) .OR. (ntime == nmdi)) THEN
      neventr (netime,ncdtpp,nobtpp) = neventr(netime,ncdtpp,nobtpp) + 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : NO OBS. DATE / TIME")')  ystid
      CYCLE read_obs
    ENDIF
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
    IF ((nobtyp == nsynop).OR.(nobtyp == ndribu)) THEN
      timlim  =  taof + MAX( wtuksue , tipmxsu )
    ELSEIF ((nobtyp == ntemp).OR.(nobtyp == npilot)) THEN
      timlim  =  taof + MAX( wtukrse , tipolmx )
    ELSE
      timlim  =  taof +        wtukare
    ENDIF
    IF (ntstep*dt/3600._ireals > timlim) THEN
      neventr (netime,ncdtpp,nobtpp) = neventr(netime,ncdtpp,nobtpp) + 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : OBSERVATION TOO OLD: ",F6.1,           &
                    &" [FORECAST HRS]")' ) ystid, taof
      CYCLE read_obs
    ENDIF
 
  ! Observation position
  !-----------------------------------------------------------------------

!   station altitude (missing)
!   --------------------------
    IF (aofbuf (naltit) == nmdi)      THEN
      IF (nobtyp == nsynop .OR. nobtyp == ntemp .OR. nobtyp == npilot)         &
        neventr (nenoal,ncdtpp,nobtpp) = neventr(nenoal,ncdtpp,nobtpp) + 1
      nelevt    =  imdi
      IF (nobtyp == nsynop) THEN
        noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        WRITE( nurej,'(" STATION ",A ," : SYNOP STA. HEIGHT MISSING",          &
                      &F6.1," HRS")' ) ystid, taof
        CYCLE read_obs
      ENDIF
      IF (nobtyp == ndribu) nelevt = 0
    ELSE
      nelevt     =  aofbuf (naltit) - nalmcn
    ENDIF

!   observation location
!   --------------------
    roblat     =  REAL ( aofbuf(nlatit) - nltmcn, ireals ) * 0.01_ireals
    roblon     =  REAL ( aofbuf(nlongt) - nlgmcn, ireals ) * 0.01_ireals

    lindenb = .FALSE.
    IF ((ystid(1:5) == '10393') .OR. (ystid(1:5) == '10394')) THEN
!     Lindenberg Met Observatory reports have to be assigned to a specific grid
!     point, that is used by Lindenberg staff for diagnostic purposes and has
!     appropriate hand-made values in the external soil / surface parameters.
!     This is done by assigning the report to pre-specified geographical coord.
!     and then assigning this to the horizontally nearest model grid pt. (by
!     setting nelevt=imdi for obs_aof_assign_gridpt)
      IF (      (roblat > 52.1_ireals) .AND. (roblat < 52.35)                  &
          .AND. (roblon > 14.1_ireals) .AND. (roblon < 14.25)) THEN
        lindenb = .TRUE.
        roblat  =  52.220_ireals
        roblon  =  14.135_ireals
        nelevtl =  nelevt
        nelevt  =  imdi
      ENDIF
    ENDIF

    rlon = rla2rlarot (roblat, roblon, pollat, pollon, polgam)
    rlat = phi2phirot (roblat, roblon, pollat, pollon)
!   observation position in grid point units (total area)
    zio_tot    = 1._ireals + (rlon - startlon_tot) /dlon
    zjo_tot    = 1._ireals + (rlat - startlat_tot) /dlat

!   If ship, set ship switch
!   ------------------------
    llship   =    .FALSE.
    IF ((ncdtyp >= nshscd).AND.(ncdtyp <= natshs)) llship = .TRUE.

!   Assign observation to grid point
!   --------------------------------
    lseaobs = llship .OR. (nobtyp == ndribu) .OR. (ncdtyp == nshtcd)           &
                     .OR. (ncdtyp == nrocsh) .OR. (ncdtyp == nshpcd)
    CALL obs_aof_assign_gridpt ( zio_tot, zjo_tot, nelevt, nobtyp , lseaobs    &
                               , iob_tot, job_tot, zsurf , lsfcob )
!   ==========================

    IF (     (ABS(fbogus) > epsy)                                              &
        .AND.((iob_tot /= ionl ) .OR. (job_tot /= jonl ))                      &
        .AND.((iob_tot /= ionl2) .OR. (job_tot /= jonl2)))  iob_tot = 0

!   reject report if observation location outside model domain
!   ----------------------------------------------------------
    IF   (iob_tot == 0) THEN
      IF (mpassiv == 0)                                                        &
        neventr (neloca,ncdtpp,nobtpp) = neventr (neloca,ncdtpp,nobtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF DOMAIN ",2F7.1    &
                   &," , ",F6.1," HRS")' ) ystid, roblon, roblat, taof
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      CYCLE read_obs
    ENDIF

!   observation location outside user-defined area
!   ----------------------------------------------
    IF (     (roblat > obnlat) .OR. (roblat < obslat)                          &
        .OR. ((obwlon >= obelon) .AND. (roblon < obwlon)                       &
                                 .AND. (roblon > obelon))                      &
        .OR. ((obwlon <  obelon) .AND. (     (roblon < obwlon)                 &
                                        .OR. (roblon > obelon)))) THEN
      IF (mpassiv == 0)                                                        &
        neventr (neloca,ncdtpp,nobtpp) = neventr (neloca,ncdtpp,nobtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : OBS. LOCATION OUT OF USER-SPECIFIED"   &
                   &," AREA ",2F7.1," , ",F6.1," HRS")' )                      &
             ystid, roblon, roblat, taof
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (.NOT. lverif) CYCLE read_obs
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvobbp )
    ENDIF

!   special Lindenberg treatment:
!   report is never rejected here, even if "model orography - station altitude"
!   is too large (which could occur if only if doromx(.) is chosen very small)
    IF (lindenb) THEN
      iob_tot = ABS( iob_tot )
      job_tot = ABS( job_tot )
      nelevt  = nelevtl
      lsfcob  = .TRUE.
    ENDIF

!   distance "model orography - station altitude" too large
!   -------------------------------------------------------
!   SYNOP / DRIBU
    IF (iob_tot <  0) THEN
      IF (mpassiv == 0)                                                        &
        neventr (nezdif,ncdtpp,nobtpp) = neventr (nezdif,ncdtpp,nobtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : HEIGHT ",I5," DIFF. TO MODEL ",F5.0    &
                   &," TOO LARGE, ",F7.1," HRS")') ystid, nelevt, zsurf, taof
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (.NOT. lverif) CYCLE read_obs
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvalbp )
      iob_tot = ABS( iob_tot )
      job_tot = ABS( job_tot )
    ENDIF
!   TEMP / PILOT
    IF (((nobtyp == ntemp) .OR. (nobtyp == npilot)) .AND. (.NOT. lsfcob)) THEN
      IF ((nelevt /= imdi) .AND. (mpassiv == 0)) THEN
        neventr (nezdif,ncdtpp,nobtpp) = neventr (nezdif,ncdtpp,nobtpp) + 1
        mrepflg = insert( mrepflg , 1 , nvalbp )
      ENDIF
      WRITE( nurej,'("STATION ",A ," : NO SURFACE-LEVEL REPORT DERIVED FROM "  &
                                     &,"TEMP / PILOT")' ) ystid
    ENDIF

  ! Observation type, code type, and black lists
  !-----------------------------------------------------------------------

!   Check if black listed ship
!   --------------------------
!   station confidence (used for blacklisted ships)
    istcon   =  IBITS (aofbuf(nstcar),nscfbp,nibits (nscfoc))

    IF ((llship) .AND. (istcon == 1))                 THEN
      IF (mpassiv == 0)                                                        &
        neventr (neblak,ncdtpp,nobtpp) = neventr(neblak,ncdtpp,nobtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : BLACKLISTED SHIP")') ystid
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (.NOT. lverif) CYCLE read_obs
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvbkbp )
    ENDIF

!   Check for bad aircraft station id
!   ---------------------------------
    IF (      (nobtyp == nairep)                                               &
        .AND. (     (ystid(1:3) == "XXX") .OR. (ystid(1:3) == "???")           &
               .OR. (ystid(1:3) == "///") .OR. (ystid(1:3) == "***"))) THEN 
      IF (mpassiv == 0)                                                        &
        neventr (neblak,ncdtpp,nobtpp) = neventr(neblak,ncdtpp,nobtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : BAD AIRCRAFT ID")')  ystid
!     WRITE( nustat,'(A20,": CAUTION: bad aircraft station ID: ",A ," ==> "    &
!                   &,"problem in data base / makeaof ???")' )  yroutine, ystid
      PRINT         '(A20,": CAUTION: bad aircraft station ID: ",A ," ==> "    &
                    &,"problem in data base / makeaof ???")' ,  yroutine, ystid
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      IF (.NOT. lverif) CYCLE read_obs
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mpassiv = 2
      mrepflg = insert( mrepflg , 1 , nvbkbp )
    ENDIF


!   Check whether observation type is to be excluded
!   ------------------------------------------------
    IF (nobtyp == nsynop .AND. .NOT.lsynop    .OR.                             &
        nobtyp == nairep .AND. .NOT.laircf    .OR.                             &
        nobtyp == nsatob .AND. .NOT.lsatob    .OR.                             &
        nobtyp == ndribu .AND. .NOT.ldribu    .OR.                             &
        nobtyp == ntemp  .AND. .NOT.ltemp     .OR.                             &
        nobtyp == npilot .AND. .NOT.lpilot    .OR.                             &
        nobtyp == nsatem .AND. .NOT.lsatem        )              THEN

      IF (      (roblat <= exnlat+atol) .AND. (roblat >= exslat-atol)          &
          .AND. (roblon >= exwlon-atol) .AND. (roblon <= exelon+atol))  THEN  
        IF (mpassiv == 0)                                                      &
          neventr (neobct,ncdtpp,nobtpp) = neventr(neobct,ncdtpp,nobtpp) + 1
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        IF (.NOT. lverif) CYCLE read_obs
        noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
        mpassiv = 2
        mrepflg = insert( mrepflg , 1 , nvexbp )
      ENDIF
    ENDIF

!   Check whether code type is to be excluded
!   -----------------------------------------
    IF (ncdtyp == nsrscd .AND. .NOT.lcd011    .OR.                             &
        ncdtyp == natscd .AND. .NOT.lcd014    .OR.                             &
        ncdtyp == nshscd .AND. .NOT.lcd021    .OR.                             &
        ncdtyp == nabscd .AND. .NOT.lcd022    .OR.                             &
        ncdtyp == nshred .AND. .NOT.lcd023    .OR.                             &
        ncdtyp == natshs .AND. .NOT.lcd024    .OR.                             &
        ncdtyp == naircd .AND. .NOT.lcd041    .OR.                             &
        ncdtyp == ncodar .AND. .NOT.lcd141    .OR.                             &
        ncdtyp == ncolba .AND. .NOT.lcd241    .OR.                             &
        ncdtyp == namdar .AND. .NOT.lcd144    .OR.                             &
        ncdtyp == nacar  .AND. .NOT.lcd244    .OR.                             &
        ncdtyp == nstbcd .AND. .NOT.lcd088    .OR.                             &
        ncdtyp == nsst   .AND. .NOT.lcd188    .OR.                             &
        ncdtyp == ndrbcd .AND. .NOT.lcd165    .OR.                             &
        ncdtyp == nbathy .AND. .NOT.lcd063    .OR.                             &
        ncdtyp == ntesac .AND. .NOT.lcd064    .OR.                             &
        ncdtyp == nldtcd .AND. .NOT.lcd035    .OR.                             &
        ncdtyp == nshtcd .AND. .NOT.lcd036    .OR.                             &
        ncdtyp == ntdrop .AND. .NOT.lcd135    .OR.                             &
        ncdtyp == nrocob .AND. .NOT.lcd039    .OR.                             &
        ncdtyp == nrocsh .AND. .NOT.lcd040    .OR.                             &
        ncdtyp == nldpcd .AND. .NOT.lcd032    .OR.                             &
        ncdtyp == nshpcd .AND. .NOT.lcd033    .OR.                             &
        ncdtyp == nwp_eu .AND. .NOT.lcd132    .OR.                             &
        ncdtyp == nra_eu .AND. .NOT.lcd133    .OR.                             &
        ncdtyp == npr_us .AND. .NOT.lcd136    .OR.                             &
        ncdtyp == nravad .AND. .NOT.lcd137    .OR.                             &
        ncdtyp == nstmcd .AND. .NOT.lcd086    .OR.                             &
        ncdtyp == nstovs .AND. .NOT.lcd186)         THEN

      IF (      (roblat <= exnlat+atol) .AND. (roblat >= exslat-atol)          &
          .AND. (roblon >= exwlon-atol) .AND. (roblon <= exelon+atol))  THEN  
        IF (mpassiv == 0)                                                      &
          neventr (neobct,ncdtpp,nobtpp) = neventr(neobct,ncdtpp,nobtpp) + 1
        IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
        IF (.NOT. lverif) noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
        IF (.NOT. lverif) CYCLE read_obs
        noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
        mpassiv = 2
        mrepflg = insert( mrepflg , 1 , nvexbp )
      ENDIF
    ENDIF


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
      IF (mpassiv == 0)                                                        &
        neventr (neloca,ncdtpp,nobtpp) = neventr(neloca,ncdtpp,nobtpp) + 1
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      WRITE( nurej,'(" STATION ",A ," : NO APPROPRIATE NODE FOUND ",           &
                   &2F7.1," , ",F6.1," HRS")' ) ystid, roblon, roblat, taof
      CYCLE read_obs
    ENDIF

!   count accepted reports
    IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) + 1

  !----------------------------------------------------------------------
  ! Section  2.3  Store this report in buffer 'nodbuf' ( --> 'iallbuf')
  !               together with the already extracted header information
  !----------------------------------------------------------------------

! passive reports are rejected, if lverpas = .FALSE.

    IF ((.NOT. lverpas) .AND. (mpassiv > 0)) THEN
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      CYCLE read_obs
    ENDIF


!   actual buffer length of array 'nodbuf'
    ipos   = ibuflen
!   report length of current observation: 
!                    header information + orig. AOF report
    irplen = nbufhed + ilen
    IF (ipos+irplen > nmxbln)       THEN
      IF (.NOT. lmulti(nobtyp))  nexceed (1)  =  nexceed(1) + 1
      IF (      lmulti(nobtyp))  nexceed (2)  =  nexceed(2) + 1
      WRITE( nustat,'(A20,": CAUTION: maxmlo, maxsgo TOO SMALL ==> REPORT RE"  &
                    &,"JECTED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:" &
                    &,I7," / ",I7)' )                  yroutine, ipos, nmxbln
      PRINT         '(A20,": CAUTION: maxmlo, maxsgo TOO SMALL ==> REPORT RE"  &
                    &,"JECTED", / ,26X,"SINCE ACTUAL / MAXIMUM BUFFER LENGTH:" &
                    &,I7," / ",I7)' ,                  yroutine, ipos, nmxbln
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!     IF (mpassiv == 0)                                                        &
        neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
      IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
      IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
      CYCLE read_obs
    ENDIF
    isurfob = 0
    IF (lsfcob) isurfob = 1
    IF (lseaobs)      mrepflg = insert( mrepflg , 1 , nvsebp )
    IF (mpassiv == 2) mrepflg = insert( mrepflg , 1 , nvpsbp )
    nodbuf (ipos+ 1) = irplen
    nodbuf (ipos+ 2) = nobtyp
    nodbuf (ipos+ 3) = ncdtyp
    nodbuf (ipos+ 4) = nvrfwr
    nodbuf (ipos+ 5) = nflag
    nodbuf (ipos+ 6) = ntime
    nodbuf (ipos+ 7) = INT( imindif , iintegers )
    nodbuf (ipos+ 8) = nelevt
    nodbuf (ipos+ 9) = NINT( zsurf   *1000._ireals )
    nodbuf (ipos+10) = NINT( zio_tot *1000._ireals )
    nodbuf (ipos+11) = NINT( zjo_tot *1000._ireals )
    nodbuf (ipos+12) = NINT( roblon  *1000._ireals )
    nodbuf (ipos+13) = NINT( roblat  *1000._ireals )
    nodbuf (ipos+14) = iob_tot
    nodbuf (ipos+15) = job_tot
    nodbuf (ipos+16) = isurfob
    nodbuf (ipos+17) = mrepflg
    nodbuf (ipos+18) = iobdat
    DO   i = 1 , ilen
      nodbuf (ipos+nbufhed+i) = aofbuf(i)
    ENDDO
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

!   Check if observation belongs to next timebox (moved to start of loop due to
!   --------------------------------------------  'CYCLE read_obs' (ODR size) )
!   IF (taof >= naoffbx+naofbox-epsy)           EXIT read_obs

    ENDDO   read_obs
!   ~~~~~~~~~~~~~~~~
    IF (MAXVAL(nexceed(1:2)) > 0)                                              &
      PRINT '("taof2",2F6.1,I4)' , taof, naoffbx, naofbox
 
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

    IF (MAXVAL(nexceed(1:2)) > 0) THEN
      OPEN (nucautn, FILE=yucautn, FORM='FORMATTED', STATUS='UNKNOWN'          &
                                 , POSITION='APPEND', IOSTAT=istat)
      IF (istat /= 0) yerr = 'OPENING OF FILE yucautn FAILED'
      IF (istat /= 0) CALL model_abort (my_cart_id, 1408, yerr, yroutine)
    ENDIF
    IF (nexceed(1) > 0) THEN
      nexceed (1)  =  NINT( nexceed(1) *MAX( c1, naofbox /MAX( wtuksua+wtuksue &
                                                             , 2*tipmxsu, c1 )))
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER FOR"          &
                     &," READING SINGLE-LEVEL OBS:",I9)' ) ntstep, nmxbln
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' ) nexceed(1)
      WRITE( nucautn,'("        (OR INCREASE NAMELIST VARIABLE maxmlo")' )
      PRINT          '("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER FOR"          &
                     &," READING SINGLE-LEVEL OBS:",I9)' , ntstep, nmxbln
      PRINT          '("   ==>  INCREASE NAMELIST VARIABLE maxsgo BY AT LEAST" &
                     &,I6)' , nexceed(1)
      PRINT          '("        (OR INCREASE NAMELIST VARIABLE maxmlo")'
    ENDIF
    IF (nexceed(2) > 0) THEN
      nexceed (2)  =  NINT( nexceed(2) *MAX( c1, naofbox /MAX( wtukrsa+wtukrse &
                                                             , 2*tipolmx, c1 )))
      WRITE( nucautn,'("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER FOR"          &
                     &," READING MULTI-LEVEL OBS:",I9)' ) ntstep, nmxbln
      WRITE( nucautn,'("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' ) nexceed(2)
      WRITE( nucautn,'("        (OR INCREASE NAMELIST VARIABLE maxsgo")' )
      PRINT          '("CAUTION !!!!! t=",I5,": TOO SMALL BUFFER FOR"          &
                     &," READING MULTI-LEVEL OBS:",I9)' , ntstep, nmxbln
      PRINT          '("   ==>  INCREASE NAMELIST VARIABLE maxmlo BY AT LEAST" &
                     &,I6)' , nexceed(2)
      PRINT          '("        (OR INCREASE NAMELIST VARIABLE maxsgo")'
    ENDIF
    IF (MAXVAL(nexceed(1:2)) > 0)  CLOSE (nucautn)

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
                           , imp_integers, icomm_cart, yerrmsl, izmplcode )
!     ====================

      IF (izmplcode /= 0)                                                      &
        CALL model_abort (my_cart_id, 11178, yerrmsl, yroutine, izmplcode)

    ENDIF

  ENDIF  ! num_compute > 1


  IF (lwonl)                                                                   &
    WRITE(nupr,'(" OBS_READ_DISTRIBUTE:  length of first report received",     &
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

!-------------------------------------------------------------------------------
! End Subroutine obs_read_distribute
!-------------------------------------------------------------------------------
END SUBROUTINE obs_read_distribute


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for saving the observation header
!-------------------------------------------------------------------------------

SUBROUTINE obs_save_head ( kpos , laddrep , nexceed )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" saves observation header
!   data in the header arrays of the observation data record ODR.
!
! Method:
!   Header data extraced in module 'obs_read_distribute' and temporarily
!   stored in the observation buffer are saved in header arrays for
!   multi- and single level reports respectivly.
!   12.11.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version  Date       Name
! -------- ---------- ----
! 1.13     1998/10/22 Michael Buchhold  : Initial release
! 1.19     1998/12/11 Christoph Schraff : Processing also of passive reports 
!                                         (for verific. mode). Revised messages.
! 1.27     1999/03/29 Christoph Schraff : Revised station characteristic format
!                                         and missing data indicator in the ODR.
! 1.31     1999/07/01 Christoph Schraff : 'mpassiv' moved from outside.
! 1.36     2000/02/24 Michael Buchhold  : Introduction of ACAR aircraft reports.
!                                         Adaptation to 32 bit AOF word length.
! 1.38     2000/04/06 Christoph Schraff : ODR station characteristics extended
!                                         (e.g. by aircraft roll angle).
! 2.5      2001/06/01 Christoph Schraff : Special value of instrument specific.
!                                         indicator for ship observations.
! 2.13     2002/01/18 Christoph Schraff : Type of instrument (bit) removed.
! 3.3      2003/04/22 Christoph Schraff : Flag for bias correction of Vaisala
!                                         RS80 humidity data.
! 3.6      2003/12/11 Christoph Schraff : Introduction of obs type SATOB.
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

  INTEGER (KIND=iintegers), INTENT(IN)     ::  &
   kpos                  ! index of last datum of previous observation
  LOGICAL, INTENT(OUT)                     ::  &
   laddrep               ! if true then add present report to ODR
  INTEGER (KIND=iintegers), INTENT(INOUT)  ::  &
   nexceed(4)            ! number of reports in excess of array size

! Local parameters: None
! ----------------

! Local variables
! ---------------

  INTEGER (KIND=iintegers) ::  &
   i, nc ,             & ! loop indices
   ist1  , ist2 ,      & ! integer equivalent of station ID (part 1 and 2)
   is    ,             & ! position of a character in ASCII collating seq.
   nstachr,            & ! station characteristics
   iobdat,             & ! observation date as yymmdd or yyyymodd
   iob, job,           & ! grid point asigned to observation (sub-domain)
   io_tot,jo_tot,      & ! grid point asigned to observation (total domain)
   nelevt,             & ! station altitude
   isurfob,            & ! = 1, if surface data are to be used
   mrepflg,            & ! report flag word as part of station charcteristics
   mpassiv,            & ! passive flag assigned to report
   nbufbit,            & ! buffer of bit pattern
   ilen,               & ! length of original AOF report
   insert,             & ! statement function to insert any bit pattern
   invar,              & ! input variable
   inval ,             & ! bit pattern to be inserted
   ibit,               & ! bit position
   ibits,              & ! statement function to unpack any bit pattern
   ibp,                & ! bit position of bit struct. to be unpacked /replaced
   icovr                 ! no. of bit occ. by the bit structure

  REAL (KIND=ireals)       :: &
   zsurf,              & ! height of model grid pt. to which obs. is assigned
   zio_tot, zjo_tot      ! station loc. in grid point units (total domain)

  CHARACTER (LEN=20)                    :: &
  yroutine               ! name of this subroutine
  CHARACTER (LEN=ilstid)                :: &
  ystid                  ! station identity

!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ibp ,icovr) = IAND (ISHFT(invar,-ibp),icovr)

! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar , ISHFT(inval,ibit))

!-------------------------------------------------------------------------------
! Begin Subroutine obs_save_header
!-------------------------------------------------------------------------------

  yroutine = 'obs_save_header'
  laddrep  = .TRUE.

  !-----------------------------------------------------------------------------
  ! Section 1.0: Get Observation header
  !-----------------------------------------------------------------------------

  nobtyp  = nodbuf(kpos + 2)
  ncdtyp  = nodbuf(kpos + 3)
  nvrfwr  = nodbuf(kpos + 4)
  nflag   = nodbuf(kpos + 5)
  ntime   = nodbuf(kpos + 6)
  taof    = REAL ( nodbuf(kpos + 7), ireals ) /60._ireals
  nelevt  = nodbuf(kpos + 8)
  zsurf   = REAL ( nodbuf(kpos + 9), ireals ) /1000._ireals
  zio_tot = REAL ( nodbuf(kpos +10), ireals ) /1000._ireals
  zjo_tot = REAL ( nodbuf(kpos +11), ireals ) /1000._ireals
  roblon  = REAL ( nodbuf(kpos +12), ireals ) /1000._ireals
  roblat  = REAL ( nodbuf(kpos +13), ireals ) /1000._ireals
  io_tot  = nodbuf(kpos +14)
  jo_tot  = nodbuf(kpos +15)
  isurfob = nodbuf(kpos +16)
  mrepflg = nodbuf(kpos +17)
  iobdat  = nodbuf(kpos +18)
  mpassiv = ibits( mrepflg , nvpsbp , nibits(nvscoc) ) * 2

! get diagnostic array position (--> nobtpp, ncdtpp)
! --------------------------------------------------
  CALL obs_pointrs ( nobtyp , ncdtyp )

  lsurfob = .FALSE.
  IF (isurfob == 1) lsurfob = .TRUE.
  ilen    = nodbuf( kpos+nbufhed+1)

! If ship, set ship switch
! ------------------------
  llship   =    .FALSE.
  IF ((ncdtyp >= nshscd).AND.(ncdtyp <= natshs)) llship = .TRUE.

! Assigned grid point in the processor domain
  iob     = io_tot - isubpos (my_cart_id,1) + 1 + nboundlines
  job     = jo_tot - isubpos (my_cart_id,2) + 1 + nboundlines

! save original AOF report in array 'aofbuf'
! -----------------------------------------
  DO i=1,ilen
    aofbuf(i) = nodbuf( kpos+nbufhed+i)
  ENDDO
  ist1=aofbuf(nstid1)
  ist2=aofbuf(nstid2)

  DO nc =3,0,-1
   is  =  IAND (ISHFT(ist1,-nc*7) , 127)
   ystid(4-nc:4-nc) =  achar(is)
   is  =  IAND (ISHFT(ist2,-nc*7) , 127)
   ystid(8-nc:8-nc) =  achar(is)
  ENDDO
  DO nc = 9, ilstid
    ystid(nc:nc) = ' '
  ENDDO

! check whether there is space in the ODR
! ---------------------------------------
  IF ( nobtyp == ntemp .OR. nobtyp == npilot ) THEN
     IF (nmlob >= maxmll) THEN
       nexceed (2)  =  nexceed(2) + 1
       IF (lwonl)                                                              &
         WRITE( nupr,'(" CAUTION !!!!! STATION ",A ,": ",I5,                   &
                      &"th MULTI-LEVEL REPORT EXCEEDS ODR SIZE MAXMLL ",I5)')  &
                ystid, nmlob, maxmll
       PRINT         '(" CAUTION !!!!! STATION ",A ,": ",I5,                   &
                      &"th MULTI-LEVEL REPORT EXCEEDS ODR SIZE MAXMLL")' ,     &
                ystid, nmlob
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!      IF (mpassiv == 0)                                                       &
         neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
       IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
       IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
       noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
       laddrep = .FALSE.
       RETURN                
     ENDIF
     IF ( nsgob >= maxsgl .AND. lsurfob ) THEN
       nexceed (1)  =  nexceed(1) + 1
       IF (lwonl)                                                              &
         WRITE( nupr,'(" CAUTION !!!!! STATION ",A ,": ",I5,                   &
                      &"th SINGLE-LEVEL REPORT EXCEEDS ODR SIZE MAXSGL ",I5)') &
                ystid, nsgob, maxsgl
       PRINT         '(" CAUTION !!!! STATION ",A ,": ",I5,                    &
                      &"th SINGLE-LEVEL REPORT EXCEEDS ODR SIZE MAXSGL")' ,    &
                ystid, nsgob
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!      IF (mpassiv == 0)                                                       &
         neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
       lsurfob = .FALSE.
     ENDIF
  ELSEIF ( nobtyp == nsynop .OR. nobtyp == ndribu .OR. nobtyp == nairep        &
                                                  .OR. nobtyp == nsatob) THEN
     IF (nsgob >= maxsgl) THEN
       nexceed (1)  =  nexceed(1) + 1
       IF (lwonl)                                                              &
         WRITE( nupr,'(" CAUTION !!!!! STATION ",A ,": ",I5,                   &
                      &"th SINGLE-LEVEL REPORT EXCEEDS ODR SIZE MAXSGL ",I5)') &
                ystid, nsgob, maxsgl
       PRINT         '(" CAUTION !!!! STATION ",A ,": ",I5,                    &
                      &"th SINGLE-LEVEL REPORT EXCEEDS ODR SIZE MAXSGL")' ,    &
                ystid, nsgob
! insufficient ODR size is the only report event, which is updated even for
! reports set to passive previously
!      IF (mpassiv == 0)                                                       &
         neventr (nesodr,ncdtpp,nobtpp) = neventr(nesodr,ncdtpp,nobtpp) + 1
       IF (mpassiv == 0) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
       IF (mpassiv >  0) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
       noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
       laddrep = .FALSE.
       RETURN                 
     ENDIF
  ELSE
     nexceed (3)  =  nexceed(3) + 1
     nexceed (4)  =  nobtyp
     IF (lwonl)                                                                &
       WRITE( nupr,'(" CAUTION !!!!! UNKNOWN OBSERVATION TYPE, nobtyp=",i2)')  &
              nobtyp
     PRINT         '(" CAUTION !!!!! UNKNOWN OBSERVATION TYPE, nobtyp=",i2)',  &
              nobtyp
     laddrep = .FALSE.
     RETURN
  ENDIF

! create ODR station characteristic word, incl. flags and instr. specification
! ----------------------------------------------------------------------------
  nstachr = mrepflg
  nbufbit = ibits ( aofbuf(nstcar) , nstcbp , nibits(nstcoc) )
  nstachr = insert( nstachr , nbufbit , nvscbp )
  nbufbit = ibits ( aofbuf(nstcar) , nscfbp , nibits(nscfoc) )
  nstachr = insert( nstachr , nbufbit , nvssbp )
  nbufbit = ibits ( aofbuf(nstcar) , nstibp , nibits(nstioc) )
  nstachr = insert( nstachr , nbufbit , nvsibp )
  nbufbit = ibits ( aofbuf(ninstr) , ninsbp , nibits(ninsoc) )
  SELECT CASE (nbufbit)
  CASE (127)
    nbufbit = 0
    IF (llship)  nbufbit = 6
  CASE ( 32)
    nbufbit = 2
    IF (llship)  nbufbit = 7
  CASE ( 21)
    nbufbit = 4
  CASE ( 22)
    nbufbit = 5
  CASE ( 18)
    nbufbit = 18
  CASE DEFAULT
    nbufbit = 1
    IF (llship)  nbufbit = 6
  END SELECT
  nstachr = insert( nstachr , nbufbit , nvinbp )
  nbufbit = ibits ( aofbuf(nstcar) , narabp , nibits(naraoc) )
  nbufbit = MIN( nbufbit, nibits(nvaaoc) )
  nstachr = insert( nstachr , nbufbit , nvaabp )
  nbufbit = ibits ( aofbuf(nstcar) , npafbp , nibits(npafoc) )
  nstachr = insert( nstachr , nbufbit , nvapbp )

 !------------------------------------------------------------------------------
 !  Section  2.0   Fill observation header arrays
 !------------------------------------------------------------------------------

  IF ( nobtyp == ntemp .OR. nobtyp == npilot ) THEN

!   multi-level reports
!   -------------------
    nmlob = nmlob + 1
    omlhed (nmlob,nhilon) = roblon
    omlhed (nmlob,nhjlat) = roblat
    IF (nelevt == imdi) THEN
      omlhed (nmlob,nhalt ) = rmdi
    ELSE
      omlhed (nmlob,nhalt ) = REAL ( nelevt, ireals )
    ENDIF
    omlhed (nmlob,nhtime) = taof
    omlhed (nmlob,nhsurf) = zsurf
    omlhed (nmlob,nhzio ) = zio_tot
    omlhed (nmlob,nhzjo ) = zjo_tot
    omlhed (nmlob,nhvcbu) = c1
    omlhed (nmlob,nhvcbt) = c1
    omlhed (nmlob,nhvcbq) = c1
    omlhed (nmlob,nhvctu) = c1
    omlhed (nmlob,nhvctt) = c1
    omlhed (nmlob,nhvctq) = c1
    momlhd (nmlob,nhio  ) = iob
    momlhd (nmlob,nhjo  ) = job
    momlhd (nmlob,nhitot) = io_tot
    momlhd (nmlob,nhjtot) = jo_tot
    momlhd (nmlob,nhobtp) = nobtyp
    momlhd (nmlob,nhcode) = ncdtyp
    momlhd (nmlob,nhschr) = nstachr
    momlhd (nmlob,nhqofl) = aofbuf(nqllta)
    momlhd (nmlob,nhpass) = mpassiv
    momlhd (nmlob,nhqcfw) = - 1
    momlhd (nmlob,nhhrmn) = ntime
    momlhd (nmlob,nhdate) = iobdat
!   READ (ydate_ini,'(I8)' )   momlhd (nmlob,nhdate)
!CS: temporary: set all radiosondes to Vaisala RS92, except the ones specified
    IF (      (nobtyp == ntemp)                                                &
        .AND. (ystid(1:3) /= '066'  ) .AND. (ystid(1:2) /= '07')               &
        .AND. (ystid(1:5) /= '10954') .AND. (ystid(1:5) /= '11747'))           &
      momlhd (nmlob,nhrtyp) = 80
!CS
    momlhd (nmlob,nhuexi) = 0
    momlhd (nmlob,nhaexi) = 0
    momlhd (nmlob,nhtexi) = 0
    momlhd (nmlob,nhqexi) = 0
    yomlhd (nmlob)        = ystid
  ENDIF
  IF (      (nobtyp == nsynop) .OR. (nobtyp == ndribu) .OR. (nobtyp == nairep) &
       .OR. (nobtyp == nsatob)                                                 &
       .OR. (lsurfob .AND. ((nobtyp == ntemp) .OR. (nobtyp == npilot)))) THEN

!   single-level reports
!   --------------------
    nsgob = nsgob + 1
    osghed (nsgob,nhilon) = roblon
    osghed (nsgob,nhjlat) = roblat
    IF (nelevt == imdi) THEN
      osghed (nsgob,nhalt ) = rmdi
    ELSE
      osghed (nsgob,nhalt ) = REAL ( nelevt, ireals )
    ENDIF
    osghed (nsgob,nhtime) = taof
    osghed (nsgob,nhsurf) = zsurf
    osghed (nsgob,nhzio ) = zio_tot
    osghed (nsgob,nhzjo ) = zjo_tot
    osghed (nsgob,nhvcbu) = c1
    osghed (nsgob,nhvcbt) = c1
    osghed (nsgob,nhvcbq) = c1
    osghed (nsgob,nhvctu) = c1
    osghed (nsgob,nhvctt) = c1
    osghed (nsgob,nhvctq) = c1
    mosghd (nsgob,nhio  ) = iob
    mosghd (nsgob,nhjo  ) = job
    mosghd (nsgob,nhitot) = io_tot
    mosghd (nsgob,nhjtot) = jo_tot
    mosghd (nsgob,nhobtp) = nobtyp
    mosghd (nsgob,nhcode) = ncdtyp
    mosghd (nsgob,nhschr) = nstachr
    mosghd (nsgob,nhqofl) = aofbuf(nqllta)
    mosghd (nsgob,nhpass) = mpassiv
    mosghd (nsgob,nhqcfw) = - 1
    mosghd (nsgob,nhhrmn) = ntime
    READ (ydate_ini,'(I8)' )   mosghd (nsgob,nhdate)
    yosghd (nsgob)        = ystid
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_save_head
!-------------------------------------------------------------------------------
END SUBROUTINE obs_save_head


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for multi-level observation processing
!-------------------------------------------------------------------------------

SUBROUTINE obs_multi_level ( lvprof )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" performs the
!   multi-level report processing
!
! Method:
!   Multi-level observation (report) data extraction. Also, an observation
!   error is assigned to each datum.
!   In addition, surface analysis data are extracted.
!   26.11.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version  Date       Name
! -------- ---------- ----
! 1.13     1998/10/22 Michael Buchhold
!  Initial release
! 1.19     1998/12/11 Christoph Schraff
! 1.19     1998/12/11 Christoph Schraff : Processing also of passive reports 
!                                        (for verif. mode). Updated data events.
!                                        2nd part of flags replaced by threshold
!                                        QC flags. ANSI violations removed.
! 1.20     1999/01/07 Guenther Doms     : Renaming of some global variables.
! 1.27     1999/03/29 Christoph Schraff : ODR flag formats revised. Station id
!                                       6-bit holleriths replaced by characters.
!                                       Revised (32-bit) missing data indicator.
! 1.31     1999/07/01 Christoph Schraff : Selection of levels made independent 
!                                         from 'lverif' to ensure consistency.
! 1.39     2000/05/03 Ulrich Schaettler : Changed variable names for lat/lon.
! 1.42     2000/06/19 Christoph Schraff : Flagged (blacklisted) geopotential
!                                         obs. stored to ODR and VOF.
! 2.5      2001/06/01 Christoph Schraff : Gross error limits revised.
! 2.13     2002/01/18 Michael Buchhold  : Use SODAR/RASS temperature profiles.
!                                         Adjustment of treatment of humidity 
!                                         obs to cloud ice scheme option.
! 3.6      2003/12/11 Christoph Schraff : No further bias correction of Vaisala 
!                                         RS80 bias-corrected humidity data.
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
  LOGICAL, INTENT(INOUT)  ::  &
   lvprof               ! .TRUE. if number of levels > 0, mult-level obs.

! Local parameters: None
! ----------------

! Local scalars:
! -------------
  LOGICAL                  ::  &
   lsurfc,            & ! .true. if active surface-level data encountered
   lsurfp,            & ! .true. if surface-level data encountered
   lbogus,            & ! .true. if bogus data are produced
   lplist,            & ! .true. if level heights given in pressure units
   lzlist,            & ! .true. if level heights given in m
   lsignu,            & ! .true. if significant level for wind data    
   lsignt,            & ! .true. if significant level for temperature
   lyes,lno,          & !  switches
   lseaobs,           & ! .TRUE. if observation platform on sea
   lseaonly,          & ! .TRUE. if observations from sea platforms
                        ! are to be used only.
   lrs80,             & ! .TRUE. if Vaisala RS80 with bias-corrected humidity
   lactrep              ! .TRUE. if the report is active


  INTEGER (KIND=iintegers) ::  &
   insert,            & ! statement function to insert any bit pattern
   invar ,            & ! word to be unpacked / to be partly replaced
   inval ,            & ! bit pattern to be inserted
   ibit  ,            & ! bit position
   ibits ,            & ! statement function to unpack any bit pattern
   ipos  ,            & ! bit position of bit struct. to be unpacked /replaced
   icovr ,            & ! no. of bit occ. by the bit structure
   ireplace,          & ! statement function to replace any bit pattern
   iboc  ,            & ! no. of bit occ. by the bit structure to be replaced
   irepl ,            & ! word containing the replacing bit structure
   iposr ,            & ! bit position of the replacing bit structure
   istat,             & ! error status variable
   j,jloop, icl,      & ! loop indices
   ilen,              & ! length of output record
   iob,job,           & ! assigned model grid point
   isurfc,            & ! surface level indicator
   istart,            & ! pointer to observation body of first level
   iskip,             & ! number of  data for one level
   iend,              & ! pointer to observation body of last level
   klev,              & ! level counter
   iduvmax,           & ! max. wind level indicator
   idtropo,           & ! topopause level indicator
   idusign,           & ! significant  wind level indicator 
   idtsign,           & ! significant temperature level indicator 
   nzaexi,            & ! height data indicator 
   nvarfw,            & ! flag word for 1 variable
   nppfwr,            & ! flag word for pressure
   izzfwr,            & ! flag word for height
   ivrfwr, nuvfwr,    & ! flag words for wind data
   nttfwr,            & ! flag word for temperature
   ntdfwr,            & ! flag word for dew point
   nrdbfl,            & ! final flag word
   nflag1,            & ! flag word
   nvarib1,           & ! extracted variable
   mpassiv              ! passive flag assigned to report


  REAL    (KIND=ireals   ) ::  &
    tt,               & ! temperature dummy value in statement funct. fspw, fspe
    alogpv,           & ! dummy for log of 'water vapour pressure / b1' in ftd
    roberr,           & ! observation error
    zstalt,           & ! station altitude
    zzaltsy,          & ! SYNOP station altitude
    zsurf,            & ! model orography height
    fisd,             & ! height difference between model orography and station
                        ! altitude
    fisduv,           & ! modified height diff. for surface wind data
    fisdtt,           & ! modified height diff. for 2m temperature data
    fisdrh,           & ! modified height diff. for 2m dew point
    rzzz,             & ! observed height
    zddd, zfff,       & ! observed wind direction and speed
    rsus, rsvs,       & ! u- and v-component in the.TRUE.geograph. system
    ruuu, rvvv,       & ! u- and v-component in rotated grid
    rttt, rtdt,       & ! observed temperature and dew point
    rheo,             & ! observed rel. humidity
    eobs,             & ! observed vapour pressure
    eswo,             & ! saturation vapour pressure over water
    eseo,             & ! saturation vapour pressure over ice
    fspw, fspe,       & ! statement functions for comput. of saturation humidity
    ftd                 ! statement function for dewpoint temperature 'td'
!   obs_model_pressure , z2p ! statement function for calculating model pressure


  CHARACTER (LEN=20)                    :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid)   :: &
    ystid,            & ! station identity
    yyes , yno             

!
!------------ End of header ----------------------------------------------------

! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ipos,icovr) = IAND (ISHFT(invar,-ipos),icovr)


! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar , ISHFT(inval,ibit))


! Statement function to replace any bit pattern by another bit pattern
! --------------------------------------------------------------------
  ireplace ( invar, ipos, iboc, irepl, iposr ) =                               &
               IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )          &
                  , ISHFT( IAND( ISHFT( irepl,-iposr ), nibits(iboc) ), ipos ) )


! Statement - functions for computation of saturation humidity
! ------------------------------------------------------------

! magnus formela for water
  fspw(tt)        = b1 * exp( b2w*(tt - b3)/(tt - b4w) )
! magnus formela for ice
  fspe(tt)        = b1 * exp( b2i*(tt - b3)/(tt - b4i) )

! statement - functions for computation of dewpoint temperature td
! ----------------------------------------------------------------

! identify: vapour pressure at t = saturation vapour pressure at td
! inverse of magnus formula  (input is log of 'water vapour pressure / b1')

  ftd (alogpv)    = ( b3*b2w -alogpv*b4w ) / ( b2w -alogpv )


!-------------------------------------------------------------------------------
! Begin Subroutine obs_multi_level
!-------------------------------------------------------------------------------

  yroutine = "obs_multi_level"

  !-----------------------------------------------------------------------------
  ! Section 1.0  Multi-level report processing (start up)
  !-----------------------------------------------------------------------------

! get assigned model grid point
  iob           = momlhd (nmlob,nhio  )
  job           = momlhd (nmlob,nhjo  )

! model height
  zsurf         = omlhed(nmlob,nhsurf) 

! get wmo station identity
  ystid         = yomlhd (nmlob)

! set level id to 0
  nlevin        =  0

! preset pressure to missing indicator
  rppp          =  rmdi

! set pressure flag to 0
  nppfwr        =  0

! set observation error to missing indicator
  roberr        =  rmdi
 
! preset surface level indicators
  isurfc        =  0
  lsurfc        = .FALSE.
  lsurfp        = .FALSE.

! station altitude
  zstalt        = omlhed (nmlob,nhalt )

! report active
  lactrep       = (momlhd(nmlob,nhpass) == 0)
  nzaexi        = 0

  !-----------------------------------------------------------------------------
  ! Section 2.0  Multi-level observation (mass and wind and humidity
  !              analysis data extraction)
  !-----------------------------------------------------------------------------

! Sort report in vertical
! -----------------------
 
  IF (nobtyp /= npilot) CALL obs_sort_levels ( nvar( 1) , iob , job , zsurf )
!                       ====================             
  IF (nobtyp == npilot) CALL obs_sort_levels ( nvar( 7) , iob , job , zsurf )
!                       ====================             

! this is not done any more !!!
! ( Convert SODAR/RASS vitual temperatures into normal air temperatures )
! IF (nobtyp == npilot .AND. ncdtyp == nra_eu)  CALL obs_RASS_tv2t (iob,job)
!                                               ==================             
 
! determine whether data are modified
! -----------------------------------
 
  lbogus =       (ABS(fbogus) >  epsy)                                         &
           .AND. (     (iob == ionl  .AND. job == jonl )                       &
                  .OR. (iob == ionl2 .AND. job == jonl2))
 
! set up  loop parameters
! -----------------------
 
  istart   =    nhdob + nlev0 (nobtyp)
  iskip    =    nlevn (nobtyp)
  iend     =    aofbuf (nrecln) - nlevn (nobtyp)

! Fill ODR of surface obs. with missing or dummy values
 
  IF (lsurfob) THEN
    osgbdy (nsgob,1:mxsbdy) = rmdi
    mosgbd (nsgob,1:mxsbdf) = imdi
  ENDIF
 
! loop over report
! ----------------
 
  klev = 0
  Report_level: DO   jloop   =   istart,     iend,   iskip
! ================
  klev = klev + 1
 
  mpassiv  =  0
 
! check if number of level exceeds ODR size
! -----------------------------------------
 
  IF (klev > maxrsl) THEN
! insufficient ODR size is the only report and data events, which is updated
! even for reports set to passive previously
    neventd (nelodr,ncdtpp,nobtpp) = neventd (nelodr,ncdtpp,nobtpp) + 1
    IF (jloop > iend-iskip) THEN
      ilen = 3 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt4
        outbuf(nacout+3) = klev
        DO icl = 1 , istrej
          outbuf(nacout+3+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF
    ENDIF
    klev = klev - 1
                                                             CYCLE  Report_level
  ENDIF
 
! Fill current level of ODR with missing or dummy values
 
  omlbdy (nmlob,klev,1:mxrbdy) = rmdi
  momlbd (nmlob,klev,1:mxrbdf) = imdi
  momlbd (nmlob,klev,  nbtflg) = 0
 
! Set logical switches for pressure/height lists
 
  lplist   =    .FALSE.
  lzlist   =    .FALSE.
 
! get pressure
! ------------

  CALL obs_extract(jloop,nvar( 1),aofbuf,nmxlob)
 
  IF ((nobtyp == npilot) .AND. ((nvarib == nmdi) .OR. (nflag >= 2))) THEN
     rppp     =    rmdi
     lzlist   =   .TRUE.
  ELSE
     IF ((nvarib == nmdi) .AND. (lactrep)) neventd (nelmis,ncdtpp,nobtpp) =    &
                                           neventd (nelmis,ncdtpp,nobtpp) + 1
     IF ((nflag  >= 2   ) .AND. (lactrep)) neventd (nelflg,ncdtpp,nobtpp) =    &
                                           neventd (nelflg,ncdtpp,nobtpp) + 1
     IF ((nvarib == nmdi) .OR. (nflag >= 2)) THEN
        klev = klev - 1
        CYCLE Report_level
     ENDIF
     rppp     =    nvarib * rfact(nvar( 1))
     lplist   =    .TRUE.

!    use no observation above 'rpplim'
!    ---------------------------------
     IF ((rppp  <  rpplim) .AND. (lactrep))  THEN
        neventd (nelext,ncdtpp,nobtpp) = neventd (nelext,ncdtpp,nobtpp) + 1
        klev = klev - 1
        CYCLE Report_level
     ENDIF
  ENDIF
  nppfwr  =  nvrfwr
 
! Extract level id
! ----------------

  IF (nwrvrf(nvar(1),nobtyp) == 0) THEN 
    IF (lactrep) neventd(nelflg,ncdtpp,nobtpp) = neventd(nelflg,ncdtpp,nobtpp)+1
    klev = klev - 1
    CYCLE  Report_level
  ENDIF
  nlevin   =  0
  nlevin   =   ibits ( aofbuf (jloop + nwrvrf(nvar ( 1),nobtyp)),              &
                                     nlinbp(nobtyp),nibits (nlinoc) )
! Extract surface level indic.
  isurfc   =  0
  isurfc   =   ibits (nlevin,nlidbp( 7),nibits (nlidoc))
  IF ( isurfc == 1 .AND. lsurfc ) THEN
!   several surface levels
    IF (lactrep) neventd(nelsfc,ncdtpp,nobtpp) = neventd(nelsfc,ncdtpp,nobtpp)+1
    ilen = 2 + istrej
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+1) = ilen
      outbuf(nacout+2) = nfmt5
      DO icl = 1 , istrej
        outbuf(nacout+2+icl) = ICHAR( yomlhd(nmlob) (icl:icl) )
      ENDDO
      nacout  = nacout + ilen
    ENDIF
    klev = klev - 1
    CYCLE Report_level
  ENDIF
 
  iduvmax  = ibits (nlevin,nlidbp( 1),nibits (nlidoc))
  idtropo  = ibits (nlevin,nlidbp( 2),nibits (nlidoc))
  idusign  = ibits (nlevin,nlidbp( 8),nibits (nlidoc))
  idtsign  = ibits (nlevin,nlidbp( 9),nibits (nlidoc))
  lsignu   = iduvmax == 1 .OR. idusign == 1
  lsignt   = idtropo == 1 .OR. idtsign == 1
 
! height
! ------
 
  Height:  DO
  rzzz     =    rmdi
  IF ((nzex(nobtyp) == 0) .AND. (.NOT.lzlist))                       EXIT Height
  CALL obs_extract (jloop,nvar  ( 7),aofbuf,nmxlob)

  IF (lzlist) THEN

!    height of observation level given in m
     IF ((nvarib == nmdi) .OR. (nflag >= 2))  THEN
        IF ((nvarib == nmdi) .AND. (lactrep)) neventd (nelmis,ncdtpp,nobtpp) = &
                                              neventd (nelmis,ncdtpp,nobtpp) + 1
        IF ((nflag  >= 2   ) .AND. (lactrep)) neventd (nelflg,ncdtpp,nobtpp) = &
                                              neventd (nelflg,ncdtpp,nobtpp) + 1
        klev = klev - 1
        CYCLE Report_level
     ENDIF
     nppfwr =  nvrfwr
     rzzz   =  nvarib * rfact (nvar  ( 7)) - nalmcn*1._ireals
     IF (rzzz < zsurf) THEN
        IF (lactrep)                                                           &
          neventd (nelnop,ncdtpp,nobtpp) = neventd (nelnop,ncdtpp,nobtpp) + 1
        klev = klev - 1
        CYCLE Report_level
     ENDIF
!    no pressure, take model pressure at same height
     rppp   = z2p ( rzzz , iob , job )
     IF (rppp < rmdich)   THEN
!    level is above top model level (levels below orography are excluded above)
        IF (lactrep)                                                           &
          neventd (nelnop,ncdtpp,nobtpp) = neventd (nelnop,ncdtpp,nobtpp) + 1
        klev = klev - 1
        CYCLE Report_level
     ENDIF
! ipolbit = ibits (nppfwr,nf1bps(2),nibits(nf1boc(2)))
  
  ELSE

!    height of observation level given in pressure units
     IF ((nvarib == nmdi) .OR. (nflag >= 2))  THEN
        IF ((nvarib == nmdi) .AND. (lactrep)) neventd (nepmis,ncdtpp,nobtpp) = &
                                              neventd (nepmis,ncdtpp,nobtpp) + 1
        IF ((nflag >= 2)     .AND. (lactrep)) neventd (nepflg,ncdtpp,nobtpp) = &
                                              neventd (nepflg,ncdtpp,nobtpp) + 1
        IF (nvarib == nmdi)                                          EXIT Height
        mpassiv = 1
     ENDIF
     rzzz  =  nvarib * rfact (nvar  ( 7)) - nalmcn*1._ireals
  ENDIF

  IF ((mpassiv > 0) .OR. (nobtyp == npilot) .OR. (zstalt <= rmdich)) THEN
    roberr = rmdi
  ELSE
    roberr = obs_error ( nvrpoe(3) )
  ENDIF

! check if level below surface and save height
! --------------------------------------------
  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nppfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nppfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nppfwr, nf1bps(3) )

  IF (      (isurfc == 1) .AND. (lsurfob)                                      &
      .AND. ((nobtyp == ntemp) .OR. (nobtyp == npilot))) THEN
     lsurfp = .TRUE.
     osgbdy (nsgob     ,nbsz  ) = rzzz
     osgbdy (nsgob     ,nbszer) = roberr
     mosgbd (nsgob     ,nbsflg) = insert( 0, nvarfw, nvfzbp )
     IF ((mpassiv == 0) .AND. (zstalt > rmdich)) lsurfc = .TRUE.
  ENDIF
  IF ((zstalt > rmdich) .AND. (nobtyp == ntemp .OR. nobtyp == npilot)) THEN
     IF ((rzzz > zstalt-MAX( rprlim*12,40.0_ireals )) .OR. (isurfc == 1)) THEN
       omlbdy (nmlob,klev,nbtz  ) = rzzz
       IF (nzaexi == 0) nzaexi = -1
!      do not use radiosonde surface pressure for surface pressure nudging
!      IF ((rzzz <= zstalt) .OR. (isurfc == 1)) THEN
       IF (rzzz <= zstalt-c2) THEN
         omlbdy (nmlob,klev,nbtzer) = rmdi
!      observed level below surface
         IF ((isurfc /= 1) .AND. (rzzz > rmdich)) THEN
           IF ((mpassiv == 0) .AND. (lactrep)) neventd(nelext,ncdtpp,nobtpp) = &
                                               neventd(nelext,ncdtpp,nobtpp) + 1
           mpassiv = 2
         ENDIF
       ELSE
         omlbdy (nmlob,klev,nbtzer) = roberr  
         IF (mpassiv == 0) nzaexi = 1
       ENDIF
       momlbd (nmlob,klev,nbtflg) = insert( 0, nvarfw, nvfzbp )
       IF (mpassiv == 2) momlbd (nmlob,klev,nbtflg) =                          &
                 insert( momlbd (nmlob,klev,nbtflg) , 1 , nvflbp )
     ELSEIF (rzzz > rmdich) THEN
!      observed level below surface
       IF (lactrep)                                                            &
         neventd (nelext,ncdtpp,nobtpp) = neventd(nelext,ncdtpp,nobtpp) + 1
       klev = klev - 1
       CYCLE Report_level
     ENDIF
  ELSEIF (nobtyp == ntemp .OR. nobtyp == npilot) THEN
     omlbdy (nmlob,klev,nbtz  ) = rzzz
     omlbdy (nmlob,klev,nbtzer) = rmdi
     momlbd (nmlob,klev,nbtflg) = insert( 0, nvarfw, nvfzbp )
     IF (nzaexi == 0) nzaexi = -1
  ENDIF

  EXIT Height
  ENDDO Height
 
! Write out pressure.
! -------------------
! Check if pressure < 'pminsig' and significant level
  
  IF (rppp <  pminsig-epsy) THEN
!   mandla  =  ibits (nlevin,nlidbp( 6),nibits (nlidoc))
    IF (     (NINT(rppp-70000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-50000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-40000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-30000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-25000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-20000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-15000._ireals) /= 0)                                   &
        .AND.(NINT(rppp-10000._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 7000._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 5000._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 3000._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 2500._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 2000._ireals) /= 0)                                   &
        .AND.(NINT(rppp- 1000._ireals) /= 0)) THEN
       IF (lactrep)                                                            &
         neventd (nelsig,ncdtpp,nobtpp) = neventd(nelsig,ncdtpp,nobtpp) + 1
!      re-set previous event counters for pressure, if current level is not used
       IF ((nvarib == nmdi) .AND. (lactrep)) neventd (nepmis,ncdtpp,nobtpp) =  &
                                             neventd (nepmis,ncdtpp,nobtpp) - 1
       IF ((nflag >= 2)     .AND. (lactrep)) neventd (nepflg,ncdtpp,nobtpp) =  &
                                             neventd (nepflg,ncdtpp,nobtpp) - 1
       klev = klev - 1
       CYCLE Report_level
    ENDIF
  ENDIF
  IF ((isurfc == 1) .AND. (lsurfob)) THEN
     osgbdy (nsgob,nbsp  ) = rppp
     mosgbd (nsgob,nbslid) = nlevin
     mosgbd (nsgob,nbsqcf) = 0
  ENDIF
  omlbdy (nmlob,klev,nbtp  ) = rppp
  momlbd (nmlob,klev,nbtlid) = nlevin
  IF (lbogus)                                                                  &
    omlbdy (nmlob,klev,nbtp) = omlbdy(nmlob,klev,nbtp) + 100._ireals*fbogus
  omlbdy (nmlob,klev,nbtlop) = LOG (omlbdy (nmlob,klev,nbtp) )

! Check if report is assigned to a sea grid pt.,
! and get height [m] of surface observation
! and height [m] difference to model orography
 
  lseaobs = (ibits(momlhd(nmlob,nhschr),nvsebp,nvscoc)  ==  1)
  IF ((isurfc == 1) .AND. lsurfob ) THEN
     IF (osgbdy(nsgob,nbsz) >  rmdich) THEN
       zzaltsy = osgbdy(nsgob,nbsz)
     ELSE
       zzaltsy = osghed(nsgob,nhalt)
     ENDIF
!    fisd      = osghed(nsgob,nhsurf) - osghed(nsgob,nhalt)
     fisd      = osghed(nsgob,nhalt) - osghed(nsgob,nhsurf)
  ENDIF
 
! Write out latitude / longit. [grid pt. units] and threshold Q.C. flag of level
! ------------------------------------------------------------------------------
  omlbdy (nmlob,klev,nbtzio) = omlhed (nmlob,nhzio)
  omlbdy (nmlob,klev,nbtzjo) = omlhed (nmlob,nhzjo)
  momlbd (nmlob,klev,nbtqcf) = 0

! Temperature and humidity
! ------------------------

  lrs80  = (ibits( momlhd(nmlob,nhschr),nvinbp,nibits(nvinoc) )  ==  18)

  Temp_Hum:  DO
! =============

  IF (mpassiv <= 1) mpassiv = 0

! Temperature
! -----------
  rttt     =    rmdi
  IF (ntex (nobtyp) == 0)                                          EXIT Temp_Hum
  CALL obs_extract (jloop,nvar( 4),aofbuf,nmxlob)
  IF ((mpassiv == 0) .AND. (lactrep)) THEN
    IF ((nvarib == nmdi) .AND. ( .NOT.lsignu .OR. lsignt .OR.(isurfc == 1)))   &
                        neventd (netmis,ncdtpp,nobtpp) =                       &
                        neventd (netmis,ncdtpp,nobtpp) + 1
    IF (nflag  >= 2   ) neventd (netflg,ncdtpp,nobtpp) =                       &
                        neventd (netflg,ncdtpp,nobtpp) + 1
  ENDIF
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif)))   EXIT Temp_Hum
  IF (nflag >= 2) mpassiv = MAX( mpassiv , i1 )
  nttfwr   =    nvrfwr
  rttt     =    nvarib * rfact(nvar  ( 4))
  IF (ntot (nobtyp) == 0)                                          EXIT Temp_Hum
  IF (      (rttt < t0_melt-90._ireals) .OR. (rttt > t0_melt+60._ireals)       &
      .OR. ((rttt > t0_melt+20) .AND. (rppp < 70000._ireals))                  &
      .OR. ((rttt > t0_melt+ 5) .AND. (rppp < 50000._ireals))                  &
      .OR. ((rttt > t0_melt- 5) .AND. (rppp < 40000._ireals))) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (netext,ncdtpp,nobtpp) = neventd(netext,ncdtpp,nobtpp) + 1
                                                                   EXIT Temp_Hum
  ENDIF
 
! Set switch for surface (2m) / upper air temperature
! ---------------------------------------------------

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nttfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nttfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nttfwr, nf1bps(3) )
! Use SODAR/RASS temperature data  provided as a PILOT report
! IF ((mpassiv > 0) .OR. (nobtyp == npilot)) THEN
  IF  (mpassiv > 0)                          THEN
    roberr = rmdi
  ELSE
    roberr = obs_error ( nvrpoe(6) )
  ENDIF
  IF ((isurfc == 1) .AND. (lsurfob)) THEN
     lyes = .FALSE.
     lno  = .FALSE.
     fisdtt = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd )) * fisd
     lseaonly = (ABS(altopsu(3))  <   epsy)
     IF (      (.NOT. lyes)                                                    &
         .AND. (     (.NOT. lseaonly .AND. (zzaltsy >  altopsu(3)))            &
                .OR. (      lseaonly .AND. .NOT. lseaobs)                      &
                .OR. (fisdtt >  doromx(3))))   lno = .TRUE.
     IF ((lno) .AND. (mpassiv == 0) .AND. (lactrep))                           &
       neventd (netalt,ncdtpp,nobtpp) = neventd(netalt,ncdtpp,nobtpp) + 1
     IF ((.NOT. lno) .OR. (lverif)) THEN
       IF (lno) roberr = rmdi
       IF (lno) nvarfw = insert( nvarfw, 1, nvfbps(4) )
       IF (lbogus) rttt = MAX( 100._ireals, rttt - 1._ireals*fbogus )
       osgbdy (nsgob,nbst  )  =  rttt
       osgbdy (nsgob,nbster)  =  roberr
       mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvftbp )
       lsurfp = .TRUE.
       IF ((mpassiv == 0) .AND. (.NOT. lno))  lsurfc = .TRUE.
     ENDIF
  ENDIF
  IF ((lbogus).AND.(rppp >= 66999._ireals).AND.(rppp <= 86000._ireals))        &
     rttt = MAX( 100._ireals, rttt + 1._ireals*fbogus )
  IF (isurfc == 1) roberr = rmdi
  omlbdy (nmlob,klev,nbtt  ) = rttt
  omlbdy (nmlob,klev,nbtter) = roberr
  momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg),nvarfw,nvftbp )
  IF ((mpassiv == 0) .AND. (isurfc /= 1)) THEN
    momlhd (nmlob,nhtexi) = 1
  ELSEIF (momlhd(nmlob,nhtexi) == 0) THEN
    momlhd (nmlob,nhtexi) = -1
  ENDIF
 
! Dew-point (only if there is temperature)
! ----------------------------------------

  rtdt     =    rmdi
  rheo     =    rmdi
  IF (ntdex (nobtyp) == 0)                                         EXIT Temp_Hum
  CALL obs_extract (jloop,nvar( 5),aofbuf,nmxlob)
  IF ((mpassiv == 0) .AND. (lactrep)) THEN
    IF (      (nvarib == nmdi) .AND. (rppp >= pqmin )                          &
        .AND. ((.NOT. lsignu) .OR. (lsignt) .OR. (isurfc == 1)))               &
                        neventd (neqmis,ncdtpp,nobtpp) =                       &
                        neventd (neqmis,ncdtpp,nobtpp) + 1
    IF (nflag  >= 2   ) neventd (neqflg,ncdtpp,nobtpp) =                       &
                        neventd (neqflg,ncdtpp,nobtpp) + 1
  ENDIF
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif)))   EXIT Temp_Hum
  IF (nflag >= 2)  mpassiv = MAX( mpassiv , i1 )
  ntdfwr   =    nvrfwr
  rtdt     =    nvarib * rfact (nvar  ( 5))
  IF (ntdot (nobtyp) == 0)                                         EXIT Temp_Hum
  IF (     (rtdt < t0_melt-150._ireals) .OR. (rtdt > t0_melt+40._ireals)       &
      .OR. ((rtdt < t0_melt-90._ireals) .AND. (isurfc == 1) .AND. (lsurfob))) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (neqlow,ncdtpp,nobtpp) = neventd(neqlow,ncdtpp,nobtpp) + 1
                                                                   EXIT Temp_Hum
  ENDIF
  IF ((rppp <  pqmin )         .AND. (mpassiv == 0) .AND. (lactrep))           &
     neventd (neq300,ncdtpp,nobtpp) = neventd (neq300,ncdtpp,nobtpp) + 1
  IF ((rppp <  pqmin )         .AND. (.NOT. lverif))               EXIT Temp_Hum
  IF  (rppp <  pqmin )  mpassiv = MAX( mpassiv , i1 )
 
! Set switch for surface (2m) / upper air dewpoint.
! Convert to model compatible relative humidity, and write to ODR
! ---------------------------------------------------------------

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), ntdfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), ntdfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), ntdfwr, nf1bps(3) )
  IF (rppp  <  pqmin )  nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
 
! eobs : observed vapour press. (measurement is inherently over ice)
!     (but td as stored in data file is over water (wmo convention),
!     which is why magnus formula over water instead of ice is used)
  eobs = fspw( rtdt )

! eswo : saturation vapour pressure over water at observed temperature
  eswo = fspw( rttt )

  rheo = eobs / eswo
  IF ((rheo >  rtshlm) .AND. (mpassiv == 0) .AND. (lactrep))                   &
     neventd(neqbig,ncdtpp,nobtpp) = neventd(neqbig,ncdtpp,nobtpp) + 1
  IF  (rheo >  rtshlm)                                             EXIT Temp_Hum

  IF ((rttt <  t0_melt) .AND. (itype_gscp <= 2)) THEN
!   eseo : saturation vapour pressure over ice at observed temperature
    eseo = fspe( rttt )

!   rheo : observed relative humidity over ice   (per def.  <= 100%)
    rheo = eobs / eseo
!        this is assumed equivalent to model relative humidity over
!        water (per def.  <= 100%; the model does not know ice) !!
!    ==> eobs-corr / eswo = eobs / eseo
    eobs = eswo * rheo
  ENDIF

  rheo = eobs / eswo
  rtdt = ftd( LOG( eobs / b1 ) )
  IF ((rheo >= rhtsat) .AND. (rheo < c1-epsy) .AND. (.NOT. lrs80)) THEN
     IF ((mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt <  t0_melt) neventd (neqsam,ncdtpp,nobtpp) =                        &
                       neventd (neqsam,ncdtpp,nobtpp) + 1 
       IF (rttt >= t0_melt) neventd (neqsap,ncdtpp,nobtpp) =                        &
                       neventd (neqsap,ncdtpp,nobtpp) + 1 
     ENDIF
     nvarfw  =  insert( nvarfw, 1, nvfbps(5) )
     rheo    =  c1
     rtdt    =  rttt
  ENDIF
  IF (rheo > c1+epsy) THEN
     IF ((mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt <  t0_melt) neventd (neqclm,ncdtpp,nobtpp) =                        &
                       neventd (neqclm,ncdtpp,nobtpp) + 1 
       IF (rttt >= t0_melt) neventd (neqclp,ncdtpp,nobtpp) =                        &
                       neventd (neqclp,ncdtpp,nobtpp) + 1 
     ENDIF
     nvarfw  =  insert( nvarfw, 1, nvfbps(6) )
  ENDIF
  IF (rheo > c1-epsy) THEN
     rheo    =  c1
     rtdt    =  rttt
  ENDIF

  roberr                        =   rherr1 / 100._ireals
  IF (rheo  <  rrhlim)  roberr  =   rherr2 / 100._ireals
  IF (rttt  <  rttlim)  roberr  =   rherr3 / 100._ireals
  IF (isurfc == 1)      roberr  = 4*rherr1 / 100._ireals
  IF (mpassiv > 0)      roberr  =   rmdi
  IF ((isurfc == 1).AND. lsurfob ) THEN
     lyes = .FALSE.
     lno  = .FALSE.
     DO   j  = 1, mxda
       WRITE(yyes ,'(i5,3x)') mrhyes( j )
       WRITE(yno  ,'(i5,3x)') mrhno ( j )
       IF (ystid == yyes) lyes = .TRUE.
       IF (ystid == yno ) lno  = .TRUE.
       IF (lwonl .AND. (lyes .OR. lno ))                                       &
         WRITE( nupr ,'(" STATION: ",a,"  RH YES/NO :",2(2x,a5))')             &
                ystid, yyes, yno
     ENDDO
     fisdrh = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd )) * fisd
     lseaonly = (ABS(altopsu(4))  <   epsy)
     IF (      (.NOT. lyes)                                                    &
         .AND. (     (.NOT. lseaonly .AND. (zzaltsy >  altopsu(4)))            &
                .OR. (      lseaonly .AND. .NOT. lseaobs)                      &
                .OR. (fisdrh >  doromx(4))))   lno = .TRUE.
     IF ((lno) .AND. (mpassiv == 0) .AND. (lactrep))                           &
       neventd (neqalt,ncdtpp,nobtpp) = neventd(neqalt,ncdtpp,nobtpp) + 1
     IF ((.NOT. lno) .OR. (lverif)) THEN
       IF (lno) roberr = rmdi
       IF (lno) nvarfw = insert( nvarfw, 1, nvfbps(4) )
       IF (lbogus) rheo = MAX( 0.01_ireals , MIN( c1, rheo +.1_ireals*fbogus ) )
       osgbdy (nsgob,nbsrh )  =  rheo
       osgbdy (nsgob,nbsqer)  =  roberr
       mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfqbp )
       lsurfp = .TRUE.
       IF ((mpassiv == 0) .AND. (.NOT. lno))  lsurfc = .TRUE.
     ENDIF
  ENDIF
  IF ((lbogus) .AND. (rppp >= 66999._ireals) .AND. (rppp <= 86000._ireals))    &
     rheo = MAX( 0.01_ireals , MIN( c1, rheo +0.1_ireals*fbogus ) )
  IF (isurfc == 1) roberr = rmdi
  omlbdy (nmlob,klev,nbtrh ) = rheo
  omlbdy (nmlob,klev,nbtqer) = roberr
  momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg),nvarfw,nvfqbp )
  IF ((mpassiv == 0) .AND. (isurfc /= 1)) THEN
    momlhd (nmlob,nhqexi) = 1
  ELSEIF (momlhd(nmlob,nhqexi) == 0) THEN
    momlhd (nmlob,nhqexi) = -1
  ENDIF

  EXIT Temp_Hum
  ENDDO Temp_Hum

! u & v components
! ----------------

  IF (mpassiv <= 1) mpassiv = 0
 
  Wind:  DO
! =========
  IF ((nuex(nobtyp) == 0) .OR. (nvex(nobtyp) == 0))                    EXIT Wind
  CALL obs_extract (jloop,nvar( 2),aofbuf,nmxlob)
  zddd     =  REAL ( nvarib, ireals )
  ivrfwr   =  nvrfwr
  nvarib1  =  nvarib
  nflag1   =  nflag
  CALL obs_extract (jloop,nvar( 3),aofbuf,nmxlob)
  IF ((mpassiv == 0) .AND. (lactrep)) THEN
    IF ((.NOT.lsignt) .OR. (lsignu) .OR. (isurfc == 1)) THEN
      IF ((nvarib1 == nmdi) .AND. (nvarib /= 0))                               &
                            neventd (nedmis,ncdtpp,nobtpp) =                   &
                            neventd (nedmis,ncdtpp,nobtpp) + 1
      IF  (nvarib  == nmdi) neventd (nefmis,ncdtpp,nobtpp) =                   &
                            neventd (nefmis,ncdtpp,nobtpp) + 1
    ENDIF
    IF ((nflag1 >= 2  ) .OR. ((ABS(nvarib1) > 360) .AND. (nvarib1 /= nmdi)     &
                                                   .AND. (nvarib /= 0)))       &
                        neventd (nedflg,ncdtpp,nobtpp) =                       &
                        neventd (nedflg,ncdtpp,nobtpp) + 1
    IF  (nflag  >= 2  ) neventd (nefflg,ncdtpp,nobtpp) =                       &
                        neventd (nefflg,ncdtpp,nobtpp) + 1
  ENDIF
  IF ((((nvarib1 == nmdi) .OR. (ABS(nvarib1) > 360)) .AND. (nvarib /= 0))      &
                          .OR. ((nflag1 >= 2) .AND. (.NOT. lverif)))   EXIT Wind
  IF   ((nvarib  == nmdi) .OR. ((nflag  >= 2) .AND. (.NOT. lverif)))   EXIT Wind
  IF ((nflag1 >= 2) .OR. (nflag >= 2))  mpassiv = MAX( mpassiv , i1 )
  zfff     =    REAL ( nvarib, ireals )
  IF  (      (zfff < -epsy) .OR. (zfff > 150.0_ireals)                         &
       .OR. ((zfff > 90.0_ireals) .AND. (rppp >= 70000._ireals))) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (nefneg,ncdtpp,nobtpp) = neventd(nefneg,ncdtpp,nobtpp) + 1
                                                                       EXIT Wind
  ENDIF
  zfff     =    MAX( zfff , c0 )
  IF   (zfff < epsy) zddd = c0
  nuvfwr   =    IOR (nvrfwr,ivrfwr)
  IF ((nuot(nobtyp) == 0) .OR. (nvot(nobtyp) == 0))                    EXIT Wind

! Transformation of wind speed and direction to wind components
! in the rotated model coordinate system
! -------------------------------------------------------------
  zfff  = zfff * rfact(nvar( 3))
  rsus  = zfff * (-SIN( zddd*degrad ))
  rsvs  = zfff * (-COS( zddd*degrad ))

  CALL uv2uvrot ( rsus,rsvs, roblat,roblon, pollat,pollon, ruuu,rvvv )
! =============

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nuvfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nuvfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nuvfwr, nf1bps(3) )

  roberr   =    obs_error ( nvrpoe( 1) )
  IF (mpassiv > 0)  roberr = rmdi
 
! Set switch for surface/non-surface wind
 
  IF ((isurfc == 1) .AND. (lsurfob)) THEN
     lyes = .FALSE.
     lno  = .FALSE.
     DO   j  = 1, mxda
       WRITE(yyes ,'(i5,3x)') muvyes( j )
       WRITE(yno  ,'(i5,3x)') muvno ( j )
       IF (ystid == yyes) lyes = .TRUE.
       IF (ystid == yno ) lno  = .TRUE.
     ENDDO
     fisduv = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd )) * fisd
     lseaonly = (ABS(altopsu(1))  <   epsy)
     IF (      (.NOT. lyes)                                                    &
         .AND. (     (.NOT. lseaonly .AND. (zzaltsy >  altopsu(1)))            &
                .OR. (      lseaonly .AND. .NOT. lseaobs)                      &
                .OR. (fisduv >  doromx(1))))   lno = .TRUE.
     IF ((lno) .AND. (mpassiv == 0) .AND. (lactrep))                           &
       neventd (nevalt,ncdtpp,nobtpp) = neventd(nevalt,ncdtpp,nobtpp) + 1
     IF ((.NOT. lno) .OR. (lverif)) THEN
       IF (lno) roberr = rmdi
       IF (lno) nvarfw = insert( nvarfw, 1, nvfbps(4) )
       IF (lbogus) THEN
          ruuu = ruuu - 1._ireals*fbogus
          rvvv = rvvv + 1._ireals*fbogus
       ENDIF
       osgbdy (nsgob,nbsu  )  =  ruuu
       osgbdy (nsgob,nbsv  )  =  rvvv
       osgbdy (nsgob,nbsuer)  =  roberr
       mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfubp )
       lsurfp = .TRUE.
       IF ((mpassiv == 0) .AND. (.NOT. lno))  lsurfc = .TRUE.
     ENDIF
  ENDIF
  IF (lbogus) THEN
     IF ((rppp >= 55999._ireals) .AND. (rppp <= 86000._ireals)) THEN
        ruuu = ruuu + 1._ireals*fbogus
        rvvv = rvvv + 1._ireals*fbogus
     ELSEIF ((rppp >= 26000._ireals).AND.(rppp <= 54000._ireals)) THEN
        ruuu = ruuu + 1._ireals*fbogus
     ELSEIF (rppp <= 24000._ireals) THEN
        ruuu = ruuu - 1._ireals*fbogus
        rvvv = rvvv - 1._ireals*fbogus
     ENDIF
  ENDIF
  IF (isurfc == 1) roberr = rmdi
  omlbdy (nmlob,klev,nbtu  ) = ruuu
  omlbdy (nmlob,klev,nbtv  ) = rvvv
  omlbdy (nmlob,klev,nbtuer) = roberr
  momlbd (nmlob,klev,nbtflg) = insert( momlbd(nmlob,klev,nbtflg),nvarfw,nvfubp )
  IF ((mpassiv == 0) .AND. (isurfc /= 1)) THEN
    momlhd (nmlob,nhuexi) = 1
  ELSEIF (momlhd(nmlob,nhuexi) == 0) THEN
    momlhd (nmlob,nhuexi) = -1
  ENDIF

  EXIT Wind
  ENDDO Wind

! Close the loop / take next level
! --------------------------------
 
  ENDDO Report_level
 
  !-----------------------------------------------------------------------------
  ! Section 3.0  Final processing of the multi-level report
  !-----------------------------------------------------------------------------

! Setting to missing value of the actual ODR
! ------------------------------------------
 
  IF ((lsurfob) .AND. (.NOT. lsurfp)) THEN
    osgbdy (nsgob,1:mxsbdy) = rmdi
    mosgbd (nsgob,1:mxsbdf) = imdi
    osghed (nsgob,1:mxshed) = rmdi
    mosghd (nsgob,1:mxshdf) = imdi
    nsgob    =  nsgob - 1
    lsurfob  =  .FALSE.
  ELSEIF ((lsurfob) .AND. (.NOT. lsurfc)) THEN
    mosghd(nsgob,nhpass) = 2
  ENDIF
  momlhd (nmlob,nhnlev) = klev
  IF (klev > 0) THEN
    IF (     (momlhd(nmlob,nhuexi) == 1)                                       &
        .OR. (momlhd(nmlob,nhtexi) == 1)                                       &
        .OR. (momlhd(nmlob,nhqexi) == 1)) THEN
      momlhd(nmlob,nhaexi) = 1
    ELSEIF (     (momlhd(nmlob,nhuexi) == -1)                                  &
            .OR. (momlhd(nmlob,nhtexi) == -1)                                  &
            .OR. (momlhd(nmlob,nhqexi) == -1)) THEN
      momlhd(nmlob,nhaexi) = -1
    ENDIF
    IF (      (MAX( momlhd(nmlob,nhaexi) , nzaexi ) <   1)                     &
        .AND. (MIN( momlhd(nmlob,nhaexi) , nzaexi ) == -1)) THEN
      IF ((lactrep) .AND. (lverif)) THEN
        neventr (nenoda,ncdtpp,nobtpp) = neventr(nenoda,ncdtpp,nobtpp) + 1
        noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
        noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      ENDIF
      momlhd(nmlob,nhpass) =  2
      IF (.NOT. lverif) THEN
        omlbdy (nmlob,1:klev,1:mxrbdy) = rmdi
        momlbd (nmlob,1:klev,1:mxrbdf) = imdi
        klev = 0
      ENDIF
    ENDIF
  ENDIF
  IF (klev <= 0) THEN
    omlhed (nmlob,1:mxrhed) = rmdi
    momlhd (nmlob,1:mxrhdf) = imdi
    nmlob = nmlob - 1
    IF (      lactrep)                                                         &
      neventr (nenoda,ncdtpp,nobtpp) = neventr(nenoda,ncdtpp,nobtpp) + 1
    IF (      lactrep) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
    IF (.NOT. lactrep) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
    noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
  ENDIF
  lvprof  =  klev > 0

!-------------------------------------------------------------------------------
! End Subroutine obs_multi_level
!-------------------------------------------------------------------------------
END SUBROUTINE obs_multi_level


!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for single-level observation processing
!-------------------------------------------------------------------------------

SUBROUTINE obs_single_level

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" performs the
!   single-level report processing (SYNOP/AIREP/DRIBU/SATOB).
!
! Method:
!   From a single-level observation (report), data are extracted and
!   observation errors are assigned (only for wind&mass and relative
!   humidity analysis).
!   Also, surface analysis ( / optional ) groups data are extracted.
!   14.11.97
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version  Date       Name
! -------- ---------- ----
! 1.13     1998/10/22 Michael Buchhold  : Initial release
! 1.19     1998/12/11 Christoph Schraff : Introduction of a verification mode
!                             (reports set passive instead of being rejected).
!                             Revised data events and flags. Allowed orography
!                             difference for T2m data independent from height.
! 1.20     1999/01/07 Guenther Doms     : Renaming of some global variables.
! 1.27     1999/03/29 Christoph Schraff : Processing also of optional groups.
!                       ODR: flag formats revised, station id 6-bit holleriths
!                       replaced by characters, 32-bit missing data indicator.
! 1.28     1999/04/19 Christoph Schraff : Processing of precip obs with missing
!                                         missing time range at synoptic times.
!                                         Corr. of index of 'nradex' (radiation)
! 1.31     1999/07/01 Christoph Schraff : 'mpassiv' moved from src_obs_processg.
! 1.36     2000/02/24 Michael Buchhold  : Get model pressure if AIREP reports
!                                         only height. Adapt to 32-bit AOF.
! 1.38     2000/04/06 Christoph Schraff : Extract pressure tendency for surface
!                                         pressure quality control threshold.
! 1.39     2000/05/03 Ulrich Schaettler : Changed variable names for lat/lon.
! 1.43     2000/07/18 Jean-Marie Bettems: For use on NEC SX5: call of 'atmhp'
!                                         using 'irealgrib'-kind reals.
! 2.5      2001/06/01 Christoph Schraff : Gross error limits revised. Data event
!                       and message for pressure tendency. Aircraft reports set
!                       to passive if less than 50m or 6hPa above the orography.
! 2.13     2002/01/18 Christoph Schraff
!  Adjustment of treatment of humidity observations to cloud ice scheme option.
!  Correction of unit conversion of cloud base height, and for 32-bit machines.
!  Discard precipitation observation if time range is not properly defined.
!  Verification of precip obs with any well-defined observation time range.
! 3.6      2003/12/11 Christoph Schraff : Introduction of obs type SATOB.
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
  LOGICAL                  ::  &
   lpprep,            & ! switch for observed/first guess pressure
   lpract,  lpractn,  & ! switch for bad reporting prctice
   lrej,              & ! switch for rejected observation
   ldat,              & ! .true. if report contains usable data
   ldatpsv,           & ! .true. if report contains any data
   lyes,lno,          & !  switches 
   lseaobs,           & ! .true. if observation platform on sea
   landsy ,           & ! .true. if land SYNOP
   lseaonly,          & ! .true. if observations from sea platforms
                        ! are to be used only.
   lactrep,           & ! .true. if the report is active
   lfog                 ! .true. if fog is derived from weather group

  INTEGER (KIND=iintegers) ::  &
   insert,            & ! statement function to insert any bit pattern
   invar ,            & ! word to be unpacked / to be partly replaced
   inval ,            & ! bit pattern to be inserted
   ibit  ,            & ! bit position
   ibits ,            & ! statement function to unpack any bit pattern
   ipos  ,            & ! bit position of bit struct. to be unpacked /replaced
   icovr ,            & ! no. of bit occ. by the bit structure
   ireplace,          & ! statement function to replace any bit pattern
   iboc  ,            & ! no. of bit occ. by the bit structure to be replaced
   irepl ,            & ! word containing the replacing bit structure
   iposr ,            & ! bit position of the replacing bit structure
   j,ivar, icl,       & ! loop indices
   ilen,              & ! length of output record
   iob,job,           & ! assigned model grid point
   istart,            & ! pointer to observation body 
   npcode, ipcode,    & ! pressure code
   nppfwr,            & ! flag word for pressure
   ntrfl,             & ! flag word for time range
   nrrfl,             & ! flag word for rain
   irrflo,            & ! original flag word for rain
   nvarfw,            & ! flag word for 1 variable   
   izzfwr,            & ! flag word for height       
   ivrfwr,nuvfwr,     & ! flag words for wind data
   nttfwr,            & ! flag word for temperature
   ntdfwr,            & ! flag word for dew point 
   nvisfl,            & ! flag for visibility
   npwwfl,            & ! flag for present weather
   ncctfl,            & ! flag for total cloud cover
   ncclfl,            & ! flag for low   cloud cover
   ncclqf,            & ! refined flag for low cloud cover
   irm                  ! status variable

  INTEGER (KIND=iintegers) ::  &
   npww,              & ! extracted present weather
   ncct,              & ! extracted total cloud cover
   nccl,              & ! derived   low   cloud cover
   nhkey,             & ! class (key) of cloud base height
   nclwgr,            & ! combined cloud and weather group (for verificat.):
                        !   (cloud base height: VUB WMO Code table 1600)
   ncwgfw,            & ! combined cloud, weather, ground, precip. flag word
   nrdbfl,            & ! final flag word
   nflag1,            & ! flag word
   nvarib1,           & ! extracted variable
   iexv,              & ! extracted variable bit pattern
   nrtr,              & ! class of precipitation time range
   nbufcl,            & ! short cloud group
   mpassiv              ! passive flag assigned to report


  REAL    (KIND=ireals   ) ::  &
    tt,               & ! temperature dummy value in statement funct. fspw, fspe
    alogpv,           & ! dummy for log of 'water vapour pressure / b1' in ftd
    zsurf ,           & ! model surface topography height
    roberr,           & ! observation error
    zadder,           & ! addition to observation error
    zstalt,           & ! station altitude
    rr,               & ! precipitation amount
    rtr,              & ! precipitation time range
    fisd,             & ! height difference between model orography and station
                        ! altitude
    fisdzz,           & ! modified height diff. for height or press. data
    fisduv,           & ! modified height diff. for surface wind data
    fisdtt,           & ! modified height diff. for 2m temperature data
    fisdrh,           & ! modified height diff. for 2m dew point
    rzzz,             & ! observed height          
    zzzz,             & ! observed height or pressure
    zddd, zfff,       & ! observed wind direction and speed
    rsus, rsvs,       & ! u- and v-component in the true geograph. system
    ruuu, rvvv,       & ! u- and v-component in rotated grid                
    rttt, rtdt,       & ! observed temperature and dew point
    rheo,             & ! observed rel. humidity
    eobs,             & ! observed vapour pressure
    eswo,             & ! saturation vapour pressure over water
    eseo,             & ! saturation vapour pressure over ice
    rttt2, rtdt2,     & ! observed 2m temperature and dew point
    fspw, fspe,       & ! statement functions for comput. of saturation humidity
    ftd                 ! statement function for dewpoint temperature 'td'


  REAL    (KIND=irealgrib) ::  &
    rzzz_grb,          & ! 
    rppp_grb             ! 


  CHARACTER (LEN=20)       :: &
    yroutine            ! name of this subroutine
  CHARACTER (LEN=ilstid)   :: &
    ystid,            & ! station identity
    yyes ,            & !
    yno                 !


!
!----------- End of header -----------------------------------------------------
 
! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ipos,icovr) = IAND (ISHFT(invar,-ipos),icovr)

! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar,ISHFT(inval,ibit) )

! Statement function to replace any bit pattern by another bit pattern
! --------------------------------------------------------------------
  ireplace ( invar, ipos, iboc, irepl, iposr ) =                               &
               IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )          &
                  , ISHFT( IAND( ISHFT( irepl,-iposr ), nibits(iboc) ), ipos ) )

 
! Statement - functions for computation of saturation humidity
! ------------------------------------------------------------

! magnus formela for water
  fspw(tt)        = b1 * exp( b2w*(tt - b3)/(tt - b4w) )
! magnus formela for ice
  fspe(tt)        = b1 * exp( b2i*(tt - b3)/(tt - b4i) )

! statement - functions for computation of dewpoint temperature td
! ----------------------------------------------------------------

! identify: vapour pressure at t = saturation vapour pressure at td
! inverse of magnus formula  (input is log of 'water vapour pressure / b1')

  ftd (alogpv)    = ( b3*b2w -alogpv*b4w ) / ( b2w -alogpv )
 
 
!-------------------------------------------------------------------------------
! Begin Subroutine obs_single_level
!-------------------------------------------------------------------------------

  yroutine = "obs_single_level"

  !-----------------------------------------------------------------------------
  ! Section 1.0  Single-level report processing: start up
  !-----------------------------------------------------------------------------

! get assigned model grid point
  iob           = mosghd (nsgob,nhio  )
  job           = mosghd (nsgob,nhjo  )

! model height
  zsurf         = osghed(nsgob,nhsurf)
 
! get wmo station identity
  ystid = yosghd (nsgob)

! single level body pointer
  istart        =  nhdob  + nlev0 (nobtyp)
 
! pressure code set to 0
  npcode        =  0

! level id set to 0
  nlevin        =  0
 
! set observation error to missing indicator
  roberr        =  rmdi
  zadder        =  0._ireals

! logical switch for observed/first guess pressure
  lpprep        = .TRUE.
  
! logical switch for bad reporting prctice (SYNOP)
  lpract        = .FALSE.
  lpractn       = .FALSE.
 
! logical switch for rejected report
  lrej          = .FALSE. 
  ldat          = .FALSE.
  ldatpsv       = .FALSE.
 
! preset pressure to missing indic.
  rppp          = rmdi
  nppfwr        = 0
 
! preset station altitude
  zstalt        = rmdi

! report active
  lactrep       = (mosghd(nsgob,nhpass) == 0)


  Extract_data: DO
! ================

  !-----------------------------------------------------------------------------
  ! Section 2.0  Pressure
  !-----------------------------------------------------------------------------
 
! Extract pressure and check missing value
! ----------------------------------------
 
  CALL obs_extract (istart,nvar( 1),aofbuf,nmxlob)
  IF ((nvarib == nmdi) .OR. (nvarib == 0) .OR. (nflag >= 2)) THEN
    lpprep     = .FALSE.
    rppp       =  rmdi
    nppfwr     =  0
    IF (((nvarib == nmdi) .OR. (nvarib == 0)) .AND. (lactrep)) THEN
      neventd (nepmis,ncdtpp,nobtpp) = neventd(nepmis,ncdtpp,nobtpp) + 1
    ELSEIF ((nflag >= 2) .AND. (lactrep)) THEN
      neventd (nepflg,ncdtpp,nobtpp) = neventd(nepflg,ncdtpp,nobtpp) + 1
    ENDIF

  ELSE
    rppp       = REAL (nvarib, ireals)
    nppfwr     = nvrfwr
 
! Pressure code information
! -------------------------
 
    IF (nobtyp == nsynop)              THEN
       npcode  = ibits (aofbuf(istart+nwrvrf(nvar( 1),nobtyp)),                &
                             npcdbp,nibits (npcdoc))
    ELSE
       rppp    = rppp * rfact (nvar  ( 1))
    ENDIF
  ENDIF
 
  !-----------------------------------------------------------------------------
  ! Section 3.0  Optional groups (surface analysis and verification)
  !-----------------------------------------------------------------------------

  CALL obs_option_groups
! ======================

! Weather and cloud groups
! ------------------------

! pre-set combined cloud and weather group, and combined flag word
  nclwgr  =  imdi
  ncwgfw  =  0

! Unpack wanted optional groups and write out
! -------------------------------------------
 
! rain group
! ----------

  Rain: DO
    rtr         =  rmdi
    ntrfl       =  0
    nrrfl       =  0
    IF (nraing  == nmdi)                                               EXIT Rain
 
! time range of precipitation observation (tr)
! ------------------------------------------
 
    iexv        =  ibits (nraing,ntrbp ,nibits( ntroc))
    IF ((nraine(1) /= 0) .AND. (iexv /= nibits( ntroc))) THEN
      ntrfl       =  ibits (nranfw,ntrfbp,nibits(nf2boc))
      IF (ntrfl  >=  2)                                                EXIT Rain
      rtr         =  iexv * rfact (16)
      IF (rtr <  rraint(1) .OR. rtr > rraint(2))                       EXIT Rain
    ENDIF

! precipitation amount (rr)
! -------------------------
 
    iexv        =  ibits (nraing,nrrbp ,nibits( nrroc))
    IF ((nraine(2) == 0) .OR. (iexv == nibits( nrroc)))                EXIT Rain
    IF (rtr     ==  rmdi .AND. iexv == 0)           THEN
      IF (ntime == 0 .OR. ntime == 1200)        THEN
         rtr =  6.0_ireals
      ELSE IF (ntime == 600 .OR. ntime == 1800) THEN
         rtr = 12.0_ireals
      ELSE
                                                                       EXIT Rain
      ENDIF
    ENDIF
    IF (rtr     ==  rmdi)                                              EXIT Rain
    irrflo      =  ibits (nranfw,nrrfbp,nibits(nf2boc))
    nrrfl       =  MAX(irrflo,ntrfl )
    IF (((nrrfl >=  2) .OR. (rtr  /=  raintp)) .AND. (.NOT. lverif))   EXIT Rain
    rr          =  iexv  *rfact (15)
    IF (iexv   == 9999) rr = 0.06_ireals
    roberr      =  1.0_ireals
!   IF((nrrfl  >=  2) .OR. (ABS( rtr-raintp ) > epsy))  roberr = rmdi
    IF (nrrfl  >=  2)  roberr = rmdi
    IF (rr      >   MIN( rrainl*rtr/6._ireals , 200._ireals ))  THEN
!      precipitation amount exceeds check limit
       IF (lactrep)                                                            &
         neventd (nerlim,ncdtpp,nobtpp) = neventd(nerlim,ncdtpp,nobtpp) + 1
       ilen = 4 + istrej
       IF (nacout+ilen <= nmxoln) THEN
         outbuf(nacout+1) = ilen
         outbuf(nacout+2) = nfmt2
         outbuf(nacout+3) = INT(rr*100._ireals)
         outbuf(nacout+4) = rtr
         DO icl = 1 , istrej
           outbuf(nacout+4+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
         ENDDO
         nacout  = nacout + ilen
       ENDIF
                                                                       EXIT Rain
    ENDIF
    ldatpsv     = .TRUE.
    IF (nrrfl < 2)  ldat = .TRUE.
!   IF((nrrfl < 2) .AND. (ABS( rtr-raintp ) <= epsy))  ldat = .TRUE.
 
! write precipitation to ODR
! --------------------------
 
    osgbdy (nsgob , nbspr ) = rr
    osgbdy (nsgob , nbsper) = roberr
    ncwgfw  =  insert ( ncwgfw , nrrfl , nvrfbp )
    nrtr    =  0
    IF (NINT( rtr ) ==  6)  nrtr = 1
    IF (NINT( rtr ) == 12)  nrtr = 2
    IF (NINT( rtr ) == 18)  nrtr = 3
    IF (NINT( rtr ) == 24)  nrtr = 4
    IF (NINT( rtr ) ==  1)  nrtr = 5
    IF (NINT( rtr ) ==  2)  nrtr = 6
    IF (NINT( rtr ) ==  3)  nrtr = 7
    ncwgfw  =  insert ( ncwgfw , nrtr  , nvtrbp )
 
    EXIT Rain
  ENDDO Rain

! Weather group
! -------------

! visibility
! ----------

  lfog       =  .FALSE.
  iexv       =  ibits (nwwgr ,nvibp ,nibits(nvioc))
  IF ((nwwgr /= nmdi) .AND. (nwwex(3) > 0) .AND. (iexv /= nibits(nvioc))) THEN
    nvisfl     =  ibits (nwwfw ,nvfbp ,nibits(nf2boc))
    IF ((nvisfl <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nvisfl <  2)  ldat = .TRUE.
 
! write visibility to ODR
      osgbdy (nsgob , nbsvis) = REAL ( MIN( iexv , 99999 ), ireals )
      ncwgfw  =  insert ( ncwgfw , nvisfl , nvvfbp )
    ENDIF

! pre-set flag for fog
    lfog  =  (nvisfl < 2) .AND. (REAL( iexv, ireals ) < vfoglim) .AND. (iexv > 0)
  ENDIF

! present weather --> fog
! -----------------------

  npww       =  ibits (nwwgr ,nwwbp  ,nibits(nwwoc))
  IF (npww  >=  100)  npww  =  nibits(nwwoc)
  IF ((nwwgr /= nmdi) .AND. (nwwex(2) > 0) .AND. (npww /= nibits(nwwoc))) THEN
    npwwfl     =  ibits (nwwfw ,nwwfbp,nibits(nf2boc))
    IF ((npwwfl <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (npwwfl <  2)  ldat = .TRUE.
 
! write present weather to buffer of ODR
      nclwgr  =  0
      nclwgr  =  insert ( nclwgr , npww   , nvwwbp )
      ncwgfw  =  insert ( ncwgfw , npwwfl , nvwfbp )
    ENDIF

! set flag for fog
!   lfog  =        (npwwfl < 2)                                                &
!            .AND. (     (npww == 28) .OR. (npww == 43) .OR. (npww == 45)      &
!                                     .OR. (npww == 47) .OR. (npww == 49)      &
!                   .OR. ((lfog) .AND. ((npww/10 == 2) .OR. (npww >= 50))))
    IF ((npwwfl < 2) .AND. (lfog) .AND. (npww/10 /= 4)) THEN
      lfog  =        (npww >= 50) .OR. (npww == 28) .OR. (npww == 17)          &
               .AND. (     (      ((npww <= 71) .OR. (npww >= 76))             &
                            .AND. (npww /= 82) .AND. (npww /= 84)              &
                            .AND. (npww /= 86) .AND. (npww /= 88)              &
                            .AND. (npww /= 90) .AND. (npww /= 94)              &
                            .AND. (npww <= 97))                                &
                      .OR. (osgbdy(nsgob,nbsvis) < 100.1_ireals))
      ilen = 6 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt20
        outbuf(nacout+3) = npww
        outbuf(nacout+4) = NINT( osgbdy(nsgob,nbsvis) )
        outbuf(nacout+5) = npwwfl
        outbuf(nacout+6) = nvisfl
        DO icl = 1 , istrej
          outbuf(nacout+6+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF
    ELSE
      lfog  =  (npwwfl < 2) .AND. (     (npww == 43) .OR. (npww == 45)         &
                                   .OR. (npww == 47) .OR. (npww == 49))
    ENDIF
  ELSE
    npwwfl = 0
    lfog   =  .FALSE.
  ENDIF

 
! General cloud group
! -------------------

! total cloud cover (only the complete general cloud group is written to ODR)
! -----------------
 
  ncctfl     =  0
  ncct       =  ibits (ngclgr,nbp   ,nibits(noc))
  IF ((ncct >= 10) .OR. (ngclex(6) == 0) .OR. (ngclgr == nmdi)) THEN
    ncct  =  nibits(noc)
  ELSE
    ncctfl     =  ibits (ngclfw,nnfbp ,nibits(nf2boc))
    IF ((ncctfl <  2) .OR. (lverif))  ldatpsv = .TRUE.
    IF  (ncctfl <  2)  ldat = .TRUE.
  ENDIF 
 
! cloud base height
! -----------------
 
! get h (height of cloud base)
  iexv       =  ibits (ngclgr,nhbp  ,nibits(nhoc))
  iexv       = iexv * 10

! convert unit of cloud base height from metres to key (VUB WMO Code table 1600)
  nhkey      =  nibits(nvhoc)
  IF (ngclgr /= nmdi) THEN
    SELECT CASE (iexv)
    CASE (   0:  49)
      nhkey =  0
    CASE (  50:  99)
      nhkey =  1
    CASE ( 100: 199)
      nhkey =  2
    CASE ( 200: 299)
      nhkey =  3
    CASE ( 300: 599)
      nhkey =  4
    CASE ( 600: 999)
      nhkey =  5
    CASE (1000:1499)
      nhkey =  6
    CASE (1500:1999)
      nhkey =  7
    CASE (2000:2499)
      nhkey =  8
    CASE (2500:16390)
      nhkey =  9
    CASE (16391)
      nhkey =  9
    CASE (16392)
      nhkey = 10
    END SELECT
  ENDIF

! low cloud cover
! ---------------
 
! get Nh (low or (,if low cloud not present,) middle cloud cover)
  nccl       =  ibits (ngclgr,nnbp  ,nibits(nnoc))
  IF (nccl  >=  10)  nccl  =  nibits(nnoc)

! IF (      (ngclgr /= nmdi) .AND. (ngclex(5) > 0) .AND. (ngclex(3) > 0))      &
!     .AND. (     ((nccl  /= nibits(nnoc)) .AND. (nhkey /= nibits(nvhoc)))     &
!            .OR. (ncct == 0) .OR. (ncct == 9))) THEN

! derive low cloud cover
  IF ((ngclgr /= nmdi) .AND. (ngclex(5) > 0) .AND. (ngclex(3) > 0)) THEN
    ncclfl  =  imdi
    ncclqf  =  imdi
! if 'N'=0 then 'Nh' and 'h' are not reported: cloudless
    IF (ncct == 0) THEN
      nccl    =  0
      ncclfl  =  ncctfl
      ncclqf  =  ncclfl
! if 'N'=9 then 'Nh' and 'h' are not reported: sky invisible, + fog flag
    ELSEIF ((ncct == 9) .AND. (lfog)) THEN
      nccl    =  8
      ncclfl  =  MAX( ncctfl , npwwfl )
      ncclqf  =  ncclfl + 1
      ilen = 5 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+1) = ilen
        outbuf(nacout+2) = nfmt21
        outbuf(nacout+3) = nccl
        outbuf(nacout+4) = ncclfl
        outbuf(nacout+5) = ncclqf
        DO icl = 1 , istrej
          outbuf(nacout+5+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF
! if 'N'=15=/ then 'Nh' and 'h' are not reported: fog flag only
!   ELSEIF ((ncct == nibits(noc)) .AND. (lfog)) THEN
!     nccl    =  8
!     ncclfl  =  MAX( 1 , npwwfl )
!     ncclqf  =  npwwfl + 2
! if 'N'>0, 'N'<9 then 'Nh' and 'h' are expected: if 'Nh' and 'h' are reported
    ELSEIF ((nccl  /= nibits(nnoc)) .AND. (nhkey /= nibits(nvhoc))) THEN
      nflag     =  ibits (ngclfw,nnnfbp,nibits(nf2boc))
      nflag1    =  ibits (ngclfw,nnhfbp,nibits(nf2boc))
! cloud base height below 2000 m
      IF      (nhkey <= 7) THEN
        ncclfl  =  MAX( nflag , nflag1 )
        ncclqf  =  ncclfl
! cloud base height above 2000 m
      ELSEIF ((nhkey == 8) .OR. (nhkey == 9)) THEN
        nccl    =  0
        ncclfl  =  nflag1
        ncclqf  =  ncclfl
! cloud base below station or invisible, + fog flag
      ELSEIF ((nhkey == 10) .AND. (lfog)) THEN
        IF (nccl == 9)  nccl  =  8
        ncclfl  =  MAX( nflag , nflag1 , npwwfl )
        ncclqf  =  ncclfl + 1
        ilen = 5 + istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt22
          outbuf(nacout+3) = nccl
          outbuf(nacout+4) = ncclfl
          outbuf(nacout+5) = ncclqf
          DO icl = 1 , istrej
            outbuf(nacout+5+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
! cloud base below station or invisible, no fog flag, 'nh' determined
      ELSEIF ((nhkey == 10) .AND. (nccl <= 8) .AND. (nclwgr /= imdi)) THEN
        ncclfl  =  MAX( nflag , nflag1 , npwwfl , 1 )
!       ncclqf  =  MAX( nflag , nflag1 , npwwfl ) + 2
        ncclqf  =  MAX( nflag , nflag1 , npwwfl ) + 1
        ilen = 5 + istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt23
          outbuf(nacout+3) = nccl
          outbuf(nacout+4) = ncclfl
          outbuf(nacout+5) = ncclqf
          DO icl = 1 , istrej
            outbuf(nacout+5+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
    ENDIF
    IF ((ncclfl <  2) .OR. ((lverif) .AND. (ncclfl /= imdi))) THEN
      roberr  =  1.0_ireals
      IF (ncclqf == 1)  roberr = 2.0_ireals
      IF (ncclqf == 2)  roberr = 3.0_ireals
      IF (ncclqf == 3)  roberr = 4.0_ireals
      IF (ncclfl >= 2)  roberr = rmdi
      ldatpsv    = .TRUE.
      IF (ncclfl <  2)  ldat = .TRUE.

! write low cloud cover to ODR, and include low cloud flag in combined flag word
      osgbdy (nsgob , nbscl ) = REAL ( nccl, ireals )
      osgbdy (nsgob , nbscer) = roberr
      ncwgfw  =  insert ( ncwgfw , ncclfl , nvcfbp )
      ncwgfw  =  insert ( ncwgfw , ncclqf , nvqfbp )
    ENDIF
  ENDIF 

! combined cloud and weather group written to ODR
! -----------------------------------------------

  IF ((ngclgr /= nmdi) .AND. (MAX( ngclex(6), ngclex(5), ngclex(3) ) > 0)) THEN

! replace cloud base height in metres by cloud base height key in combined group
    IF (nclwgr == imdi)  nclwgr  =  0
!   nbufcl  =  ibits ( ngclgr , nchbp , nchoc + ncmoc )
!   nclwgr  =  insert ( nclwgr , nbufcl , nvchbp )
!   nclwgr  =  insert ( nclwgr , nhkey  , nvhbp  )
!   nbufcl  =  ibits ( ngclgr , nclbp , ncloc + nnoc  + noc )
!   nclwgr  =  insert ( nclwgr , nbufcl , nvclbp )
    nclwgr  =  ireplace( nclwgr , nvchbp , nchoc + ncmoc , ngclgr , nchbp )
    nclwgr  =  insert ( nclwgr , nhkey  , nvhbp  )
    nclwgr  =  ireplace( nclwgr , nvclbp , ncloc + nnoc + noc , ngclgr , nclbp )

! include general cloud group flags in combined flag word
!   nbufcl  =  ibits ( ngclfw , nchfbp, 6* nv2foc )
!   ncwgfw  =  insert ( ncwgfw , nbufcl , nvufbp )
    ncwgfw  =  ireplace( ncwgfw , nvufbp , 6* nv2foc , ngclfw , nchfbp )

! write combined cloud and weather group to ODR
    mosgbd (nsgob , nbscwg) = nclwgr
  ENDIF


! Ground group
! ------------

! ground temperature (minimum temperature during past 12 hrs 5 cm above surface)
! ------------------
 
  iexv       =  ibits (ngrngr,ntgbp ,nibits(ntgoc))
  IF ((ngrngr /= nmdi) .AND. (ngrnex(3) > 0) .AND. (iexv /= nibits(ntgoc))) THEN
    nflag      =  ibits (ngrnfw,ntgfbp,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write ground temperature to ODR
      osgbdy (nsgob , nbstg ) = iexv * rfact(nvar(4))
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvgfbp )
    ENDIF
  ENDIF


! Extreme temperature group
! -------------------------

! minimum temperature (during past 12 hrs 2 m above surface)
! -------------------
 
  iexv       =  ibits (nextgr,ntnbp ,nibits(ntnoc))
  IF ((nextgr /= nmdi) .AND. (nextex(1) > 0) .AND. (iexv /= nibits(ntnoc))) THEN
    nflag      =  ibits (nextfw,nttnbp,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write minimum temperature to ODR
      osgbdy (nsgob , nbstn ) = iexv * rfact(nvar(4))
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvtfbp )
    ENDIF
  ENDIF

! maximum temperature (during past 12 hrs 2 m above surface)
! -------------------
 
  iexv       =  ibits (nextgr,ntxbp ,nibits(ntxoc))
  IF ((nextgr /= nmdi) .AND. (nextex(2) > 0) .AND. (iexv /= nibits(ntxoc))) THEN
    nflag      =  ibits (nextfw,nttxbp,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write maximum temperature to ODR
      osgbdy (nsgob , nbstx ) = iexv * rfact(nvar(4))
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvxfbp )
    ENDIF
  ENDIF


! Special phenomena group
! -----------------------

! maximum wind speed of 10 minute mean wind (during past weather period)
! -----------------------------------------
 
  iexv       =  ibits (nsppgr(1),nlspbp,nibits(nlspoc))
  IF (      (nsppgr(1) /= nmdi) .AND. (nsppex(1) > 0)                          &
      .AND. (iexv /= nibits(nlspoc))) THEN
    nflag      =  ibits (nsppfw(1),nlspfb,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write maximum 10' mean wind speed to ODR
      osgbdy (nsgob , nbsfme) = REAL ( iexv, ireals )
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvefbp )
    ENDIF
  ENDIF

! maximum wind speed of gusts (during past weather period)
! ---------------------------
 
  iexv       =  ibits (nsppgr(1),nuspbp,nibits(nuspoc))
  IF (      (nsppgr(1) /= nmdi) .AND. (nsppex(2) > 0)                          &
      .AND. (iexv /= nibits(nuspoc))) THEN
    nflag      =  ibits (nsppfw(1),nuspfb,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write maximum gust wind speed to ODR
      osgbdy (nsgob , nbsfgu) = REAL ( iexv, ireals )
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvffbp )
    ENDIF
  ENDIF


! Radiation group
! ---------------

! global radiation, sum over past hour
! ------------------------------------
 
  iexv       =  ibits (nradgr,nffbp ,nibits(nffoc))
  IF ((nradgr /= nmdi) .AND. (nradex(1) > 0) .AND. (iexv /= nibits(nffoc))) THEN
    nflag      =  ibits (nradfw,nfffbp,nibits(nf2boc))
    IF ((nflag  <  2) .OR. (lverif)) THEN
      ldatpsv    = .TRUE.
      IF (nflag  <  2)  ldat = .TRUE.
 
! write global radiation to ODR
      osgbdy (nsgob , nbsrad) = iexv * 10000.0_ireals
      nflag   =  nflag / 2
      ncwgfw  =  insert ( ncwgfw , nflag , nvsfbp )
    ENDIF
  ENDIF


! write combined cloud, weather, precipitation and temperature flag word to ODR
! -----------------------------------------------------------------------------

  mosgbd (nsgob , nbscfw) = ncwgfw


  !-----------------------------------------------------------------------------
  ! Section 4.0  Single level report (wind&mass/rel. hum. analysis)
  !-----------------------------------------------------------------------------

!  Check if SYNOP which is assigned to a land grid pt.,
!  and get height [m] difference to model orography
 
  lseaobs = (ibits(mosghd(nsgob,nhschr),nvsebp,nvscoc) == 1)
  landsy  = (nobtyp == nsynop) .AND. (.NOT.lseaobs)
! fisd    = osghed(nsgob,nhsurf) - osghed(nsgob,nhalt)
  fisd    = osghed(nsgob,nhalt) - osghed(nsgob,nhsurf)
 
! height: is derived directly from station height
! =======          or from pressure report (SYNOP)
!                  or from height report (AIREP)
 
  rzzz        =   rmdi
  zadder      =   0._ireals
  izzfwr      =   0
 
! SYNOP case:
! -----------
! prepare pressure or height
 
  IF (nobtyp == nsynop) THEN
     zstalt   = osghed (nsgob,nhalt )
 
! assign zzzz (and change rppp) according to pressure code
 
     IF (lpprep) THEN
       ipcode =  npcode + 1
       zzzz   =  rppp          * rfact (nvar  ( 7)) * nzfact(ipcode) +         &
                 rzset(ipcode) * rfact (nvar  ( 7))
       rppp   =  rppp          * rfact (nvar  ( 1)) * npfact(ipcode) +         &
                 rpset(ipcode) * rfact (nvar  ( 1))
       izzfwr =  nppfwr                             * nzflfc(ipcode)
       nppfwr =  nppfwr                             * npflfc(ipcode)
     ENDIF
 
! possible error in pmsl (if rppp  >   1060 hpa)
     IF (rppp  > rerrpp) THEN
       lpract = .TRUE.
     ENDIF

! if pressure is reported at station height
     IF (npcode == 1) THEN
        rzzz       =  zstalt
     ELSEIF (lpprep)  THEN

! if pressure is not reported at station height
        rzzz       =  zzzz

! check station reporting practice
! --------------------------------
!        bad reporting practice is reckoned to be if:
!         -station altitude above 800m and reports pmsl
!         -station altitude above 800m and reports 1000mb level
!         -station altitude above 1700m and reports 900mb level
!         -station altitude below 300m and reports 900mb level
!         -station altitude above 2300m and reports 850mb level
!         -station altitude below 800m and reports 850mb level
!         -station altitude above 3700m and reports 700mb level
!         -station altitude below 2300m and reports 700mb level
!         -station altitude below 3700m and reports 500mb level

        IF(npcode ==  0 .AND. zstalt > rh800 )  lpract = .TRUE.
        IF(npcode == 10 .AND. zstalt > rh800 )  lpract = .TRUE.
        IF(npcode ==  9 .AND. zstalt > rh1700)  lpract = .TRUE.
        IF(npcode ==  9 .AND. zstalt < rh300 )  lpract = .TRUE.
        IF(npcode ==  2 .AND. zstalt > rh2300)  lpract = .TRUE.
        IF(npcode ==  2 .AND. zstalt < rh800 )  lpract = .TRUE.
        IF(npcode ==  3 .AND. zstalt > rh3700)  lpract = .TRUE.
        IF(npcode ==  3 .AND. zstalt < rh2300)  lpract = .TRUE.
        IF(npcode == 11 .AND. zstalt < rh3700)  lpract = .TRUE.
        IF(npcode ==  4 .AND. ABS(zstalt- 500._ireals) > rh800) lpract=.TRUE.
        IF(npcode ==  5 .AND. ABS(zstalt-1000._ireals) > rh800) lpract=.TRUE.
        IF(npcode ==  6 .AND. ABS(zstalt-2000._ireals) > rh800) lpract=.TRUE.
        IF(npcode ==  7 .AND. ABS(zstalt-3000._ireals) > rh800) lpract=.TRUE.
        IF(npcode ==  8 .AND. ABS(zstalt-4000._ireals) > rh800) lpract=.TRUE.
 
! neglect pressure obs., if extrapolation too large (cf. *ps_obs_sing*)
        IF (ABS(zstalt-rzzz) > doromx(2)) lpractn = .TRUE.
 
! additional observation error due to extrapolation
        zadder = 0.04_ireals * ABS( zstalt - rzzz )
! convey height error to pressure error (rough approximation)
        zadder = zadder * rppp / (r_d * t0_melt)
     ENDIF
  ENDIF
 
! DRIBU case
! ----------
 
  IF (nobtyp == ndribu) rzzz  =  0._ireals
 
! AIREP / SATOB case: 
! ------------------

  IF ((nobtyp == nairep) .OR. (nobtyp == nsatob)) THEN
    IF (.NOT. lpprep) THEN
!   in case of no pressure get observed height
      CALL obs_extract (istart,nvar( 7),aofbuf,nmxlob)
!     PRINT *,ystid," nvarib_zzz:",nvarib
      IF ((nvarib == nmdi) .OR. (nvarib == 0) .OR. (nflag >= 2)) THEN
!       reject AIREP without pressure and height
        IF (lactrep)                                                           &
          neventr(nenops,ncdtpp,nobtpp) = neventr(nenops,ncdtpp,nobtpp) + 1
        ilen = 2 + istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt1
          DO icl = 1 , istrej
            outbuf(nacout+2+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
        lrej = .TRUE.
!       (set (ldat) so that the IF-loop after 'ENDDO Extract_data' is omitted)
        ldat = .TRUE.
        EXIT  Extract_data
      ENDIF
      rzzz = REAL(nvarib, ireals)
      IF (rzzz < zsurf) THEN
        neventr (nenops,ncdtpp,nobtpp) = neventr(nenops,ncdtpp,nobtpp) + 1
        lrej = .TRUE.
        ldat = .TRUE.
        EXIT  Extract_data
      ENDIF
 
!   AIREP:  calculate pressure from height according to standard atmosphere
      IF (nobtyp == nairep) THEN
!       rzzz_grb = REAL(rzzz, irealgrib)
!       CALL atmhp ('single', irm, 'US1976', rzzz_grb, rppp_grb )
!       rppp = REAL(rppp_grb,ireals)

        CALL std_atmosphere ( 'z', rzzz, 'p', rppp, g, r_d, irm )
!       ===================

!   SATOB:  calculate pressure from height according to model atmosphere
      ELSEIF (nobtyp == nsatob) THEN
        rppp  =  z2p ( rzzz , iob , job )
!                ===
      ENDIF
!     print *,ystid," rzzz=",rzzz,"  rppp=",rppp
      IF (irm /= 0)   THEN
        IF (lactrep)                                                           &
          neventr (nenops,ncdtpp,nobtpp) = neventr (nenops,ncdtpp,nobtpp) + 1
        lrej = .TRUE.
        ldat=.TRUE.
        EXIT  Extract_data
      ENDIF
      rppp = rppp * 100._ireals
      lpprep = .true.
      nppfwr = 0
    ENDIF
    rzzz = rmdi
  ENDIF

! set aircraft rep. to passive if less than  50m or  6hPa  \  above (model)
! set  SATOB   rep. to passive if less than 200m or 20hPa  /    orography
! -------------------------------------------------------------------------
 
  IF ((nobtyp == nairep) .OR. (nobtyp == nsatob)) THEN
    mpassiv  =  0
    IF (nobtyp == nairep) THEN
      IF ((rzzz > rmdich) .AND. (rzzz < zsurf+50._ireals))  mpassiv = 2
      IF ((rzzz < rmdich) .AND. (rppp > psurf(iob,job) -600._ireals))          &
          mpassiv = 2
    ELSEIF (nobtyp == nsatob) THEN
      IF ((rzzz > rmdich) .AND. (rzzz < zsurf+200._ireals))  mpassiv = 2
      IF ((rzzz < rmdich) .AND. (rppp > psurf(iob,job) -2000._ireals))         &
          mpassiv = 2
    ENDIF
    IF ((lactrep) .AND. (mpassiv == 2)) THEN
      neventr (nezdif,ncdtpp,nobtpp) = neventr (nezdif,ncdtpp,nobtpp) + 1
      noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
      mosghd (nsgob,nhpass)  =  2
      lactrep                = .FALSE.
    ENDIF
    rzzz = rmdi
  ENDIF

  IF ((rppp > rmdich) .AND. (rzzz > rmdich)) ldat = .TRUE.

! use no report above 'rpplim'
! ----------------------------
 
  IF ((rppp < rpplim) .AND. (rppp > rmdich)) THEN
    IF (lactrep)                                                               &
      neventr (nenops,ncdtpp,nobtpp) = neventr(nenops,ncdtpp,nobtpp) + 1
    lrej = .TRUE.
    EXIT  Extract_data
  ENDIF
 
! height observation error
! ------------------------
 
  lno  =  .FALSE.
  IF ((nobtyp == nairep) .OR. (nobtyp == nsatob)) THEN
    roberr = rmdi
  ELSE 
    lseaonly = (ABS(altopsu(2)) < epsy)
    fisdzz = ((fdoro(2)-c1)/c2 + SIGN( (fdoro(2)+c1)/c2 , fisd )) * fisd
    IF (rzzz > rmdich)                                                         &
      fisdzz = fisdzz + ABS( osghed(nsgob,nhalt) -rzzz )
    IF ((.NOT. lpprep) .OR. (nzot(nobtyp) == 0) .OR. (nzex(nobtyp) == 0)) THEN
       roberr = rmdi
    ELSEIF (lpract) THEN
       roberr = rmdi
       IF (lactrep)                                                            &
         neventd (neprac,ncdtpp,nobtpp) = neventd(neprac,ncdtpp,nobtpp) + 1
    ELSEIF ((lpractn) .OR. ((.NOT.lseaonly)   .AND. (rzzz >  altopsu(2)))      &
                      .OR. ((     lseaonly)   .AND. (landsy))                  &
                      .OR. ((nobtyp == nsynop).AND. (fisdzz > doromx(2)))) THEN
       lno    = .TRUE.
       roberr = rmdi
       IF (lactrep)                                                            &
         neventd (nepalt,ncdtpp,nobtpp) = neventd(nepalt,ncdtpp,nobtpp) + 1
    ELSE
       roberr = obs_error ( nvrpoe(3) )
       roberr = roberr + zadder
    ENDIF
  ENDIF
 
! write pressure and height to ODR
! --------------------------------
 
  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nppfwr+izzfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nppfwr+izzfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nppfwr+izzfwr, nf1bps(3) )
  IF  (lno)     nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
  IF  (lpract)  nvarfw  =  insert( nvarfw, 1, nvfbps(5) )
  osgbdy (nsgob,nbsp  ) = rppp
  osgbdy (nsgob,nbsz  ) = rzzz
  osgbdy (nsgob,nbszer) = roberr
! mosgbd (nsgob,nbspfl) = nppfwr + izzfwr
  mosgbd (nsgob,nbsflg) = insert( 0, nvarfw, nvfzbp )
  mosgbd (nsgob,nbslid) = insert( 0, 1, nvlidp(7) )
  IF (nobtyp == nairep) mosgbd (nsgob,nbslid) = 0
  IF (nobtyp == nsatob) mosgbd (nsgob,nbslid) = 0
  IF (nobtyp == nsynop) mosgbd (nsgob,nbslid) = npcode

  mosgbd (nsgob,nbsqcf) = 0
 
! pressure tendency (for surface pressure quality control threshold only)
! -----------------
 
  IF (nobtyp == nsynop) THEN
    CALL obs_extract (istart,nvar( 6),aofbuf,nmxlob)
    IF (nvarib /= nmdi) THEN
      IF ((nflag < 2) .AND. (ABS(nvarib-500) < 400)) THEN
        ldatpsv    = .TRUE.
        osgbdy (nsgob,nbspst) = REAL ( (nvarib - 500) * 10, ireals )
      ELSEIF (lactrep) THEN
        neventd (nepsdt,ncdtpp,nobtpp) = neventd(nepsdt,ncdtpp,nobtpp) + 1
        ilen = 4 + istrej
        IF (nacout+ilen <= nmxoln) THEN
          outbuf(nacout+1) = ilen
          outbuf(nacout+2) = nfmt6
          outbuf(nacout+3) = nflag
          outbuf(nacout+4) = nvarib - 500
          DO icl = 1 , istrej
            outbuf(nacout+4+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
          ENDDO
          nacout  = nacout + ilen
        ENDIF
      ENDIF
    ENDIF
  ENDIF

 
! u,v and 10m u,v components
! --------------------------
 
  Wind: DO
! ========

  mpassiv  =  0
 
  IF ((nuex(nobtyp) == 0) .OR. (nvex(nobtyp) == 0))                    EXIT Wind
  CALL obs_extract (istart,nvar( 2),aofbuf,nmxlob)
  zddd     =   REAL  ( nvarib, ireals )
! if (ncdtyp==244) then
!   print *,ystid," nvarib=",nvarib," zddd=",zddd
! endif
  ivrfwr   =   nvrfwr
  nvarib1  =   nvarib
  nflag1   =   nflag
  CALL obs_extract (istart,nvar( 3),aofbuf,nmxlob)
  IF (lactrep) THEN
    IF ((nvarib1 == nmdi) .AND. (nvarib /= 0))                                 &
                         neventd (nedmis,ncdtpp,nobtpp) =                      &
                         neventd (nedmis,ncdtpp,nobtpp) + 1
    IF (nvarib  == nmdi) neventd (nefmis,ncdtpp,nobtpp) =                      &
                         neventd (nefmis,ncdtpp,nobtpp) + 1
    IF ((nflag1 >= 2   ) .OR. ((ABS(nvarib1) > 360) .AND. (nvarib1 /= nmdi)    &
                                                    .AND. (nvarib /= 0)))      &
                         neventd (nedflg,ncdtpp,nobtpp) =                      &
                         neventd (nedflg,ncdtpp,nobtpp) + 1
    IF  (nflag  >= 2   ) neventd (nefflg,ncdtpp,nobtpp) =                      &
                         neventd (nefflg,ncdtpp,nobtpp) + 1
  ENDIF
! if (ncdtyp==244) then
!   print *,ystid," nvarib_ff=",nvarib," nflag1,nflag:",nflag1,nflag
! endif

  IF ((((nvarib1 == nmdi) .OR. (ABS(nvarib1) > 360)) .AND. (nvarib /= 0))      &
                          .OR. ((nflag1 >= 2) .AND. (.NOT. lverif)))   EXIT Wind
  IF   ((nvarib  == nmdi) .OR. ((nflag  >= 2) .AND. (.NOT. lverif)))   EXIT Wind
  IF ((nflag1 >= 2) .OR. (nflag >= 2))  mpassiv = 1
  zfff     =    REAL  ( nvarib, ireals )
  IF  (      (zfff < -epsy)                                                    &
       .OR. ((zfff >  90.0_ireals) .AND. (     (nobtyp == nsynop)              &
                                          .OR. (nobtyp == ndribu)))            &
       .OR. ((zfff > 150.0_ireals) .AND. (     (nobtyp == nairep)              &
                                          .OR. (nobtyp == nsatob))) ) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        & 
      neventd (nefneg,ncdtpp,nobtpp) = neventd(nefneg,ncdtpp,nobtpp) + 1
                                                                       EXIT Wind
  ENDIF
  zfff     =    MAX( zfff , c0 )
  IF   (zfff < epsy) zddd = c0
  nuvfwr   =    IOR (nvrfwr, ivrfwr)
  IF ((nuot(nobtyp) == 0) .OR. (nvot(nobtyp) == 0))                    EXIT Wind

! transformation of wind speed and direction to wind components
! in the rotated model coordinate system
! -------------------------------------
  zfff  = zfff * rfact(nvar( 3))
  rsus  = zfff * (-SIN( zddd*degrad ))
  rsvs  = zfff * (-COS( zddd*degrad ))

  CALL uv2uvrot ( rsus,rsvs, roblat,roblon, pollat,pollon, ruuu,rvvv )
! =============

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nuvfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nuvfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nuvfwr, nf1bps(3) )
 
! check if 10m (u,v) wind is to be accepted:
! -----------------------------------------
! (observations assigned to land model grid pts. can be excluded
! by choosing 'altopsu(1) == 0.')
 
  IF ((nobtyp == nsynop) .OR. (nobtyp == ndribu))            THEN
     lyes = .FALSE.
     lno  = .FALSE.
     DO   j  = 1, mxda
       WRITE(yyes ,'(i5,3x)') muvyes( j )
       WRITE(yno  ,'(i5,3x)') muvno ( j )
       IF (ystid == yyes) lyes = .TRUE.
       IF (ystid == yno ) lno  = .TRUE.
     ENDDO
     fisduv = ((fdoro(1)-c1)/c2 + SIGN( (fdoro(1)+c1)/c2 , fisd )) * fisd
     lseaonly = (ABS(altopsu(1))  <   epsy)
     IF (      (.NOT. lyes)                                                    &
         .AND. (     (.NOT. lseaonly .AND. (zstalt  >  altopsu(1)))            &
                .OR. (      lseaonly .AND. landsy)                             &
                .OR. (fisduv >  doromx(1))))   lno = .TRUE.
 
! write 10m u,v wind to ODR
! -------------------------
 
     IF ((lno) .OR. (mpassiv > 0)) THEN
       roberr  =  rmdi
       IF ((mpassiv == 0) .AND. (lactrep))                                     &
         neventd (nevalt,ncdtpp,nobtpp) = neventd(nevalt,ncdtpp,nobtpp) + 1
       IF (lverif)  ldatpsv = .TRUE.
     ELSEIF ((nobtyp == ndribu) .AND. (zfff < epsy)) THEN
       roberr  =  rmdi
       IF ((mpassiv == 0) .AND. (lactrep))                                     &
         neventd (nefneg,ncdtpp,nobtpp) = neventd(nefneg,ncdtpp,nobtpp) + 1
       IF (lverif)  ldatpsv = .TRUE.
       nvarfw  =  insert( nvarfw, 1, nvfbps(5) )
     ELSE
       roberr  =  obs_error ( nvrpoe(1) )
       ldat    = .TRUE.
     ENDIF
     IF  (lno)  nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
     osgbdy (nsgob,nbsu  )  =  ruuu
     osgbdy (nsgob,nbsv  )  =  rvvv
     osgbdy (nsgob,nbsuer)  =  roberr
     mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfubp )
  
  ELSEIF ((nobtyp == nairep) .OR. (nobtyp == nsatob))   THEN
 
! write aircraft (upper air) u,v wind to ODR
! ------------------------------------------
 
     IF (mpassiv == 0) roberr = obs_error ( nvrpoe(1) )
     IF (mpassiv >  0) roberr = rmdi
     osgbdy (nsgob,nbsu  )  =  ruuu
     osgbdy (nsgob,nbsv  )  =  rvvv
     osgbdy (nsgob,nbsuer)  =  roberr
     mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfubp )
     mosgbd (nsgob,nbslid)  =  insert( mosgbd(nsgob,nbslid), nvlido, nvlidp(8) )
     ldatpsv   = .TRUE.
     IF (mpassiv == 0) ldat   = .TRUE.

  ENDIF
  EXIT Wind
  ENDDO Wind


! Temperature (not 2m temperature)
! --------------------------------

  tem_hum    : DO
! ===============

  mpassiv  =  0

  rttt     =    rmdi
  IF (ntex (nobtyp) == 0)                                           EXIT tem_hum
  CALL obs_extract (istart,nvar( 4),aofbuf,nmxlob)
  IF ((nvarib == nmdi) .AND. (lactrep)) neventd (netmis,ncdtpp,nobtpp) =       &
                                        neventd (netmis,ncdtpp,nobtpp) + 1
  IF ((nflag  >= 2   ) .AND. (lactrep)) neventd (netflg,ncdtpp,nobtpp) =       &
                                        neventd (netflg,ncdtpp,nobtpp) + 1
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif)))    EXIT tem_hum
  IF (nflag >= 2) mpassiv = 1
  nttfwr   =    nvrfwr
  rttt     =    nvarib * rfact(nvar  ( 4))
  IF (ntot (nobtyp) == 0)                                           EXIT tem_hum
  IF (      (rttt < t0_melt-90._ireals) .OR. (rttt > t0_melt+60._ireals)       &
      .OR. ((rttt > t0_melt+20) .AND. (rppp < 70000._ireals))                  &
      .OR. ((rttt > t0_melt+ 5) .AND. (rppp < 50000._ireals))                  &
      .OR. ((rttt > t0_melt- 5) .AND. (rppp < 40000._ireals))) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (netext,ncdtpp,nobtpp) = neventd(netext,ncdtpp,nobtpp) + 1
                                                                    EXIT tem_hum
  ENDIF
 
! write upper-air (aircraft) temperature to ODR
! ---------------------------------------------
 
  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nttfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nttfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nttfwr, nf1bps(3) )
 
  IF (mpassiv == 0) roberr = obs_error ( nvrpoe(6) )
  IF (mpassiv >  0) roberr = rmdi
  osgbdy (nsgob,nbst  )  =  rttt
  osgbdy (nsgob,nbster)  =  roberr
  mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvftbp )
  mosgbd (nsgob,nbslid)  =  insert( mosgbd(nsgob,nbslid), nvlido, nvlidp(9) )
  ldatpsv   = .TRUE.
  IF (mpassiv == 0) ldat   = .TRUE.
 
! dew-point and relative humidity (not 2m dew-point)
! --------------------------------------------------
 
  rtdt     =    rmdi
  IF (ntdex (nobtyp) == 0)                                          EXIT tem_hum
  CALL obs_extract (istart,nvar( 5),aofbuf,nmxlob)
  IF ((mpassiv == 0) .AND. (lactrep)) THEN
    IF (nvarib == nmdi) neventd (neqmis,ncdtpp,nobtpp) =                       &
                        neventd (neqmis,ncdtpp,nobtpp) + 1
    IF (nflag  >= 2   ) neventd (neqflg,ncdtpp,nobtpp) =                       &
                        neventd (neqflg,ncdtpp,nobtpp) + 1
  ENDIF
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif)))    EXIT tem_hum
  IF (nflag >= 2) mpassiv = 1
  ntdfwr   =    nvrfwr
  rtdt     =    nvarib * rfact (nvar  ( 5))
  IF (ntdot (nobtyp) == 0)                                          EXIT tem_hum
  IF ((rtdt <  t0_melt-150._ireals) .AND. (mpassiv == 0) .AND. (lactrep))      &
     neventd (neqlow,ncdtpp,nobtpp) = neventd (neqlow,ncdtpp,nobtpp) + 1
  IF  (rtdt <  t0_melt-150._ireals)                                 EXIT tem_hum
  IF ((rppp <  pqmin )         .AND. (mpassiv == 0) .AND. (lactrep))           &
     neventd (neq300,ncdtpp,nobtpp) = neventd (neq300,ncdtpp,nobtpp) + 1
  IF ((rppp <  pqmin )         .AND. (.NOT. lverif))                EXIT tem_hum
  IF  (rppp <  pqmin )  mpassiv = MAX( mpassiv , i1 )

! convert humidity to model compatible relative humidity, and write to ODR
! ------------------------------------------------------------------------

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), ntdfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), ntdfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), ntdfwr, nf1bps(3) )
  IF (rppp  <  pqmin )  nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
 
! eobs : observed vapour press. (measurement is inherently over ice)
!       (but td as stored in data file is over water (WMO convention),
!       which is why magnus formula over water instead of ice is used)
  eobs = fspw( rtdt )

! eswo : saturation vapour pressure over water at observed temperat.
  eswo = fspw( rttt )

  IF ((rttt < t0_melt) .AND. (itype_gscp <= 2))   THEN

!   eseo : saturation vapour pressure over ice at observed temperature
    eseo = fspe( rttt )

!   rheo : observed relative humidity over ice   (per def.  <= 100%)
    rheo = eobs / eseo
    !! this is assumed equivalent to model relative humidity over
!     water (per def.  <= 100%; the model does not know ice) !!
!     ==> eobs-corr / eswo = eobs / eseo

    eobs = eswo * rheo
  ENDIF

  rheo = eobs / eswo
  rtdt = ftd( LOG( eobs / b1 ) )
  IF ((rheo >  rtshlm) .AND. (mpassiv == 0) .AND. (lactrep))                   &
     neventd(neqbig,ncdtpp,nobtpp) = neventd(neqbig,ncdtpp,nobtpp) + 1
  IF  (rheo >  rtshlm)                                              EXIT tem_hum
  IF (rheo >= rhtsat) THEN
     IF ((rheo <  1.0_ireals-epsy) .AND. (mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt <  t0_melt) neventd (neqsam,ncdtpp,nobtpp) =                   &
                       neventd (neqsam,ncdtpp,nobtpp) + 1
       IF (rttt >= t0_melt) neventd (neqsap,ncdtpp,nobtpp) =                   &
                       neventd (neqsap,ncdtpp,nobtpp) + 1
     ENDIF
     IF ((rheo >  1.0_ireals+epsy) .AND. (mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt <  t0_melt) neventd (neqclm,ncdtpp,nobtpp) =                   &
                       neventd (neqclm,ncdtpp,nobtpp) + 1
       IF (rttt >= t0_melt) neventd (neqclp,ncdtpp,nobtpp) =                   &
                       neventd (neqclp,ncdtpp,nobtpp) + 1
     ENDIF
     IF (rheo < c1-epsy)  nvarfw  =  insert( nvarfw, 1, nvfbps(5) )
     IF (rheo > c1+epsy)  nvarfw  =  insert( nvarfw, 1, nvfbps(6) )
     rheo  =  1.0_ireals
     rtdt  =  rttt
  ENDIF
 
  roberr   =    rherr1 / 100._ireals
  IF (rheo  <  rrhlim)      roberr    =  rherr2 / 100._ireals
  IF (rttt  <  rttlim)      roberr    =  rherr3 / 100._ireals
  IF (mpassiv > 0)  roberr = rmdi
  osgbdy (nsgob,nbsrh )  =  rheo
  osgbdy (nsgob,nbsqer)  =  roberr
  mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfqbp )
  mosgbd (nsgob,nbslid)  =  insert( mosgbd(nsgob,nbslid), nvlido, nvlidp(6) )
  ldatpsv   = .TRUE.
  IF (mpassiv == 0) ldat   = .TRUE.

  EXIT tem_hum
  ENDDO tem_hum

 
! 2m temperature
! --------------

  tem_hum_2m:  DO
! ===============

  mpassiv  =  0

  rttt2    =    rmdi
  IF (nt2ex (nobtyp) == 0)                                       EXIT tem_hum_2m
  CALL obs_extract (istart,nvar( 4),aofbuf,nmxlob)
  IF ((nvarib == nmdi) .AND. (lactrep)) neventd (netmis,ncdtpp,nobtpp) =       &
                                        neventd (netmis,ncdtpp,nobtpp) + 1
  IF ((nflag  >= 2   ) .AND. (lactrep)) neventd (netflg,ncdtpp,nobtpp) =       &
                                        neventd (netflg,ncdtpp,nobtpp) + 1
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif))) EXIT tem_hum_2m
  IF (nflag >= 2) mpassiv = 1
  nttfwr   =    nvrfwr
  rttt2    =    nvarib * rfact(nvar  ( 4))
  IF (nt2ot (nobtyp) == 0)                                       EXIT tem_hum_2m
  IF ((rttt2 < t0_melt-90._ireals) .OR. (rttt2 > t0_melt+60._ireals)) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (netext,ncdtpp,nobtpp) = neventd(netext,ncdtpp,nobtpp) + 1
                                                                 EXIT tem_hum_2m
  ENDIF

! write 2m temperature to ODR
! ---------------------------

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), nttfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), nttfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), nttfwr, nf1bps(3) )
 
  lno  = .FALSE.

! check station model height difference

  fisdtt = ((fdoro(3)-c1)/c2 + SIGN( (fdoro(3)+c1)/c2 , fisd )) * fisd
  lseaonly = (ABS(altopsu(3))  <   epsy)
  IF    (     (.NOT. lseaonly .AND. zstalt >  altopsu(3))                      &
         .OR. (      lseaonly .AND. landsy )                                   &
         .OR. (fisdtt >      doromx(3)                    ) ) lno = .TRUE.
!        .OR. (fisdtt > MAX( doromx(3), zstalt/5._ireals) ) ) lno = .TRUE.
  IF ((lno) .AND. (mpassiv == 0) .AND. (lactrep))                              &
     neventd (netalt,ncdtpp,nobtpp) = neventd(netalt,ncdtpp,nobtpp) + 1
  IF ((.NOT. lno) .OR. (lverif)) THEN
! note: extrapolation of T2m to the model surface level is done in 'sing_local'
     roberr = 1._ireals + 0.003_ireals *ABS(fisd)
     IF ((lno) .OR. (mpassiv > 0))  roberr = rmdi
     IF  (lno)  nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
     osgbdy (nsgob,nbst  )  =  rttt2
     osgbdy (nsgob,nbster)  =  roberr
     mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvftbp )
     ldatpsv   = .TRUE.
     IF ((mpassiv == 0) .AND. (.NOT. lno))  ldat   = .TRUE.
  ENDIF
 
! 2m dewpoint and relative humidity
! ---------------------------------
 
  rtdt2    =    rmdi
  IF (ntd2ex(nobtyp) == 0)                                       EXIT tem_hum_2m
  CALL obs_extract (istart,nvar( 5),aofbuf,nmxlob)
  IF ((mpassiv == 0) .AND. (lactrep)) THEN
    IF (nvarib == nmdi) neventd (neqmis,ncdtpp,nobtpp) =                       &
                        neventd (neqmis,ncdtpp,nobtpp) + 1
    IF (nflag  >= 2   ) neventd (neqflg,ncdtpp,nobtpp) =                       &
                        neventd (neqflg,ncdtpp,nobtpp) + 1
  ENDIF
  IF ((nvarib == nmdi) .OR. ((nflag >= 2) .AND. (.NOT. lverif))) EXIT tem_hum_2m
  IF (nflag >= 2)  mpassiv = 1
  ntdfwr   =    nvrfwr
  rtdt2    =    nvarib * rfact(nvar  ( 5))
  IF (ntd2ot(nobtyp) == 0)                                       EXIT tem_hum_2m
  IF ((rtdt2 < t0_melt-90._ireals) .OR. (rtdt2 > t0_melt+40._ireals)) THEN
    IF ((mpassiv == 0) .AND. (lactrep))                                        &
      neventd (neqlow,ncdtpp,nobtpp) = neventd(neqlow,ncdtpp,nobtpp) + 1
                                                                 EXIT tem_hum_2m
  ENDIF

! convert to model compatible relative humidity, and write to ODR
! ---------------------------------------------------------------

  nvarfw  =  ireplace(      0, nvfbps(1), nvfboc(1), ntdfwr, nf1bps(4) )
  nvarfw  =  ireplace( nvarfw, nvfbps(2), nvfboc(2), ntdfwr, nf1bps(5) )
  nvarfw  =  ireplace( nvarfw, nvfbps(3), nvfboc(3), ntdfwr, nf1bps(3) )

! eobs : observed vapour press. (measurement is inherently over ice)
!     (but td as stored in data file is over water (WMO convention),
!     which is why magnus formula over water instead of ice is used)
  eobs = fspw( rtdt2 )

! eswo : saturation vapour pressure over water at observed temperat.
  eswo = fspw( rttt2 )

  IF ((rttt2 < t0_melt) .AND. (itype_gscp <= 2))  THEN
!   eseo : saturation vapour pressure over ice at observed temperature
    eseo = fspe( rttt2 )

!   rheo : observed relative humidity over ice   (per def.  <= 100%)
    rheo = eobs / eseo
    !!   this is assumed equivalent to model relative humidity over
!        water (per def.  <= 100%; the model does not know ice) !!
!    ==> eobs-corr / eswo = eobs / eseo
    eobs = eswo * rheo
  ENDIF

  rheo  = eobs / eswo
  rtdt2 = ftd( LOG( eobs / b1 ) )
  IF ((rheo >  rtshlm) .AND. (mpassiv == 0) .AND. (lactrep))                   &
     neventd(neqbig,ncdtpp,nobtpp) = neventd(neqbig,ncdtpp,nobtpp) + 1
  IF  (rheo >  rtshlm)                                           EXIT tem_hum_2m
  IF (rheo >= rhtsat) THEN
     IF ((rheo <  1.0_ireals-epsy) .AND. (mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt2 <  t0_melt) neventd (neqsam,ncdtpp,nobtpp) =                  &
                        neventd (neqsam,ncdtpp,nobtpp) + 1
       IF (rttt2 >= t0_melt) neventd (neqsap,ncdtpp,nobtpp) =                  &
                        neventd (neqsap,ncdtpp,nobtpp) + 1
     ENDIF
     IF ((rheo >  1.0_ireals+epsy) .AND. (mpassiv == 0) .AND. (lactrep)) THEN
       IF (rttt2 <  t0_melt) neventd (neqclm,ncdtpp,nobtpp) =                  &
                        neventd (neqclm,ncdtpp,nobtpp) + 1
       IF (rttt2 >= t0_melt) neventd (neqclp,ncdtpp,nobtpp) =                  &
                        neventd (neqclp,ncdtpp,nobtpp) + 1
     ENDIF
     IF (rheo < c1-epsy)  nvarfw  =  insert( nvarfw, 1, nvfbps(5) )
     IF (rheo > c1+epsy)  nvarfw  =  insert( nvarfw, 1, nvfbps(6) )
     rheo   =  1.0_ireals
     rtdt2  =  rttt2
  ENDIF
 
  lyes = .FALSE.
  lno  = .FALSE.
  DO   j  = 1, mxda
    WRITE(yyes ,'(i5,3x)') mrhyes( j )
    WRITE(yno  ,'(i5,3x)') mrhno ( j )
    IF (ystid == yyes) lyes = .TRUE.
    IF (ystid == yno ) lno  = .TRUE.
    IF (lwonl .AND. (lyes .OR. lno))                                           &
      WRITE( nupr ,'(" STATION: ",a,"  RH YES/NO :",2(2x,a))')                 &
             ystid, yyes, yno
  ENDDO 
  fisdrh = ((fdoro(4)-c1)/c2 + SIGN( (fdoro(4)+c1)/c2 , fisd )) * fisd
  lseaonly = (ABS(altopsu(4))  <   epsy)
  IF (     (.NOT. lyes)                                                        &
      .AND.(     (.NOT. lseaonly .AND. (zstalt >  altopsu(4)))                 &
            .OR. (      lseaonly .AND.  landsy )                               &
            .OR. (fisdrh >  doromx(4))))   lno = .true.
  IF ((lno) .AND. (mpassiv == 0) .AND. (lactrep))                              &
     neventd (neqalt,ncdtpp,nobtpp) = neventd(neqalt,ncdtpp,nobtpp) + 1
  IF ((.NOT. lno) .OR. (lverif)) THEN
     roberr    =  rherr1 / 100._ireals  +  0.0004_ireals *ABS(fisd)
     IF ((lno) .OR. (mpassiv > 0))  roberr = rmdi
     IF  (lno)  nvarfw  =  insert( nvarfw, 1, nvfbps(4) )
     osgbdy (nsgob,nbsrh )  =  rheo
     osgbdy (nsgob,nbsqer)  =  roberr
     mosgbd (nsgob,nbsflg)  =  insert( mosgbd(nsgob,nbsflg), nvarfw, nvfqbp )
     ldatpsv   = .TRUE.
     IF ((mpassiv == 0) .AND. (.NOT. lno))  ldat   = .TRUE.
  ENDIF

  EXIT tem_hum_2m
  ENDDO tem_hum_2m
 
  EXIT Extract_data
  ENDDO Extract_data

! Check if report contains any usable information
  IF ((.NOT. ldat) .AND. (.NOT. ldatpsv)) THEN
!    no accepted data in report
     IF (lactrep)                                                              &
       neventr(nenoda,ncdtpp,nobtpp) = neventr(nenoda,ncdtpp,nobtpp) + 1
     ilen = 2 + istrej
     IF (nacout+ilen <= nmxoln) THEN
       outbuf(nacout+1) = ilen
       outbuf(nacout+2) = nfmt3
       DO icl = 1 , istrej
         outbuf(nacout+2+icl) = ICHAR( yosghd(nsgob) (icl:icl) )
       ENDDO
       nacout  = nacout + ilen
     ENDIF
     lrej = .TRUE.
  ENDIF

  !-----------------------------------------------------------------------------
  ! Section 5.0  Reject the present report
  !-----------------------------------------------------------------------------

  IF (lrej) THEN

     DO  ivar = 1 , mxsbdy
       osgbdy (nsgob,ivar) = rmdi
     ENDDO    
     DO  ivar = 1 , mxsbdf
       mosgbd (nsgob,ivar) = imdi
     ENDDO    
     DO  ivar = 1 , mxshed
       osghed (nsgob,ivar) = rmdi
     ENDDO    
     DO  ivar = 1 , mxshdf
       mosghd (nsgob,ivar) = imdi
     ENDDO
     nsgob = nsgob - 1
     IF (      lactrep) noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
     IF (.NOT. lactrep) noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
     noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1

     lsurfob  =  .FALSE.
! IF (nobtyp == nairep) THEN
!   print *,'r',ystid,ncdtyp,nobtyp,rppp,rttt,ruuu,rvvv
! endif


  ELSEIF (.NOT. lrej) THEN

     lsurfob  =  .TRUE.
     IF ((.NOT. ldat) .AND. (lactrep)) THEN
       neventr(nenoda,ncdtpp,nobtpp) = neventr(nenoda,ncdtpp,nobtpp) + 1
       noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
       noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
     ENDIF
     IF (.NOT. ldat) mosghd(nsgob,nhpass) = 2
  ENDIF

!-------------------------------------------------------------------------------
! End Subroutine obs_single_level
!-------------------------------------------------------------------------------
END SUBROUTINE obs_single_level

 
!-------------------------------------------------------------------------------
!+ Module procedure in "src_obs_proc_aof" for rejection of redundant reports/data
!-------------------------------------------------------------------------------

SUBROUTINE obs_redundancy ( lvprof , nsgoba )

!-------------------------------------------------------------------------------
! Description:
!   This module procedure of module "src_obs_proc_aof" rejects redundant reports /
!   data.
!
! Method:
!   Separate treatment for redundancy between single-level reports (including
!   surface level of TEMPs), between multi-level reports, and between single-
!   level and multi-level aircraft (obs. type: AIREP) reports.
!   Compare with O.I. (DWD):
!   Check includes - removal of redundant information
!                  - multi level check
!                  - level selection
!   The 'first guess departure check' is done later in modified form as
!   'threshold quality control'.
!   Active reports are never replaced by passive reports.
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! S-story:
! Version  Date       Name
! -------- ---------- ----
! 1.13     1998/10/22 Michael Buchhold & Christoph Schraff : Initial release
! 1.19     1998/12/11 Christoph Schraff : Processing also of passive data.
!                                        2nd part of flags replaced by threshold
!                                        QC flags (for a verification mode).
! 1.27     1999/03/29 Christoph Schraff : METARs redundant against SYNOPs with
!                       pressure. ODR flag formats revised. Low cloud introduced.
!                       Station id 6-bit holleriths replaced by characters.
!                       Revised (32-bit) missing data indicator.
! 1.29     1999/05/11 Ulrich Schaettler : Preparations for including MPE_IO
! 1.31     1999/07/01 Christoph Schraff : Bug corr. at merging multi-level reps:
!                       above top of active report; when merging T-profile with
!                       wind profile parts; removal of approx. colocated levels.
!                       Complementation with mandatory levels from rejected rep.
! 1.36     2000/02/24 Michael Buchhold  : Introduction of ACAR aircraft reports.
! 1.38     2000/04/06 Christoph Schraff : Redundant rep. set passive for verif.
! 1.39     2000/05/03 Ulrich Schaettler : Changed variable names for lat/lon.
! 2.5      2001/06/01 Christoph Schraff : Correction to redundancy checking
!                                         between more than 2 colocated reports.
! 3.6      2003/12/11 Christoph Schraff : Inclusion of observation type SATOB.
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
  LOGICAL, INTENT(IN)      ::  &
   lvprof               ! .TRUE. if a multi-level report has been
                        !          added to the ODR, and is to be
                        !          checked for redundancy now

  INTEGER, INTENT(INOUT)   ::  &
   nsgoba               ! number of single-level reports before having
                        ! read new data at current timestep

! Local parameters
! ----------------

  REAL    (KIND=ireals)    , PARAMETER  :: &
    c0e = 0.0001_ireals ! small constant less than accuracy of values

! Local scalars
! -------------

  LOGICAL                  ::  &
   lairc, lair2,      & ! .TRUE., if aircraft report
   lairmlc,           & ! .TRUE., if aircraft report is to be checked against 
                        !  multi-level reports
   lsatbc, lsatb2,    & ! .TRUE., if satob report
   lship, lship2,     & ! .TRUE., if ship report
   lvprofc, lvprof2,  & ! .TRUE., if TEMP or PILOT report
   lvprofr,           & ! .TRUE., if rejected report is TEMP or PILOT
   lcotim,            & ! .TRUE. for temporal colocation
   lcoloc,            & ! .TRUE. for spatial colocation
   lstcor,            & ! .TRUE., if second report is a station correction
   lpassiv,           & ! .TRUE., if passive bit set
   lconven, lconven2, & ! .TRUE., if obs type is TEMP or PILOT
   lraso,   lraso2,   & ! .TRUE., if obs type is radiosonde
   lwprofl, lwprof2,  & ! .TRUE., if obs type is wind profiler, RASS or radarVAD
   ltop,              & ! .TRUE., if top level
   lrepla,            & ! .TRUE., if data from the redundant report may
                        !  replace missing data of the active report
   lreplace,          & ! .TRUE., if redundant data is replaced
   lmand,             & ! .TRUE., if rejected level is mandatory level
   lsub,              & ! .TRUE., if current report is contained in active rep.
   lchecksub,         & ! .TRUE., if only check whether lsub is true
   laddu, laddt, laddq  ! .TRUE., if current quantity is to be added

  INTEGER (KIND=iintegers) ::  &
   insert,            & ! statement function to insert any bit pattern
   invar ,            & ! word to be unpacked / to be partly replaced
   inval ,            & ! bit pattern to be inserted
   ibit  ,            & ! bit position
   ibits ,            & ! statement function to unpack any bit pattern
   ipos  ,            & ! bit position of bit struct. to be unpacked /replaced
   icovr ,            & ! no. of bit occ. by the bit structure
   ireplace,          & ! statement function to replace any bit pattern
   iboc  ,            & ! no. of bit occ. by the bit structure to be replaced
   irepl ,            & ! word containing the replacing bit structure
   iposr ,            & ! bit position of the replacing bit structure
   nstchr, nstchr2,   & ! station characteristic
   nstcor, nstcor2,   & ! station correction bit
   npassiv,npassv2,   & ! passive (bit) flag
   nsglo ,            & ! ODR index of single-level report
   nmlvo ,            & ! ODR index of multi-level report
   nobtyp2,           & ! observation type
   ncdtyp2,           & ! code type
   ncdtypa,           & ! code type of active report
   nodract,           & ! index of active ODR report
   nodrrej,           & ! index of rejected ODR report
   nreplpr,           & ! number of replaced data   
   nactpr,            & ! active report
   nobtypr, ncdtypr,  & ! observation/code type of rejected report
   igplon, igplon2,   & ! gridpoint in zonal direction assigned to observ.
   igplat, igplat2,   & ! gridpoint in meridional direction assig. to observ.
   klev, krej, kact,  & ! level indices
   krjc, kclo, klvi,  & !
   nflaga, nflagr,    & ! variable flags
   irdbfl,            & ! report data base flag
   iovrida, iovridr,  & ! override flags
   istcor,            & ! == 1, if active report is a station correction
   ilevin,            & ! level information
   ilen,              & ! length of output record
   iiv, icl             ! loop indices

  REAL (KIND=ireals)       ::  &
   rdegkm2,           & ! factor for conversion of degree to km (**2)
   rscale,            & ! scaling factor
   obtime, obtime2,   & ! observation date/time
   obbalt, obbalt2,   & ! observation altitude
   cosolat,           & ! cos (lat)
   rtmlmt,            & ! colocation threshold for time
   rhzlmt,            & ! colocation threshold for horizontal distance
   rvtlmt,            & ! colocation threshold for vertical distance
   rdplmt, rprlmt,    & ! colocation thresholds for vert. dist. in pres. units
   dprjc,  dprjca,    & ! pressure differences
   pnexlmt,           & ! pressure limit imposed by next active level above
   rhzob,             & ! horizontal distance in km
   rmod,              & ! name of statement function
   ztti,   ziv          ! dummy arguments of statement function

  CHARACTER (LEN=ilstid)   :: &
   ystid , ystid2,    & ! station identity
   ystids,            & ! station identity
   ystida, ystidr       ! station identity of active and rejected report resp.

  CHARACTER (LEN=20)    yroutine

! Local (automatic) arrays:
! -------------------------
  REAL (KIND=ireals)        ::  &
   zmlbdy(maxrsl,mxrbdy), & ! temporary real array for multi-level report body
   zmlhed(mxrhed)           ! temporary array for one multi-level report header

  INTEGER (KIND=iintegers)  ::  &
   iomlbd(maxrsl,mxrbdf), & ! temporary integer array for multi-level rep. body
   iomlhd(mxrhdf)           ! temporary array for one  multi-level rep. header

!
!------------ End of header ------------------------------------------


! statement function to unpack any bit pattern
! --------------------------------------------
  ibits (invar,ipos,icovr)=IAND (ISHFT(invar,-ipos),icovr)


! Statement function to insert any bit pattern
! --------------------------------------------
  insert(invar,inval,ibit) = IOR (invar , ishft(inval,ibit) )


! Statement function to replace any bit pattern by another bit pattern
! --------------------------------------------------------------------
  ireplace ( invar, ipos, iboc, irepl, iposr ) =                               &
               IOR( IAND( invar, NOT( ISHFT( nibits(iboc), ipos ) ) )          &
                  , ISHFT( IAND( ISHFT( irepl,-iposr ), nibits(iboc) ), ipos ) )


!---------------------------------------------------------------------------
! Begin Subroutine obs_redundancy
!---------------------------------------------------------------------------

! MOD function for positive reals
rmod (ztti,ziv) = ztti - ziv *INT( ztti/ziv + epsy ) + epsy

  yroutine = "obs_redundancy"

  !-------------------------------------------------------------------------
  ! Section 0.   Preset some variable
  !-------------------------------------------------------------------------

  lairc   = .FALSE.
  lairmlc = .FALSE.
  rdegkm2 = 111._ireals **2
  rscale  = 100._ireals*(1._ireals+epsy)

  !-------------------------------------------------------------------------
  ! Section 1.0  Redundancy check for single-level reports
  !-------------------------------------------------------------------------

IF (lsurfob)                   THEN
 
!  Get some information of current report
!  --------------------------------------
 
!  Observation and code type
 
  nobtyp   =  mosghd (nsgob,nhobtp)
  ncdtyp   =  mosghd (nsgob,nhcode)
  lship    =       (ncdtyp == nshscd) .OR. (ncdtyp == nabscd)                  &
              .OR. (ncdtyp == nshred) .OR. (ncdtyp == natshs)
  lvprofc  =  (nobtyp == ntemp) .OR. (nobtyp == npilot)
  lairc    =  (nobtyp == nairep)
  lsatbc   =  (nobtyp == nsatob)
 
!  Get WMO station identity, station correction bit, and passive flag
 
  ystid    =  yosghd (nsgob)
  nstchr   =  mosghd(nsgob,nhschr)
  nstcor   =  ibits (nstchr,nvscbp,nibits (nvscoc))
  npassiv  =  mosghd(nsgob,nhpass)
 
!  Get 4-d location, and cos(latitude) for distances
 
  obtime = osghed (nsgob,nhtime)
  igplon = mosghd (nsgob,nhitot)
  igplat = mosghd (nsgob,nhjtot)
  obbalt = osghed (nsgob,nhalt )
  IF (lairc .OR. lsatbc) obbalt = - osgbdy (nsgob,nbsp  )
  cosolat = COS( (startlat_tot+(igplat-1)*dlat) * degrad )
 
!  Preset thresholds for 4-d colocation
 
  IF ((nobtyp == nairep) .OR. (nobtyp == nsatob)) THEN
    rtmlmt = rtmlair
    rhzlmt = rhzlair
    rvtlmt = rvtlair * 100._ireals
  ELSE
    rtmlmt = rtmlim
    rhzlmt = rhzlim
    rvtlmt = rvtlim
  ENDIF
  rhzlmt = rhzlmt **2

!  Preset ODR index for check against all prior reports
 
  nsglo  = 0


!  Get next report as 'second' report
!  ----------------------------------
  Next_report: DO
 
!  Get ODR index of the second report
  nsglo  = nsglo + 1
 
!  Leave 'loop' if no redundancy with present report
 
  IF (nsglo == nsgob)                                           EXIT Next_report
 
!  Check if current / second report is (partly) redundant
!  ------------------------------------------------------
 
!  Get 4-d location and pressure of second report
 
  obtime2 = osghed (nsglo,nhtime)
  igplon2 = mosghd (nsglo,nhitot)
  igplat2 = mosghd (nsglo,nhjtot)
  obbalt2 = osghed (nsglo,nhalt )
  IF (lairc .OR. lsatbc) obbalt2 = - osgbdy (nsglo,nbsp  )
 
!  Check time difference

  lcotim = ( ABS( obtime2 - obtime )  <=  rtmlmt )
  IF (.NOT. lcotim)                                            CYCLE Next_report
 
!  Check horizontal and vertical (quasi-) colocation
 
  rhzob  = (cosolat *(igplon2-igplon) * dlon)**2 +                             &
                    ((igplat2-igplat) * dlat)**2
  rhzob  = rhzob * rdegkm2
  lcoloc =      (      rhzob               <=  rhzlmt )                        &
          .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmt )
  IF (.NOT. lcoloc)                                            CYCLE Next_report
 
!  Get observation and code type

  nobtyp2  =  mosghd (nsglo ,nhobtp)
  ncdtyp2  =  mosghd (nsglo ,nhcode)
  lship2   =       (ncdtyp2 == nshscd) .OR. (ncdtyp2 == nabscd)                &
              .OR. (ncdtyp2 == nshred) .OR. (ncdtyp2 == natshs)
  lvprof2  =  (nobtyp2 == ntemp) .OR. (nobtyp2 == npilot)
  lair2    =  (nobtyp2 == nairep)
  lsatb2   =  (nobtyp2 == nsatob)
 
!  Check if station id is identical, or if ship,
!  or if 'colocated' temp/pilot and synop/dribu
!  check also if both reports are aireps or not aireps
 
  ystid2   =  yosghd (nsglo)

  IF (.NOT. (     (ystid2 == ystid)                                            &
             .OR. ((      lvprofc) .AND. (.NOT. lvprof2))                      &
             .OR. ((.NOT. lvprofc) .AND. (      lvprof2))                      &
             .OR. ((lship) .AND. (lship2))))                   CYCLE Next_report
  IF (     ((      lairc) .AND. (.NOT. lair2))                                 &
      .OR. ((.NOT. lairc) .AND. (      lair2)))                CYCLE Next_report
  IF (     ((      lsatbc) .AND. (.NOT. lsatb2))                               &
      .OR. ((.NOT. lsatbc) .AND. (      lsatb2)))              CYCLE Next_report
 
!  Decide which report is redundant
!  --------------------------------
 
!  Get station correction bit of 'second' report

  nstchr2  =  mosghd(nsglo,nhschr)
  nstcor2  =  ibits (nstchr2,nvscbp,nibits (nvscoc))
  npassv2  =  mosghd(nsglo,nhpass)
 
!  Decide whether current or 'second' report is redundant
 
  lstcor = .FALSE.
  lrepla = .TRUE.
  IF     ((npassiv <  2) .AND. (npassv2 == 2)) THEN
    nodract  = nsgob
    nodrrej  = nsglo
    lrepla   = .FALSE.
  ELSEIF ((npassiv == 2) .AND. (npassv2 <  2)) THEN
    nodract  = nsglo
    nodrrej  = nsgob
    lrepla   = .FALSE.
  ELSEIF ((nobtyp == nsynop) .AND. (nobtyp2 /= nsynop)) THEN
    nodract  = nsgob
    nodrrej  = nsglo
  ELSEIF ((nobtyp /= nsynop) .AND. (nobtyp2 == nsynop)) THEN
    nodract  = nsglo
    nodrrej  = nsgob
  ELSEIF (nstcor >  nstcor2) THEN
    nodract  = nsgob
    nodrrej  = nsglo
    lstcor   = .TRUE.
  ELSEIF (nstcor <  nstcor2) THEN
    nodract  = nsglo
    nodrrej  = nsgob
    lstcor   = .TRUE.
  ELSEIF ((nobtyp == nsynop) .AND. (nobtyp2 == nsynop)) THEN
    IF (      (osgbdy(nsglo,nbszer) < rmdich)                                  &
        .AND. (osgbdy(nsgob,nbszer) > rmdich)) THEN
! report 1 is SYNOP, report 2 is METAR and rejected
      nodract  = nsgob
      nodrrej  = nsglo
    ELSE
      nodract  = nsglo
      nodrrej  = nsgob
    ENDIF
  ELSE
    nodract  = nsglo
    nodrrej  = nsgob
  ENDIF
 
!  Determine whether the 'old' report has the passive bit set, which means that,
!  if it's of type 'aircraft', it is part of a multi-level report.
 
! lpassiv = ibits (mosghd(nsglo,nhschr),nvpsbp,nvscoc)  ==  1
  lpassiv = mosghd(nsglo,nhpass)  ==  1
 
!  Replace / reject redundant data
!  -------------------------------
 
!  Replace missing data of active report by available data of rejected report
 
  nreplpr = 0
  IF (     (.NOT.lvprofc) .AND. (.NOT.lvprof2) .AND. (.NOT.lairc)              &
                                               .AND. (.NOT.lsatbc)             &
      .AND.(osgbdy(nodract,nbszer) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbszer) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbsp  ) = osgbdy (nodrrej,nbsp  )
    osgbdy (nodract,nbsz  ) = osgbdy (nodrrej,nbsz  )
    osgbdy (nodract,nbszer) = osgbdy (nodrrej,nbszer)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvfzbp, nvfaoc &
                                      , mosgbd(nodrrej,nbsflg), nvfzbp)
    nreplpr = nreplpr + 1
  ENDIF
  IF (     (osgbdy(nodract,nbsuer) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbsuer) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbsu  ) = osgbdy (nodrrej,nbsu  )
    osgbdy (nodract,nbsv  ) = osgbdy (nodrrej,nbsv  )
    osgbdy (nodract,nbsuer) = osgbdy (nodrrej,nbsuer)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvfubp, nvfaoc &
                                      , mosgbd(nodrrej,nbsflg), nvfubp)
    nreplpr = nreplpr + 2
  ENDIF
  IF (     (osgbdy(nodract,nbster) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbster) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbst  ) = osgbdy (nodrrej,nbst  )
    osgbdy (nodract,nbster) = osgbdy (nodrrej,nbster)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvftbp, nvfaoc &
                                      , mosgbd(nodrrej,nbsflg), nvftbp)
    nreplpr = nreplpr + 4
  ENDIF
  IF (     (osgbdy(nodract,nbsqer) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbsqer) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbsrh ) = osgbdy (nodrrej,nbsrh )
    osgbdy (nodract,nbsqer) = osgbdy (nodrrej,nbsqer)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvfqbp, nvfaoc &
                                      , mosgbd(nodrrej,nbsflg), nvfqbp)
    nreplpr = nreplpr + 8
  ENDIF
  IF (     (osgbdy(nodract,nbsper) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbsper) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbspr ) = osgbdy (nodrrej,nbspr )
    osgbdy (nodract,nbsper) = osgbdy (nodrrej,nbsper)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvrfbp, nv2foc &
                                      , mosgbd(nodrrej,nbsflg), nvrfbp)
    nreplpr = nreplpr + 16
  ENDIF
  IF (     (osgbdy(nodract,nbscer) <  rmdich)                                  &
      .AND.(osgbdy(nodrrej,nbscer) >  rmdich) .AND. (lrepla)) THEN
    osgbdy (nodract,nbscl ) = osgbdy (nodrrej,nbscl )
    osgbdy (nodract,nbscer) = osgbdy (nodrrej,nbscer)
    mosgbd (nodract,nbsflg) = ireplace( mosgbd(nodract,nbsflg), nvcfbp, nv2foc &
                                      , mosgbd(nodrrej,nbsflg), nvcfbp)
    nreplpr = nreplpr + 32
  ENDIF

!  Set rejected report to passive, if verification of passive reports is on, or
!  replace present report by active report
!  and set last (current) report to missing.

  IF ((lverif) .AND. (lverpas)) THEN
    mosghd (nodrrej,nhpass)  = 2
    mosghd (nodrrej,nhschr)  = insert( mosghd(nodrrej,nhschr) , 1 , nvrdbp )
  ELSEIF (.NOT. ((lairc) .AND. (nodrrej <= nsgoba))) THEN
    osgbdy (nsglo ,1:mxsbdy) = osgbdy (nodract,1:mxsbdy)
    osgbdy (nsgob ,1:mxsbdy) = rmdi
    mosgbd (nsglo ,1:mxsbdf) = mosgbd (nodract,1:mxsbdf)
    mosgbd (nsgob ,1:mxsbdf) = imdi
    osghed (nsglo ,1:mxshed) = osghed (nodract,1:mxshed)
    osghed (nsgob ,1:mxshed) = rmdi
    mosghd (nsglo ,1:mxshdf) = mosghd (nodract,1:mxshdf)
    mosghd (nsgob ,1:mxshdf) = imdi
    yosghd (nsglo )          = yosghd (nodract)
    yosghd (nsgob )          = '        '
  ELSE
 
!  If the current report replaces an 'old' (read at a
!  previous timestep) aircraft report: 
!  First replace that old aircraft report by the 'last
!  old' report, then replace the 'last old' report by
!  the current report, then set the last (current)
!  report to missing.
!  Also, reduce the number of 'old' reports by one.
!  (this procedure is necessary for the correct produc-
!  tion of multi-level aircraft reports, see *airporg*)
 
    osgbdy (nsglo ,1:mxsbdy) = osgbdy (nsgoba ,1:mxsbdy)
    osgbdy (nsgoba,1:mxsbdy) = osgbdy (nodract,1:mxsbdy)
    osgbdy (nsgob ,1:mxsbdy) = rmdi
    mosgbd (nsglo ,1:mxsbdf) = mosgbd (nsgoba ,1:mxsbdf)
    mosgbd (nsgoba,1:mxsbdf) = mosgbd (nodract,1:mxsbdf)
    mosgbd (nsgob ,1:mxsbdf) = imdi
    osghed (nsglo ,1:mxshed) = osghed (nsgoba ,1:mxshed)
    osghed (nsgoba,1:mxshed) = osghed (nodract,1:mxshed)
    osghed (nsgob ,1:mxshed) = rmdi
    mosghd (nsglo ,1:mxshdf) = mosghd (nsgoba ,1:mxshdf)
    mosghd (nsgoba,1:mxshdf) = mosghd (nodract,1:mxshdf)
    mosghd (nsgob ,1:mxshdf) = imdi
    yosghd (nsglo )          = yosghd (nsgoba )
    yosghd (nsgoba)          = yosghd (nodract)
    yosghd (nsgob )          = '        '
    nsglo   = nsgoba
    nodrrej = nsgoba
    nsgoba  = nsgoba - 1
  ENDIF

!  Get diagnostic array position (--> nobtpp, ncdtpp)
!  and update statistics and report events
 
  IF (nodrrej == nsglo) THEN
    nobtypr = nobtyp2
    ncdtypr = ncdtyp2
    lvprofr = lvprof2
  ELSE
    nobtypr = nobtyp
    ncdtypr = ncdtyp
    lvprofr = lvprofc
  ENDIF
  CALL obs_pointrs ( nobtypr , ncdtypr )
  IF (.NOT. lvprofr) THEN
    IF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
      noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
      neventr (neredn,ncdtpp,nobtpp) = neventr(neredn,ncdtpp,nobtpp) + 1
    ELSE
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
    ENDIF
    IF ((lverif) .AND. (lverpas)) THEN
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
    ELSE
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
    ENDIF
  ELSEIF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
    neventr (neredx,ncdtpp,nobtpp) = neventr(neredx,ncdtpp,nobtpp) + 1
  ENDIF

!  Print message

  IF (lprodr)  THEN
    IF (nodrrej == nsglo) THEN
      ncdtypa = ncdtyp
      ystida  = ystid
      ncdtypr = ncdtyp2
      ystidr  = ystid2
      nactpr  = 2
    ELSE
      ncdtypa = ncdtyp2
      ystida  = ystid2
      ncdtypr = ncdtyp
      ystidr  = ystid
      nactpr  = 1
    ENDIF
    ilen = 15 + 2*istrej
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt10
      outbuf(nacout+ 3) = ncdtypr
      DO iiv = 1,5
        IF (osgbdy(nsglo,iiv) > rmdich) THEN
           outbuf(nacout+3+iiv) = INT(osgbdy(nsglo,iiv)*rscale)
        ELSE
           outbuf(nacout+3+iiv) = imdi
        ENDIF
      ENDDO
      outbuf(nacout+ 9) = ncdtypa
      outbuf(nacout+10) = INT(osghed(nsglo,nhilon)*rscale)
      outbuf(nacout+11) = INT(osghed(nsglo,nhjlat)*rscale)
      outbuf(nacout+12) = nactpr
      istcor = 0
      IF (lstcor) istcor = 1
      outbuf(nacout+13) = istcor
      outbuf(nacout+14) = nreplpr
      outbuf(nacout+15) = INT(osghed(nsglo,nhtime)*rscale)
      DO icl = 1 , istrej
        outbuf(nacout+15       +icl) = ICHAR( ystidr (icl:icl) )
        outbuf(nacout+15+istrej+icl) = ICHAR( ystida (icl:icl) )
      ENDDO
      nacout  = nacout + ilen
    ENDIF
  ENDIF ! lprodr == .true.
 
! Adjust total number of single-level reports
 
  IF (.NOT. ((lverif) .AND. (lverpas)))  nsgob = nsgob - 1
 
!  Prepare for further processing (++++)
!  ------------------------------
 
!  An aircraft report should be checked against multi-level reports, only if it
!  does not appear as single-level report with index 'nsglo > nsgoba', and if
!  the active single-level report is part of a multi-level report (since it may
!  complement missing data of the active report / level).
!  conditions: - aircraft report
!              - current report redundant
!              - active report read at previous timestep
!              - passive bit of active present report set
!  (otherwise, all reports with current station ID will be re-processed in subr.
!   *obs_air_org_mult*.)

  lairmlc  =  (lairc) .AND. (nodract == nsglo) .AND. (lpassiv)                 &
                      .AND. (nsgoba >= nsglo) .AND. (lrepla)
 
  IF (.NOT. ((lverif) .AND. (lverpas)))                        EXIT  Next_report
  ENDDO  Next_report
 
ENDIF

  !-------------------------------------------------------------------------
  ! Section 2.0  Redundancy check for multi-level reports
  !-------------------------------------------------------------------------
 
IF (lvprof )                                        THEN      


!  (For redundancy of multi-level AIREPs (not active),
!   lowest data have to be (quasi-) colocated)
 
!  Get some information of current report
!  --------------------------------------
 
!  Observation and code type
 
  nobtyp   =  momlhd (nmlob,nhobtp)
  ncdtyp   =  momlhd (nmlob,nhcode)
  lconven  =  (nobtyp == ntemp) .OR. (nobtyp == npilot)
  lraso    =  (nobtyp == npilot) .AND. (     (ncdtyp == nldpcd)                &
                                        .OR. (ncdtyp == nshpcd))
  lraso    =  (nobtyp == ntemp) .OR. (lraso)
  lwprofl  =  (nobtyp == npilot) .AND. (     (ncdtyp == nwp_eu)                &
                                        .OR. (ncdtyp == nra_eu)                &
                                        .OR. (ncdtyp == npr_us)                &
                                        .OR. (ncdtyp == nravad))
  lairc    =  (nobtyp == nairep)
   
!  Get WMO station identity, station correction bit, and passive flag
 
  ystid    =  yomlhd (nmlob)

  nstchr   =  momlhd(nmlob,nhschr)
  nstcor   =  ibits (nstchr,nvscbp,nibits (nvscoc))
  npassiv  =  momlhd(nmlob,nhpass)

!  get 4-d location, and cos(latitude) for distances
 
  obtime = omlhed (nmlob,nhtime)
  igplon = momlhd (nmlob,nhitot)
  igplat = momlhd (nmlob,nhjtot)

  obbalt = omlhed (nmlob,nhalt )
  IF (lairc) obbalt = - osgbdy (nmlob,nbtp)
  cosolat = COS( (startlat_tot+(igplat-1)*dlat) * degrad )

!  Preset thresholds for 4-d colocation

  IF (nobtyp == nairep) THEN
    rtmlmt = rtmlair
    rhzlmt = rhzlair
    rvtlmt = rvtlair * 100._ireals
    rdplmt = rvtlair * 100._ireals
    rprlmt = rprlim  * 100._ireals
  ELSE
    rtmlmt = rtmlim
    rhzlmt = rhzlim
    rvtlmt = rvtlim
    rdplmt = rdplim * 100._ireals
    rprlmt = rprlim * 100._ireals
  ENDIF
  rhzlmt = rhzlmt **2

!  Preset ODR index for check against all prior reports:
!  at first, the report is auto-checked against itself to get rid of (approx.)
!  colocated levels within the report

  nmlvo  = nmlob + 1

!  Get next report as 'second' report
!  ----------------------------------

  Next_mlrep:  DO
! ~~~~~~~~~~~~~~~

!  Get ODR index of the second report

  nmlvo  = nmlvo - 1

!  Leave 'loop' if no redundancy with present report

  IF (nmlvo == 0)                                                EXIT Next_mlrep
 
!  Check if current / second report is (partly) redundant
!  ------------------------------------------------------
 
!  Get 4-d location and pressure of second report
 
  obtime2 = omlhed (nmlvo,nhtime)
  igplon2 = momlhd (nmlvo,nhitot)
  igplat2 = momlhd (nmlvo,nhjtot)
  obbalt2 = omlhed (nmlvo,nhalt )
  IF (lairc) obbalt2 = - osgbdy (nmlvo,nbtp)


!  Check time difference
 
  lcotim =      ( ABS( obtime2 - obtime )  <=  rtmlmt )
  IF (.NOT.lcotim)                                              CYCLE Next_mlrep

!  Check horizontal and vertical (quasi-) colocation

  rhzob  = (cosolat *(igplon2-igplon) * dlon)**2 +                             &
                    ((igplat2-igplat) * dlat)**2
  rhzob  = rhzob * rdegkm2
  lcoloc =      (      rhzob               <=  rhzlmt )                        &
          .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmt )
  IF (.NOT.lcoloc)                                              CYCLE Next_mlrep

!  Get observation and code type

  nobtyp2  =  momlhd (nmlvo ,nhobtp)
  ncdtyp2  =  momlhd (nmlvo ,nhcode)
  lconven2 =  (nobtyp2 == ntemp) .OR. (nobtyp2 == npilot)
  lraso2   =  (nobtyp2 == npilot) .AND. (     (ncdtyp2 == nldpcd)              &
                                         .OR. (ncdtyp2 == nshpcd))
  lraso2   =  (nobtyp2 == ntemp) .OR. (lraso2)
  lwprof2  =  (nobtyp2 == npilot) .AND. (     (ncdtyp2 == nwp_eu)              &
                                         .OR. (ncdtyp2 == nra_eu)              &
                                         .OR. (ncdtyp2 == npr_us)              &
                                         .OR. (ncdtyp2 == nravad))
  lair2    =  (nobtyp2 == nairep)
 
!  No redundancy between reports from the different
!  observing systems  TEMP/PILOT , AIREP and PROFILER
 
  ystid2   =  yomlhd (nmlvo)

  IF  (lraso2  .NEQV. lraso)                                    CYCLE Next_mlrep
  IF  (lwprof2 .NEQV. lwprofl)                                  CYCLE Next_mlrep
  IF  (lair2   .NEQV. lairc)                                    CYCLE Next_mlrep
  IF ((nobtyp /= nobtyp2) .AND. (.NOT. (lraso .AND. lraso2)))   CYCLE Next_mlrep
 
!  No redundancy (check), if current report is active and past report is passive
!  (in this way, the current active report is given the chance to be checked
!   against another past report which is active (in case of lverpas = .TRUE.))
!   In this case, and only if the reason for the past report to be passive is 
!   redundancy, check only whether the current report is completely contained
!   in the past report. In this case, reject the current report.

  lchecksub = .FALSE.
  npassv2  =  momlhd(nmlvo,nhpass)
! IF ((npassiv <  2) .AND. (npassv2 == 2))                      CYCLE Next_mlrep
  IF ((npassiv <  2) .AND. (npassv2 == 2)) THEN
    IF (ibits(momlhd(nmlvo,nhschr),nvrdbp,nibits(1)) == 0)      CYCLE Next_mlrep
    lchecksub = .TRUE.
  ENDIF

!  Decide which report is redundant
!  --------------------------------

!  Get station correction bit of 'second' report
  nstchr2  =  momlhd(nmlvo,nhschr)
  nstcor2  =  ibits (nstchr2,nvscbp,nibits (nvscoc))
  npassv2  =  momlhd(nmlvo,nhpass)

!  Determine indices of active resp. redundant report
!  (if airep, both reports may be used)

  lstcor   = .FALSE.
  lrepla = .TRUE.
  IF     ((npassiv <  2) .AND. (npassv2 == 2)) THEN
    nodract  = nmlob
    nodrrej  = nmlvo
    IF (lchecksub) THEN
      nodract  = nmlvo
      nodrrej  = nmlob
    ENDIF
    lrepla   = .FALSE.
  ELSEIF ((npassiv == 2) .AND. (npassv2 <  2)) THEN
    nodract  = nmlvo
    nodrrej  = nmlob
    lrepla   = .FALSE.
  ELSEIF ((nobtyp == ntemp ).AND.(nobtyp2 /= ntemp )) THEN
    nodract  = nmlob
    nodrrej  = nmlvo
  ELSEIF ((nobtyp /= ntemp ).AND.(nobtyp2 == ntemp )) THEN
    nodract  = nmlvo
    nodrrej  = nmlob
  ELSEIF ((lairc).AND.(obbalt2-obbalt >  rvtlmt)) THEN
    nodract  = nmlob
    nodrrej  = nmlvo
  ELSEIF (lairc) THEN
    nodract  = nmlvo
    nodrrej  = nmlob
  ELSEIF (nstcor >  nstcor2) THEN
    nodract  = nmlob
    nodrrej  = nmlvo
    lstcor   = .TRUE.
  ELSEIF (nstcor <  nstcor2) THEN
    nodract  = nmlvo
    nodrrej  = nmlob
    lstcor   = .TRUE.
  ELSE
    nodract  = nmlvo
    nodrrej  = nmlob
  ENDIF

!  Fill data of active report and supplementary data
!  of redundant report into temporary field
!  -------------------------------------------------
 
!  Fill header of active report into temporary field;
!  fill body of temporary field with missing values

  zmlhed (1:mxrhed) = omlhed (nodract,1:mxrhed)
  iomlhd (1:mxrhdf) = momlhd (nodract,1:mxrhdf)
  ystids            = yomlhd (nodract)
  zmlbdy (1:maxrsl,1:mxrbdy) = rmdi
  iomlbd (1:maxrsl,1:mxrbdf) = imdi
 
!  Loop over levels of body of active report
 
  lsub = .TRUE.
  ltop = .FALSE.
  klev = 0
  krej = 1

  Levels_acr: DO  kact = 1 , momlhd(nodract,nhnlev)
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! 'while' loop: get all rejected levels which are more than 'rdplmt' [pa]
!               below the active level 'kact'.
!               Of these levels, fill those levels, which are not within
!               'rdplmt' [pa] of any level of the active report,
!               into the temporary field
 
  DO WHILE (      (   omlbdy(nodrrej,krej,nbtp)                                &
                   >= omlbdy(nodract,kact,nbtp)+rdplmt)                        &
            .AND. (.NOT. ltop) .AND. (klev < maxrsl))
     IF (lrepla) THEN
       klev = klev + 1
       zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
       iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
       IF (zmlbdy(klev,nbtuer) > rmdich)  iomlhd (nhuexi) = 1
       IF (zmlbdy(klev,nbtter) > rmdich)  iomlhd (nhtexi) = 1
       IF (zmlbdy(klev,nbtqer) > rmdich)  iomlhd (nhqexi) = 1
     ENDIF
     IF (krej <  momlhd(nodrrej,nhnlev)) THEN
       krej = krej + 1
     ELSE
       ltop = .TRUE.
     ENDIF
     lsub = .FALSE.
  ENDDO

! 'while' loop: get all rejected levels which are less than 'rdplmt' [pa],
!               but more than 'rprlmt' [pa] below the active level 'kact'.
!               Of these levels, fill those observed quantities, which are not
!               present in the active level 'kact', into the temporary field.
!        (Note: This may be important if e.g. a TEMP report containing only
!               temperature and humidity and a PILOT report have been derived
!               and disseminated from one and the same radiosonde ascent.)

  DO WHILE (      (   omlbdy(nodrrej,krej,nbtp)                                &
                   >= omlbdy(nodract,kact,nbtp)+rprlmt)                        &
            .AND. (.NOT. ltop) .AND. (klev < maxrsl))
     laddu = .FALSE.
     laddt = .FALSE.
     laddq = .FALSE.
     IF (lrepla) THEN
       lmand =     (rmod( omlbdy(nodrrej,krej,nbtp),10000.0_ireals ) < 2*epsy) &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)- 5000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-15000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-25000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-85000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-92500.0_ireals ) < epsy)
       laddu =        ((omlbdy(nodract,kact,nbtuer) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtuer) >= rmdich)     
       laddt =        ((omlbdy(nodract,kact,nbtter) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtter) >= rmdich)     
       laddq =        ((omlbdy(nodract,kact,nbtqer) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtqer) >= rmdich)     
     ENDIF
     IF ((laddu) .OR. (laddt) .OR. (laddq)) THEN
       klev = klev + 1
       zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
       iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
       IF (.NOT. laddu) zmlbdy (klev,nbtu  ) = rmdi
       IF (.NOT. laddu) zmlbdy (klev,nbtv  ) = rmdi
       IF (.NOT. laddu) zmlbdy (klev,nbtuer) = rmdi
       IF (.NOT. laddu) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvfubp,nvfaoc, 0,nvfubp)
       IF (.NOT. laddt) zmlbdy (klev,nbtt  ) = rmdi
       IF (.NOT. laddt) zmlbdy (klev,nbtter) = rmdi
       IF (.NOT. laddt) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvftbp,nvfaoc, 0,nvftbp)
       IF (.NOT. laddq) zmlbdy (klev,nbtrh ) = rmdi
       IF (.NOT. laddq) zmlbdy (klev,nbtqer) = rmdi
       IF (.NOT. laddq) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvfqbp,nvfaoc, 0,nvfqbp)
       IF (zmlbdy(klev,nbtuer) > rmdich)  iomlhd (nhuexi) = 1
       IF (zmlbdy(klev,nbtter) > rmdich)  iomlhd (nhtexi) = 1
       IF (zmlbdy(klev,nbtqer) > rmdich)  iomlhd (nhqexi) = 1
     ENDIF
     IF (krej <  momlhd(nodrrej,nhnlev)) THEN
       krej = krej + 1
     ELSE
       ltop = .TRUE.
     ENDIF
     lsub = .FALSE.
  ENDDO
 
!  check if there is still space in the buffer

  IF (klev == maxrsl)                                            EXIT Levels_acr

!  (approx.) colocated levels are removed
!  (however, the active level is complemented by colocated level data
!   in the auto-checking)

  IF ((kact > 1) .AND. (  omlbdy(nodract,MAX(kact-1,1),nbtp)-rprlmt)           &
                        < omlbdy(nodract,kact,nbtp))            CYCLE Levels_acr
 
!  Fill next level of active report into temporary field
 
  klev = klev + 1
  zmlbdy (klev,1:mxrbdy) = omlbdy (nodract,kact,1:mxrbdy)
  iomlbd (klev,1:mxrbdf) = momlbd (nodract,kact,1:mxrbdf)
  krjc  = krej
  kclo  = krej
  dprjc = - rdplmt
 
! 'while' loop: get all rejected levels which are within 'rprlmt' [pa] of the
!               active level 'kact';
!               find the rejected level 'kclo' most distant to the act. level 
!               'kact' (since that level is most likely to contain additional
!               information, particularly when auto-checking).

  DO WHILE (      (   omlbdy(nodrrej,krej,nbtp)                                &
                   >  omlbdy(nodract,kact,nbtp)-rprlmt)                        &
            .AND. (.NOT. ltop))
    dprjca = ABS(  omlbdy(nodrrej,krej,nbtp)                                   &
                 - omlbdy(nodract,kact,nbtp) )
    IF (dprjca >= dprjc) THEN
      dprjc = dprjca
      kclo  = krej
    ENDIF
    IF (krej <  momlhd(nodrrej,nhnlev)) THEN
      krej = krej + 1
    ELSE
      ltop = .TRUE.
    ENDIF
  ENDDO

!  If level 'kclo' is close enough to level 'kact', complement missing data
!  and replace data if override bit is higher or flag is lower.
 
  IF ((ABS( omlbdy(nodrrej,kclo,nbtp)                                          &
           -omlbdy(nodract,kact,nbtp)) <= rprlmt) .AND. (lrepla)) THEN

!   If level of rejected report is mandatory, replace pressure and height

    lmand  =     (rmod( omlbdy(nodrrej,kclo,nbtp),10000.0_ireals ) < 2*epsy)   &
            .OR. ( ABS( omlbdy(nodrrej,kclo,nbtp)- 5000.0_ireals ) < epsy)     &
            .OR. ( ABS( omlbdy(nodrrej,kclo,nbtp)-15000.0_ireals ) < epsy)     &
            .OR. ( ABS( omlbdy(nodrrej,kclo,nbtp)-25000.0_ireals ) < epsy)     &
            .OR. ( ABS( omlbdy(nodrrej,kclo,nbtp)-85000.0_ireals ) < epsy)     &
            .OR. ( ABS( omlbdy(nodrrej,kclo,nbtp)-92500.0_ireals ) < epsy)
    IF (      (ABS( omlbdy(nodrrej,kclo,nbtp)                                  &
                   -omlbdy(nodract,kact,nbtp)) > epsy)                         &
        .AND. (     (      (omlbdy(nodract,kact,nbtzer) < rmdich)              &
                     .AND. (omlbdy(nodrrej,kclo,nbtzer) > rmdich))             &
               .OR. (      (lmand)                                             &
                     .AND. (     (omlbdy(nodract,kact,nbtzer) < rmdich)        &
                            .OR. (omlbdy(nodrrej,kclo,nbtzer) > rmdich))))) THEN
      ilen = 8 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+ 1) = ilen
        outbuf(nacout+ 2) = nfmt16
        outbuf(nacout+ 3) = INT(omlbdy(nodract,kact,nbtp  )*rscale)
        outbuf(nacout+ 4) = INT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
        outbuf(nacout+ 5) = imdi
        outbuf(nacout+ 6) = imdi
        outbuf(nacout+ 7) = imdi
        outbuf(nacout+ 8) = imdi
        IF (omlbdy(nodract,kact,nbtzer) > rmdich)                            &
          outbuf(nacout+ 5) = INT(omlbdy(nodract,kact,nbtzer)*rscale)
        IF (omlbdy(nodrrej,kclo,nbtzer) > rmdich)                            &
          outbuf(nacout+ 6) = INT(omlbdy(nodrrej,kclo,nbtzer)*rscale)
        IF (omlbdy(nodract,kact,nbtz  ) > rmdich)                            &
          outbuf(nacout+ 7) = INT(omlbdy(nodract,kact,nbtz  )*rscale)
        IF (omlbdy(nodrrej,kclo,nbtz  ) > rmdich)                            &
          outbuf(nacout+ 8) = INT(omlbdy(nodrrej,kclo,nbtz  )*rscale)
        DO icl = 1 , istrej
          outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF

      IF (num_compute == 1)                                                    &
         WRITE( nuodr,'(1x,a,'' Z: ACT/REJ: P/ZERR/Z:'',2f8.0,2f6.1,2f7.0)')   &
                ystid                                                          &
               ,omlbdy(nodract,kact,nbtp)                                      &
               ,omlbdy(nodrrej,kclo,nbtp)                                      &
               ,omlbdy(nodract,kact,nbtzer)                                    &
               ,omlbdy(nodrrej,kclo,nbtzer)                                    &
               ,omlbdy(nodract,kact,nbtz  )                                    &
               ,omlbdy(nodrrej,kclo,nbtz  )
      zmlbdy (klev,nbtp  ) = omlbdy(nodrrej,kclo,nbtp)
      zmlbdy (klev,nbtz  ) = omlbdy(nodrrej,kclo,nbtz)
      zmlbdy (klev,nbtzer) = omlbdy(nodrrej,kclo,nbtzer)
      iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg), nvfzbp, nvfaoc     &
                                     , momlbd(nodrrej,kclo,nbtflg), nvfzbp )
    ENDIF
 
!   a: wind components
 
    irdbfl  =  ibits (momlbd(nodract,kact,nbtflg), nvfubp, nibits(nvfaoc))
    iovrida =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflaga  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    irdbfl  =  ibits (momlbd(nodrrej,kclo,nbtflg), nvfubp, nibits(nvfaoc))
    iovridr =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflagr  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    IF (     (      (omlbdy(nodract,kact,nbtuer) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtuer) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtu  ) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtu  ) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtuer) > rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtuer) > rmdich)                     &
              .AND. (nflagr <  nflaga))) THEN
      ilen = 8 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+ 1) = ilen
        outbuf(nacout+ 2) = nfmt13
        outbuf(nacout+ 3) = INT(omlbdy(nodract,kact,nbtp  )*rscale)
        outbuf(nacout+ 4) = INT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
        outbuf(nacout+ 5) = imdi
        outbuf(nacout+ 6) = imdi
        IF (omlbdy(nodract,kact,nbtuer) > rmdich)                              &
          outbuf(nacout+ 5) = INT(omlbdy(nodract,kact,nbtuer)*rscale)
        IF (omlbdy(nodrrej,kclo,nbtuer) > rmdich)                              &
          outbuf(nacout+ 6) = INT(omlbdy(nodrrej,kclo,nbtuer)*rscale)
        outbuf(nacout+ 7) = nflaga
        outbuf(nacout+ 8) = nflagr
        DO icl = 1 , istrej
          outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF

      IF (num_compute == 1)                                                    &
         WRITE( nuodr,'(1x,a,'' U: ACT/REJ: P/U/FLAG:'',2f8.0,2f6.1,2i2)')     &
                ystid                                                          &
               ,omlbdy(nodract,kact,nbtp)                                      &
               ,omlbdy(nodrrej,kclo,nbtp)                                      &
               ,omlbdy(nodract,kact,nbtuer)                                    &
               ,omlbdy(nodrrej,kclo,nbtuer), nflaga, nflagr
      zmlbdy (klev,nbtu  ) = omlbdy (nodrrej,kclo,nbtu  )
      zmlbdy (klev,nbtv  ) = omlbdy (nodrrej,kclo,nbtv  )
      zmlbdy (klev,nbtuer) = omlbdy (nodrrej,kclo,nbtuer)
      iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg), nvfubp, nvfaoc     &
                                     , momlbd(nodrrej,kclo,nbtflg), nvfubp )
      iomlhd (     nhuexi) = 1
    ENDIF

!   b: temperature
 
    irdbfl  =  ibits (momlbd(nodract,kact,nbtflg), nvftbp, nibits(nvfaoc))
    iovrida =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflaga  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    irdbfl  =  ibits (momlbd(nodrrej,kclo,nbtflg), nvftbp, nibits(nvfaoc))
    iovridr =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflagr  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    IF (     (      (omlbdy(nodract,kact,nbtter) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtter) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtt  ) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtt  ) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtter) > rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtter) > rmdich)                     &
              .AND. (nflagr <  nflaga))) THEN
      ilen = 8 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+ 1) = ilen
        outbuf(nacout+ 2) = nfmt14
        outbuf(nacout+ 3) = INT(omlbdy(nodract,kact,nbtp  )*rscale)
        outbuf(nacout+ 4) = INT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
        outbuf(nacout+ 5) = imdi
        outbuf(nacout+ 6) = imdi
        IF (omlbdy(nodract,kact,nbtter) > rmdich)                              &
          outbuf(nacout+ 5) = INT(omlbdy(nodract,kact,nbtter)*rscale)
        IF (omlbdy(nodrrej,kclo,nbtter) > rmdich)                              &
          outbuf(nacout+ 6) = INT(omlbdy(nodrrej,kclo,nbtter)*rscale)
        outbuf(nacout+ 7) = nflaga
        outbuf(nacout+ 8) = nflagr
        DO icl = 1 , istrej
          outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF

      IF (num_compute == 1)                                                    &
         WRITE( nuodr,'(''T: ACT/REJ: P/T/FLAG:'',2f8.0,2f6.1,2i2)')           &
                omlbdy(nodract,kact,nbtp)                                      &
               ,omlbdy(nodrrej,kclo,nbtp)                                      &
               ,omlbdy(nodract,kact,nbtter)                                    &
               ,omlbdy(nodrrej,kclo,nbtter), nflaga, nflagr
      zmlbdy (klev,nbtt  ) = omlbdy (nodrrej,kclo,nbtt  )
      zmlbdy (klev,nbtter) = omlbdy (nodrrej,kclo,nbtter)
      iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg), nvftbp, nvfaoc     &
                                     , momlbd(nodrrej,kclo,nbtflg), nvftbp )
      iomlhd (     nhtexi) = 1
    ENDIF
 
!   c: humidity
 
    irdbfl  =  ibits (momlbd(nodract,kact,nbtflg), nvfqbp, nibits(nvfaoc))
    iovrida =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflaga  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    irdbfl  =  ibits (momlbd(nodrrej,kclo,nbtflg), nvfqbp, nibits(nvfaoc))
    iovridr =  ibits (irdbfl,nvfbps( 3),nibits (nvfboc( 3)))
    nflagr  =  ibits (irdbfl,nvfbps( 1),nibits (nvfboc( 1)))
    IF (     (      (omlbdy(nodract,kact,nbtqer) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtqer) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtrh ) < rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtrh ) > rmdich))                    &
        .OR. (      (omlbdy(nodract,kact,nbtqer) > rmdich)                     &
              .AND. (omlbdy(nodrrej,kclo,nbtqer) > rmdich)                     &
              .AND. (nflagr <  nflaga))) THEN
      ilen = 8 + istrej
      IF (nacout+ilen <= nmxoln) THEN
        outbuf(nacout+ 1) = ilen
        outbuf(nacout+ 2) = nfmt15
        outbuf(nacout+ 3) = INT(omlbdy(nodract,kact,nbtp  )*rscale)
        outbuf(nacout+ 4) = INT(omlbdy(nodrrej,kclo,nbtp  )*rscale)
        outbuf(nacout+ 5) = imdi
        outbuf(nacout+ 6) = imdi
        IF (omlbdy(nodract,kact,nbtqer) > rmdich)                              &
          outbuf(nacout+ 5) = INT(omlbdy(nodract,kact,nbtqer)*rscale)
        IF (omlbdy(nodrrej,kclo,nbtqer) > rmdich)                              &
          outbuf(nacout+ 6) = INT(omlbdy(nodrrej,kclo,nbtqer)*rscale)
        outbuf(nacout+ 7) = nflaga
        outbuf(nacout+ 8) = nflagr
        DO icl = 1 , istrej
          outbuf(nacout+8+icl) = ICHAR( ystid (icl:icl) )
        ENDDO
        nacout  = nacout + ilen
      ENDIF

      IF (num_compute == 1)                                                    &
         WRITE( nuodr,'(''Q: ACT/REJ: P/Q/FLAG:'',2f8.0,2f6.1,2i2)')           &
                omlbdy(nodract,kact,nbtp)                                      &
               ,omlbdy(nodrrej,kclo,nbtp)                                      &
               ,omlbdy(nodract,kact,nbtqer)                                    &
               ,omlbdy(nodrrej,kclo,nbtqer), nflaga, nflagr
      zmlbdy (klev,nbtrh ) = omlbdy (nodrrej,kclo,nbtrh )
      zmlbdy (klev,nbtqer) = omlbdy (nodrrej,kclo,nbtqer)
      iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg), nvfqbp, nvfaoc     &
                                     , momlbd(nodrrej,kclo,nbtflg), nvfqbp )
      iomlhd (     nhqexi) = 1
    ENDIF
  ENDIF
 
  IF ((ABS( omlbdy(nodrrej,kclo,nbtp)                                          &
           -omlbdy(nodract,kact,nbtp)) <= rprlmt)) THEN
    IF (MAX( ABS( MAX( omlbdy(nodract,kact,nbtu) ,-c3600 )                     &
                 -MAX( omlbdy(nodrrej,kclo,nbtu) ,-c3600 ) )                   &
           , ABS( MAX( omlbdy(nodract,kact,nbtv) ,-c3600 )                     &
                 -MAX( omlbdy(nodrrej,kclo,nbtv) ,-c3600 ) )                   &
           , ABS( MAX( omlbdy(nodract,kact,nbtt) , c0 )                        &
                 -MAX( omlbdy(nodrrej,kclo,nbtt) , c0 ) )                      &
           , ABS( MAX( omlbdy(nodract,kact,nbtrh),-c1 )                        &
                 -MAX( omlbdy(nodrrej,kclo,nbtrh),-c1 ) ) ) > c0e) lsub =.FALSE.
  ENDIF

! 'while' loop: get all rejected levels which are more than 'rprlmt' [pa],
!               but less than 'rdplmt' [pa] above the active level 'kact',
!               and less than 'rprlmt' [pa] below the level 'kact+1'.
!               Of these levels, fill those observed quantities, which are not
!               present in either of the 2 active levels, into the temporary
!               field.
!        (Note: This may be important if e.g. a TEMP report containing only
!               temperature and humidity and a PILOT report have been derived
!               and disseminated from one and the same radiosonde ascent.)

  pnexlmt  =  c0
  IF (kact < momlhd(nodract,nhnlev))  pnexlmt  =  omlbdy(nodract,kact+1,nbtp)  &
                                                + rprlmt
  DO WHILE (      (   omlbdy(nodrrej,krej,nbtp)                                &
                   >= omlbdy(nodract,kact,nbtp)-rdplmt)                        &
            .AND. (   omlbdy(nodrrej,krej,nbtp) >= pnexlmt)                    &
            .AND. (.NOT. ltop) .AND. (klev < maxrsl))
     laddu = .FALSE.
     laddt = .FALSE.
     laddq = .FALSE.
     IF (lrepla) THEN
       lmand =     (rmod( omlbdy(nodrrej,krej,nbtp),10000.0_ireals ) < 2*epsy) &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)- 5000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-15000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-25000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-85000.0_ireals ) < epsy)   &
              .OR. ( ABS( omlbdy(nodrrej,krej,nbtp)-92500.0_ireals ) < epsy)
       laddu =        ((omlbdy(nodract,kact,nbtuer) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtuer) >= rmdich)     
       laddt =        ((omlbdy(nodract,kact,nbtter) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtter) >= rmdich)     
       laddq =        ((omlbdy(nodract,kact,nbtqer) <  rmdich) .OR. (lmand))   &
                .AND.  (omlbdy(nodrrej,krej,nbtqer) >= rmdich)     
     ENDIF
     IF ((lrepla) .AND. (kact < momlhd(nodract,nhnlev))) THEN
       IF (omlbdy(nodrrej,krej,nbtp) <= omlbdy(nodract,kact+1,nbtp)+rdplmt) THEN
         laddu = laddu .AND. (omlbdy(nodract,kact+1,nbtuer) < rmdich .OR. lmand)
         laddt = laddt .AND. (omlbdy(nodract,kact+1,nbtter) < rmdich .OR. lmand)
         laddq = laddq .AND. (omlbdy(nodract,kact+1,nbtqer) < rmdich .OR. lmand)
       ENDIF
     ENDIF
     IF ((laddu) .OR. (laddt) .OR. (laddq)) THEN
       klev = klev + 1
       zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
       iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
       IF (.NOT. laddu) zmlbdy (klev,nbtu  ) = rmdi
       IF (.NOT. laddu) zmlbdy (klev,nbtv  ) = rmdi
       IF (.NOT. laddu) zmlbdy (klev,nbtuer) = rmdi
       IF (.NOT. laddu) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvfubp,nvfaoc, 0,nvfubp)
       IF (.NOT. laddt) zmlbdy (klev,nbtt  ) = rmdi
       IF (.NOT. laddt) zmlbdy (klev,nbtter) = rmdi
       IF (.NOT. laddt) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvftbp,nvfaoc, 0,nvftbp)
       IF (.NOT. laddq) zmlbdy (klev,nbtrh ) = rmdi
       IF (.NOT. laddq) zmlbdy (klev,nbtqer) = rmdi
       IF (.NOT. laddq) iomlbd (klev,nbtflg) = ireplace( iomlbd(klev,nbtflg)   &
                                                       ,nvfqbp,nvfaoc, 0,nvfqbp)
       IF (zmlbdy(klev,nbtuer) > rmdich)  iomlhd (nhuexi) = 1
       IF (zmlbdy(klev,nbtter) > rmdich)  iomlhd (nhtexi) = 1
       IF (zmlbdy(klev,nbtqer) > rmdich)  iomlhd (nhqexi) = 1
     ENDIF
     IF (krej <  momlhd(nodrrej,nhnlev)) THEN
       krej = krej + 1
     ELSE
       ltop = .TRUE.
     ENDIF
     lsub = .FALSE.
  ENDDO

!  Close the loop over levels of body of active report

  ENDDO Levels_acr
! ~~~~~~~~~~~~~~~~

!  Get all rejected levels which are more than 'rdplmt' [pa] above
!  the highest active level. Fill those levels into the temporary field.

  DO WHILE ((.NOT. ltop) .AND. (klev < maxrsl))
     IF (lrepla) THEN
       klev = klev + 1
       zmlbdy (klev,1:mxrbdy) = omlbdy (nodrrej,krej,1:mxrbdy)
       iomlbd (klev,1:mxrbdf) = momlbd (nodrrej,krej,1:mxrbdf)
     ENDIF
     IF (krej <  momlhd(nodrrej,nhnlev)) THEN
       krej = krej + 1
     ELSE
       ltop = .TRUE.
     ENDIF
     lsub = .FALSE.
  ENDDO

! determine number of levels

  iomlhd (nhnlev) = klev

!  Replace / reject redundant data
!  -------------------------------

!  Copy rejected report to last ('nmlob') report for output

  IF ((nodrrej /= nmlob) .AND. (nmlvo < nmlob)) THEN
    omlbdy (nmlob ,1:maxrsl,1:mxrbdy) = omlbdy (nodrrej,1:maxrsl,1:mxrbdy)
    momlbd (nmlob ,1:maxrsl,1:mxrbdf) = momlbd (nodrrej,1:maxrsl,1:mxrbdf)
    momlhd (nmlob          ,1:mxrhdf) = momlhd (nodrrej         ,1:mxrhdf)
  ENDIF

!  Copy temporary field into ODR
 
  omlbdy (nmlvo ,1:maxrsl,1:mxrbdy) = zmlbdy (1:maxrsl,1:mxrbdy)
  momlbd (nmlvo ,1:maxrsl,1:mxrbdf) = iomlbd (1:maxrsl,1:mxrbdf)
  omlhed (nmlvo ,1:mxrhed) = zmlhed (1:mxrhed)
  momlhd (nmlvo ,1:mxrhdf) = iomlhd (1:mxrhdf)
  yomlhd (nmlvo )          = ystids

!  If the report is being auto-checked, proceed to check against other reports

  IF (nmlvo == nmlob)                                           CYCLE Next_mlrep
 
!  Get diagnostic array position (--> nobtpp, ncdtpp)
!  and update statistics and report events
 
  IF ((.NOT. lchecksub) .OR. (lsub)) THEN
    IF (nodrrej == nmlvo) THEN
      ncdtyp = ncdtyp2
      CALL obs_pointrs ( nobtyp2, ncdtyp2)
    ELSE
      CALL obs_pointrs ( nobtyp , ncdtyp )
    ENDIF
    IF ((npassiv < 2) .AND. (npassv2 < 2)) THEN
      noctac (nobtpp,ncdtpp) = noctac(nobtpp,ncdtpp) - 1
      neventr (neredn,ncdtpp,nobtpp) = neventr(neredn,ncdtpp,nobtpp) + 1
    ELSE
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) - 1
    ENDIF
!   IF ((lverif) .AND. (lverpas)) THEN
!   IF ((lverif) .AND. (lverpas) .AND. (.NOT. ((lsub) .AND. (lwprofl)))) THEN
    IF ((lverif) .AND. (lverpas) .AND. (.NOT. lsub)) THEN
      noctps (nobtpp,ncdtpp) = noctps(nobtpp,ncdtpp) + 1
    ELSE
      noctrj (nobtpp,ncdtpp) = noctrj(nobtpp,ncdtpp) + 1
    ENDIF
  ENDIF

!  Print message
 
  IF (((.NOT. lchecksub) .OR. (lsub)) .AND. (lprodr)) THEN
! IF (lprodr)                             THEN
    IF (nodrrej == nmlvo) THEN
      ncdtypa = ncdtyp
      ystida  = ystid
      ncdtypr = ncdtyp2
      ystidr  = ystid2
    ELSE
      ncdtypa = ncdtyp2
      ystida  = ystid2
      ncdtypr = ncdtyp
      ystidr  = ystid
    ENDIF
    ilen = 12 + 2*istrej + (momlhd(nmlvo,nhnlev) + momlhd(nmlob,nhnlev)) *11
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt11
      outbuf(nacout+ 3) = ncdtypr
      outbuf(nacout+ 4) = momlhd(nmlvo,nhitot)
      outbuf(nacout+ 5) = momlhd(nmlvo,nhjtot)
      outbuf(nacout+ 6) = NINT(omlhed(nmlvo,nhilon)*rscale)
      outbuf(nacout+ 7) = NINT(omlhed(nmlvo,nhjlat)*rscale)
      IF (omlhed(nmlvo,nhalt) > rmdich) THEN
        outbuf(nacout+8) = NINT(omlhed(nmlvo,nhalt))
      ELSE
        outbuf(nacout+8) = imdi
      ENDIF
      outbuf(nacout+ 9) = NINT(omlhed(nmlvo,nhtime)*rscale)
      outbuf(nacout+10) = ncdtypa
      outbuf(nacout+11) = momlhd(nmlvo,nhnlev)
      outbuf(nacout+12) = momlhd(nmlob,nhnlev)
      DO icl = 1 , istrej
        outbuf(nacout+12       +icl) = ICHAR( ystidr (icl:icl) )
        outbuf(nacout+12+istrej+icl) = ICHAR( ystida (icl:icl) )
      ENDDO
      nacout = nacout + 12 + 2*istrej
      DO klvi = 1, momlhd(nmlvo,nhnlev) + momlhd(nmlob,nhnlev)
        DO iiv = 1, 10
          outbuf (nacout+ iiv) = imdi
        ENDDO
        nactpr = nmlob
        klev   = klvi  
        IF (klvi > momlhd(nmlob,nhnlev)) THEN
          nactpr = nmlvo
          klev   = klvi - momlhd(nmlob,nhnlev)
        ENDIF
        IF (omlbdy(nactpr,klev,nbtu  ) > rmdich)                               &
          outbuf (nacout+ 1) = NINT(omlbdy(nactpr,klev,nbtu  )*rscale)
        IF (omlbdy(nactpr,klev,nbsv  ) > rmdich)                               &
          outbuf (nacout+ 2) = NINT(omlbdy(nactpr,klev,nbtv  )*rscale)
        IF (omlbdy(nactpr,klev,nbtt  ) > rmdich)                               &
          outbuf (nacout+ 3) = NINT(omlbdy(nactpr,klev,nbtt  )*rscale)
        IF (omlbdy(nactpr,klev,nbtrh ) > rmdich)                               &
          outbuf (nacout+ 4) = NINT(omlbdy(nactpr,klev,nbtrh )*rscale)
        IF (omlbdy(nactpr,klev,nbtp  ) > rmdich)                               &
          outbuf (nacout+ 5) = NINT(omlbdy(nactpr,klev,nbtp  )*rscale)
        IF (omlbdy(nactpr,klev,nbtz  ) > rmdich)                               &
          outbuf (nacout+ 6) = NINT(omlbdy(nactpr,klev,nbtz  )*rscale)
        IF (omlbdy(nactpr,klev,nbtuer) > rmdich)                               &
          outbuf (nacout+ 7) = NINT(omlbdy(nactpr,klev,nbtuer)*rscale)
        IF (omlbdy(nactpr,klev,nbtter) > rmdich)                               &
          outbuf (nacout+ 8) = NINT(omlbdy(nactpr,klev,nbtter)*rscale)
        IF (omlbdy(nactpr,klev,nbtqer) > rmdich)                               &
          outbuf (nacout+ 9) = NINT(omlbdy(nactpr,klev,nbtqer)*rscale)
        IF (omlbdy(nactpr,klev,nbtzer) > rmdich)                               &
          outbuf (nacout+10) = NINT(omlbdy(nactpr,klev,nbtzer)*rscale)
        outbuf(nacout+11) = momlbd(nactpr,klev,nbtlid)
        nacout  = nacout + 11
      ENDDO
    ENDIF 

  ENDIF ! lprodr

!  Reset last report to missing or passive (for verification of passive reports)

! IF ((lverif) .AND. (lverpas)) THEN
! IF ((lverif) .AND. (lverpas) .AND. (.NOT. ((lsub) .AND. (lwprofl)))) THEN
  IF ((lverif) .AND. (lverpas) .AND. (.NOT. lsub)) THEN
    IF (.NOT. lchecksub) THEN
      momlhd (nmlob ,nhpass)   = 2
      momlhd (nmlob ,nhschr)   = insert( momlhd(nmlob,nhschr) , 1 , nvrdbp )

      EXIT  Next_mlrep
    ENDIF
  ELSE
    omlbdy (nmlob ,1:maxrsl,1:mxrbdy) = rmdi
    momlbd (nmlob ,1:maxrsl,1:mxrbdf) = imdi
    omlhed (nmlob ,1:mxrhed) = rmdi
    momlhd (nmlob ,1:mxrhdf) = imdi
    yomlhd (nmlob )          = '        '
!  Adjust total number of vertical profiles
    nmlob = nmlob - 1

    EXIT  Next_mlrep
  ENDIF

! EXIT  Next_mlrep
  ENDDO Next_mlrep
! ~~~~~~~~~~~~~~~~

!-----------------------------------------------------------------------
!     nqlltar = momlhd (nmlob,nhqofl)
!     nredund = 0
!     IF (lredund) nredund = 1
!     nqlltar = insert( nqlltar , nredund , nf3rbp )
!     momlhd (nmlob,nhqofl) = nqlltar
!-----------------------------------------------------------------------
 
ELSE IF ( lairmlc)     THEN

  !-------------------------------------------------------------------------
  ! Section 3.0  Redundancy check for single-level aireps
  !              versus multi-level airep reports.
  !-------------------------------------------------------------------------

! (This is only meaningful, if the current report was found to be redundant,
!  but may contain data that are missing in a multi-level report)


!  Preset ODR index for check against all prior reports
 
  nmlvo  = 0
 
!  Get next multi-level airep report
!  ---------------------------------
 
  Next_airep : DO
 
!  Get ODR index of the multi-level report
 
  nmlvo  = nmlvo + 1
 
!  Leave 'loop' if no redundancy with any report
 
  IF (nmlvo == nmlob+1)                                          EXIT Next_airep
 
!  Leave 'loop' if only one of the two reports is fixed passive

  IF (momlhd(nmlvo,nhpass)/2 /= mosghd(nsglo,nhpass)/2)          EXIT Next_airep
 
!  Check if observation type is airep
 
  ncdtyp2  =  momlhd (nmlvo ,nhcode)
  lair2    =      (ncdtyp2 == naircd).OR.(ncdtyp2 == ncodar)                   &
              .OR.(ncdtyp2 == ncolba).OR.(ncdtyp2 == namdar)                   &
              .OR.(ncdtyp2 == nacar)
  IF (.NOT.lair2)                                               CYCLE Next_airep
  nobtyp2  =  nairep
 
!  Check 'station' id, i.e. flight number
!  --------------------------------------
 
  ystid2   =  yomlhd (nmlvo)

  IF (ystid2 /= ystid)                                          CYCLE Next_airep
 
!  Check for redundancy of current report with one level
!  of multi-level airep report
!  -----------------------------------------------------
 
!  Initialize level index
 
  klev = 0
 
!  Get level index of next level
 
  Lev_index  : DO

  klev = klev + 1
 
!  Leave 'loop' if no redundancy with present report
 
  IF (klev == momlhd(nmlvo,nhnlev)+1)                           CYCLE Next_airep
 
!  Get 4-d location of present level
 
  obtime2 = omlhed (nmlvo,nhtime)
  igplon2 = momlhd (nmlvo,nhitot)
  igplat2 = momlhd (nmlvo,nhjtot)

  obbalt2 = - omlbdy (nmlvo,klev,nbtp)
 
!  Check time difference
 
  lcotim =      ( ABS( obtime2 - obtime )  <=  rtmlmt )
  IF (.NOT.lcotim)                                               CYCLE Lev_index
 
!  Check horizontal and vertical (quasi-) colocation
 
  rhzob  = (cosolat *(igplon2-igplon) * dlon)**2 +                             &
                    ((igplat2-igplat) * dlat)**2
  rhzob  = rhzob * rdegkm2
  lcoloc =      (      rhzob               <=  rhzlmt )                        &
          .AND. ( ABS( obbalt2 - obbalt )  <=  rvtlmt )
  IF (.NOT.lcoloc)                                               CYCLE Lev_index

  EXIT Lev_index
  ENDDO Lev_index


!  Decide which report is redundant
!  --------------------------------
 
!  Get station correction bit of multi-level report
 
  nstchr2  =  momlhd(nmlvo,nhschr)
  nstcor2  =  ibits (nstchr2,nvscbp,nibits (nvscoc))

!  Determine which report is redundant
!  (note: current report is always redundant, c.f. ++++ )
 
  lreplace = .FALSE.
  IF (nstcor >  nstcor2) lreplace = .TRUE.
 
!  Replace / reject redundant data
!  -------------------------------
 
!  Replace data if current report is active (sta. corr.)
!                            (this is in fact excluded)
!                   or if level of multi-level report is not
!                         station correction and lacks data

  ilevin = momlbd (nmlob,klev,nbtlid)
  IF (lreplace) THEN
    omlbdy (nmlvo,klev,nbtp  ) = osgbdy (nsglo,nbsp  )
    omlbdy (nmlvo,klev,nbtz  ) = osgbdy (nsglo,nbsz  )
    omlbdy (nmlvo,klev,nbtzer) = osgbdy (nsglo,nbszer)
    momlbd (nmlvo,klev,nbtflg) = ireplace( momlbd(nmlvo,klev,nbtflg), nvfzbp   &
                                         , nvfaoc, mosgbd(nsglo,nbsflg), nvfzbp)
    omlbdy (nmlvo,klev,nbtzio) = osghed (nsglo,nhzio )
    omlbdy (nmlvo,klev,nbtzjo) = osghed (nsglo,nhzjo )
  ENDIF
  IF (     ((lreplace) .AND. (     (omlbdy(nmlvo,klev,nbtuer) <  rmdich)       &
                              .OR. (osgbdy(nsglo,nbsuer) >  rmdich)))          &
      .OR. (      (omlbdy(nmlvo,klev,nbtuer) <  rmdich)                        &
            .AND. (osgbdy(nsglo,nbsuer) >  rmdich))) THEN
    omlbdy (nmlvo,klev,nbtu  ) = osgbdy (nsglo,nbsu  )
    omlbdy (nmlvo,klev,nbtv  ) = osgbdy (nsglo,nbsv  )
    omlbdy (nmlvo,klev,nbtuer) = osgbdy (nsglo,nbsuer)
    momlbd (nmlvo,klev,nbtflg) = ireplace( momlbd(nmlvo,klev,nbtflg), nvfubp   &
                                         , nvfaoc, mosgbd(nsglo,nbsflg), nvfubp)
    IF (osgbdy(nsglo,nbsuer) >  rmdich) THEN
      momlhd (nmlvo,nhuexi) = 1
      momlhd (nmlvo,nhaexi) = 1
      ilevin = insert( ilevin , 1 , nlidbp(8) )
    ELSE
      ilevin = insert( ilevin , 0 , nlidbp(8) )
    ENDIF
  ENDIF
  IF (     ((lreplace) .AND. (     (omlbdy(nmlvo,klev,nbtter) <  rmdich)       &
                              .OR. (osgbdy(nsglo,nbster) >  rmdich)))          &
      .OR. (      (omlbdy(nmlvo,klev,nbtter) <  rmdich)                        &
            .AND. (osgbdy(nsglo,nbster) >  rmdich))) THEN
    omlbdy (nmlvo,klev,nbtt  ) = osgbdy (nsglo,nbst  )
    omlbdy (nmlvo,klev,nbtter) = osgbdy (nsglo,nbster)
    momlbd (nmlvo,klev,nbtflg) = ireplace( momlbd(nmlvo,klev,nbtflg), nvftbp   &
                                         , nvfaoc, mosgbd(nsglo,nbsflg), nvftbp)
    IF (osgbdy(nsglo,nbster) >  rmdich) THEN
      momlhd (nmlvo,nhtexi) = 1
      momlhd (nmlvo,nhaexi) = 1
      ilevin = insert( ilevin , 1 , nlidbp(9) )
    ELSE
      ilevin = insert( ilevin , 0 , nlidbp(9) )
    ENDIF
  ENDIF
  IF (     ((lreplace) .AND. (     (omlbdy(nmlvo,klev,nbtqer) <  rmdich)       &
                              .OR. (osgbdy(nsglo,nbsqer) >  rmdich)))          &
      .OR. (      (omlbdy(nmlvo,klev,nbtqer) <  rmdich)                        &
            .AND. (osgbdy(nsglo,nbsqer) >  rmdich))) THEN
    omlbdy (nmlvo,klev,nbtrh ) = osgbdy (nsglo,nbsrh )
    omlbdy (nmlvo,klev,nbtqer) = osgbdy (nsglo,nbsqer)
    momlbd (nmlvo,klev,nbtflg) = ireplace( momlbd(nmlvo,klev,nbtflg), nvfqbp   &
                                         , nvfaoc, mosgbd(nsglo,nbsflg), nvfqbp)
    IF (osgbdy(nsglo,nbsqer) >  rmdich) THEN
      momlhd (nmlvo,nhqexi) = 1
      momlhd (nmlvo,nhaexi) = 1
      ilevin = insert( ilevin , 1 , nlidbp(6) )
    ELSE
      ilevin = insert( ilevin , 0 , nlidbp(6) )
    ENDIF
  ENDIF
  momlbd (nmlob,klev,nbtlid) = ilevin
  IF ((lreplace).AND.(klev == 1)) THEN
    omlhed (nmlvo,nhilon) = osghed (nsglo,nhilon)
    omlhed (nmlvo,nhjlat) = osghed (nsglo,nhjlat)
    omlhed (nmlvo,nhalt ) = osghed (nsglo,nhalt )
    omlhed (nmlvo,nhtime) = osghed (nsglo,nhtime)
    omlhed (nmlvo,nhzio ) = osghed (nsglo,nhzio )
    omlhed (nmlvo,nhzjo ) = osghed (nsglo,nhzjo )
    momlhd (nmlvo,nhio  ) = mosghd (nsglo,nhio  )
    momlhd (nmlvo,nhjo  ) = mosghd (nsglo,nhjo  )
    momlhd (nmlvo,nhobtp) = mosghd (nsglo,nhobtp)
    momlhd (nmlvo,nhcode) = mosghd (nsglo,nhcode)
    momlhd (nmlvo,nhschr) = mosghd (nsglo,nhschr)
    momlhd (nmlvo,nhqofl) = mosghd (nsglo,nhqofl)
    yomlhd (nmlvo)        = yosghd (nsglo)
  ENDIF
 
!  If current single-level report was redundant,
!  replace it by last single-level report
!  and set 'last' single-level report to missing
!  
!  note: this is already done in section 1
 
!  Get diagnostic array position (--> nobtpp, ncdtpp)
!  and update statistics and report events
!
!  note: this is already done in section 1

!  Print message
 
  IF (lprodr) THEN
    ilen = 13 + 2*istrej
    IF (nacout+ilen <= nmxoln) THEN
      outbuf(nacout+ 1) = ilen
      outbuf(nacout+ 2) = nfmt12
      DO iiv=1,5
      IF (omlbdy(nmlvo,klev,iiv) > rmdich) THEN
        outbuf(nacout+ 2+iiv) = INT(omlbdy(nmlvo,klev,iiv)*rscale)
      ELSE
        outbuf(nacout+ 2+iiv) = imdi
      ENDIF
      ENDDO
      outbuf(nacout+ 8) = momlhd(nmlvo,nhcode)
      outbuf(nacout+ 9) = INT(omlbdy(nmlvo,klev,nbtzio)*rscale)
      outbuf(nacout+10) = INT(omlbdy(nmlvo,klev,nbtzjo)*rscale)
      outbuf(nacout+11) = momlhd(nmlvo,nhio)
      outbuf(nacout+12) = momlhd(nmlvo,nhjo)
      outbuf(nacout+13) = INT(omlhed(nmlvo,nhtime)*rscale)
      DO icl = 1 , istrej
        outbuf(nacout+13+icl) = ICHAR( ystid (icl:icl) )
      ENDDO
      nacout  = nacout + ilen
    ENDIF

!   IF (num_compute == 1)                                                      &
!     WRITE( nuodr,'('' --> REDUNDANT AIREP: ACTIVE:'',3f7.1,f5.2,f8.0         &
!                  &,1x, a5, i4, 2f5.1, 2i4, f5.1)' )                          &
!            (omlbdy(nmlvo,klev,iiv), iiv=1,5)  , ystid                        &
!          , momlhd(nmlvo,nhcode)                                              &
!          ,(omlbdy(nmlvo,klev,iiv), iiv=nbtzio,nbtzjo)                        &
!          ,(momlhd(nmlvo,iiv),      iiv=nhio  ,nhjo  )                        &
!          , omlhed(nmlvo,nhtime)
!   ENDIF
  ENDIF
!  Adjust total number of single level reports
!  if current single-level report was rejected
!  note: this is already done in section 1.
 
  EXIT  Next_airep
  ENDDO Next_airep
 
ENDIF     ! lvprof/lairml

!-------------------------------------------------------------------------------
! End Subroutine obs_redundancy
!-------------------------------------------------------------------------------
END SUBROUTINE obs_redundancy

!-------------------------------------------------------------------------------

END MODULE src_obs_proc_aof
