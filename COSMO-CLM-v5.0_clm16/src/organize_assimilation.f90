!+ External procedure for organizing the assimilation
!-------------------------------------------------------------------------------

SUBROUTINE organize_assimilation (yaction, ierror, yerrmsg)

!-------------------------------------------------------------------------------
!
! Description:
!   This procedure is the driving routine for calling the different
!   modules of the assimilation. It is called first at the beginning of the 
!   program just to initialize the NAMELIST variables. Later it is 
!   called during the time-stepping.
!
! Method:
!
! Current Code Owner: DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  christoph.schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.39       2000/05/03 Ulrich Schaettler
!  Initial release
! 2.4        2001/01/29 Christoph Schraff
!  Namelist variables 'khumbal' and 'qcfpst' added, 'qctf' and 'qctfsu' removed.
!  Most of the defaults of the namelist parameters set to operational values.
! 2.5        2001/06/01 Christoph Schraff
!  New namelist variable 'yaofpath' for the path(-name) of input obs. file AOF.
! 2.12       2001/11/07 Christoph Schraff
!  Bug correction at check of AOF path name.
! 2.13       2002/01/18 Michael Buchhold
!  Additional namelist parameters defining active / passive status of
!  wind profiler and SODAR/RASS reports
! 2.19       2002/10/24 Michael Buchhold + Christoph Schraff
!  - Additional namelist parameters defining active / passive status of
!    Radar VAD wind reports.
!  - New namelist parameters defining analysis times of surface wind speed
!  - Reading new namelist variables: 'qgeotop', 'vpblsu', 'loiqv2m', 'lqfqv2m'
! 3.3        2003/04/22 Christoph Schraff
!  New namelist variables 'gnudggp', 'maxgpo', 'lgps', 'lcd096', and 'lgpsbias'
!  for assimilation of ground-based GPS-derived integrated water vapour data.
! 3.6        2003/12/11 Christoph Schraff
!  New namelist variable 'rhfgps' for horizontal correlation scale for GPS data.
!  Default for 'tconbox' doubled.
!  Checked the IOSTAT-value for reading the namelist group /nudging/
! 3.7        2004/02/18 Ulrich Schaettler
!  Introduced handling for unit numbers of ascii files
! 3.13       2004/12/03 Klaus Stephan, Stefan Klink
!  Introduced latent heat nudging
! 3.14       2005/01/25 Klaus Stephan
!  Introduced namelist variable rlhn_scale_dp
! 3.16       2005/07/22 Ulrich Schaettler
!  Bug Correction: rlhn_scale_dp was erroneously defined as Integer
! 3.18       2006/03/03 Christoph Schraff
!  New namelist variables 'kmultw' for multiple observations, 'lseparw' for
!  multiple observing systems, 'qcciq' and 'qcsiq' for QC thresholds of IWV,
!  and 'rhfpsdd' for scaling horizontal correlation scale for surface prssure
!  depending on data density.
!   and from Klaus Stephan
!  LHN namelist parameter moved to data_lheat_nudge to avoid to many dependencies
!  New namelist parameter for integrated precipitation flux and radar blacklist
! 3.19       2006/04/25 Klaus Stephan
!  Some adaptations for Latent Heat Nudging
! 3.21       2006/12/04 Klaus Stephan
!  Adaptations and new Namelist variables for Latent Heat Nudging
! V3_23        2007/03/30 Klaus Stephan
!  Modifications for Namelist Input for Latent Heat Nudging
! V3_24        2007/04/26 Klaus Stephan
!  Eliminated Namelist variables lhn_diagprec, rlhn_scale_dp for Latent Heat Nudging
! V3_28        2007/07/18 Klaus Stephan
!  New Namelist variables for LHN:
!  - data input: noobs_date,n_noobs
!  - reference precipitaion: rqrsgflux
!  - vertical restriction: ktop_temp
! V4_1         2007/12/04 Klaus Stephan
!  New Namelist variables for LHN:
!  - use radar beam height: lhn_height, height_file
!  - bright band detection algorithm: lhn_bright
! V4_2         2007/12/05 Klaus Stephan
!  Changed default value for file with radar height
! V4_5         2008/09/10 Christoph Schraff
!  To allow for reading from NetCDF observation input files instead of an AOF
!  file: new namelist variables 'itype_obfile' and 'ycdfdir'.
! V4_12        2010/05/11 Daniel Leuenberger
!  - removed namelist switch lhn_radar_dx
!  - added namelist switch lhn_spqual       (default value: false)
!  - added namelist lhn_dt_obs              (default vaule: 5min)
!  - added namelist nradar                  (default value: 32)
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_15        2010/11/19 Klaus Stephan
!  New namelist variables:
!  - 'mqcorr92': switch for bias correction of Vaisala RS92 radiosonde humidity
!  - set 'nradar_d' to 33, considering the new radar stattion MEM
! V4_19        2011/08/01 Ulrich Schaettler
!  Introduced conditional compilation for GRIBDWD: check settings of llhn, lsurfa
! V4_22        2012/01/31 Christoph Schraff
!  New namelist variables:
!  - 'mveripr' : switches (code) for writing VOF and/or feedobs files
!  - 'mpsgcor' : switch for applying geostrophic pressure increments
!  - 'qgeops'  : reduction factor to the geostrophic pressure increments
!  - 'lscatt'  : for active use of observation type 'scatterometer'
!  - 'maxmlv'  : max. number of observation levels in multi-level report
!  - 'mxfrep'  : max. number of reports in NetCDF feedobs file
!  - 'mxfobs'  : max. number of observations in NetCDF feedobs file
!  - for use of new observation code types: 'lcd037' (TEMP MOBILE),
!    'lcd038' (PILOT MOBILE), 'lcd140' (METAR, not thoroughly tested),
!    'lcd123' (scatterometer ASCAT), and 'lcd122' (QuickScat)
!  - 'igpscen' : GPS processing centres assigned as active obs code types,
!                ordered according to their preference in the redundancy check.
!                Default means that no GPS processing centre is used actively.
!  - for new weighting of multiple observations and observation systems:
!    'nwtyp'   : if > 1 then compute net obs. increments for 'nwtyp' different
!                sets of observing systems separately
!    'niwtyp'  : number of observing systems for each set of obs systems
!    'iwtyp'   : defines these sets of observation systems
!    'kwtyp'   : defines the multiple weighting for each of these sets
!  Other changes to the namelist:
!  - Namelist parameters 'lseparw' and 'kmultw' are removed.
!  - Modified values for (and, due to modified 'fdoro', meaning of) 'doromx'.
!  Preparation for 1DVAR and the processing of satellite data and production of
!    retrievals. Several new variables are defined which will become namelist
!    variables only when 1DVAR will be included in the official version:
!    - 'l1dvar' : general on/off switch for 1DVar to derive satellite retrievals
!    - 'maxtvo' : ODR size for satellite retrievals
!    - 'gnudgtv' : nudging coefficients for satellite retrievals
!    - 'rhfrtv'  : factor for horizont. weight function for satellite retrievals
!    - 'mcdmsg1', 'mcdmsg2', 'mcdno15', 'mcdno16, 'mcdno17', 'mcdno18' : for use
!                  of new obs types MSG-1 or -2 resp. NOAA-15, -16, -17, or -18.
!    Note: An experimental version for use of satellite data, based on V4_18,
!          is available from christoph.schraff_at_dwd.de .
! V4_26        2012/12/06 Ulrich Schaettler
!  Increased icdfdirlen to 250 (length of directory name ycdfdir)
! V4_28        2013/07/12 Christoph Schraff, Klaus Stephan, Ulrich Schaettler
!  New namelist variables (by CS):
!  - 'qcflbcp' : factor to threshold used for check of ps-obs against LBC fields
!  - 'irun_osse': if feedback file 'fof' present: model run used as obs for OSSE
!  - 'losse_fg' : first guess check flag from 'fof' converted to 'dataset flag'
!  - 'fperturb' : factor to perturb (simulated or original) obs from 'fof' file
!  - 'iseed'    : external seed for random number generator.
!  Default value of nradar set to 35. (KS)
!  New namelist variables for near surface analysis (US)
!  - 'ydir_lansfc' : to specify directory where to write the files
!  - 'yform_lansfc': to specify format of the files ('grb1','api1','api2')
! V4_29        2013/10/04 Davide Cesari
!  Corrected a defined-statement for GRIBDWD/GRIBAPI
! V5_00_clm9   2016/05/11 H.-J. Panitz, IMK/KIT
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!   adapted from COSMO_5.1
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================

USE data_parameters,    ONLY:  ireals, iintegers
USE data_parallel,      ONLY:  my_cart_id, lcompute_pe, nproc, imp_reals,      &
                               imp_integers, imp_logical, imp_character,       &
                               icomm_cart, my_world_id, icomm_world,           &
                               intbuf, realbuf, logbuf, charbuf
USE data_modelconfig,   ONLY:  dt, ie_tot, je_tot, ke_tot
USE data_runcontrol,    ONLY:  nuspecif, nstop, lout_anai,                     &
    luseobs      ,& ! on - off switch for using observational data for
    luse_rttov      ! on - off switch for producing synthetic satellite images

!-------------------------------------------------------------------------------

USE data_nudge_all  , ONLY :   &

! 1. parameters and related variables
! -----------------------------------

    icdfdirlen   ,& ! max. length of name of directory where
                    !   NetCDF observation input files reside
    mxda         ,& ! max. number of data pts 'always' / never to be nudged
    mxgpc        ,& ! max. number of GPS processing centres used
    mxwt         ,& ! max. number of sets of obs systems for which
                    !      net obs. increments are computed separately
    mxtyw        ,& ! max. number of obs or code types defining sets of obs syst

! 2. namelist variables controlling the data assimilation
! -------------------------------------------------------

    lnudge       ,& ! .f. : on - off switch for nudging
    lverif       ,& ! .f. : on - off switch for verification
    lverpas      ,& ! .t. : on - off switch for verif. also of passive reports
!   lobdens      ,& ! .f. : if TRUE then compute obs. density at obs. locations
    nwtyp        ,& ! 1   : if > 1 then compute net obs. increments for 'nwtyp'
                    !       different sets of observing systems separately
    niwtyp       ,& ! 1,0,0,..: number of obs or code types which belong to
                    !           the sets of observing systems
    iwtyp        ,& ! 0,0,0,..: obs or code types belonging to a set of obs
                    !           system, specified successively for each set
    kwtyp        ,& ! 1   : function for weights W for multiple observations
    mruntyp      ,& ! -1  : type of current model run used for increments in VOF
    mveripr      ,& ! 3   : type of verification/observation file(s) written
    nudgsta      ,& ! 0   : start of nudging period in timesteps
    nudgend      ,& ! 0   : end of nudging period in timesteps
    nversta      ,& ! 0   : start of verification period in timesteps
    nverend      ,& ! 0   : end of verification period in timesteps
    hversta      ,& ! 0   : start of verification period in 'model integ. hours'
    hverend      ,& ! 0   : end of verification period in 'model integr. hours'
    tconbox         ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)

USE data_nudge_all  , ONLY :   &

    luvgcor      ,& ! .t. : .t. ==> geostrophic wind correction applied
    mpsgcor      ,& ! 1   : mode to apply geostrophic pressure correction
    ntpscor      ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (i.e. p. at lowest model level) is nudged
    khumbal      ,& ! 100 : range around convectively precipitating grid pts, at
                    !       which specified (not relative) humidity is preserved
    ptpstop      ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    qgeo         ,& ! .3  : factor to the geostrophic wind increments at 1000hPa
    qgeotop      ,& ! .5  : factor to the geostrophic wind increments at ptpstop
    qgeops       ,& ! .9  : factor to the geostrophic pressure increments
    gnudg        ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar      ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgsu      ,& ! 6,12,0,6*10^-4: nudging coef. for surface-level data [1/s]
    gnudggp      ,& !        0*10^-4: nudging coef. for GPS-derived IWV [1/s]
    ltipol       ,& ! .t. : .t. ==> linear interpolation in time of upper-air
                    !               data which are less than 'tipolmx' hrs apart
    ltipsu       ,& ! .t. : .t. ==> linear interpolation in time of surface-lev.
                    !               data which are less than 'tipmxsu' hrs apart
    tipolmx      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       for upper-air data (set tipolmx = 0, if .NOT ltipol)
    tipmxsu      ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       of surface + GPS data  (tipmxsu = 0, if .NOT ltipsu)
    wtukrsa      ,& ! 3.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukrse      ,& ! 1.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukara      ,& ! 1.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for AIRCRAFT data
    wtukare      ,& ! 0.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for AIRCRAFT data
    wtuksua      ,& ! 1.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       rel. to the obs. time   for surface-level + GPS data
    wtuksue         ! 0.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       rel. to the obs. time   for surface-level + GPS data

USE data_nudge_all  , ONLY :   &

    msprpar      ,& ! 1   : switch specifying the surfaces along which observat.
                    !       increments of upper-air data are (primarily) spread
    msprpsu      ,& ! 0   : switch specifying the surface along which surface-
                    !       level data increments are primarily spread
    vcorls       ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
    vcutof       ,& ! 2* .75,: cut-off of the vertical correlation.       Units:
                    ! 2*1.     value of correlation at cut-off is [exp(-vcutof)]
    wablua       ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    vcorlsu      ,& ! 2*.013,: square of the vertical correlation scale,
                    ! 2*.002   i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpsu == 1) or
                    !          in potential temperature if (msprpsu == 2)
    vcutosu      ,& ! 2* .75 : cut-off of the vertical correlation.       Units:
                    ! 2*4.     value of correlation at cut-off : [exp(-vcutosu)]
    vpblsu       ,& ! 2*99. ,: Gaussian vertical influence 'radius' of potential
                    ! 2*99.    temperature differences between obs. level and
                    !          model level, to render vertical weights depend on
                    !          the stability (in the PBL) even if (msprpsu <= 1)
    wablsu       ,& ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
    lsvcorl      ,& ! .t. : .t. ==> adjustment of vertical correlation scales
                    !               in the presence of close observations
    rhinfl       ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0., 0.   autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac       ,& ! 1., 0.,: multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac       ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    rhfgps       ,& ! 0.45   : scaling (reduction) factor of the total 'COSAC'
                    !          for humidity derived from GPS IWV
    rhfpsdd      ,& ! 1.0    : minimum scaling (reduction) factor of the total
                    !          'COSAC' for surface pressure dep. on data density
    cutofr       ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni        ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which obs. increments are spread laterally
    rhiflsu      ,& ! 2 *70.,: constant part of the 'correlation scale of the
                    !   100.,  autoregressive horiz. correlation function'
                    !    70.   (='COSAC') [km]  (at the obs. time)
    rhtfsu       ,& ! 1.,    : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.,      relative to the 'COSAC' at the obs. time as given
                    ! 1.       by 'rhiflsu')
    cutofsu      ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    vcsnisu      ,& ! 2*2500.: square of Gaussian vertical influence 'radius'
                    ! 2*   9.  in potential temperature (if msprpsu <= 1) or
                    !          log( pressure ) (if msprpsu == 2) on surfaces
                    !          along which obs. increments are spread laterally
    cnondiv      ,& ! .1  : constant part of the factor to the non-divergence
                    !       correction ('fanodicor') in 2-dim wind correlations
    fnondiv      ,& ! .8  : multiplication factor to the vertically varying part
                    !       of the factor to the non-divergence correction
                    !       (as defined in 'data_nudge_all')
    tnondiv      ,& ! 1.1 : temporal multiplication factor to the factor to the
                    !       non-divergence correction for the beginning and end
                    !       of the nudging period for 1 obs relative to that
                    !       given by 'cnondiv', 'fnondiv' for the obs. time
    lscadj       ,& ! .T.,   : .F. ==> linear vertical interpolation (in log(p))
                    ! .T.,     instead of vertical scale adjustment (by vertical
                    ! .T.,     averaging over the model layer) for conveying the
                    ! .F.      observational information to the model levels
    topobs       ,& !  849., : threshold [hPa]: above this level (p < topobs),
                    ! 1099.,   only obs. increments at model levels are used,
                    !  799.,699.  i.e. obs. incr. at obs. levels are not used
    botmod       ,& ! 3*1099.: threshold [hPa]: below this level (p > botmod),
                    !    899.  only obs. increments at obs. levels are used,
                    !          i.e. obs. incr. at model levels are not computed
    loiqv2m      ,& ! .f. : .t. ==> 2-m humidity observation increments as
                    !               differences of specific humidity
    lqfqv2m      ,& ! .f. : .t. ==> quality factor for 2-m humidity observations
                    !               dependent on T-2m differences
    dtqc         ,& ! 720.: timestep (in [s]) for the threshold quality control
    qcc          ,& !  0.,500: constant parts of the quality control thresholds
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf         ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (as def. in 'data_nudge_local',
                    !          not available for pressure ps and humidity RH)
    qccsu        ,& ! 12.,500: constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcciq        ,& ! 1.  : constant part of QC threshold for IWV
    qcsiq        ,& ! .15 : IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile
    qcflbcp      ,& ! 1.4 : enhancement factor to the threshold used for the
                    !       check of ps against lateral boundary fields
    qcfpst          ! 1.5 : maximum enhancement of the weight for surface
                    !       pressure observations due to the pressure tendency

USE data_nudge_all  , ONLY :   &

    itype_obfile ,& !   1 : type of observation input file(s): 1: AOF, 2: NetCDF
    irun_osse    ,& !   0 : model run to derive obs values from file yfofin='fof'
    losse_fg     ,& ! .f. : f.g.check flag from 'fof converted to 'dataset flag'
    ycdfdir      ,& ! '/.'   : directory with NetCDF obs input + blacklist files
    yaofpath     ,& ! 'aof'  : path(name) of input observation file 'AOF'
    obnlat       ,& !   90.  : northern boundary of observation area
    obslat       ,& !  -90.  : southern boundary of observation area
    obwlon       ,& ! -180.  : western boundary of observation area
    obelon       ,& !  180.  : eastern boundary of observation area
    exnlat       ,& !   90.  : northern boundary for exclusion area
    exslat       ,& !  -90.  : southern boundary for exclusion area
    exwlon       ,& ! -180.  : western boundary for exclusion area
    exelon       ,& !  180.  : eastern boundary for exclusion area
    doromx       ,& !  100., : vertical extrapolation cut-off and gaussian
                    !  150.,   radius of height differences between model
                    !  150.,   orography and surface station height for a factor
                    !  150.,   contributing to the quality weight factor as part
                    !          of the nudging weights
    altopsu      ,& !  100., : SYNOP obs. above height 'altopsu' are not assimi-
                    ! 3*5000.  lated. If (altopsu == 0.) then SYNOP / surf. TEMP
                    !          assigned to land grid pts. are not assimilated
    thairh       ,& !   20.  : maximum horizontal distance [km] between the
                    !          lowest report and any single level report that
                    !          is added to a multi-level AIRCRAFT report
    fperturb     ,& !    0.  : factor to the obs error variances to define the
                    !          size of random perturbations added to the obs
                    !          (only for data from feedback file yfofin='fof')
    lgpsbias     ,& !  .f.   : .t. ==> bias correction to GPS IWV applied
    maxmlo       ,& !  300   : max. number of multi-level reports in the ODR
    maxsgo       ,& ! 3000   : max. number of (surface-level and upper-air)
                    !                         single-level reports in the ODR
    maxuso       ,& !  900   : max. number of upper-air single-level rep. in ODR
    maxgpo       ,& ! 3000   : max. number of GPS reports in ODR on total domain
    maxmlv       ,& !  100   : max. number of obs levels in multi-level reports
    mxfrep       ,& !   -1   : max. number of reports in NetCDF feedobs file
    mxfobs       ,& !   -1   : max. number of observations in feedobs file
    nolbc        ,& !    5   : number of grid rows at lateral boundaries
                    !          where obs are neglected
    mqcorr92     ,& !    0   : switch for bias correction for Vaisala RS92
                    !            radiosonde humidity
                    !          = 0 : no correction for humidity
                    !          = 1 : correct only solar radiation bias
                    !          = 2 : correct total bias (incl. nighttime bias)
    iseed        ,& !    0   : external seed for random number generator
    mrhyes       ,& ! mxda*0 : station ID's of SURFACE humidity data that are
                    !          always used (irrespective of 'doromx', 'altopsy)
    mrhno        ,& ! mxda*0 : station ID's of SURFACE humidity data never used
    muvyes       ,& ! mxda*0 : station ID's of SURFACE wind data always used
    muvno        ,& ! mxda*0 : station ID's of SURFACE wind data never used
    igpscen         ! X* -1  : array of processing centres of GPS data

USE data_nudge_all  , ONLY :   &

    lsynop       ,& ! .t.    : .t. if SYNOP data is used
    laircf       ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob       ,& ! .false.: .t. if SATOB data is used
    ldribu       ,& ! .t.    : .t. if BUOY  data is used (drifting buoy)
    ltemp        ,& ! .t.    : .t. if TEMP  data is used
    lpilot       ,& ! .t.    : .t. if PILOT data is used
    lsatem       ,& ! .false.: .t. if SATEM data is used
    lgps         ,& ! .false.: .t. if GPS   data is used
    lscatt       ,& ! .t.    : .t. if SCATT data is used (scatterometer)
    lcd011       ,& ! .t.    : .t. if synop code  11 data is used (land synop)
    lcd014       ,& ! .t.    : .t. if synop code  14 data is used (automatic)
    lcd021       ,& ! .t.    : .t. if synop code  21 data is used (ship)
    lcd022       ,& ! .t.    : .t. if synop code  22 data is used (ship abbrev.)
    lcd023       ,& ! .t.    : .t. if synop code  23 data is used (shred)
    lcd024       ,& ! .t.    : .t. if synop code  24 data is used (autom. ship)
    lcd140       ,& ! .t.    : .t. if synop code 140 data is used (metar)
    lcd041       ,& ! .t.    : .t. if airep code  41 data is used (codar)
    lcd141       ,& ! .t.    : .t. if airep code 141 data is used (airep)
    lcd241       ,& ! .t.    : .t. if airep code 241 data is used (colba)
    lcd144       ,& ! .t.    : .t. if airep code 144 data is used (amdar)
    lcd244       ,& ! .t.    : .t. if airep code 244 data is used (acars)
    lcd088       ,& ! .t.    : .t. if satob code  88 data is used (satob)
    lcd188       ,& ! .false.: .t. if satob code 188 data is used (sst)
    lcd063       ,& ! .t.    : .t. if dribu code  63 data is used (bathy)
    lcd064       ,& ! .t.    : .t. if dribu code  64 data is used (tesac)
    lcd165          ! .t.    : .t. if dribu code 165 data is used (drift. buoy)

USE data_nudge_all  , ONLY :   &

    lcd035       ,& ! .t.    : .t. if temp  code  35 data is used (land temp)
    lcd036       ,& ! .t.    : .t. if temp  code  36 data is used (temp ship)
    lcd037       ,& ! .t.    : .t. if temp  code  37 data is used (mobile)
    lcd135       ,& ! .t.    : .t. if temp  code 135 data is used (dropsonde)
    lcd039       ,& ! .t.    : .t. if temp  code  39 data is used (rocob)
    lcd040       ,& ! .t.    : .t. if temp  code  40 data is used (rocob ship)
    lcd032       ,& ! .t.    : .t. if pilot code  32 data is used (land pilot)
    lcd033       ,& ! .t.    : .t. if pilot code  33 data is used (pilot ship)
    lcd038       ,& ! .t.    : .t. if pilot code  38 data is used (mobile)
    lcd132       ,& ! .t.    : .t. if pilot code 132 data is used (win-prof eu)
    lcd133       ,& ! .t.    : .t. if pilot code 133 data is used (sod/rass eu)
    lcd136       ,& ! .t.    : .t. if pilot code 136 data is used (pro/rass us)
    lcd137       ,& ! .t.    : .t. if pilot code 137 data is used (Radar VAD)
    lcd086       ,& ! .t.    : .t. if satem code  86 data is used (satem)
    lcd186       ,& ! .t.    : .t. if atovs code 186 data is used (hi-res ATOVS)
    lcd122       ,& ! .t.    : .t. if scatt code 122 data is used (QuickScat)
    lcd123       ,& ! .t.    : .t. if scatt code 123 data is used (ASCAT)
    lcd096          ! .t.    : .t. if gps data from COST ASCII file is used

USE data_nudge_all  , ONLY :   &

    lsurfa       ,& ! .f.    : .t. if surface fields are analysed
    lt2m         ,& ! .f.    : .t. if 2m temperat. field is analysed
    lrh2m        ,& ! .f.    : .t. if 2m rel. hum. field is analysed
    lprecp       ,& ! .f.    : .t. if precipitation is analysed
    lff10m       ,& ! .f.    : .t. if 10m wind speed is analysed
    ht2a         ,& ! 999.   : time of 1. T2m-ana in hours since model start
    ht2i         ,& ! 999.   : time increment to next T2m analysis
    hh2a         ,& ! 999.   : time of 1. RH2m-ana in hours since model start
    hh2i         ,& ! 999.   : time increment to next RH2m analysis
    hffa         ,& ! 999.   : time of 1. 10m wind-ana in hrs since model start
    hffi         ,& ! 999.   : time increment to next wind analysis
    hprc         ,& ! 999.   : time of prec-ana in hours since model start
    raintp       ,& ! 12.    : time period of precipitation analysis
    ydir_lansfc  ,& ! './'   : directory where to write the 2-D analyses
    yform_lansfc ,& ! 'grb1' : format for the 2-D analyses files
    lpraof       ,& ! .f.    : .t. for diagnostic print of analysis obs file AOF
    lprodr       ,& ! .t.    : .t. for diagnostic print of obs data records ODR
    ldiasa       ,& ! .f.    : .t. for diagnostics of surface analysis
    ionl         ,& ! 167    : / grid point coordinates
    jonl         ,& ! 103    : \ for standard output on nudging
    ionl2        ,& ! 167    : / 2nd grid pt coordinates
    jonl2        ,& ! 103    : \ for other standard output on nudging
    noctrq       ,& !   9    : observation/code type of diagnostic print
                    !          (if (noctrq == 9) then print all obs./code types)
    dinlat       ,& !   0.   : northern boundary of diagnostic print area
    dislat       ,& !   0.   : southern boundary of diagnostic print area
    diwlon       ,& !   0.   : western boundary of diagnostic print area
    dielon          !   0.   : eastern boundary of diagnostic print area

USE data_nudge_all  , ONLY :   &

    l1dvar       ,& ! .f. : on - off switch for 1DVar
    gnudgtv      ,& ! 0, 0,6,6*10^-4: nudging coeffic. for sat retrivals [1/s]
    rhfrtv       ,& ! 1., 1.,: scaling (reduction) factor of the total 'COSAC'
                    ! 1., 0.5  for temperature / humidity satellite retrievals
    maxtvo       ,& !    1   : max. number of sat retrievals within total domain
    mcdmsg1      ,& !  0     : processing / use of MSG1   code  71 data
    mcdmsg2      ,& !  0     : processing / use of MSG2   code  72 data
    mcdno15      ,& !  0     : processing / use of NOAA15 code 206 data
    mcdno16      ,& !  0     : processing / use of NOAA16 code 207 data
    mcdno17      ,& !  0     : processing / use of NOAA17 code 208 data
    mcdno18      ,& !  0     : processing / use of NOAA18 code 209 data

    nsfc         ,& ! information about the surface analysis
    nuaof        ,& ! AOF (input): analysis observation file
    nugps        ,& ! GPS observation file unit
    nuaofex      ,& ! (output:) expanded AOF

    isetyp0      ,& ! index in 'niwtyp' which points to the first index with
                    ! observation type '0' in 'iwtyp' (denoting the set of obs
                    ! systems with contains all ('remaining') obs types that are
                    ! not specified explicitly in 'iwtyp')
    isetyp          ! defines for each observation or code type in 'iwtyp' which
                    ! set of observing systems (index of 'niwtyp') it belongs to

! end of data_nudge_all

USE data_obs_lib_cosmo, ONLY :   &

    nusatin      ,& ! satellite meta data input files (bias correction,
                    !                                  channel selection)
    nucautn      ,& ! caution messages if too many obs for current ODR size
    nuqc         ,& ! data rejected by threshold quality control at obs time
    nustat       ,& ! statistics of processed reports
    nurej        ,& ! direct reporting of rejected obs. reports
    nuodr        ,& ! observations stored in the observation data record
    nuverif      ,& ! VOF (output): verification observation file (observations
                    !   incl. quality control flag for verification)
    nupr         ,& ! all the remaining information
    epsy         ,& ! = 1.E-8_ireals : commonly used very small value > 0
    mxobtp          ! number of observation types

! end of data_obs_lib_cosmo

USE data_lheat_nudge  , ONLY :   &

! Namelist variables controlling the latent heat nudging
! ---------------------------------------------------------

    llhn           ,&  ! on/off switch for latent heat nudging (lhn)
    llhnverif      ,&  ! .f. : on - off switch for verification
    lhn_search     ,&  ! search for appropriate nearby model heating profile
    lhn_filt       ,&  ! horizontal filtering of lhn t-increments
    lhn_relax      ,&  ! vertical filtering of lhn t-increments
    lhn_limit      ,&  ! impose absolute limit on lhn t-increments (abs_lhn_lim)
    lhn_hum_adj    ,&  ! apply a humidity adjustment along with t-increments
    lhn_spqual     ,&  ! switch for the use of a spatial quality function
    lhn_black      ,&  ! use blacklist for radar data
    lhn_incloud    ,&  ! apply the LHN-scaling in cloudy layers only
    lhn_diag       ,&  ! produce more detailed diagnostic output during lhn
    lhn_qrs        ,&  ! calculate the integraed precipitation flux
    lhn_logscale   ,&  ! apply logarithmic scaling factors
    lhn_wweight    ,&  ! apply a weighting with respect to the mean horizontal wind
    nlhn_start     ,&  ! start of latent heat nudging period in timesteps
    nlhn_end       ,&  ! end of latent heat nudging period in timesteps
    nlhnverif_start,&  ! start of latent heat nudging period in timesteps
    nlhnverif_end  ,&  ! end of latent heat nudging period in timesteps
    noobs_date     ,&  ! date (yyyymmddhhxx) of missing observations
    n_noobs        ,&  ! max number of missing observations (12 per hour)
    rlhn_search    ,&  ! radius (gridpoints) for profiles search (if lhn_search)
    ktop_lhn       ,&  ! index for uppest model layer for which lhn is performed
    ktop_temp      ,&  ! temperature of  uppest model layer for which lhn is performed
    kbot_lhn       ,&  ! index for lowest model layer for which lhn is performed
    lhn_dt_obs     ,&  ! time step of input data in minutes
    nradar         ,&  ! max. number of radar stations within input data
    lhn_coef       ,&  ! factor for reduction of lhn t-increments
    thres_lhn      ,&  ! threshold of rain rates to be consinderd within lhn approach
    rqrsgmax       ,&  ! ratio of maximum of qrsgflux, needed for reference precipitation
    nlhn_relax     ,&  ! number of interations of horizontal filtering
    abs_lhn_lim    ,&  ! absolute limit for lhn t-increm. (imposed if lhn_limit)
    fac_lhn_search ,&  ! factor for searching nearby profiles
    fac_lhn_up     ,&  ! limiting factor for upscaling of model heating profile
    fac_lhn_down   ,&  ! limiting factor for downscaling model heating profile
    rad_wobs_lhn   ,&  ! max. distance to radar for full observation weight
    radar_in       ,&  ! directory for reading radar-files
    blacklist_file     ! filename for blacklist

USE data_lheat_nudge  , ONLY :   &
    lhn_height     ,&  ! use height infos for radar data
    lhn_bright     ,&  ! apply bright band detection
    height_file    ,&  ! dxheight_file_name
    nulhn

! end of data_lheat_nudge
!-------------------------------------------------------------------------------

USE parallel_utilities, ONLY:  distribute_values
USE environment,        ONLY:  get_free_unit

USE src_obs_use_org,    ONLY:  organize_nudging
USE src_lheat_nudge,    ONLY:  organize_lhn
USE src_sfcana,         ONLY:  organize_sfcana, init_sfcana

!===============================================================================

IMPLICIT NONE

!===============================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
!                                                           
! Parameter list:
CHARACTER (LEN= *),       INTENT(IN)            ::                      &
  yaction      ! action to be performed

INTEGER (KIND=iintegers), INTENT(OUT)           ::                      &
  ierror       ! error status 

CHARACTER (LEN=  *)                             ::                      &
  yerrmsg      ! error message

!-------------------------------------------------------------------------------
!
! Local variables: 
INTEGER (KIND=iintegers)   ::   izmemstat, izerrstat, nuin
CHARACTER (LEN= 9)         ::   yinput
CHARACTER (LEN=15)         ::   yzroutine

!-------------------------------------------------------------------------------
!- End of header
!-------------------------------------------------------------------------------
 
!-------------------------------------------------------------------------------
!- Begin Subroutine organize_assimilation
!-------------------------------------------------------------------------------

  yzroutine = 'organize_assimilation'
  izerrstat = 0_iintegers
  ierror    = 0_iintegers
  yerrmsg   = '    '

!-------------------------------------------------------------------------------
! Section 1: Input of the NAMELIST Input
!-------------------------------------------------------------------------------

  IF (yaction == 'input') THEN

    ! Open NAMELIST-INPUT file
    IF (my_world_id == 0) THEN
      PRINT *,'    INPUT OF THE NAMELISTS FOR ASSIMILATION'
      yinput   = 'INPUT_ASS'
      nuin     =  1
      OPEN(nuin   , FILE=yinput  , FORM=  'FORMATTED', STATUS='UNKNOWN',  &
           IOSTAT=izerrstat)
      IF(izerrstat /= 0) THEN
        yerrmsg  = ' ERROR    *** Error while opening file INPUT_ASS *** '
        ierror   = 2
        RETURN
      ENDIF
    ENDIF

    ! read the NAMELIST group nudging:
    CALL input_nudging (nuspecif, nuin, izerrstat)

    IF (my_world_id == 0) THEN
      ! Close file for input of the NAMELISTS
      CLOSE (nuin    , STATUS='KEEP')
    ENDIF

    IF (izerrstat < 0) THEN
      yerrmsg = 'ERROR *** while reading NAMELIST Group /NUDGING/ ***'
      ierror  = 3
    ELSEIF (izerrstat > 0) THEN
      yerrmsg = 'ERROR *** Wrong values occured in NAMELIST INPUT_ASS ***'
      ierror  = 4
    ENDIF

!-------------------------------------------------------------------------------
! Section 2: 
!-------------------------------------------------------------------------------

  ELSEIF (yaction == 'init') THEN

    ! initialize unit numbers for the ascii output files
    CALL get_free_unit (nsfc)
    CALL get_free_unit (nupr)
    CALL get_free_unit (nuqc)
    CALL get_free_unit (nuverif)
    CALL get_free_unit (nustat)
    CALL get_free_unit (nuaof)
    CALL get_free_unit (nuaofex)
    CALL get_free_unit (nurej)
    CALL get_free_unit (nuodr)
    CALL get_free_unit (nugps)
    CALL get_free_unit (nusatin)
    CALL get_free_unit (nucautn)
    CALL get_free_unit (nulhn)

    ! initialize meta data for output of near-surface analysis fields
    CALL init_sfcana

!-------------------------------------------------------------------------------
! Section 3:
!-------------------------------------------------------------------------------

  ELSEIF (yaction == 'nudge') THEN

    CALL organize_nudging
 
!-------------------------------------------------------------------------------
! Section 4:
!-------------------------------------------------------------------------------

  ELSEIF (yaction == 'lhn') THEN

    CALL organize_lhn

!-------------------------------------------------------------------------------
! Section 5:
!-------------------------------------------------------------------------------

  ELSEIF (yaction == 'surface') THEN

    CALL organize_sfcana

!-------------------------------------------------------------------------------
! Section 6: All other actions are wrong
!-------------------------------------------------------------------------------

  ELSE

    ierror  = 1
    yerrmsg = 'ERROR *** No valid action for the assimilation ***'

  ENDIF

!-------------------------------------------------------------------------------
! Internal procedures
!-------------------------------------------------------------------------------

CONTAINS

!===============================================================================
!+ Internal procedure in "organize_assimilation" for NAMELIST input
!-------------------------------------------------------------------------------

SUBROUTINE input_nudging (nuspecif, nuin, ierrstat)

!-------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST-group nudge. 
!   The group nudge contains variables for controlling the data 
!   assimilation
!     - general parameters, and parameters for upper-air data
!     - parameters for nudging data from surface stations
!     - parameters for printing output : general and related to nudging
!       vertical profiles
!     - various parameters, including some related to lateral spreading
!       of obs. increments
!     - parameters for printing output related to nudging data from 
!       surface stations
!     - parameters for controlling observation processing
!
! Method:
!   All variables are initialized with default values and then read in
!   from the file INPUT. The input values are checked for errors and for
!   consistency. If wrong input values are detected the program prints
!   an error message. The program is not stopped in this routine but an
!   error code is returned to the calling routine that aborts the
!   program after reading in all other namelists.
!   In parallel mode, the variables are distributed to all nodes with
!   the environment-routine distribute_namelists.
!   Both, default and input values are written to the file YUSPECIF
!   (specification of the run).
!
!-------------------------------------------------------------------------------

! Subroutine arguments:
! --------------------

  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Local parameters:
! ----------------

  REAL    (KIND=ireals)    , PARAMETER  :: &
    c0     = 0.0_ireals     ,& ! standard real constant 0.0
    c1     = 1.0_ireals        ! standard real constant 1.0

! Local variables:
! ---------------

! CHARACTER (LEN= *) yuspecif    ! file for protocolling the task
  CHARACTER (LEN=34) yrformat    ! format for printing real    namelist variable
  CHARACTER (LEN=34) yiformat    ! format for printing integer namelist variable
  CHARACTER (LEN=34) ylformat    ! format for printing logical namelist variable
  CHARACTER (LEN=34) yr1format   ! format for printing real    namelist variable
  CHARACTER (LEN=18) ygpformat   ! format for printing 'igpscen'
  CHARACTER (LEN=18) ywtformat   ! format for printing 'niwtyp', 'kwtyp'
  CHARACTER (LEN=27) yotformat   ! format for printing 'iwtyp'
  CHARACTER (LEN=20)       :: &
    yroutine     ,& ! name of this subroutine
    yerrmsg         ! error message

  INTEGER   (KIND=iintegers)    ::        &
!   nvarint      ,& ! number of integer   namelist parameter values
!   nvarreal     ,& ! number of real      namelist parameter values
!   nvarlog      ,& ! number of logical   namelist parameter values
!   nvarchar     ,& ! number of character namelist parameter values
    nnl      (4) ,& ! numbers of namelist parameter values of the differ. types
    nactgpc      ,& ! number of active GPS processing centres
    ierr, istat, k, n, noffset, ivar, iwtmin, iwtmax, iiwt, iiw2, nzylen,      &
    nlhnv_start_dum, nlhn_start_dum, ndelobs_lhn, nn, iz_err

  LOGICAL                    ::       &
    liwident, liwrest

  REAL (KIND=ireals)         ::       &
    dtdoh           ! dt / 3600 seconds

! Variables for default values
! ----------------------------

! For all arrays of length 4: 1: horiz. wind, 2: 'surface' pressure,
!                             3: temperature, 4: (relative) humidity

!      0.    General steering switches
!      -------------------------------

  LOGICAL                    ::       &
    lnudge_d     ,& ! .f. : on - off switch for nudging
    lverif_d     ,& ! .f. : on - off switch for verification
    lverpas_d       ! .t. : on - off switch for verif. also of passive reports

!      1.    General variables controlling the nudging
!      -----------------------------------------------

! LOGICAL                    ::       &
!   lobdens_d       ! .f. : if TRUE then compute obs. density at obs. locations
!   lseparw_d       ! .f. : if TRUE then compute net obs. increments for
!                   !               different observing systems separately

  INTEGER (KIND=iintegers)   ::       &
    nudgsta_d    ,& ! 0   : start of nudging period in timesteps
    nudgend_d    ,& ! 0   : end of nudging period in timesteps
    nversta_d    ,& ! 0   : start of verification period in timesteps
    nverend_d    ,& ! 0   : end of verification period in timesteps
    mruntyp_d    ,& ! -1  : type of current model run used for increments in VOF
    mveripr_d       ! 3   : type of verification/observation file(s) written

  INTEGER (KIND=iintegers)   ::       &
    nwtyp_d        ,& ! 1 : if > 1 then compute net obs. increments for 'nwtyp'
                      !     different sets of observing systems separately
    niwtyp_d(mxwt) ,& ! 1,0,0,..: for each of these sets observing systems:
                      !     number of observation or code types which belong
                      !     to that set of observing systems
    iwtyp_d (mxtyw),& ! 0,0,0,..: observation types (for values > 0) or code 
                      !     types (for values < 0) belonging to a set of obs.
                      !     systems, specified successively for each of these
                      !     sets in array 'iwtyp'
    kwtyp_d (mxwt+2)  ! 1 : mode of weights W for multiple observations,
                      !     specified for each set of obseration systems:
                      !     1 : W=  w**2 /sum(w) ;  w= weight of 1 single obs
                      !     2 : W= (w**2 +w) /(1+sum(w))

  REAL (KIND=ireals)         ::       &
    hnudgsta_d   ,& ! 0   : start of nudging period in 'model integration hours'
    hnudgsta     ,& ! 0   : start of nudging period in 'model integration hours'
    hnudgend_d   ,& ! 0   : end of nudging period in 'model integration hours'
    hnudgend     ,& ! 0   : end of nudging period in 'model integration hours'
    hnudgend_dum ,& ! 0   : end of nudging period in 'model integration hours'
    hversta_d    ,& ! 0   : start of verification period in 'model integ. hours'
    hverend_d    ,& ! 0   : end of verification period in 'model integr. hours'
    hverend_dum  ,& ! -3. : end of verification period in 'model integr. hours'
    tconbox_d       ! 6*dt: timestep [s] for computing analysis increments
                    !       (i.e. time box of constant analysis increments)

!      2.    Corrections to balance the analysis increments
!      ----------------------------------------------------

  LOGICAL                    ::       &
    luvgcor_d       ! .t. : .t. ==> geostrophic wind correction applied

  INTEGER (KIND=iintegers)   ::       &
    mpsgcor_d    ,& ! 1   : mode to apply geostrophic pressure correction
                    !       (balances near-surface wind analysis incements
                    !        obtained by isotropic lateral spreading):
                    !       = 0 : no pressure correction
                    !       = 1 : corr. balacing wind from scatterometer only
                    !       = 2 : corr. balacing scatt. + in-situ 10-m wind obs
    ntpscor_d    ,& ! 1   : switch for temperature correction when the 'surface'
                    !       pressure (i.e. p. at lowest model level) is nudged,
                    !       so that the final geopotential change due to surface
                    !       pressure nudging is zero above the level 'ptpstop'
                    !       0 : no temperature correction
                    !       1 : cfi(k) = ck(k)^2 * EXP( (1.-ck(k)^3) /4.)
                    !       2 : cfi(k) = .5 *ck(k) *(1.+ck(k))
                    !           where  ck(k) =  (p(k) - ptpstop)
                    !                         / (ps - ptpstop)   for p > ptpstop
                    !                  ck(k) = 0.                for p < ptpstop
                    !           and cfi is the correlation of the geopotential
                    !                   change at level k to that at the surface
    khumbal_d       ! 100 : radius (in grid pts. units) of the area around a
                    !       convectively precipitating grid point, at which
                    !       specified instead of relative humidity is preserved
                    !       when nudging temperature.
                    !       If 'khumbal' <= -1 or >= 99, then relative resp. 
                    !       specific humidity is preserved everywhere.
                    !       At the temperature correction (ps-nudging), relative
                    !       humidity is preserved except if 'khumbal' >= 100

  REAL (KIND=ireals)         ::       &
    ptpstop_d    ,& ! 400.: pressure [hPa] at upper boundary of the temperature
                    !       correction for 'surface' pressure nudging
    qgeo_d       ,& ! .3  : factor to the geostrophic wind increments at 1000hPa
    qgeotop_d    ,& ! .5  : factor to the geostrophic wind increments at ptpstop
    qgeops_d     ,& ! 1.  : factor to the geostrophic pressure increments

!      3.    Nudging coefficients
!      --------------------------

    gnudg_d  (4) ,& ! 6,12,6,6*10^-4: nudging coefficients for TEMP / PILOT data
    gnudgar_d(4) ,& ! 6, 0,6,0*10^-4: nudging coeffic. for AIRCRAFT data [1/s]
    gnudgsu_d(4) ,& ! 6,12,0,6*10^-4: nudging coef. for surface-level data [1/s]
    gnudggp_d       !        0*10^-4: nudging coef. for GPS-derived IWV [1/s]

!      4.    Temporal weights
!      ----------------------

  LOGICAL                    ::       &
    ltipol_d     ,& ! .t. : .t. ==> linear interpolation in time of upper-air
                    !               data which are less than 'tipolmx' hrs apart
                    !           ==> at most 2 reports of same type per station
                    !               used at each timestep
    ltipsu_d        ! .t. : .t. ==> linear interpolation in time of surface-lev.
                    !        data which are less than 'tipmxsu' hrs apart

  REAL (KIND=ireals)         ::       &
    tipolmx_d    ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       for upper-air data (set tipolmx = 0, if .NOT ltipol)
    tipmxsu_d    ,& ! 1.0 : max. time span (hrs) allowed for linear interpolat.
                    !       of surface + GPS data  (tipmxsu = 0, if .NOT ltipsu)
    wtukrsa_d    ,& ! 3.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukrse_d    ,& ! 1.0 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for TEMP / PILOT
    wtukara_d    ,& ! 1.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       relative to the obs. time .   for AIRCRAFT data
    wtukare_d    ,& ! 0.5 : for saw-tooth-shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       relative to the obs. time .   for AIRCRAFT data
    wtuksua_d    ,& ! 1.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the past
                    !       rel. to the obs. time   for surface-level + GPS data
    wtuksue_d       ! 0.5 : for saw-tooth shaped temporal weights:
                    !       temporal radius of influence towards the future
                    !       rel. to the obs. time   for surface-level + GPS data

!      5.    Spatial weights : Spreading of observational information
!      --------------------------------------------------------------

!      5.1   Mode of spreading
!      -----------------------

  INTEGER (KIND=iintegers)   ::       &
    msprpar_d    ,& ! 1   : switch specifying the surfaces along which obs.
                    !       increments of upper-air data are (primarily) spread
                    !       0: spreading along model levels, vertical correl-
                    !          ations depend approx. on 'ln(p)' differences
                    !       1: spreading along horizontal surfaces, vertical
                    !          correlations depend approx. on 'ln(p)' differenc.
                    !       2: spreading along isentropic surfaces, vertical
                    !          correlations depend potential temperature differ.
    msprpsu_d       ! 0   : parameter specifying the surface along which
                    !       surface-level data increments are primarily
                    !       spreaded, i.e. the parameter which scales the
                    !       vertical distance between the grid pt. assigned
                    !       to the obs. and the target grid pt.
                    !       1 : vertical weights dependent on (scaled) height
                    !           (approx. log( pressure ) ) differences
                    !          ==> spreading primarily horizontal
                    !       2 : vertical weights dependent on potential
                    !           temperature differences
                    !          ==> spreading primarily along isentropes

!      5.2   Vertical weights
!      ----------------------

  REAL (KIND=ireals)         ::       &
!      for upper-air data
    vcorls_d (4) ,& ! 2*.333,: square of the vertical correlation scale,
                    ! 2*.04    i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpar <= 1), or
                    !          in potential temperature if (msprpar == 2)
                    !          (reasonable values in latter case: 2*275., 2*33.)
    vcutof_d (4) ,& ! 2*.75, : cut-off of the vertical correlation.       Units:
                    ! 2*1.     value of correlation at cut-off is [exp(-vcutof)]
                    !          (c.f UKMO: cut-off at: |ln(p/pobs))| = 0.5
                    !                     correl. fn:  exp(-3 *(ln(p/pobs))**2 )
                    !               ==> value of correl. at cut-off: exp(-.75) )
                    !          (i.e. for (msprpar <= 1) and default values for
                    !                'vcorls', 'vcutof', the cut-off 'radius'
                    !                for wind is 0.5 (density) scale heights.)
    wablua_d (4) ,& ! 4*  1. : factor to weights within the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)
!      for surface-level data
    vcorlsu_d(4) ,& ! 2*.013,: square of the vertical correlation scale,
                    ! 2*.002   i.e. of the Gaussian vertical influence 'radius'
                    !          in log( pressure ) if (msprpsu == 1) or
                    !          in potential temperature if (msprpsu == 2)
                    !          (e-folding decay height of ca. 300m for T, RH, if
                    !          default, or 2*11.1, 2*1.33 for (msprpsu == 2))
    vcutosu_d(4) ,& ! 2* .75 : cut-off of the vertical correlation.       Units:
                    ! 2*4.     value of correlation at cut-off : [exp(-vcutosu)]
    vpblsu_d (4) ,& ! 2*99. ,: Gaussian vertical influence 'radius' of potential
                    ! 2*99.    temperature differences between obs. level and
                    !          model level, to render vertical weights depend on
                    !          the stability (in the PBL) even if (msprpsu <= 1)
    wablsu_d (4)    ! 4*  1. : factor to weights above the ABL
                    !          (atmospheric boundary layer, i.e. mixed layer)

  LOGICAL                    ::       &
    lsvcorl_d       ! .t. : .t. ==> adjustment of vertical correlation scales
                    !               in the presence of close observations

!      5.3   Horizontal weights
!      ------------------------

  REAL (KIND=ireals)         ::       &
!      for upper-air data
    rhinfl_d (4) ,& ! 0.,70.,: constant part of the 'correlation scale of the
                    ! 0., 0.   autoregressive horiz. correlation function'
                    !          (='COSAC') [km]
    rhvfac_d (4) ,& ! 1., 0. : multiplication factor to the vertically varying
                    ! 2* .83   part of the 'COSAC' (as def. in 'data_nudge_all')
    rhtfac_d (4) ,& ! 1.3,   : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.3,     relative to the 'COSAC' at the obs. time as given
                    ! 1.3      by 'rhinfl', 'rhvfac')
    rhfgps_d     ,& ! 0.45   : scaling (reduction) factor of the total 'COSAC'
                    !          for humidity derived from GPS IWV
    rhfpsdd_d    ,& ! 1.0    : minimum scaling (reduction) factor of the total
                    !          'COSAC' for surface pressure dep. on data density
    cutofr_d (4) ,& ! 4* 3.5 : cut-off in 'COSAC' units of the horizontal
                    !          correlation function
    vcsni_d  (4) ,& ! 4*2500.: square of Gaussian vertical influence 'radius'
                    !          in potential temperature (if msprpar <= 1) or
                    !          log( pressure ) (if msprpar == 2) on surfaces
                    !          along which the obs. increments are spread
                    !          laterally. This determines the non-isotropy of
                    !          the (1-d) horizontal correlation functions.
!      for surface-level data
    rhiflsu_d(4) ,& ! 2 *70.,: constant part of the 'correlation scale of the
                    !   100.,  autoregressive horiz. correlation function'
                    !    70.   (='COSAC') [km]  (at the obs. time)
    rhtfsu_d (4) ,& ! 1.,    : scaling factor of the total 'COSAC' for the
                    ! 1.43,    beginning and end of the nudging period for 1 obs
                    ! 1.,      relative to the 'COSAC' at the obs. time as given
                    ! 1.       by 'rhiflsu')
    cutofsu_d(4) ,& ! 2., 3.5: cut-off in 'COSAC' units of the horizontal
                    ! 2., 2.   correlation function
    vcsnisu_d(4) ,& ! 2*2500.: square of Gaussian vertical influence 'radius'
                    ! 2*   9.  in potential temperature (if msprpsu <= 1) or
                    !          log( pressure ) (if msprpsu == 2) on surfaces
                    !          along which the obs. increments are spread
                    !          laterally. This determines the non-isotropy of
                    !          the (1-d) horizontal correlation functions.
!      degree of non-divergence of 2-dim. horizontal wind correlation functions
    cnondiv_d    ,& ! .1  : constant part of the factor to the non-divergence
                    !       correction ('fanodicor') in 2-dim wind correlations
    fnondiv_d    ,& ! .8  : multiplication factor to the vertically varying part
                    !       of the factor to the non-divergence correction
                    !       (as defined in 'data_nudge_all')
    tnondiv_d       ! 1.1 : temporal multiplication factor to the factor to the
                    !       non-divergence correction for the beginning and end
                    !       of the nudging period for 1 obs relative to that
                    !       given by 'cnondiv', 'fnondiv' for the obs. time

!      6.    Computation of observation increments
!      -------------------------------------------------

!      6.1.  Vertical profiles of observation increments
!      -------------------------------------------------

  LOGICAL                    ::       &
    lscadj_d (4)    ! .T.,   : .F. ==> linear vertical interpolation (in log(p))
                    ! .T.,     instead of vertical scale adjustment (by vertical
                    ! .T.,     averaging over the model layer) for conveying the
                    ! .F.      observational information to the model levels
                    !          (for comput. obs. incr. at model levels)

  REAL (KIND=ireals)         ::       &
    topobs_d (4) ,& !  849., : threshold [hPa]: above this level (p < topobs),
                    ! 1099.,   only obs. increments at model levels are used,
                    !  799.,   i.e. obs. incr. at obs. levels are not used
                    !  699.    ('topobs' is fixed at 1099. if (msprpar == 0))
    botmod_d (4)    ! 1099., : threshold [hPa]: below this level (p > botmod),
                    ! 1099.,   only obs. increments at obs. levels are used,
                    ! 1099.,   i.e. obs. incr. at model levels are not computed
                    !  899.    ((botmod >= topobs), and
                    !           'botmod' is fixed at 1099. if (msprpar == 0))

!      6.2   Surface-level observation increments
!      ------------------------------------------

  LOGICAL                    ::       &
    loiqv2m_d    ,& ! .f. : .t. ==> 2-m humidity observation increments as
                    !               differences of specific humidity instead of
                    !               relative humidity
    lqfqv2m_d       ! .f. : .t. ==> quality factor for 2-m humidity observations
                    !               dependent on T-2m differences

!      7.    Quality control and quality weights
!      -----------------------------------------

!      7.1   (Threshold) Quality control
!      ---------------------------------

  REAL (KIND=ireals)         ::       &
    dtqc_d       ,& ! 720.: timestep (in [s]) for the threshold quality control
                    !       (the quality ctrl is always applied to observations
                    !        at the first time when they are used)
!      for upper-air data
    qcc_d    (4) ,& !  0.,500: constant parts of the quality control thresholds
                    !  0., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
    qcvf_d   (4) ,& !  5., 1.: multiplication factor to the vertically varying
                    ! 10., 0.  part of the QCT (as def. in 'data_nudge_local')
!      for surface-level data
    qccsu_d  (4) ,& ! 12.,500: constant parts of the quality control thresholds
                    ! 12., .7  (='QCT'). (u,v):[m/s], ps: [Pa], T: [k], RH: [ ]
!      for integrated water vapour (IWV) derived from GPS or radiosonde data
    qcciq_d      ,& ! 1.  : constant part of QC threshold for IWV
    qcsiq_d      ,& ! .15 : IWV QC threshold, as a fraction of IWV of saturated
                    !       model temperature profile
    qcflbcp_d    ,& ! 1.4 : enhancement factor to the threshold used for the
                    !       check of ps against lateral boundary fields
                    !       (if (qcflbcp <= 0) then no LBC check for ps)

!      7.2   Quality weights
!      ---------------------

    qcfpst_d        ! 1.5 : maximum enhancement of the weight for surface
                    !       pressure observations due to the pressure tendency

!      8.    Observation processing
!      ----------------------------

!      8.0   Reading of observation reports
!      ------------------------------------

  INTEGER (KIND=iintegers)   ::       &
    irun_osse_d  ,& !   0 : model run from which are simulated obs from file
                    !       yfofin='fof' are used to derive the obs values
                    !       (if = 0 then obs values = obs values from 'fof')
    itype_obfile_d  !   1 : type of observation input file(s): 1: AOF, 2: NetCDF

  LOGICAL                    ::       &
    losse_fg_d      ! .f. : if true then first guess check flag from 'fof' is
                    !               converted into 'dataset' pre-processing flag
                    !       if false then first guess check flag is discarded
                    !               and obs may be used actively

  CHARACTER (LEN=icdfdirlen) :: ycdfdir_d   ! directory with NetCDF input files
  CHARACTER (LEN=70)         :: yaofpath_d  ! full path name of AOF file

!      8.1   Use of stations / reports
!      -------------------------------

  REAL (KIND=ireals)         ::       &
    obnlat_d     ,& !   90.  : northern boundary of observation area
    obslat_d     ,& !  -90.  : southern boundary of observation area
    obwlon_d     ,& ! -180.  : western boundary of observation area
    obelon_d     ,& !  180.  : eastern boundary of observation area
    exnlat_d     ,& !   90.  : northern boundary for exclusion area
    exslat_d     ,& !  -90.  : southern boundary for exclusion area
    exwlon_d     ,& ! -180.  : western boundary for exclusion area
    exelon_d     ,& !  180.  : eastern boundary for exclusion area
    doromx_d (4) ,& !  100., : vertical extrapolation cut-off and gaussian
                    !  150.,   radius of height differences between model
                    !  150.,   orography and surface station height for a factor
                    !  150.,   contributing to the quality weight factor as part
                    !          of the nudging weights
                    !          (height diff of SYNOP/ GPS with (z-obs > z-model,
                    !           --> interpolation instead of extrapolation) are
                    !           divided by 4 (fdoro) for surf pressure/ IWV obs)
    altopsu_d(4) ,& !   100.,: SYNOP obs. above height 'altopsu' are not assimi-
                    ! 3*5000.  lated. If (altopsu == 0.) then SYNOP / surf. TEMP
                    !          assigned to land grid pts. are not assimilated
    thairh_d     ,& !   20.  : maximum horizontal distance [km] between the
                    !          lowest report and any single level report that
                    !          is added to a multi-level AIRCRAFT report
    fperturb_d      !    0.  : factor to the obs error variances to define the
                    !          size of random perturbations added to the obs
                    !          (only for data from feedback file yfofin='fof')

  LOGICAL                    ::       &
    lgpsbias_d      !  .f.   : .t. ==> bias correction to GPS IWV applied

  INTEGER (KIND=iintegers)   ::       &
!   size def. of the 'ODR' (obs. data record) for internal storage of all obs.
!   within the max. time window given by 'wtuk??a', wtuk??e', 'tip??mx'
    maxmlo_d     ,& !  300   : max. number of multi-level reports in the ODR
    maxsgo_d     ,& ! 3000   : max. number of (surface-level and upper-air)
                    !                         single-level reports in the ODR
    maxuso_d     ,& !  900   : max. number of upper-air single-level rep. in ODR
    maxgpo_d     ,& ! 3000   : max. number of GPS reports in ODR on total domain
    maxmlv_d     ,& !  100   : max. number of obs levels in multi-level reports
    mxfrep_d     ,& !   -1   : max. number of reports in NetCDF feedobs file
    mxfobs_d     ,& !   -1   : max. number of observations in feedobs file
    nolbc_d      ,& !    5   : number of grid rows at lateral boundaries
                    !          where obs are neglected
    mqcorr92_d   ,& !    0   : switch for bias correction for Vaisala RS92
                    !            radiosonde humidity
    iseed_d      ,& !    0   : external seed for random number generator
    mrhyes_d(mxda),& ! mxda*0 : station ID's of SURFACE humidity data that are
                     !          always used (irrespective of 'doromx', 'altopsu)
    mrhno_d (mxda),& ! mxda*0 : station ID's of SURFACE humidity data never used
    muvyes_d(mxda),& ! mxda*0 : station ID's of SURFACE wind data always used
    muvno_d (mxda),& ! mxda*0 : station ID's of SURFACE wind data never used
    igpscen_d(mxgpc) ! X* -1  : array of used GPS processing centres 

!      8.2   Use of observation and code types
!      ---------------------------------------

  LOGICAL                    ::       &
    lsynop_d     ,& ! .t.    : .t. if SYNOP data is used
    laircf_d     ,& ! .t.    : .t. if AIREP data is used (aircraft)
    lsatob_d     ,& ! .false.: .t. if SATOB data is used
    ldribu_d     ,& ! .t.    : .t. if DRIBU data is used (drifting buoy)
    ltemp_d      ,& ! .t.    : .t. if TEMP  data is used
    lpilot_d     ,& ! .t.    : .t. if PILOT data is used
    lsatem_d     ,& ! .false.: .t. if SATEM data is used
    lgps_d       ,& ! .false.: .t. if GPS   data is used
    lscatt_d        ! .t.    : .t. if SCATT data is used (scatterometer)

  LOGICAL                    ::       &
    lcd011_d     ,& ! .t.    : .t. if synop code  11 data is used (land synop)
    lcd014_d     ,& ! .t.    : .t. if synop code  14 data is used (automatic)
    lcd021_d     ,& ! .t.    : .t. if synop code  21 data is used (ship)
    lcd022_d     ,& ! .t.    : .t. if synop code  22 data is used (ship abbrev.)
    lcd023_d     ,& ! .t.    : .t. if synop code  23 data is used (shred)
    lcd024_d     ,& ! .t.    : .t. if synop code  24 data is used (autom. ship)
    lcd140_d     ,& ! .t.    : .t. if synop code 140 data is used (metar)
    lcd041_d     ,& ! .t.    : .t. if airep code  41 data is used (codar)
    lcd141_d     ,& ! .t.    : .t. if airep code 141 data is used (airep)
    lcd241_d     ,& ! .t.    : .t. if airep code 241 data is used (colba)
    lcd144_d     ,& ! .t.    : .t. if airep code 144 data is used (amdar)
    lcd244_d     ,& ! .t.    : .t. if airep code 244 data is used (acars)
    lcd088_d     ,& ! .t.    : .t. if satob code  88 data is used (satob)
    lcd188_d     ,& ! .false.: .t. if satob code 188 data is used (sst)
    lcd063_d     ,& ! .t.    : .t. if dribu code  63 data is used (bathy)
    lcd064_d     ,& ! .t.    : .t. if dribu code  64 data is used (tesac)
    lcd165_d     ,& ! .t.    : .t. if dribu code 165 data is used (drift. buoy)
    lcd035_d     ,& ! .t.    : .t. if temp  code  35 data is used (land temp)
    lcd036_d     ,& ! .t.    : .t. if temp  code  36 data is used (temp ship)
    lcd037_d     ,& ! .t.    : .t. if temp  code  37 data is used (mobile)
    lcd135_d     ,& ! .t.    : .t. if temp  code 135 data is used (dropsonde)
    lcd039_d     ,& ! .t.    : .t. if temp  code  39 data is used (rocob)
    lcd040_d     ,& ! .t.    : .t. if temp  code  40 data is used (rocob ship)
    lcd032_d     ,& ! .t.    : .t. if pilot code  32 data is used (land pilot)
    lcd033_d     ,& ! .t.    : .t. if pilot code  33 data is used (pilot ship)
    lcd038_d     ,& ! .t.    : .t. if pilot code  38 data is used (mobile)
    lcd132_d     ,& ! .t.    : .t. if pilot code 132 data is used (win-prof eu)
    lcd133_d     ,& ! .t.    : .t. if pilot code 133 data is used (sod/rass eu)
    lcd136_d     ,& ! .t.    : .t. if pilot code 136 data is used (pro/rass us)
    lcd137_d     ,& ! .t.    : .t. if pilot code 137 data is used (Radar VAD)
    lcd086_d     ,& ! .t.    : .t. if satem code  86 data is used (satem)
    lcd186_d     ,& ! .t.    : .t. if atovs code 186 data is used (hi-res ATOVS)
    lcd122_d     ,& ! .t.    : .t. if scatt code 122 data is used (QuickScat)
    lcd123_d     ,& ! .t.    : .t. if scatt code 123 data is used (ASCAT)
    lcd096_d        ! .t.    : .t. if gps data from COST ASCII file is used

!      9.    2-D analyses
!      ------------------

  LOGICAL                    ::       &
    lsurfa_d     ,& ! .f.    : .t. if surface fields are analysed
    lt2m_d       ,& ! .f.    : .t. if 2m temperat. field is analysed
    lrh2m_d      ,& ! .f.    : .t. if 2m rel. hum. field is analysed
    lprecp_d     ,& ! .f.    : .t. if precipitation is analysed
    lff10m_d        ! .f.    : .t. if 10m wind speed is analysed

  REAL (KIND=ireals)         ::       &
    ht2a_d       ,& ! 999.   : time of 1. T2m-ana in hours since model start
    ht2i_d       ,& ! 999.   : time increment to next T2m analysis
    hh2a_d       ,& ! 999.   : time of 1. RH2m-ana in hours since model start
    hh2i_d       ,& ! 999.   : time increment to next RH2m analysis
    hffa_d       ,& ! 999.   : time of 1. 10m wind-ana in hours since model start
    hffi_d       ,& ! 999.   : time increment to next wind analysis
    hprc_d       ,& ! 999.   : time of prec-ana in hours since model start
    raintp_d        ! 12.    : time period of precipitation analysis

  CHARACTER (LEN=250)        ::       &
    ydir_lansfc_d   ! './'   : directory where to write the 2-D analyses

  CHARACTER (LEN=  4)        ::       &
    yform_lansfc_d  ! 'grb1' : format for the 2-D analyses files

!      10.   Diagnostic output
!      -----------------------

  LOGICAL                    ::       &
    lpraof_d     ,& ! .f.    : .t. for diagnostic print of analysis obs file AOF
    lprodr_d     ,& ! .t.    : .t. for diagnostic print of obs data records ODR
    ldiasa_d        ! .f.    : .t. for diagnostics of surface analysis

  INTEGER (KIND=iintegers)   ::       &
    ionl_d       ,& ! 167    : / grid point coordinates
    jonl_d       ,& ! 103    : \ for standard output on nudging
    ionl2_d      ,& ! 167    : / 2nd grid pt coordinates
    jonl2_d      ,& ! 103    : \ for other standard output on nudging
    noctrq_d        !   9    : observation/code type of diagnostic print
                    !          (if (noctrq == 9) then print all obs./code types)

  REAL (KIND=ireals)         ::       &
    dinlat_d     ,& !   0.   : northern boundary of diagnostic print area
    dislat_d     ,& !   0.   : southern boundary of diagnostic print area
    diwlon_d     ,& !   0.   : western boundary of diagnostic print area
    dielon_d        !   0.   : eastern boundary of diagnostic print area

!      11.   Latent heat nudging
!      -------------------------

  LOGICAL                          ::           &
    llhn_d           ,& ! on/off switch for latent heat nudging (lhn)
    llhnverif_d      ,& ! on/off switch for verification against radar
    lhn_search_d     ,& ! search for appropriate nearby model heating profile
    lhn_filt_d       ,& ! vertical filtering of lhn t-increments
    lhn_relax_d      ,& ! horizontal filtering of lhn t-increments
    lhn_limit_d      ,& ! impose absolute limit on lhn t-increm. (abs_lhn_lim)
    lhn_hum_adj_d    ,& ! apply a humidity adjustment along with t-increments
    lhn_incloud_d    ,& ! apply the LHN-scaling in cloudy layers only
    lhn_spqual_d     ,& ! switch for use of a spatial quality function
    lhn_black_d      ,& ! use blacklist for radar data
    lhn_qrs_d        ,& ! calculate the integrated precipitation flux
    lhn_logscale_d   ,& ! apply logarithmic scaling factors
    lhn_wweight_d    ,& ! apply a weighting with respect to the mean horizontal wind
    lhn_height_d     ,& ! use height infos for radar data
    lhn_bright_d     ,& ! apply bright band detection
    lhn_diag_d          ! produce more detailed diagnostic output during lhn

  INTEGER (KIND=iintegers) ::  &
    nlhn_start_d     ,& ! start of latent heat nudging period in timesteps
    nlhn_end_d       ,& ! end of latent heat nudging period in timesteps
    nlhnverif_start_d,& ! start of latent heat nudging period in timesteps
    nlhnverif_end_d  ,& ! end of latent heat nudging period in timesteps
    rlhn_search_d    ,& ! radius (grid pts) for profiles search (if lhn_search)
    nlhn_relax_d     ,& ! number of interations of horizontal filtering
    ktop_lhn_d       ,& ! index for uppest model layer at which lhn is performed
    kbot_lhn_d       ,& ! index for lowest model layer at which lhn is performed
    nradar_d            ! max. number of radar stations within input data

  REAL (KIND=ireals)         ::       &
    hlhn_start            ,&
    hlhn_start_d          ,&
    hlhn_end              ,&
    hlhn_end_d            ,&
    hlhn_end_dum          ,&
    hlhnverif_start       ,&
    hlhnverif_start_d     ,&
    hlhnverif_end         ,&
    hlhnverif_end_d       ,&
    hlhnverif_end_dum

  REAL (KIND=ireals)               ::           &
    lhn_coef_d        ,& ! factor for reduction of lhn t-increments
    lhn_dt_obs_d      ,& ! time step of input data in minutes
    abs_lhn_lim_d     ,& ! absolute limit for lhn t-increments (used if lhn_limit)
    fac_lhn_search_in ,& ! factor when searching will applied incoming
    fac_lhn_search_d  ,& ! factor when searching will applied
    fac_lhn_up_in     ,& ! limiting factor for upscaling of model heating profile incoming
    fac_lhn_up_d      ,& ! limiting factor for upscaling of model heating profile
    fac_lhn_down_in   ,& ! limiting factor for downscaling model heating profile incoming
    fac_lhn_down_d    ,& ! limiting factor for downscaling model heating profile
    thres_lhn_d       ,& ! threshold of rain rates to be consinderd within lhn approach
    rad_wobs_lhn_d    ,& ! max. distance to radar for full observation weight
    rqrsgmax_d        ,& ! ratio of maximum of qrsgflux, needed for reference precipitation
    ktop_temp_d          ! temperature of uppest model layer for which lhn is performed

 CHARACTER (LEN=100)              ::           &
    radar_in_d           ,& ! directory for reading radar-files
    blacklist_file_d     ,& ! filename of blacklist for radar data
    height_file_d           ! dxheight_file_name

 CHARACTER (LEN=12)               ::           &
    noobs_date_d

!      12.   1DVAR: these variables will become namelist variables
!                   only when 1DVAR is included in the official version
!      ----------------------------------------------------------------

  LOGICAL                    ::       &
    l1dvar_d        ! .f. : on - off switch for 1DVar

  REAL    (KIND=ireals)      ::       &
    gnudgtv_d(4) ,& ! 0, 0,6,6*10^-4: nudging coeffic. for sat retrivals [1/s]
    rhfrtv_d (4)    ! 1., 1.,: scaling (reduction) factor of the total 'COSAC'
                    ! 1., 0.5  for temperature / humidity satellite retrievals

  INTEGER (KIND=iintegers)   ::       &
                    !max. # reports currently active (i.e within the time window
                    !given by 'wtuk??a','wtuk??e','tip*mx*') in the total domain
    maxtvo_d        !    1   : max. number of sat retrievals within total domain

  INTEGER (KIND=iintegers)   ::       &
    ! Note: Satellite data cannot yet be produced with current official
    !       version, i.e. the following namelist variables have no effect.
    !       An experimental version for use of satellite data, based on
    !       V4_18, is available from christoph.schraff_at_dwd.de .
                    ! mcdxxxx = 0 : xxxx is not processed
                    ! mcdxxxx = 1 : xxxx is processed, but not used actively
                    ! mcdxxxx = 2 : xxxx is used actively (and processed)
                    ! --> if (mcdxxxx >= 1) then for satellite xxxx, 3 input
                    !     files are stringently expected (for bias correction,
                    !     channel selection, and RTTOV coefficients)
    mcdmsg1_d    ,& !  0     : processing / use of MSG1   code  71 data
    mcdmsg2_d    ,& !  0     : processing / use of MSG2   code  72 data
    mcdno15_d    ,& !  0     : processing / use of NOAA15 code 206 data
    mcdno16_d    ,& !  0     : processing / use of NOAA16 code 207 data
    mcdno17_d    ,& !  0     : processing / use of NOAA17 code 208 data
    mcdno18_d       !  0     : processing / use of NOAA18 code 209 data

!HJP Begin 2016-05-11
  CHARACTER(LEN=250)         :: iomsg_str
!HJP END   2016-05-11

! Define the namelist group
! -------------------------

  NAMELIST /nudging/  lnudge   ,lverif   ,lverpas                     ,        &
                      nudgsta  ,nudgend  ,tconbox  ,hnudgsta ,hnudgend,        &
                      nversta  ,nverend  ,mruntyp  ,hversta  ,hverend ,        &
                      mveripr  ,nwtyp    ,niwtyp   ,iwtyp    ,kwtyp   ,        &
                      ntpscor  ,ptpstop  ,khumbal  ,luvgcor  ,qgeo    ,        &
                      gnudggp            ,mpsgcor  ,qgeops   ,qgeotop ,        &
                      gnudg    ,ltipol   ,tipolmx  ,wtukrsa  ,wtukrse ,        &
                      gnudgar                      ,wtukara  ,wtukare ,        &
                      gnudgsu  ,ltipsu   ,tipmxsu  ,wtuksua  ,wtuksue ,        &
                      msprpar  ,vcorls   ,vcutof   ,wablua   ,lsvcorl ,        &
                      msprpsu  ,vcorlsu  ,vcutosu  ,wablsu   ,vpblsu  ,        &
                      rhinfl   ,rhvfac   ,rhtfac   ,cutofr   ,vcsni   ,        &
                      rhiflsu  ,rhfgps   ,rhtfsu   ,cutofsu  ,vcsnisu ,        &
                      rhfpsdd  ,cnondiv  ,fnondiv  ,tnondiv           ,        &
                      lscadj   ,topobs   ,botmod   ,loiqv2m  ,lqfqv2m ,        &
                      dtqc     ,qcc      ,qcvf     ,qccsu    ,qcfpst  ,        &
                      qcciq    ,qcsiq    ,qcflbcp

  NAMELIST /nudging/  obnlat   ,obslat   ,obwlon   ,obelon   ,ycdfdir ,        &
                      exnlat   ,exslat   ,exwlon   ,exelon   ,yaofpath,        &
                      doromx   ,altopsu  ,thairh   ,lgpsbias ,itype_obfile,    &
                      mqcorr92 ,mxfrep   ,mxfobs                      ,        &
                      irun_osse,losse_fg ,fperturb ,iseed             ,        &
                      maxmlo   ,maxsgo   ,maxuso   ,maxgpo            ,maxmlv ,&
                      nolbc    ,mrhyes   ,mrhno    ,muvyes   ,muvno   ,igpscen,&
                      lsynop   ,laircf   ,lsatob   ,ldribu   ,ltemp   ,        &
                      lpilot   ,lsatem   ,lgps     ,lscatt            ,        &
                      lcd011   ,lcd014   ,lcd021   ,lcd022   ,lcd023  ,lcd024 ,&
                      lcd140   ,lcd041   ,lcd141   ,lcd241   ,lcd144  ,        &
                      lcd088   ,lcd188   ,lcd063   ,lcd064   ,lcd165  ,        &
                      lcd035   ,lcd036   ,lcd037   ,lcd135   ,lcd039  ,lcd040 ,&
                      lcd032   ,lcd033   ,lcd038   ,lcd132   ,lcd133  ,lcd136 ,&
                      lcd137   ,lcd086   ,lcd186   ,lcd244   ,lcd096  ,        &
                      lcd122   ,lcd123                                ,        &
                      lsurfa   ,lt2m     ,lrh2m    ,lprecp   ,raintp  ,        &
                      ht2a     ,ht2i     ,hh2a     ,hh2i     ,hprc    ,        &
                      hffa     ,hffi     ,lff10m   ,ydir_lansfc, yform_lansfc, &
                      lpraof   ,lprodr   ,ldiasa   ,noctrq            ,        &
                      ionl     ,jonl     ,ionl2    ,jonl2             ,        &
                      dinlat   ,dislat   ,diwlon   ,dielon

! NAMELIST /nudging/  l1dvar   ,maxtvo   ,gnudgtv  ,rhfrtv            ,        &
!                     mcdmsg1  ,mcdmsg2  ,mcdno15  ,mcdno16  ,mcdno17 ,        &
!                     mcdno18

  NAMELIST /nudging/  llhn         ,llhnverif                  ,           &
                      nlhn_start   ,nlhn_end                   ,           &
                      hlhn_start   ,hlhn_end                   ,           &
                      nlhnverif_start ,nlhnverif_end           ,           &
                      hlhnverif_start ,hlhnverif_end           ,           &
                      lhn_coef, fac_lhn_up  ,fac_lhn_down      ,           &
                      thres_lhn    ,noobs_date                 ,           &
                      ktop_temp    ,rqrsgmax                   ,           &
                      ktop_lhn     ,kbot_lhn                   ,           &
                      radar_in     ,rad_wobs_lhn               ,           &
                      lhn_black    ,blacklist_file             ,           &
                      lhn_search   ,rlhn_search, fac_lhn_search,           &
                      lhn_filt     ,lhn_hum_adj                ,           &
                      lhn_limit    ,abs_lhn_lim                ,           &
                      lhn_relax    ,nlhn_relax                 ,           &
                      lhn_incloud  ,lhn_diag, lhn_qrs          ,           &
                      lhn_logscale ,lhn_wweight                ,           &
                      lhn_bright   ,lhn_height                 ,           &
                      height_file  ,lhn_spqual                 ,           &
                      lhn_dt_obs   ,nradar

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_nudging
!-------------------------------------------------------------------------------

  yroutine    = 'input_nudging'
  ierrstat    = 0_iintegers
  iz_err      = 0_iintegers

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 1: Initialize the default variables
!-------------------------------------------------------------------------------
!          1a: NUDGING
!-------------------------------------------------------------------------------

  lnudge_d    = .FALSE.
  lverif_d    = .FALSE.
  lverpas_d   = .TRUE.
  l1dvar_d    = .FALSE.
! lobdens_d   = .FALSE.
  nudgsta_d   = 0
  nudgend_d   = 0
  nversta_d   = 0
  nverend_d   = 0
  mruntyp_d   = -1
  mveripr_d   = 3
  nwtyp_d     = 1
  niwtyp_d    = 0
  niwtyp_d(1) = 1
  iwtyp_d     = 0
  kwtyp_d     = 1
  hnudgsta_d  = 0.0_ireals
  hnudgend_d  = 0.0_ireals
  hnudgend_dum=-3.0_ireals
  hversta_d   = 0.0_ireals
  hverend_d   = 0.0_ireals
  hverend_dum =-3.0_ireals
  tconbox_d   = 6 * dt
  ntpscor_d   = 1
  khumbal_d   = 100
  ptpstop_d   = 400.0_ireals
  luvgcor_d   = .TRUE.
  mpsgcor_d   = 1
  qgeo_d      = 0.3_ireals
  qgeotop_d   = 0.5_ireals
  qgeops_d    = 0.9_ireals
  gnudg_d     = (/0.0006_ireals, 0.0012_ireals, 0.0006_ireals, 0.0006_ireals/)
  gnudgar_d   = (/0.0006_ireals, 0.0000_ireals, 0.0006_ireals, 0.0000_ireals/)
  gnudgsu_d   = (/0.0006_ireals, 0.0012_ireals, 0.0000_ireals, 0.0006_ireals/)
  gnudgtv_d   = (/0.0000_ireals, 0.0000_ireals, 0.0006_ireals, 0.0006_ireals/)
  gnudggp_d   =   0.0000_ireals
  ltipol_d    = .TRUE.
  ltipsu_d    = .TRUE. 
  tipolmx_d   = 1.0_ireals
  tipmxsu_d   = 1.0_ireals
  wtukrsa_d   = 3.0_ireals
  wtukrse_d   = 1.0_ireals
  wtukara_d   = 1.5_ireals
  wtukare_d   = 0.5_ireals
  wtuksua_d   = 1.5_ireals
  wtuksue_d   = 0.5_ireals
  msprpar_d   = 1
  msprpsu_d   = 0
  vcorls_d    = (/0.333_ireals, 0.333_ireals, 0.040_ireals, 0.040_ireals/)
  vcorlsu_d   = (/0.013_ireals, 0.013_ireals, 0.002_ireals, 0.002_ireals/)
  vcutof_d    = (/0.75_ireals, 0.75_ireals, 1.00_ireals, 1.00_ireals/)
  vcutosu_d   = (/0.75_ireals, 0.75_ireals, 4.00_ireals, 4.00_ireals/)
! vpblsu_d    = (/2.0_ireals, 2.0_ireals, 1.0_ireals, 1.0_ireals/)
  vpblsu_d    = (/99.0_ireals, 99.0_ireals, 99.0_ireals, 99.0_ireals/)
  wablua_d    = (/1.0_ireals, 1.0_ireals, 1.0_ireals, 1.0_ireals/)
  wablsu_d    = (/1.0_ireals, 1.0_ireals, 1.0_ireals, 1.0_ireals/)
  lsvcorl_d   = .TRUE.
  rhinfl_d    = (/  0.0_ireals,  70.0_ireals,   0.0_ireals,   0.0_ireals/)
  rhiflsu_d   = (/ 70.0_ireals,  70.0_ireals, 100.0_ireals,  70.0_ireals/)
  rhvfac_d    = (/1.00_ireals, 0.00_ireals, 0.83_ireals, 0.83_ireals/)
  rhtfac_d    = (/1.30_ireals, 1.43_ireals, 1.30_ireals, 1.30_ireals/)
  rhtfsu_d    = (/1.00_ireals, 1.43_ireals, 1.00_ireals, 1.00_ireals/)
  cutofr_d    = (/3.5_ireals, 3.5_ireals, 3.5_ireals, 3.5_ireals/)
  cutofsu_d   = (/2.0_ireals, 3.5_ireals, 2.0_ireals, 2.0_ireals/)
  vcsni_d     = (/2500.0_ireals, 2500.0_ireals, 2500.0_ireals, 2500.0_ireals/)
  vcsnisu_d   = (/2500.0_ireals, 2500.0_ireals,    9.0_ireals,    9.0_ireals/)
  rhfrtv_d    = (/1.0_ireals, 1.0_ireals, 1.0_ireals, 0.5_ireals/)
  rhfgps_d    = 0.45_ireals
  rhfpsdd_d   = 1.0_ireals
  cnondiv_d   = 0.1_ireals
  fnondiv_d   = 0.8_ireals
  tnondiv_d   = 1.1_ireals
  lscadj_d    = (/.TRUE., .TRUE.,.TRUE., .FALSE./)
  topobs_d    = (/ 849.0_ireals, 1099.0_ireals,  799.0_ireals,  699.0_ireals/)
  botmod_d    = (/1099.0_ireals, 1099.0_ireals, 1099.0_ireals,  899.0_ireals/)
  loiqv2m_d   = .FALSE.
  lqfqv2m_d   = .FALSE.
  dtqc_d      = 720.0_ireals
  qcc_d       = (/ 0.0_ireals, 500.0_ireals,  0.0_ireals, 0.7_ireals/)
  qccsu_d     = (/12.0_ireals, 500.0_ireals, 12.0_ireals, 0.7_ireals/)
  qcvf_d      = (/ 5.0_ireals,   1.0_ireals, 10.0_ireals, 0.0_ireals/)
  qcciq_d     = 1.0_ireals
  qcsiq_d     = .15_ireals
  qcflbcp_d   = 1.4_ireals
  qcfpst_d    = 1.5_ireals

  itype_obfile_d = 1
  irun_osse_d = 0
  ycdfdir_d   = './'
  yaofpath_d  = 'aof'
  obnlat_d    = 90.0_ireals
  obslat_d    = -90.0_ireals
  obwlon_d    = -180.0_ireals
  obelon_d    = 180.0_ireals
  exnlat_d    = 90.0_ireals
  exslat_d    = -90.0_ireals
  exwlon_d    = -180.0_ireals
  exelon_d    = 180.0_ireals
  doromx_d    = (/100.0_ireals, 150.0_ireals, 150.0_ireals, 150.0_ireals/)
  altopsu_d   = (/100.0_ireals, 5000.0_ireals, 5000.0_ireals, 5000.0_ireals/)
  thairh_d    = 20.0_ireals
  fperturb_d  = 0.0_ireals
  losse_fg_d  = .FALSE.
  lgpsbias_d  = .FALSE.
  maxmlo_d    = 300
  maxsgo_d    = 3000
  maxuso_d    = 900
  maxgpo_d    = 3000
  maxtvo_d    = 1
  maxmlv_d    = 100
  mxfrep_d    = -1
  mxfobs_d    = -1
  nolbc_d     = 5
  mqcorr92_d  = 0
  iseed_d     = 0
  mrhyes_d    = 0
  mrhno_d     = 0
  muvyes_d    = 0
  muvno_d     = 0
  igpscen_d   = -1
  lsynop_d    = .TRUE.
  laircf_d    = .TRUE.
  lsatob_d    = .FALSE.
  ldribu_d    = .TRUE.
  ltemp_d     = .TRUE.
  lpilot_d    = .TRUE.
  lsatem_d    = .FALSE.
  lgps_d      = .FALSE.
  lscatt_d    = .TRUE.
  lcd011_d    = .TRUE.
  lcd014_d    = .TRUE.
  lcd021_d    = .TRUE.
  lcd022_d    = .TRUE.
  lcd023_d    = .TRUE.
  lcd024_d    = .TRUE.
  lcd140_d    = .TRUE.
  lcd041_d    = .TRUE.
  lcd141_d    = .TRUE.
  lcd241_d    = .TRUE.
  lcd144_d    = .TRUE.
  lcd244_d    = .TRUE.
  lcd088_d    = .TRUE.
  lcd188_d    = .TRUE.
  lcd063_d    = .TRUE.
  lcd064_d    = .TRUE.
  lcd165_d    = .TRUE.
  lcd035_d    = .TRUE.
  lcd036_d    = .TRUE.
  lcd037_d    = .TRUE.
  lcd135_d    = .TRUE.
  lcd039_d    = .TRUE.
  lcd040_d    = .TRUE.
  lcd032_d    = .TRUE.
  lcd033_d    = .TRUE.
  lcd038_d    = .TRUE.
  lcd132_d    = .TRUE.
  lcd133_d    = .TRUE.
  lcd136_d    = .FALSE.
  lcd137_d    = .FALSE.
  lcd086_d    = .FALSE.
  lcd186_d    = .FALSE.
  lcd122_d    = .TRUE.
  lcd123_d    = .TRUE.
  lcd096_d    = .TRUE.
  mcdmsg1_d   = 0
  mcdmsg2_d   = 0
  mcdno15_d   = 0
  mcdno16_d   = 0
  mcdno17_d   = 0
  mcdno18_d   = 0
  lsurfa_d    = .FALSE.
  lt2m_d      = .FALSE.
  lrh2m_d     = .FALSE.
  lprecp_d    = .FALSE.
  lff10m_d    = .FALSE.
  ht2a_d      = 999.0_ireals
  ht2i_d      = 999.0_ireals
  hh2a_d      = 999.0_ireals
  hh2i_d      = 999.0_ireals
  hffa_d      = 999.0_ireals
  hffi_d      = 999.0_ireals
  hprc_d      = 999.0_ireals
  raintp_d    = 12.0_ireals
  ydir_lansfc_d  = './'
  yform_lansfc_d = 'grb1'
  lpraof_d    = .FALSE.
  lprodr_d    = .TRUE.
  ldiasa_d    = .FALSE.
  ionl_d      = 167
  jonl_d      = 103
  ionl2_d     = 167
  jonl2_d     = 103
  noctrq_d    = 9
  dinlat_d    = 0.0_ireals
  dislat_d    = 0.0_ireals
  diwlon_d    = 0.0_ireals
  dielon_d    = 0.0_ireals

!-------------------------------------------------------------------------------
!          1b: LHN
!-------------------------------------------------------------------------------

  llhn_d             = .FALSE.
  llhnverif_d        = .FALSE.
  lhn_search_d       = .TRUE.
  lhn_filt_d         = .TRUE.
  lhn_relax_d        = .TRUE.
  lhn_limit_d        = .TRUE.
  lhn_hum_adj_d      = .TRUE.
  lhn_spqual_d       = .FALSE.
  lhn_black_d        = .TRUE.
  lhn_incloud_d      = .TRUE.
  lhn_qrs_d          = .TRUE.
  lhn_logscale_d     = .TRUE.
  lhn_wweight_d      = .FALSE.
  lhn_diag_d         = .TRUE.
  lhn_height_d       = .TRUE.
  lhn_bright_d       = .TRUE.
  nlhn_start_d       = 0
  nlhn_end_d         = 360
  nlhnverif_start_d  = 0
  nlhnverif_end_d    = 360
  hlhn_start_d       = 0.0_ireals
  hlhn_end_d         = 3.0_ireals
  hlhn_end_dum       =-3.0_ireals
  hlhnverif_start_d  = 0.0_ireals
  hlhnverif_end_d    = 3.0_ireals
  hlhnverif_end_dum  =-3.0_ireals
  rlhn_search_d      = 10
  ktop_lhn_d         = 1
  kbot_lhn_d         = ke_tot
  nradar_d           = 35
  nlhn_relax_d       = 2_iintegers
  lhn_dt_obs_d       = 5.0_ireals
  lhn_coef_d         = 1.0_ireals
  abs_lhn_lim_d      = 50._ireals / 3600._ireals    ! max. change in heating: 4K/h
  fac_lhn_search_d   = 5._ireals
  fac_lhn_up_d       = 2.0_ireals
  fac_lhn_down_d     = 1.0_ireals / 2.0_ireals
  rad_wobs_lhn_d     = 100._ireals
  thres_lhn_d        = 0.1_ireals / 3600._ireals
  radar_in_d         = '.'
  blacklist_file_d   = 'blacklist_dx.grib1'
  height_file_d      = 'height_dx.grib1'
  noobs_date_d(:)    = '            '
  rqrsgmax_d         = 0.0
  ktop_temp_d        = -999.9

!-------------------------------------------------------------------------------
!- Section 2: Initialize variables with defaults
!-------------------------------------------------------------------------------
!         2a: NUDGING
!-------------------------------------------------------------------------------

  lnudge      = lnudge_d
  lverif      = lverif_d
  lverpas     = lverpas_d
  l1dvar      = l1dvar_d
! lobdens     = lobdens_d
  nudgsta     = nudgsta_d
  nudgend     = nudgend_d
  nversta     = nversta_d
  nverend     = nverend_d
  mruntyp     = mruntyp_d
  mveripr     = mveripr_d
  nwtyp       = nwtyp_d
  niwtyp      = niwtyp_d
  iwtyp       = iwtyp_d
  kwtyp       = kwtyp_d
  hnudgsta    = hnudgsta_d
  hnudgend    = hnudgend_dum
  hversta     = hversta_d
  hverend     = hverend_dum
  tconbox     = tconbox_d
  ntpscor     = ntpscor_d
  khumbal     = khumbal_d
  ptpstop     = ptpstop_d
  luvgcor     = luvgcor_d
  mpsgcor     = mpsgcor_d
  qgeo        = qgeo_d 
  qgeotop     = qgeotop_d
  qgeops      = qgeops_d
  gnudg       = gnudg_d
  gnudgar     = gnudgar_d
  gnudgsu     = gnudgsu_d
  gnudgtv     = gnudgtv_d
  gnudggp     = gnudggp_d
  ltipol      = ltipol_d
  ltipsu      = ltipsu_d
  tipolmx     = tipolmx_d
  tipmxsu     = tipmxsu_d
  wtukrsa     = wtukrsa_d
  wtukrse     = wtukrse_d
  wtukara     = wtukara_d
  wtukare     = wtukare_d
  wtuksua     = wtuksua_d
  wtuksue     = wtuksue_d
  msprpar     = msprpar_d
  msprpsu     = msprpsu_d
  vcorls      = vcorls_d
  vcorlsu     = vcorlsu_d
  vcutof      = vcutof_d
  vcutosu     = vcutosu_d
  vpblsu      = vpblsu_d
  wablua      = wablua_d
  wablsu      = wablsu_d
  lsvcorl     = lsvcorl_d
  rhinfl      = rhinfl_d
  rhiflsu     = rhiflsu_d
  rhvfac      = rhvfac_d 
  rhtfac      = rhtfac_d
  rhtfsu      = rhtfsu_d
  cutofr      = cutofr_d
  cutofsu     = cutofsu_d
  vcsni       = vcsni_d
  vcsnisu     = vcsnisu_d
  rhfrtv      = rhfrtv_d 
  rhfgps      = rhfgps_d 
  rhfpsdd     = rhfpsdd_d 
  cnondiv     = cnondiv_d
  fnondiv     = fnondiv_d
  tnondiv     = tnondiv_d
  lscadj      = lscadj_d
  topobs      = topobs_d
  botmod      = botmod_d
  loiqv2m     = loiqv2m_d
  lqfqv2m     = lqfqv2m_d
  dtqc        = dtqc_d
  qcc         = qcc_d
  qccsu       = qccsu_d 
  qcvf        = qcvf_d
  qcciq       = qcciq_d
  qcsiq       = qcsiq_d
  qcflbcp     = qcflbcp_d
  qcfpst      = qcfpst_d

  itype_obfile= itype_obfile_d
  irun_osse   = irun_osse_d
  ycdfdir     = ycdfdir_d
  yaofpath    = yaofpath_d
  obnlat      = obnlat_d
  obslat      = obslat_d
  obwlon      = obwlon_d
  obelon      = obelon_d
  exnlat      = exnlat_d
  exslat      = exslat_d
  exwlon      = exwlon_d
  exelon      = exelon_d
  doromx      = doromx_d
  altopsu     = altopsu_d
  thairh      = thairh_d
  fperturb    = fperturb_d
  losse_fg    = losse_fg_d
  lgpsbias    = lgpsbias_d
  maxmlo      = maxmlo_d
  maxsgo      = maxsgo_d
  maxuso      = maxuso_d
  maxgpo      = maxgpo_d
  maxtvo      = maxtvo_d
  maxmlv      = maxmlv_d
  mxfrep      = mxfrep_d
  mxfobs      = mxfobs_d
  nolbc       = nolbc_d
  mqcorr92    = mqcorr92_d
  iseed       = iseed_d
  mrhyes      = mrhyes_d
  mrhno       = mrhno_d
  muvyes      = muvyes_d 
  muvno       = muvno_d  
  lsynop      = lsynop_d
  laircf      = laircf_d
  lsatob      = lsatob_d
  ldribu      = ldribu_d
  ltemp       = ltemp_d
  lpilot      = lpilot_d
  lsatem      = lsatem_d
  lgps        = lgps_d
  lscatt      = lscatt_d
  lcd011      = lcd011_d
  lcd014      = lcd014_d
  lcd021      = lcd021_d
  lcd022      = lcd022_d
  lcd023      = lcd023_d
  lcd024      = lcd024_d
  lcd140      = lcd140_d
  lcd041      = lcd041_d
  lcd141      = lcd141_d
  lcd241      = lcd241_d
  lcd144      = lcd144_d
  lcd244      = lcd244_d
  lcd088      = lcd088_d
  lcd188      = lcd188_d
  lcd063      = lcd063_d
  lcd064      = lcd064_d
  lcd165      = lcd165_d
  lcd035      = lcd035_d
  lcd036      = lcd036_d
  lcd037      = lcd037_d
  lcd135      = lcd135_d
  lcd039      = lcd039_d
  lcd040      = lcd040_d
  lcd032      = lcd032_d
  lcd033      = lcd033_d
  lcd038      = lcd038_d
  lcd132      = lcd132_d
  lcd133      = lcd133_d
  lcd136      = lcd136_d
  lcd137      = lcd137_d
  lcd086      = lcd086_d
  lcd186      = lcd186_d
  lcd122      = lcd122_d
  lcd123      = lcd123_d
  lcd096      = lcd096_d
  mcdmsg1     = mcdmsg1_d
  mcdmsg2     = mcdmsg2_d
  mcdno15     = mcdno15_d
  mcdno16     = mcdno16_d
  mcdno17     = mcdno17_d
  mcdno18     = mcdno18_d
  lsurfa      = lsurfa_d
  lt2m        = lt2m_d
  lrh2m       = lrh2m_d
  lprecp      = lprecp_d
  lff10m      = lff10m_d
  ht2a        = ht2a_d
  ht2i        = ht2i_d
  hh2a        = hh2a_d
  hh2i        = hh2i_d
  hffa        = hffa_d
  hffi        = hffi_d
  hprc        = hprc_d
  raintp      = raintp_d
  ydir_lansfc = ydir_lansfc_d
  yform_lansfc= yform_lansfc_d
  lpraof      = lpraof_d
  lprodr      = lprodr_d
  ldiasa      = ldiasa_d
  ionl        = ionl_d
  jonl        = jonl_d
  ionl2       = ionl2_d
  jonl2       = jonl2_d
  noctrq      = noctrq_d
  dinlat      = dinlat_d
  dislat      = dislat_d
  diwlon      = diwlon_d
  dielon      = dielon_d
  igpscen     = igpscen_d

!-------------------------------------------------------------------------------
!         2a: LHN
!-------------------------------------------------------------------------------

  llhn            = llhn_d
  llhnverif       = llhnverif_d
  lhn_search      = lhn_search_d
  lhn_filt        = lhn_filt_d
  lhn_relax       = lhn_relax_d
  lhn_limit       = lhn_limit_d
  lhn_hum_adj     = lhn_hum_adj_d
  lhn_diag        = lhn_diag_d
  lhn_qrs         = lhn_qrs_d
  lhn_logscale    = lhn_logscale_d
  lhn_wweight     = lhn_wweight_d
  lhn_spqual      = lhn_spqual_d
  lhn_black       = lhn_black_d
  lhn_incloud     = lhn_incloud_d
  lhn_height      = lhn_height_d
  lhn_bright      = lhn_bright_d
  nlhn_start      = nlhn_start_d
  nlhn_end        = nlhn_end_d
  nlhnverif_start = nlhnverif_start_d
  nlhnverif_end   = nlhnverif_end_d
  hlhn_start      = hlhn_start_d
  hlhn_end        = hlhn_end_dum
  hlhnverif_start = hlhnverif_start_d
  hlhnverif_end   = hlhnverif_end_dum
  rlhn_search     = rlhn_search_d
  ktop_lhn        = ktop_lhn_d
  kbot_lhn        = kbot_lhn_d
  nradar          = nradar_d
  nlhn_relax      = nlhn_relax_d
  lhn_dt_obs      = lhn_dt_obs_d
  lhn_coef        = lhn_coef_d
  abs_lhn_lim     = abs_lhn_lim_d
  fac_lhn_search  = fac_lhn_search_d
  fac_lhn_up      = fac_lhn_up_d
  fac_lhn_down    = fac_lhn_down_d
  rad_wobs_lhn    = rad_wobs_lhn_d
  thres_lhn       = thres_lhn_d
  rqrsgmax        = rqrsgmax_d
  ktop_temp       = ktop_temp_d
  radar_in        = radar_in_d
  blacklist_file  = blacklist_file_d
  height_file     = height_file_d
  noobs_date(:)   = noobs_date_d(:)

!-------------------------------------------------------------------------------
!- Section 3: Input of the namelist values
!-------------------------------------------------------------------------------

!HJP Begin 2016-05-11
! READ (nuin, nudging, IOSTAT=iz_err)
  iomsg_str(:) = ' '
  READ (nuin, nudging, IOSTAT=iz_err, IOMSG=iomsg_str)

  IF (iz_err /= 0) WRITE (*,'(A,A)') 'Namelist-ERROR NUDGING: ', TRIM(iomsg_str)
!HJP END   2016-05-11
ENDIF

IF (nproc > 1) THEN
  ! distribute error status to all processors
  CALL distribute_values  (iz_err, 1, 0, imp_integers,  icomm_world, ierr)
ENDIF

IF (iz_err /= 0) THEN
  ierrstat = -1
  RETURN
ENDIF

IF (my_world_id == 0) THEN

!-------------------------------------------------------------------------------
!- Section 4: Check values for errors and consistency
!-------------------------------------------------------------------------------

  dtdoh = dt / 3600.0_ireals

  IF (.NOT. luseobs) lnudge    = .FALSE.
  IF (.NOT. luseobs) lverif    = .FALSE.
  IF (.NOT. luseobs) l1dvar    = .FALSE.
  IF (.NOT. luseobs) lsurfa    = .FALSE.
  IF (.NOT. luseobs) llhn      = .FALSE.
  IF (.NOT. luseobs) llhnverif = .FALSE.
  IF (.NOT. lverif)  lverpas   = .FALSE.

! Check lsurfa and llhn: these must be false, unless GRIBDWD or GRIBAPI is defined
#if !(defined GRIBDWD) && !defined(GRIBAPI)
  IF (llhn) THEN
    PRINT *, ' ERROR  *** llhn = .TRUE., but model is not compiled to use GRIB library ***'
    ierrstat = 9003
  ENDIF

  IF (lsurfa) THEN
    PRINT *, ' ERROR  *** lsurfa = .TRUE., but model is not compiled to use GRIB library ***'
    ierrstat = 9004
  ENDIF
#endif

  IF ((l1dvar) .AND. (luse_rttov)) THEN
    WRITE( nuspecif,'("luse_rttov =.TRUE. ==> set l1dvar =.FALSE. !!!")' )
    l1dvar = .FALSE.
  ENDIF

! hnudgend = c1

  IF (MAX( ABS(hnudgsta-hnudgsta_d) , ABS(hnudgend-hnudgend_dum) ) > epsy) THEN
    nudgsta  =  INT( (hnudgsta - epsy) /dtdoh + c1 )
    nudgend  =  INT( (hnudgend + epsy) /dtdoh      )
  ELSEIF (MAX( ABS(nudgsta-nudgsta_d) , ABS(nudgend-nudgend_d) ) > 0) THEN
    hnudgsta =  nudgsta * dtdoh
    hnudgend =  nudgend * dtdoh
  ELSE
    hnudgend =  hnudgend_d
    nudgsta  =  INT( (hnudgsta - epsy) /dtdoh + c1 )
    nudgend  =  INT( (hnudgend + epsy) /dtdoh      ,iintegers)
  ENDIF
  IF (nudgsta > nudgend) lnudge = .FALSE.
  IF (.NOT. lnudge) THEN
    nudgsta  =  0
    nudgend  = -1
    hnudgsta = 0.0_ireals
    hnudgend = 0.0_ireals
  ENDIF
  IF (MAX( ABS( hversta-hversta_d ) , ABS( hverend-hverend_d ) ) > epsy) THEN
    nversta  =  INT( (hversta - epsy) /dtdoh + c1 )
    nverend  =  INT( (hverend + epsy) /dtdoh      ,iintegers )
  ELSEIF (MAX( ABS( nversta-nversta_d ) , ABS( nverend-nverend_d ) ) > 0) THEN
    hversta  =  nversta * dtdoh
    hverend  =  nverend * dtdoh
  ELSE
    hverend  =  hverend_d
    nversta  =  INT( (hversta - epsy) /dtdoh + c1 )
    nverend  =  INT( (hverend + epsy) /dtdoh      ,iintegers )
  ENDIF
  IF (nversta < 0) THEN
    nversta  =  0
    hversta  =  0.0_ireals
  ENDIF
  IF (nverend > nstop) THEN
    nverend  =  nstop
    hverend  =  nverend * dtdoh
  ENDIF
  IF (nversta == 1) hversta  =  MAX( hversta , 0.0001_ireals )
  IF (nversta > nverend) lverif = .FALSE.
  IF (.NOT. lverif) THEN
    nversta  =  0
    nverend  = -999
    hversta  =  0.0_ireals
    hverend  = -1.0_ireals
  ENDIF

  IF (.NOT. lverif)  mruntyp = -1
  IF (.NOT. lverif)  mveripr =  0
  IF (      lverif)  mveripr = MAX( 1, MIN( 3, mveripr ) )
  ! if obs are read from AOF file then writing NetCDF feedobs file is prohibited
  IF (itype_obfile == 1)  mveripr = (mveripr /2 ) * 2
  irun_osse  =  ABS( irun_osse )

  IF (nwtyp >= 1)  niwtyp (nwtyp+1:mxwt) = 0
  IF (     (nwtyp <= 0) .OR. (nwtyp > mxwt)                                    &
      .OR. (MINVAL( niwtyp ) <= -1) .OR. (SUM( niwtyp ) > mxtyw)) THEN
    IF ((nwtyp <= 0) .OR. (nwtyp > mxwt)) THEN
      WRITE( nuspecif,'("nwtyp not in [1,mxwt] ==> nwtyp set to 1 ,")' )
    ELSEIF (MINVAL( niwtyp ) <= -1) THEN
      WRITE( nuspecif,'("negative values present but not allowed in niwtyp "   &
                      &,"==> nwtyp set to 1 ,")' )
    ELSE
      WRITE( nuspecif,'("sum of values in niwtyp > mxtyw (not allowed) "       &
                      &,"==> nwtyp set to 1 ,")' )
    ENDIF
    WRITE( nuspecif,'("  niwtyp set to (1,0,0,...) a"                          &
                    &,"nd iwtyp set to (0,0,0,...)  !!")')
    nwtyp     =  1
    niwtyp    =  0
    niwtyp(1) =  1
    iwtyp     =  0
    kwtyp     =  1
  ENDIF
! in array 'iwtyp': - values < -999 or > mxobtp must not occur
!                   - all values must not occur more than once
!                   - if '0' does not occur, write CAUTION
  iwtmin   =  1000
  iwtmax   = -1000
  liwident = .FALSE.
  liwrest  = .FALSE.
  DO iiwt = SUM( niwtyp ) , 1 , -1
    IF (iwtyp(iiwt) < iwtmin)  iwtmin  = iwtyp(iiwt)
    IF (iwtyp(iiwt) > iwtmax)  iwtmax  = iwtyp(iiwt)
    IF (iwtyp(iiwt) ==  0   )  liwrest = .TRUE.
    DO iiw2 = iiwt-1 , 1 , -1
      IF (iwtyp(iiw2) == iwtyp(iiwt)) liwident = .TRUE.
    ENDDO
  ENDDO
  IF ((iwtmin < -999) .OR. (iwtmax > mxobtp) .OR. (liwident)) THEN
    IF ((iwtmin < -999) .OR. (iwtmax > mxobtp)) THEN
      WRITE( nuspecif,'("values of iwtyp must be within [-999,mxobtp] "        &
                      &,"==> nwtyp set to 1 ,")' )
    ELSE
      WRITE( nuspecif,'("identical values in iwtyp are not allowed "           &
                      &,"==> nwtyp set to 1 ,")' )
    ENDIF
    WRITE( nuspecif,'("  niwtyp set to (1,0,0,...) a"                          &
                    &,"nd iwtyp set to (0,0,0,...)  !!")')
    nwtyp     =  1
    niwtyp    =  0
    niwtyp(1) =  1
    iwtyp     =  0
    kwtyp     =  1
  ELSEIF (.NOT. liwrest) THEN
    IF (SUM( niwtyp ) < mxtyw) THEN
      DO iiwt = SUM( niwtyp ) + 1 , mxtyw 
        iwtyp (iiwt) = 0
      ENDDO
      nwtyp          =  nwtyp + 1
      niwtyp (nwtyp) = 1
    ELSE
      WRITE( nuspecif,'("Caution: be sure that iwtyp really contains all "     &
                      &,"observation (code) types !")' )
    ENDIF
  ENDIF
  DO iiwt = 1 , nwtyp+2
    kwtyp (iiwt) = MAX( 1 , MIN( 2 , kwtyp(iiwt) ) )
  ENDDO

  IF (tconbox <= dt     +2*epsy) tconbox = dt
  IF (dtqc    <= tconbox+2*epsy) dtqc    = tconbox

  IF ((mpsgcor > 2) .OR. (mpsgcor < 0)) THEN
    mpsgcor = mpsgcor_d
    WRITE( nuspecif,'(''mpsgcor not in [0,2] ==> mpsgc'' &
                  &,''or set to'',I2,'' !!'')' ) mpsgcor
  ENDIF
  IF ((ntpscor > 2) .OR. (ntpscor < 0)) THEN
    ntpscor = ntpscor_d
    WRITE( nuspecif,'(''ntpscor not in [0,2] ==> ntpsc'' &
                  &,''or set to'',I2,'' !!'')' ) ntpscor
  ENDIF
! use dfi-dfis instead of dp-dps -correlation for temperature correction
  IF (ntpscor == 2)  ntpscor = 4
  IF ((msprpar > 2) .OR. (msprpar < 0)) THEN
    msprpar = msprpar_d
    WRITE( nuspecif,'(''msprpar not in [0,2] ==> msprp'' &
                  &,''ar set to'',I2,'' !!'')' ) msprpar
  ENDIF
  IF ((msprpsu > 2) .OR. (msprpsu < 0)) THEN
    msprpsu = msprpsu_d
    WRITE( nuspecif,'(''msprpsu not in [0,2] ==> msprp'' &
                  &,''su set to'',I2,'' !!'')' ) msprpsu
  ENDIF
! checking of ptpstop: needs values for pp,p0 !!!  ==> check in 'gather_varia.h'
  IF ((qgeo    < c0) .OR. (qgeo    > c1)) qgeo    = qgeo_d
  IF ((qgeotop < c0) .OR. (qgeotop > c1)) qgeotop = qgeotop_d
  IF ((qgeops  < c0) .OR. (qgeops  > c1)) qgeops  = qgeops_d

  IF (.NOT. ltipol) tipolmx = c0
  IF (.NOT. ltipsu) tipmxsu = c0
  IF (lnudge) THEN
    tipolmx = MAX( tipolmx , c0 )
    tipmxsu = MAX( tipmxsu , c0 )
    wtukrsa = MAX( wtukrsa , dtdoh )
    wtukrse = MAX( wtukrse , dtdoh )
    wtukara = MAX( wtukara , dtdoh )
    wtukare = MAX( wtukare , dtdoh )
    wtuksua = MAX( wtuksua , dtdoh )
    wtuksue = MAX( wtuksue , dtdoh )
  ELSE
    tipolmx = c0
    tipmxsu = c0
    wtukrsa = dtdoh
    wtukrse = dtdoh
    wtukara = dtdoh
    wtukare = dtdoh
    wtuksua = dtdoh
    wtuksue = dtdoh
  ENDIF
  IF (tipolmx < epsy) ltipol = .FALSE.
  IF (tipmxsu < epsy) ltipsu = .FALSE.
  rhfgps  = MIN( MAX( rhfgps  , epsy ) , c1 )
  rhfpsdd = MIN( MAX( rhfpsdd , epsy ) , c1 )
  cnondiv = MAX( cnondiv , c0 )
  fnondiv = MAX( fnondiv , c0 )
  tnondiv = MAX( tnondiv , c0 )
  qcsiq   = MIN( qcsiq   , c1 )
  qcfpst  = MAX( MIN( qcfpst , 2.0_ireals ) , c1 )

  DO ivar = 1 , 4
    IF (lnudge) THEN
      gnudg   (ivar) = MAX( gnudg  (ivar) , c0 )
      gnudgar (ivar) = MAX( gnudgar(ivar) , c0 )
      gnudgsu (ivar) = MAX( gnudgsu(ivar) , c0 )
      gnudgtv (ivar) = MAX( gnudgtv(ivar) , c0 )
      gnudggp        = MAX( gnudggp       , c0 )
    ELSE
      gnudg   (ivar) = c0
      gnudgar (ivar) = c0
      gnudgsu (ivar) = c0
      gnudgtv (ivar) = c0
      gnudggp        = c0
    ENDIF
    vcorls  (ivar) = MAX( vcorls (ivar) , 0.00001_ireals )
    vcorlsu (ivar) = MAX( vcorlsu(ivar) , 0.00001_ireals )
    vpblsu  (ivar) = MAX( vpblsu (ivar) , 0.001_ireals )
    vcutof  (ivar) = MAX( vcutof (ivar) , c0 )
    vcutosu (ivar) = MAX( vcutosu(ivar) , c0 )
    wablua  (ivar) = MAX( wablua (ivar) , c0 )
    wablsu  (ivar) = MAX( wablsu (ivar) , c0 )
    rhinfl  (ivar) = MAX( rhinfl (ivar) , c0 )
    rhiflsu (ivar) = MAX( rhiflsu(ivar) , c0 )
    rhvfac  (ivar) = MAX( rhvfac (ivar) , c0 )
    rhtfac  (ivar) = MAX( rhtfac (ivar) , c0 )
    rhtfsu  (ivar) = MAX( rhtfsu (ivar) , c0 )
    cutofr  (ivar) = MAX( cutofr (ivar) , c0 )
    cutofsu (ivar) = MAX( cutofsu(ivar) , c0 )
    vcsni   (ivar) = MAX( vcsni  (ivar) , 0.001_ireals )
    vcsnisu (ivar) = MAX( vcsnisu(ivar) , 0.001_ireals )
    rhfrtv  (ivar) = MIN( MAX( rhfrtv(ivar)  , epsy ) , c1 )
    IF (msprpar == 0) THEN
      topobs (ivar) = 1099.0_ireals
      botmod (ivar) = 1099.0_ireals
    ENDIF
    topobs  (ivar) = MAX( topobs (ivar) , c1 )
    botmod  (ivar) = MAX( botmod (ivar) , topobs(ivar) )
    qcc     (ivar) = MAX( qcc    (ivar) , c0 )
    qccsu   (ivar) = MAX( qccsu  (ivar) , c0 )
    qcvf    (ivar) = MAX( qcvf   (ivar) , c0 )
    doromx  (ivar) = MAX( doromx (ivar) , 0.001_ireals )
  ENDDO
  maxmlo   = MAX( maxmlo  , 1_iintegers )
  maxsgo   = MAX( maxsgo  , 1_iintegers )
  maxuso   = MAX( maxuso  , 1_iintegers )
  maxgpo   = MAX( maxgpo  , 1_iintegers )
  maxtvo   = MAX( maxtvo  , 1_iintegers )
  maxmlv   = MAX( maxmlv  , 1_iintegers )
  nolbc    = MAX( nolbc   , 2_iintegers )
  mqcorr92 = MAX( MIN( mqcorr92, 2_iintegers ) , 0_iintegers )

  DO ivar = mxgpc, 1, -1
    IF ((igpscen(ivar) < 0) .OR. (igpscen(ivar) > 37)) THEN
      igpscen(ivar) = igpscen_d(ivar)
    ENDIF
    IF (ivar > 1) THEN
      DO k = 1 , ivar-1
        IF (igpscen(ivar) == igpscen(k)) igpscen(ivar) = -1
      ENDDO
    ENDIF
  ENDDO
  nactgpc = 0
  noffset = 0
  DO ivar = 1, mxgpc
    IF (igpscen(ivar) == -1) THEN
      noffset = noffset + 1
    ELSE
      nactgpc = nactgpc + 1
      IF (noffset >= 1) THEN
        igpscen (ivar-noffset) = igpscen(ivar)
        igpscen (ivar)         = -1
      ENDIF
    ENDIF
  ENDDO

  ! Check length of the directory name for NetCDF obs input files
  nzylen  = LEN_TRIM(ycdfdir)
  IF (nzylen >= LEN( ycdfdir )) THEN
    PRINT *,' WARNING  *** max. length of "ycdfdir" is ',icdfdirlen,' letters ***'
  ENDIF

  ! Check length of the path name for AOF
  nzylen  = LEN_TRIM(yaofpath)
  IF (nzylen == 0) THEN
    PRINT *,' ERROR    *** yaofpath is empty    *** '
    ierrstat = 9002
  ELSEIF (nzylen >= LEN( yaofpath )) THEN
    PRINT *,' WARNING  *** max. length of "yaofpath" is 70 letters *** '
  ENDIF

!         Check for valid observation area boundary
  IF ( obnlat.GT.90.1_ireals    .OR.    obnlat.LT.-90.1_ireals .OR.            &
       obslat.GT.90.1_ireals    .OR.    obslat.LT.-90.1_ireals .OR.            &
       obwlon.LT.-180.1_ireals  .OR.    obwlon.GT.360.1_ireals .OR.            &
       obelon.LT.-180.1_ireals  .OR.    obelon.GT.360.1_ireals) THEN
      WRITE(nuspecif,                                                          &
      & '(///,"***   ***   ***   ***   ***   ***   ***   ***",///)')
      WRITE(nuspecif,'(" INVALID OBSERVATION AREA BOUNDARY")')
      WRITE(nuspecif,'(" OBNLAT = ",F20.6)') obnlat
      WRITE(nuspecif,'(" OBSLAT = ",F20.6)') obslat
      WRITE(nuspecif,'(" OBWLON = ",F20.6)') obwlon
      WRITE(nuspecif,'(" OBELON = ",F20.6)') obelon
!US   yerrmsg = 'invalid obs. area'
!US   CALL model_abort (my_world_id, 13211, yerrmsg, yroutine)
      PRINT *, ' *** invalid observation area ***'
      ierrstat = 9001
  ENDIF

!         Check for empty observation area
  IF (obnlat.LE.obslat.OR. obwlon.EQ.obelon) THEN
      WRITE(nuspecif,                                                          &
      & '(///,"***   ***   ***   ***   ***   ***   ***   ***",///)')
      WRITE(nuspecif,'(" OBSERVATION AREA EMPTY",/,                            &
               & " NO OBSERVATIONS WILL BE USED IN THE ANALYSIS")')
      WRITE(nuspecif,'(" OBNLAT = ",F20.6)') obnlat
      WRITE(nuspecif,'(" OBSLAT = ",F20.6)') obslat
      WRITE(nuspecif,'(" OBWLON = ",F20.6)') obwlon
      WRITE(nuspecif,'(" OBELON = ",F20.6)') obelon
      WRITE(nuspecif,'(" OBELON = ",F20.6)') obelon
      WRITE(nuspecif,                                                          &
      & '(///,"***   ***   ***   ***   ***   ***   ***   ***",///)')
  ENDIF
!
  IF (obwlon.GT.180.0_ireals) obwlon=obwlon-360.0_ireals
  IF (obelon.GT.180.0_ireals) obelon=obelon-360.0_ireals

!         If lsurfa is FALSE reset submode switches and vice verse
!         Surface analyses possible only in nudging run mode
! IF (.NOT. lnudge)    lsurfa=.FALSE.
  IF (.NOT. lsurfa)                              THEN
     lt2m   = .FALSE.
     lrh2m  = .FALSE.
     lprecp = .FALSE.
     lff10m = .FALSE.
  ELSEIF (.NOT.lt2m .AND. .NOT.lrh2m .AND. .NOT.lff10m .AND. .NOT.lprecp ) THEN
     lsurfa = .FALSE.
  ENDIF

  IF ((TRIM(yform_lansfc) /= 'grb1') .AND.                                      &
      (TRIM(yform_lansfc) /= 'api1') .AND. (TRIM(yform_lansfc) /= 'api2') ) THEN
    PRINT *,' ERROR    *** yform_lansfc not valid ', TRIM(yform_lansfc)
    ierrstat = 1002
#ifndef GRIBDWD
  ELSEIF (yform_lansfc == 'grb1') THEN
    PRINT *, ' ERROR  *** yform_lansfc = grb1, but model is not compiled to use DWD GRIB library ***'
    ierrstat = 1002
#endif
#ifndef GRIBAPI
  ELSEIF (TRIM(yform_lansfc(1:3)) == 'api') THEN
    PRINT *,' ERROR    *** yform_lansfc = api, but model is not compiled to use grib-api library *** '
    ierrstat = 1002
#endif
  ENDIF

! LHN

  IF (MAX( ABS(hlhn_start-hlhn_start_d) ,                 &
           ABS(hlhn_end-hlhn_end_dum) ) > epsy) THEN
    nlhn_start  =  INT( (hlhn_start - epsy) /dtdoh + c1 )
    nlhn_end    =  INT( (hlhn_end + epsy) /dtdoh      )
  ELSEIF (MAX( ABS(nlhn_start-nlhn_start_d) ,             &
               ABS(nlhn_end-nlhn_end_d) ) > 0) THEN
    hlhn_start  =  nlhn_start * dtdoh
    hlhn_end    =  nlhn_end * dtdoh
  ELSE
    hlhn_end    =  hlhn_end_d
    nlhn_start  =  INT( (hlhn_start - epsy) /dtdoh + c1 )
    nlhn_end    =  INT( (hlhn_end + epsy) /dtdoh      ,iintegers)
  ENDIF
  IF (nlhn_start == 0) nlhn_start = 1
  IF (nlhn_start > nlhn_end) llhn = .FALSE.
  IF (.NOT. llhn) THEN
    nlhn_start  =  0
    nlhn_end    = -1
    hlhn_start  = 0.0_ireals
    hlhn_end    = 0.0_ireals
  ENDIF

  IF (MAX( ABS(hlhnverif_start-hlhnverif_start_d) ,       &
           ABS(hlhnverif_end-hlhnverif_end_dum) ) > epsy) THEN
    nlhnverif_start  =  INT( (hlhnverif_start - epsy) /dtdoh + c1 )
    nlhnverif_end    =  INT( (hlhnverif_end + epsy) /dtdoh      )
  ELSEIF (MAX( ABS(nlhnverif_start-nlhnverif_start_d) ,   &
               ABS(nlhnverif_end-nlhnverif_end_d) ) > 0) THEN
    hlhnverif_start  =  nlhnverif_start * dtdoh
    hlhnverif_end    =  nlhnverif_end * dtdoh
  ELSE
    hlhnverif_end    =  hlhnverif_end_d
    nlhnverif_start  =  INT( (hlhnverif_start - epsy) /dtdoh + c1 )
    nlhnverif_end    =  INT( (hlhnverif_end + epsy) /dtdoh      ,iintegers)
  ENDIF
  IF (nlhnverif_start == 0) nlhnverif_start = 1
  IF (nlhnverif_start > nlhnverif_end) llhnverif = .FALSE.
  IF (.NOT. llhnverif) THEN
    nlhnverif_start  =  0
    nlhnverif_end    = -1
    hlhnverif_start  = 0.0_ireals
    hlhnverif_end    = 0.0_ireals
  ENDIF

! check latent heat nudging parameters, impose limits
  IF (rlhn_search > ie_tot/2.) THEN
    rlhn_search = INT(MIN (REAL(rlhn_search,ireals),    &
                           REAL(ie_tot,ireals)/2.0_ireals))
  ENDIF
  IF (rlhn_search > je_tot/2.) THEN
    rlhn_search = INT(MIN (REAL(rlhn_search,ireals),    &
                           REAL(je_tot,ireals)/2.0_ireals))
  ENDIF
  ktop_lhn      = MAX(ktop_lhn,1)
  kbot_lhn      = MIN(kbot_lhn,ke_tot)
! the following parameters are to be definitively specified later
   nlhn_relax    = MAX(nlhn_relax,1_iintegers)
   abs_lhn_lim   = MIN(abs_lhn_lim,10.0_ireals)
   fac_lhn_up    = MIN(fac_lhn_up,10.0_ireals)
   fac_lhn_down  = MAX(fac_lhn_down,1.0_ireals/10.0_ireals)
   fac_lhn_search = MAX(fac_lhn_search,1.5_ireals*fac_lhn_up)
   rad_wobs_lhn  = MIN(rad_wobs_lhn,200.0_ireals)
   thres_lhn     = MIN(thres_lhn,0.5_ireals/3600._ireals)
   if (thres_lhn < 0.0_ireals) thres_lhn = 0.0_ireals
   fac_lhn_down_in   = fac_lhn_down
   fac_lhn_up_in     = fac_lhn_up
   fac_lhn_search_in = fac_lhn_search

  IF (lhn_logscale) THEN
     fac_lhn_down   = 1_ireals + log(fac_lhn_down)
     fac_lhn_up     = 1_ireals + log(fac_lhn_up)
     fac_lhn_search = 1_ireals + log(fac_lhn_search)
  ENDIF

  ndelobs_lhn  = INT(lhn_dt_obs*60/dt,iintegers)

  IF ( nlhn_start > 1 .AND. nlhnverif_start > 1 ) then
   nlhn_start_dum  = ndelobs_lhn - mod(nlhn_start,ndelobs_lhn)
   nlhnv_start_dum = ndelobs_lhn - mod(nlhnverif_start,ndelobs_lhn)
   IF (nlhn_start <= nlhnverif_start .AND. nlhn_start_dum /= ndelobs_lhn)  &
      nlhn_start = nlhn_start + nlhn_start_dum
   IF (nlhnverif_start <= nlhn_start .AND. nlhnv_start_dum /= ndelobs_lhn) &
      nlhnverif_start = nlhnverif_start + nlhnv_start_dum
   IF (nlhn_start_dum /= ndelobs_lhn .OR. nlhnv_start_dum /= ndelobs_lhn) then
    PRINT *,'WARNING (LHN): '
    PRINT *,'Radarinformation not available for lhn_start or lhnverif_start.'
    PRINT *,'The chosen time does not match the interval of observation'
    PRINT *,'Timestep of start is set to the next observational time: ',   &
               nlhn_start, nlhnverif_start
   ENDIF
  ENDIF

  IF (.NOT.llhn .AND. .NOT.llhnverif) THEN
   lhn_qrs      =.false.
  ENDIF
  IF ( lhn_bright ) lhn_height = .TRUE.
ENDIF  ! (my_world_id == 0)

!-------------------------------------------------------------------------------
!- Section 5: Distribute variables to all nodes
!-------------------------------------------------------------------------------

IF (nproc > 1) THEN
  nnl = 0

  IF (my_world_id == 0) THEN

!   logbuf for nudging
    logbuf (  1) =  lnudge
    logbuf (  2) =  lverif
    logbuf (  3) =  lverpas
    logbuf (  4) =  l1dvar
    logbuf (  5) =  luvgcor
    logbuf (  6) =  lsvcorl
    logbuf (  7) =  ltipol
    logbuf (  8) =  ltipsu  
    logbuf (  9) =  loiqv2m
    logbuf ( 10) =  lqfqv2m
    logbuf ( 11) =  losse_fg
    logbuf ( 12) =  lgpsbias
!   logbuf ( 13) =  lobdens 
    noffset  =  12
    DO    k = 1,4
      logbuf ( noffset+k) =  lscadj(k)
    ENDDO
    noffset  =  noffset + 4
    logbuf (noffset+ 1) = lsynop
    logbuf (noffset+ 2) = laircf
    logbuf (noffset+ 3) = lsatob
    logbuf (noffset+ 4) = ldribu
    logbuf (noffset+ 5) = ltemp 
    logbuf (noffset+ 6) = lpilot
    logbuf (noffset+ 7) = lsatem
    logbuf (noffset+ 8) = lgps
    logbuf (noffset+ 9) = lscatt
    logbuf (noffset+10) = lsurfa
    logbuf (noffset+11) = lt2m  
    logbuf (noffset+12) = lrh2m 
    logbuf (noffset+13) = lprecp
    logbuf (noffset+14) = lff10m
    noffset  =  noffset + 14
    logbuf (noffset+ 1) = lcd011
    logbuf (noffset+ 2) = lcd014
    logbuf (noffset+ 3) = lcd021
    logbuf (noffset+ 4) = lcd022
    logbuf (noffset+ 5) = lcd023
    logbuf (noffset+ 6) = lcd024
    logbuf (noffset+ 7) = lcd140
    logbuf (noffset+ 8) = lcd041
    logbuf (noffset+ 9) = lcd141
    logbuf (noffset+10) = lcd241
    logbuf (noffset+11) = lcd144
    logbuf (noffset+12) = lcd244
    logbuf (noffset+13) = lcd088
    logbuf (noffset+14) = lcd188
    logbuf (noffset+15) = lcd063
    logbuf (noffset+16) = lcd064
    logbuf (noffset+17) = lcd165
    logbuf (noffset+18) = lcd035
    logbuf (noffset+19) = lcd036
    logbuf (noffset+20) = lcd037
    logbuf (noffset+21) = lcd135
    logbuf (noffset+22) = lcd039
    logbuf (noffset+23) = lcd040
    logbuf (noffset+24) = lcd032
    logbuf (noffset+25) = lcd033
    logbuf (noffset+26) = lcd038
    logbuf (noffset+27) = lcd132
    logbuf (noffset+28) = lcd133
    logbuf (noffset+29) = lcd136
    logbuf (noffset+30) = lcd137
    logbuf (noffset+31) = lcd086
    logbuf (noffset+32) = lcd186
    logbuf (noffset+33) = lcd122
    logbuf (noffset+34) = lcd123
    logbuf (noffset+35) = lcd096
    logbuf (noffset+36) = lpraof
    logbuf (noffset+37) = lprodr
    logbuf (noffset+38) = ldiasa

!   logbuf for LHN
    noffset  =  noffset + 38
    logbuf (noffset+1) = llhn
    logbuf (noffset+2) = llhnverif
    logbuf (noffset+3) = lhn_search
    logbuf (noffset+4) = lhn_filt
    logbuf (noffset+5) = lhn_relax
    logbuf (noffset+6) = lhn_limit
    logbuf (noffset+7) = lhn_hum_adj
    logbuf (noffset+8) = lhn_spqual
    logbuf (noffset+9) = lhn_incloud
    logbuf (noffset+10) = lhn_diag
    logbuf (noffset+11) = lhn_qrs
    logbuf (noffset+12) = lhn_logscale
    logbuf (noffset+13) = lhn_black
    logbuf (noffset+14) = lhn_wweight
    logbuf (noffset+15) = lhn_bright
    logbuf (noffset+16) = lhn_height
    noffset = noffset + 16
    nnl (1) = noffset
!   nnl (1) = 12 + 4 + 14 + 38 + 16 = 84

!   intbuf for nudging
    intbuf (  1) =  nudgsta
    intbuf (  2) =  nudgend
    intbuf (  3) =  nversta
    intbuf (  4) =  nverend
    intbuf (  5) =  mruntyp
    intbuf (  6) =  mveripr
    intbuf (  7) =  nwtyp
    intbuf (  8) =  ntpscor
    intbuf (  9) =  mpsgcor
    intbuf ( 10) =  khumbal
    intbuf ( 11) =  msprpar
    intbuf ( 12) =  msprpsu
    intbuf ( 13) =  maxmlo
    intbuf ( 14) =  maxsgo
    intbuf ( 15) =  maxuso
    intbuf ( 16) =  maxgpo
    intbuf ( 17) =  maxtvo
    intbuf ( 18) =  maxmlv
    intbuf ( 19) =  mxfrep
    intbuf ( 20) =  mxfobs
    intbuf ( 21) =  nolbc
    intbuf ( 22) =  mqcorr92
    intbuf ( 23) =  itype_obfile
    intbuf ( 24) =  irun_osse
    intbuf ( 25) =  iseed
    intbuf ( 26) =  mcdmsg1
    intbuf ( 27) =  mcdmsg2
    intbuf ( 28) =  mcdno15
    intbuf ( 29) =  mcdno16
    intbuf ( 30) =  mcdno17
    intbuf ( 31) =  mcdno18
    intbuf ( 32) =  ionl
    intbuf ( 33) =  jonl
    intbuf ( 34) =  ionl2
    intbuf ( 35) =  jonl2
    intbuf ( 36) =  noctrq
    noffset = 36
    DO k = 1 , mxgpc
      intbuf (noffset+k) =  igpscen(k)
    ENDDO
    noffset = noffset + mxgpc
    DO k = 1 , mxwt
      intbuf (noffset     +k) =  niwtyp(k)
      intbuf (noffset+mxwt+k) =  kwtyp (k)
    ENDDO
    intbuf (noffset+2*mxwt+1) =  kwtyp (mxwt+1)
    intbuf (noffset+2*mxwt+2) =  kwtyp (mxwt+2)
    noffset = noffset + 2*mxwt + 2
    DO k = 1 , mxtyw
      intbuf (noffset+k) =  iwtyp(k)
    ENDDO
    noffset = noffset + mxtyw
    DO    k = 1 , mxda
      intbuf (noffset       +k) =  mrhyes(k)
      intbuf (noffset+  mxda+k) =  mrhno(k)
      intbuf (noffset+2*mxda+k) =  muvyes(k)
      intbuf (noffset+3*mxda+k) =  muvno(k)
    ENDDO

    noffset = noffset + 4*mxda
!   intbuf for LHN
    intbuf (noffset+1) = nlhn_start
    intbuf (noffset+2) = nlhn_end
    intbuf (noffset+3) = nlhnverif_start
    intbuf (noffset+4) = nlhnverif_end
    intbuf (noffset+5) = rlhn_search
    intbuf (noffset+6) = ktop_lhn
    intbuf (noffset+7) = kbot_lhn
    intbuf (noffset+8) = nlhn_relax
    intbuf (noffset+9) = nradar
    noffset = noffset + 9
    nnl (2) = noffset
!   nnl (2) = 36 + mxgpc + 2*mxwt+2 + mxtyw + 4*mxda + 9
!           = 36 + 20 + 42 + 50 + 80 + 9 = 237

!   realbuf for nudging
    realbuf (  1) = tconbox
    realbuf (  2) = ptpstop
    realbuf (  3) = qgeo
    realbuf (  4) = tipolmx
    realbuf (  5) = tipmxsu
    realbuf (  6) = wtukrsa
    realbuf (  7) = wtukrse
    realbuf (  8) = wtukara
    realbuf (  9) = wtukare
    realbuf ( 10) = wtuksua
    realbuf ( 11) = wtuksue
    realbuf ( 12) = rhfgps
    realbuf ( 13) = rhfpsdd
    realbuf ( 14) = cnondiv
    realbuf ( 15) = fnondiv
    realbuf ( 16) = tnondiv
    realbuf ( 17) = dtqc
    realbuf ( 18) = thairh
    realbuf ( 19) = fperturb
    realbuf ( 20) = obnlat
    realbuf ( 21) = obslat
    realbuf ( 22) = obwlon
    realbuf ( 23) = obelon
    realbuf ( 24) = exnlat
    realbuf ( 25) = exslat
    realbuf ( 26) = exwlon
    realbuf ( 27) = exelon
    realbuf ( 28) = ht2a
    realbuf ( 29) = ht2i
    realbuf ( 30) = hh2a
    realbuf ( 31) = hh2i
    realbuf ( 32) = hffa
    realbuf ( 33) = hffi
    realbuf ( 34) = hprc
    realbuf ( 35) = raintp
    realbuf ( 36) = dinlat
    realbuf ( 37) = dislat
    realbuf ( 38) = diwlon
    realbuf ( 39) = dielon
    realbuf ( 40) = hversta
    realbuf ( 41) = hverend
    realbuf ( 42) = qcfpst
    realbuf ( 43) = qgeotop
    realbuf ( 44) = qgeops
    realbuf ( 45) = gnudggp
    realbuf ( 46) = qcciq
    realbuf ( 47) = qcsiq
    realbuf ( 48) = qcflbcp
    noffset = 48
    DO   k = 1,4
      realbuf (noffset     +k) = gnudg(k)
      realbuf (noffset+   4+k) = gnudgar(k)
      realbuf (noffset+ 2*4+k) = gnudgsu(k)
      realbuf (noffset+ 3*4+k) = gnudgtv(k)
      realbuf (noffset+ 4*4+k) = vcorls(k)
      realbuf (noffset+ 5*4+k) = vcorlsu(k)
      realbuf (noffset+ 6*4+k) = vcutof(k)
      realbuf (noffset+ 7*4+k) = vcutosu(k)
      realbuf (noffset+ 8*4+k) = vpblsu(k)
      realbuf (noffset+ 9*4+k) = wablua(k)
      realbuf (noffset+10*4+k) = wablsu(k)
      realbuf (noffset+11*4+k) = rhinfl(k)
      realbuf (noffset+12*4+k) = rhiflsu(k)
      realbuf (noffset+13*4+k) = rhvfac(k)
      realbuf (noffset+14*4+k) = rhtfac(k)
      realbuf (noffset+15*4+k) = rhtfsu(k)
      realbuf (noffset+16*4+k) = rhfrtv(k)
      realbuf (noffset+17*4+k) = cutofr(k)
      realbuf (noffset+18*4+k) = cutofsu(k)
      realbuf (noffset+19*4+k) = vcsni(k)
      realbuf (noffset+20*4+k) = vcsnisu(k)
      realbuf (noffset+21*4+k) = topobs(k)
      realbuf (noffset+22*4+k) = botmod(k)
      realbuf (noffset+23*4+k) = qcc(k)
      realbuf (noffset+24*4+k) = qccsu(k)
      realbuf (noffset+25*4+k) = qcvf(k)
      realbuf (noffset+26*4+k) = doromx(k)
      realbuf (noffset+27*4+k) = altopsu(k)
    ENDDO

!   realbuf for LHN
    noffset = noffset + 28*4
    realbuf (noffset+1) = lhn_coef
    realbuf (noffset+2) = abs_lhn_lim
    realbuf (noffset+3) = fac_lhn_up
    realbuf (noffset+4) = fac_lhn_down
    realbuf (noffset+5) = rad_wobs_lhn
    realbuf (noffset+6) = thres_lhn
    realbuf (noffset+7) = lhn_dt_obs
    realbuf (noffset+8) = fac_lhn_search
    realbuf (noffset+9) = ktop_temp
    realbuf (noffset+10) = rqrsgmax
    noffset = noffset + 10
    nnl (3) = noffset
!   nnl (3) = 48 + 4*28 + 10 = 170

!   charbuf for nudging and LHN
    charbuf ( 1)(  1:100) = ycdfdir(  1:100)
    charbuf ( 2)(  1:100) = ycdfdir(101:200)
    charbuf ( 3)(  1: 50) = ycdfdir(201:250)
    charbuf ( 4) = yaofpath
    charbuf ( 5) = radar_in
    charbuf ( 6) = blacklist_file
    charbuf ( 7) = height_file
    charbuf ( 8) = yform_lansfc
    charbuf ( 9)(  1:100) = ydir_lansfc(  1:100)
    charbuf (10)(  1:100) = ydir_lansfc(101:200)
    charbuf (11)(  1: 50) = ydir_lansfc(201:250)
    noffset = 11
    ! put noobs_dates in charbuf
    DO nn = 1 , n_noobs
       charbuf (noffset+nn) = noobs_date(nn)
    ENDDO
    noffset = noffset + n_noobs
    nnl (4) = noffset

  ENDIF ! my_world_id == 0

  CALL distribute_values (nnl    ,   4   , 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (logbuf , nnl(1), 0, imp_logical,  icomm_world, ierr)
  CALL distribute_values (intbuf , nnl(2), 0, imp_integers, icomm_world, ierr)
  CALL distribute_values (realbuf, nnl(3), 0, imp_reals,    icomm_world, ierr)
  CALL distribute_values (charbuf, nnl(4), 0, imp_character,icomm_world, ierr)

! CALL distribute_values (logbuf ,  83, 0, imp_logical,  icomm_world, ierr)
! CALL distribute_values (intbuf , 235, 0, imp_integers, icomm_world, ierr)
! CALL distribute_values (realbuf, 168, 0, imp_reals,    icomm_world, ierr)
! CALL distribute_values (charbuf, 5+n_noobs, 0, imp_character,icomm_world,ierr)

  IF (my_world_id /= 0) THEN

!   logbuf for nudging
    lnudge   = logbuf ( 1)
    lverif   = logbuf ( 2)
    lverpas  = logbuf ( 3)
    l1dvar   = logbuf ( 4)
    luvgcor  = logbuf ( 5)
    lsvcorl  = logbuf ( 6)
    ltipol   = logbuf ( 7)
    ltipsu   = logbuf ( 8)
    loiqv2m  = logbuf ( 9)
    lqfqv2m  = logbuf (10)
    losse_fg = logbuf (11)
    lgpsbias = logbuf (12)
!   lobdens  = logbuf (13)
    noffset  = 12
    DO    k = 1,4
      lscadj(k) = logbuf ( noffset+k)
    ENDDO
    noffset  = noffset + 4
    lsynop  = logbuf (noffset+ 1)
    laircf  = logbuf (noffset+ 2)
    lsatob  = logbuf (noffset+ 3)
    ldribu  = logbuf (noffset+ 4)
    ltemp   = logbuf (noffset+ 5)
    lpilot  = logbuf (noffset+ 6)
    lsatem  = logbuf (noffset+ 7)
    lgps    = logbuf (noffset+ 8)
    lscatt  = logbuf (noffset+ 9)
    lsurfa  = logbuf (noffset+10)
    lt2m    = logbuf (noffset+11)
    lrh2m   = logbuf (noffset+12)
    lprecp  = logbuf (noffset+13)
    lff10m  = logbuf (noffset+14)
    noffset  = noffset + 14
    lcd011  = logbuf (noffset+ 1)
    lcd014  = logbuf (noffset+ 2)
    lcd021  = logbuf (noffset+ 3)
    lcd022  = logbuf (noffset+ 4)
    lcd023  = logbuf (noffset+ 5)
    lcd024  = logbuf (noffset+ 6)
    lcd140  = logbuf (noffset+ 7)
    lcd041  = logbuf (noffset+ 8)
    lcd141  = logbuf (noffset+ 9)
    lcd241  = logbuf (noffset+10)
    lcd144  = logbuf (noffset+11)
    lcd244  = logbuf (noffset+12)
    lcd088  = logbuf (noffset+13)
    lcd188  = logbuf (noffset+14)
    lcd063  = logbuf (noffset+15)
    lcd064  = logbuf (noffset+16)
    lcd165  = logbuf (noffset+17)
    lcd035  = logbuf (noffset+18)
    lcd036  = logbuf (noffset+19)
    lcd037  = logbuf (noffset+20)
    lcd135  = logbuf (noffset+21)
    lcd039  = logbuf (noffset+22)
    lcd040  = logbuf (noffset+23)
    lcd032  = logbuf (noffset+24)
    lcd033  = logbuf (noffset+25)
    lcd038  = logbuf (noffset+26)
    lcd132  = logbuf (noffset+27)
    lcd133  = logbuf (noffset+28)
    lcd136  = logbuf (noffset+29)
    lcd137  = logbuf (noffset+30)
    lcd086  = logbuf (noffset+31)
    lcd186  = logbuf (noffset+32)
    lcd122  = logbuf (noffset+33)
    lcd123  = logbuf (noffset+34)
    lcd096  = logbuf (noffset+35)
    lpraof  = logbuf (noffset+36)
    lprodr  = logbuf (noffset+37)
    ldiasa  = logbuf (noffset+38)

!   logbuf for LHN
    noffset      = noffset + 38
    llhn         = logbuf (noffset+1)
    llhnverif    = logbuf (noffset+2)
    lhn_search   = logbuf (noffset+3)
    lhn_filt     = logbuf (noffset+4)
    lhn_relax    = logbuf (noffset+5)
    lhn_limit    = logbuf (noffset+6)
    lhn_hum_adj  = logbuf (noffset+7)
    lhn_spqual   = logbuf (noffset+8)
    lhn_incloud  = logbuf (noffset+9)
    lhn_diag     = logbuf (noffset+10)
    lhn_qrs      = logbuf (noffset+11)
    lhn_logscale = logbuf (noffset+12)
    lhn_black    = logbuf (noffset+13)
    lhn_wweight  = logbuf (noffset+14)
    lhn_bright   = logbuf (noffset+15)
    lhn_height   = logbuf (noffset+16)

!   intbuf for nudging
    nudgsta      = intbuf (  1)
    nudgend      = intbuf (  2)
    nversta      = intbuf (  3)
    nverend      = intbuf (  4)
    mruntyp      = intbuf (  5)
    mveripr      = intbuf (  6)
    nwtyp        = intbuf (  7)
    ntpscor      = intbuf (  8)
    mpsgcor      = intbuf (  9)
    khumbal      = intbuf ( 10)
    msprpar      = intbuf ( 11)
    msprpsu      = intbuf ( 12)
    maxmlo       = intbuf ( 13)
    maxsgo       = intbuf ( 14)
    maxuso       = intbuf ( 15)
    maxgpo       = intbuf ( 16)
    maxtvo       = intbuf ( 17)
    maxmlv       = intbuf ( 18)
    mxfrep       = intbuf ( 19)
    mxfobs       = intbuf ( 20)
    nolbc        = intbuf ( 21)
    mqcorr92     = intbuf ( 22)
    itype_obfile = intbuf ( 23)
    irun_osse    = intbuf ( 24)
    iseed        = intbuf ( 25)
    mcdmsg1      = intbuf ( 26)
    mcdmsg2      = intbuf ( 27)
    mcdno15      = intbuf ( 28)
    mcdno16      = intbuf ( 29)
    mcdno17      = intbuf ( 30)
    mcdno18      = intbuf ( 31)
    ionl         = intbuf ( 32)
    jonl         = intbuf ( 33)
    ionl2        = intbuf ( 34)
    jonl2        = intbuf ( 35)
    noctrq       = intbuf ( 36)
    noffset = 36
    DO k = 1 , mxgpc
      igpscen(k)  =  intbuf (noffset+k)
    ENDDO
    noffset = noffset + mxgpc
    DO k = 1 , mxwt
      niwtyp (k)  =  intbuf (noffset     +k)
      kwtyp  (k)  =  intbuf (noffset+mxwt+k)
    ENDDO
    kwtyp (mxwt+1)  =  intbuf (noffset+2*mxwt+1)
    kwtyp (mxwt+2)  =  intbuf (noffset+2*mxwt+2)
    noffset = noffset + 2*mxwt + 2
    DO k = 1 , mxtyw
      iwtyp  (k)  =  intbuf (noffset+k)
    ENDDO
    noffset = noffset + mxtyw
    DO    k = 1 , mxda
      mrhyes(k)  = intbuf (noffset       +k)
      mrhno(k)   = intbuf (noffset+  mxda+k)
      muvyes(k)  = intbuf (noffset+2*mxda+k)
      muvno(k)   = intbuf (noffset+3*mxda+k)
    ENDDO
 
!   intbuf for LHN
    noffset        = noffset + 4*mxda
    nlhn_start     = intbuf (noffset+1)
    nlhn_end       = intbuf (noffset+2)
    nlhnverif_start= intbuf (noffset+3)
    nlhnverif_end  = intbuf (noffset+4)
    rlhn_search    = intbuf (noffset+5)
    ktop_lhn       = intbuf (noffset+6)
    kbot_lhn       = intbuf (noffset+7)
    nlhn_relax     = intbuf (noffset+8)
    nradar         = intbuf (noffset+9)

!   realbuf for nudging
    tconbox  = realbuf (  1)
    ptpstop  = realbuf (  2)
    qgeo     = realbuf (  3)
    tipolmx  = realbuf (  4)
    tipmxsu  = realbuf (  5)
    wtukrsa  = realbuf (  6)
    wtukrse  = realbuf (  7)
    wtukara  = realbuf (  8)
    wtukare  = realbuf (  9)
    wtuksua  = realbuf ( 10)
    wtuksue  = realbuf ( 11)
    rhfgps   = realbuf ( 12)
    rhfpsdd  = realbuf ( 13)
    cnondiv  = realbuf ( 14)
    fnondiv  = realbuf ( 15)
    tnondiv  = realbuf ( 16)
    dtqc     = realbuf ( 17)
    thairh   = realbuf ( 18)
    fperturb = realbuf ( 19)
    obnlat   = realbuf ( 20)
    obslat   = realbuf ( 21)
    obwlon   = realbuf ( 22)
    obelon   = realbuf ( 23)
    exnlat   = realbuf ( 24)
    exslat   = realbuf ( 25)
    exwlon   = realbuf ( 26)
    exelon   = realbuf ( 27)
    ht2a     = realbuf ( 28)
    ht2i     = realbuf ( 29)
    hh2a     = realbuf ( 30)
    hh2i     = realbuf ( 31)
    hffa     = realbuf ( 32)
    hffi     = realbuf ( 33)
    hprc     = realbuf ( 34)
    raintp   = realbuf ( 35)
    dinlat   = realbuf ( 36)
    dislat   = realbuf ( 37)
    diwlon   = realbuf ( 38)
    dielon   = realbuf ( 39)
    hversta  = realbuf ( 40)
    hverend  = realbuf ( 41)
    qcfpst   = realbuf ( 42)
    qgeotop  = realbuf ( 43)
    qgeops   = realbuf ( 44)
    gnudggp  = realbuf ( 45)
    qcciq    = realbuf ( 46)
    qcsiq    = realbuf ( 47)
    qcflbcp  = realbuf ( 48)
    noffset = 48
    DO   k = 1,4
      gnudg  (k) = realbuf (noffset     +k)
      gnudgar(k) = realbuf (noffset+   4+k)
      gnudgsu(k) = realbuf (noffset+ 2*4+k)
      gnudgtv(k) = realbuf (noffset+ 3*4+k)
      vcorls (k) = realbuf (noffset+ 4*4+k)
      vcorlsu(k) = realbuf (noffset+ 5*4+k)
      vcutof (k) = realbuf (noffset+ 6*4+k)
      vcutosu(k) = realbuf (noffset+ 7*4+k)
      vpblsu (k) = realbuf (noffset+ 8*4+k)
      wablua (k) = realbuf (noffset+ 9*4+k)
      wablsu (k) = realbuf (noffset+10*4+k)
      rhinfl (k) = realbuf (noffset+11*4+k)
      rhiflsu(k) = realbuf (noffset+12*4+k)
      rhvfac (k) = realbuf (noffset+13*4+k)
      rhtfac (k) = realbuf (noffset+14*4+k)
      rhtfsu (k) = realbuf (noffset+15*4+k)
      rhfrtv (k) = realbuf (noffset+16*4+k)
      cutofr (k) = realbuf (noffset+17*4+k)
      cutofsu(k) = realbuf (noffset+18*4+k)
      vcsni  (k) = realbuf (noffset+19*4+k)
      vcsnisu(k) = realbuf (noffset+20*4+k)
      topobs (k) = realbuf (noffset+21*4+k)
      botmod (k) = realbuf (noffset+22*4+k)
      qcc    (k) = realbuf (noffset+23*4+k)
      qccsu  (k) = realbuf (noffset+24*4+k)
      qcvf   (k) = realbuf (noffset+25*4+k)
      doromx (k) = realbuf (noffset+26*4+k)
      altopsu(k) = realbuf (noffset+27*4+k)
    ENDDO

!   realbuf for LHN
    noffset         = noffset + 28*4
    lhn_coef        = realbuf (noffset+1)
    abs_lhn_lim     = realbuf (noffset+2)
    fac_lhn_up      = realbuf (noffset+3)
    fac_lhn_down    = realbuf (noffset+4)
    rad_wobs_lhn    = realbuf (noffset+5)
    thres_lhn       = realbuf (noffset+6)
    lhn_dt_obs      = realbuf (noffset+7)
    fac_lhn_search  = realbuf (noffset+8)
    ktop_temp       = realbuf (noffset+9)
    rqrsgmax        = realbuf (noffset+10)

!   charbuf for nudging and LHN
    ycdfdir(  1:100)  = charbuf ( 1)(1:100)
    ycdfdir(101:200)  = charbuf ( 2)(1:100)
    ycdfdir(201:250)  = charbuf ( 3)(1: 50)
    yaofpath          = charbuf ( 4)
    radar_in          = charbuf ( 5)
    blacklist_file    = charbuf ( 6)
    height_file       = charbuf ( 7)
    yform_lansfc      = charbuf ( 8) 
    ydir_lansfc(  1:100) = charbuf ( 9)(  1:100)
    ydir_lansfc(101:200) = charbuf (10)(  1:100)
    ydir_lansfc(201:250) = charbuf (11)(  1: 50)
    DO nn=1,n_noobs
       noobs_date(nn) = charbuf ( 11+nn)
    ENDDO

  ENDIF ! my_world_id /= 0

ENDIF ! nproc > 0

!-------------------------------------------------------------------------------
!- Section 6: Correction of namelist variables 'luseobs', if required
!-------------------------------------------------------------------------------

luseobs   =  (lnudge) .OR. (lsurfa) .OR. (lverif) .OR. (llhn) .OR. (llhnverif)
lout_anai =  (lnudge)

IF (my_world_id == 0) THEN
  WRITE (nuspecif, '(A)') '  '
  WRITE (nuspecif, '(A,L12)')                                                  &
           '  CHANGED luseobs IN NAMELIST runctl:  luseobs = ',luseobs
  WRITE (nuspecif, '(A)') '  '
ENDIF

! definition of 'isetyp' and 'isetyp0'
  noffset = 0
  DO iiwt = 1 , nwtyp
    isetyp (noffset+1:noffset+niwtyp(iiwt)) = iiwt  
    noffset = noffset + niwtyp(iiwt)
  ENDDO
  DO iiwt = SUM( niwtyp ) , 1 , -1
    IF (iwtyp(iiwt) == 0)  isetyp0 = isetyp(iiwt)
  ENDDO

!-------------------------------------------------------------------------------
!- Section 7: Output of the namelist variables and their default values
!-------------------------------------------------------------------------------

IF (my_world_id == 0) THEN

  WRITE( yrformat, '(''(T8,A,T21,F12.2,T40,F12.2,T60,A)'')' )
  ! format with more positions after decimal point
  WRITE( yr1format, '(''(T8,A,T21,F12.6,T40,F12.6,T60,A)'')' )
  WRITE( yiformat, '(''(T8,A,T21,I12  ,T40,I12  ,T60,A)'')' )
  WRITE( ylformat, '(''(T8,A,T21,L12  ,T40,L12  ,T60,A)'')' )
! WRITE( ygpformat, '(''(T8,A,T21,'',I3.3,''I3)'')' )  mxgpc
! WRITE( ywtformat, '(''(T8,A,T21,'',I3.3,''I3)'')' )  mxwt

  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A24)') '0     NAMELIST:  nudging'
  WRITE (nuspecif, '(A24)') '      ------------------'
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value'      &
                                             , 'Default Value', 'Format'

  WRITE (nuspecif, '(T11,A)') '-->  General variables controlling the nudging'
  WRITE (nuspecif, ylformat) 'lnudge  ',lnudge  ,lnudge_d  ,'L'
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                             'nudgsta ',nudgsta ,nudgsta_d ,'I'                &
                                   ,' : hnudgsta ',hnudgsta
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                             'nudgend ',nudgend ,nudgend_d ,'I'                &
                                   ,' : hnudgend ',hnudgend
  WRITE (nuspecif, ylformat) 'lverif  ',lverif  ,lverif_d  ,'L'
  WRITE (nuspecif, ylformat) 'lverpas ',lverpas ,lverpas_d ,'L'
  WRITE (nuspecif, ylformat) 'l1dvar  ',l1dvar  ,l1dvar_d  ,'L'
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                             'nversta ',nversta ,nversta_d ,'I'                &
                                   ,' : hversta  ',hversta
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                             'nverend ',nverend ,nverend_d ,'I'                &
                                   ,' : hverend  ',hverend
  WRITE (nuspecif, yiformat) 'mruntyp ',mruntyp ,mruntyp_d ,'I'
  WRITE (nuspecif, yiformat) 'mveripr ',mveripr ,mveripr_d ,'I'
  WRITE (nuspecif, yrformat) 'tconbox ',tconbox ,tconbox_d ,'R'
  WRITE (nuspecif, yiformat) 'nwtyp   ',nwtyp   ,nwtyp_d   ,'I'
  DO ivar = 1 , 4
    IF (ivar <= 2)  iwtmax = (MIN( mxwt   , nwtyp  +1 ) - 19 + 19-1) / 19
    IF (ivar >= 3)  iwtmax = (MIN( mxwt+2 , nwtyp+2+1 ) - 19 + 19-1) / 19
    IF (ivar <= 2)  iiwt   = MIN( 19 , mxwt   )
    IF (ivar >= 3)  iiwt   = MIN( 19 , mxwt+2 )
    WRITE( ywtformat, '(''(T8,A,T24,'',I3.3,''I3)'')' )  iiwt
!   WRITE( 0,*) 'ywtformat 1  ', ivar, ywtformat
    IF (ivar == 1) WRITE(nuspecif,ywtformat) 'niwtyp (actual) ',niwtyp  (1:iiwt)
    IF (ivar == 2) WRITE(nuspecif,ywtformat) 'niwtyp (default)',niwtyp_d(1:iiwt)
    IF (ivar == 3) WRITE(nuspecif,ywtformat) 'kwtyp  (actual) ',kwtyp   (1:iiwt)
    IF (ivar == 4) WRITE(nuspecif,ywtformat) 'kwtyp  (default)',kwtyp_d (1:iiwt)
    DO iiw2 = 1 , iwtmax
      iwtmin = 19 + (iiw2-1)*19
      IF (ivar <= 2)  iiwt   = MIN( mxwt   - iwtmin , 19 )
      IF (ivar >= 3)  iiwt   = MIN( mxwt+2 - iwtmin , 19 )
      WRITE( ywtformat, '(''(     T24,'',I3.3,''I3)'')' )  iiwt
!     WRITE( 0,*) 'ywtformat 2  ', ivar, iiw2, ywtformat
      IF (ivar == 1) WRITE (nuspecif, ywtformat)  niwtyp  (iwtmin+1:iwtmin+iiwt)
      IF (ivar == 2) WRITE (nuspecif, ywtformat)  niwtyp_d(iwtmin+1:iwtmin+iiwt)
      IF (ivar == 3) WRITE (nuspecif, ywtformat)  kwtyp   (iwtmin+1:iwtmin+iiwt)
      IF (ivar == 4) WRITE (nuspecif, ywtformat)  kwtyp_d (iwtmin+1:iwtmin+iiwt)
    ENDDO
  ENDDO
  iwtmax = (MIN( mxtyw, SUM( niwtyp )+1 ) - 8 + 16-1) / 16
  DO ivar = 1 , 2
    iiwt = MIN( 8 , mxtyw )
    WRITE( yotformat, '(''(T8,A,T28,I2,A,T49,'',I3.3,''I4)'')' )  iiwt
!   WRITE( 0,*) 'yotformat 1  ', ivar, yotformat
    IF (ivar == 1) WRITE (nuspecif, yotformat) 'iwtyp (array length=', mxtyw,  &
                                            '), actual values: ',iwtyp  (1:iiwt)
    IF (ivar == 2) WRITE (nuspecif, yotformat) 'iwtyp (array length=', mxtyw,  &
                                            '), default values:',iwtyp_d(1:iiwt)
    DO iiw2 = 1 , iwtmax
      iwtmin = 8 + (iiw2-1)*16
      iiwt   = MIN( mxtyw-iwtmin , 16 )
      WRITE( yotformat, '(''(              T17,'',I3.3,''I4)'')' )  iiwt
!     WRITE( 0,*) 'yotformat 2  ', ivar, iiw2, yotformat
      IF (ivar == 1) WRITE (nuspecif, yotformat)  iwtyp  (iwtmin+1:iwtmin+iiwt)
      IF (ivar == 2) WRITE (nuspecif, yotformat)  iwtyp_d(iwtmin+1:iwtmin+iiwt)
    ENDDO
  ENDDO

  WRITE (nuspecif, '(T11,2A)') '-->  Corrections to balance the analysis '     &
                              ,'increments'
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value'      &
                                             , 'Default Value', 'Format'
  WRITE (nuspecif, yiformat) 'khumbal ',khumbal ,khumbal_d ,'I'
  WRITE (nuspecif, yiformat) 'ntpscor ',ntpscor ,ntpscor_d ,'I'
  WRITE (nuspecif, yrformat) 'ptpstop ',ptpstop ,ptpstop_d ,'R'
  WRITE (nuspecif, ylformat) 'luvgcor ',luvgcor ,luvgcor_d ,'L'
  WRITE (nuspecif, yiformat) 'mpsgcor ',mpsgcor ,mpsgcor_d ,'I'
  WRITE (nuspecif, yrformat) 'qgeo    ',qgeo    ,qgeo_d    ,'R'
  WRITE (nuspecif, yrformat) 'qgeotop ',qgeotop ,qgeotop_d ,'R'
  WRITE (nuspecif, yrformat) 'qgeops  ',qgeops  ,qgeops_d  ,'R'

  WRITE (nuspecif, '(T11,A)') '-->  Nudging coefficients'
  WRITE (nuspecif, '(T7,A,T18,A,T50,A)') 'Variable'                            &
                                      , ' <--   Actual Values (R)   --> '      &
                                      , ' <--   Default Values (R)   -->'
  WRITE (nuspecif, '(T8,A,T16,4F8.5,T49,4F8.5)') 'gnudg   ',gnudg   ,gnudg_d
  WRITE (nuspecif, '(T8,A,T16,4F8.5,T49,4F8.5)') 'gnudgar ',gnudgar ,gnudgar_d
  WRITE (nuspecif, '(T8,A,T16,4F8.5,T49,4F8.5)') 'gnudgsu ',gnudgsu ,gnudgsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.5,T49,4F8.5)') 'gnudgtv ',gnudgtv ,gnudgtv_d
  WRITE (nuspecif, '(T8,A,T40, F8.5,T73, F8.5)') 'gnudggp ',gnudggp ,gnudggp_d

  WRITE (nuspecif, '(T11,A)') '-->  Temporal weights'
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value'      &
                                             , 'Default Value', 'Format'
  WRITE (nuspecif, ylformat) 'ltipol  ',ltipol  ,ltipol_d  ,'L'
  WRITE (nuspecif, ylformat) 'ltipsu  ',ltipsu  ,ltipsu_d  ,'L'
  WRITE (nuspecif, yrformat) 'tipolmx ',tipolmx ,tipolmx_d ,'R'
  WRITE (nuspecif, yrformat) 'tipmxsu ',tipmxsu ,tipmxsu_d ,'R'
  WRITE (nuspecif, yrformat) 'wtukrsa ',wtukrsa ,wtukrsa_d ,'R'
  WRITE (nuspecif, yrformat) 'wtukrse ',wtukrse ,wtukrse_d ,'R'
  WRITE (nuspecif, yrformat) 'wtukara ',wtukara ,wtukara_d ,'R'
  WRITE (nuspecif, yrformat) 'wtukare ',wtukare ,wtukare_d ,'R'
  WRITE (nuspecif, yrformat) 'wtuksua ',wtuksua ,wtuksua_d ,'R'
  WRITE (nuspecif, yrformat) 'wtuksue ',wtuksue ,wtuksue_d ,'R'

  WRITE (nuspecif, '(T11,A)') '-->  Spatial weights'
  WRITE (nuspecif, yiformat) 'msprpar ',msprpar ,msprpar_d ,'I'
  WRITE (nuspecif, yiformat) 'msprpsu ',msprpsu ,msprpsu_d ,'I'
  WRITE (nuspecif, yrformat) 'rhfgps  ',rhfgps  ,rhfgps_d  ,'R'
  WRITE (nuspecif, yrformat) 'rhfpsdd ',rhfpsdd ,rhfpsdd_d ,'R'
  WRITE (nuspecif, yrformat) 'fnondiv ',fnondiv ,fnondiv_d ,'R'
  WRITE (nuspecif, yrformat) 'cnondiv ',cnondiv ,cnondiv_d ,'R'
  WRITE (nuspecif, yrformat) 'tnondiv ',tnondiv ,tnondiv_d ,'R'
  WRITE (nuspecif, ylformat) 'lsvcorl ',lsvcorl ,lsvcorl_d ,'L'
  WRITE (nuspecif, '(T7,A,T18,A,T50,A)') 'Variable'                            &
                                      , ' <--   Actual Values (R)   --> '      &
                                      , ' <--   Default Values (R)   -->'
  WRITE (nuspecif, '(T8,A,T16,4F8.4,T49,4F8.4)') 'vcorls  ',vcorls  ,vcorls_d
  WRITE (nuspecif, '(T8,A,T16,4F8.4,T49,4F8.4)') 'vcorlsu ',vcorlsu ,vcorlsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.4,T49,4F8.4)') 'vcutof  ',vcutof  ,vcutof_d
  WRITE (nuspecif, '(T8,A,T16,4F8.4,T49,4F8.4)') 'vcutosu ',vcutosu ,vcutosu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.4,T49,4F8.4)') 'vpblsu  ',vpblsu  ,vpblsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'wablua  ',wablua  ,wablua_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'wablsu  ',wablsu  ,wablsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhvfac  ',rhvfac  ,rhvfac_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhinfl  ',rhinfl  ,rhinfl_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhiflsu ',rhiflsu ,rhiflsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhtfac  ',rhtfac  ,rhtfac_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhtfsu  ',rhtfsu  ,rhtfsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'rhfrtv  ',rhfrtv  ,rhfrtv_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'cutofr  ',cutofr  ,cutofr_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'cutofsu ',cutofsu ,cutofsu_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'vcsni   ',vcsni   ,vcsni_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'vcsnisu ',vcsnisu ,vcsnisu_d

  WRITE (nuspecif, '(T11,A)') '-->  Computation of observation increments'
  WRITE (nuspecif, '(T8,A,T16,4L8  ,T49,4L8  )') 'lscadj  ',lscadj  ,lscadj_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'topobs  ',topobs  ,topobs_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'botmod  ',botmod  ,botmod_d
  WRITE (nuspecif, '(T8,A,T16, L8  ,T49, L8  )') 'loiqv2m ',loiqv2m ,loiqv2m_d
  WRITE (nuspecif, '(T8,A,T16, L8  ,T49, L8  )') 'lqfqv2m ',lqfqv2m ,lqfqv2m_d

  WRITE (nuspecif, '(T11,A)') '-->  Quality control and quality weights'
  WRITE (nuspecif, '(T8,A,T16, F8.2,T49, F8.2)') 'dtqc    ',dtqc    ,dtqc_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'qcvf    ',qcvf    ,qcvf_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'qcc     ',qcc     ,qcc_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'qccsu   ',qccsu   ,qccsu_d
  WRITE (nuspecif, '(T8,A,T16, F8.2,T49, F8.2)') 'qcciq   ',qcciq   ,qcciq_d
  WRITE (nuspecif, '(T8,A,T16, F8.2,T49, F8.2)') 'qcsiq   ',qcsiq   ,qcsiq_d
  WRITE (nuspecif, '(T8,A,T16, F8.2,T49, F8.2)') 'qcflbcp ',qcflbcp ,qcflbcp_d
  WRITE (nuspecif, '(T8,A,T16, F8.2,T49, F8.2)') 'qcfpst  ',qcfpst  ,qcfpst_d

  WRITE (nuspecif, '(T11,A)') '-->  Observation processing'
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'doromx  ',doromx  ,doromx_d
  WRITE (nuspecif, '(T8,A,T16,4F8.2,T49,4F8.2)') 'altopsu ',altopsu ,altopsu_d
  WRITE (nuspecif, '(T8,A,T21,A,T70,A)')  'ycdfdir  ', TRIM(ycdfdir) , 'C*250'
  WRITE (nuspecif, '(T8,A,T21,A,T70,A)')  'yaofpath ', yaofpath      , 'C*70'
! WRITE (nuspecif, '(T8,A,T25,20I3)') 'igpscen (actual) ',igpscen 
! WRITE (nuspecif, '(T8,A,T25,20I3)') 'igpscen (default)',igpscen_d                 
! WRITE (nuspecif, ygpformat) 'igpscen (actual) ',igpscen 
! WRITE (nuspecif, ygpformat) 'igpscen (default)',igpscen_d                 
  DO ivar = 1 , 2
    iwtmax = (MIN( mxgpc , nactgpc +3 ) - 19 + 19-1) / 19
    iiwt   = MIN( 19 , mxgpc  )
    WRITE( ygpformat, '(''(T8,A,T24,'',I3.3,''I3)'')' )  iiwt
!   WRITE( 0,*) 'ygpformat ', ygpformat
    IF (ivar==1) WRITE(nuspecif,ygpformat) 'igpscen (actual) ',igpscen  (1:iiwt)
    IF (ivar==2) WRITE(nuspecif,ygpformat) 'igpscen (default)',igpscen_d(1:iiwt)
    DO iiw2 = 1 , iwtmax
      iwtmin = 19 + (iiw2-1)*19
      iiwt   = MIN( mxgpc  - iwtmin , 19 )
      WRITE( ygpformat, '(''(     T24,'',I3.3,''I3)'')' )  iiwt
!     WRITE( 0,*) 'ygpformat ', ygpformat
      IF (ivar == 1) WRITE (nuspecif, ygpformat) igpscen  (iwtmin+1:iwtmin+iiwt)
      IF (ivar == 2) WRITE (nuspecif, ygpformat) igpscen_d(iwtmin+1:iwtmin+iiwt)
    ENDDO
  ENDDO
  WRITE (nuspecif, '(T7,A,T21,A,T39,A,T58,A)') 'Variable', 'Actual Value'      &
                                             , 'Default Value', 'Format'
  WRITE (nuspecif, yiformat) 'itype_obfile',itype_obfile,itype_obfile_d,'I'
  WRITE (nuspecif, yiformat) 'irun_osse'   ,irun_osse   ,irun_osse_d   ,'I'
  WRITE (nuspecif, ylformat) 'losse_fg',losse_fg,losse_fg_d,'L'
  WRITE (nuspecif, ylformat) 'lgpsbias',lgpsbias,lgpsbias_d,'L'
  WRITE (nuspecif, yiformat) 'mqcorr92',mqcorr92,mqcorr92_d,'I'
  WRITE (nuspecif, yiformat) 'iseed   ',iseed   ,iseed_d   ,'I'
  WRITE (nuspecif, yrformat) 'fperturb',fperturb,fperturb_d,'R'
  WRITE (nuspecif, yrformat) 'thairh  ',thairh  ,thairh_d  ,'R'
  WRITE (nuspecif, yrformat) 'obnlat  ',obnlat  ,obnlat_d  ,'R'
  WRITE (nuspecif, yrformat) 'obslat  ',obslat  ,obslat_d  ,'R'
  WRITE (nuspecif, yrformat) 'obwlon  ',obwlon  ,obwlon_d  ,'R'
  WRITE (nuspecif, yrformat) 'obelon  ',obelon  ,obelon_d  ,'R'
  WRITE (nuspecif, yrformat) 'exnlat  ',exnlat  ,exnlat_d  ,'R'
  WRITE (nuspecif, yrformat) 'exslat  ',exslat  ,exslat_d  ,'R'
  WRITE (nuspecif, yrformat) 'exwlon  ',exwlon  ,exwlon_d  ,'R'
  WRITE (nuspecif, yrformat) 'exelon  ',exelon  ,exelon_d  ,'R'
  WRITE (nuspecif, yiformat) 'maxmlo  ',maxmlo  ,maxmlo_d  ,'I'
  WRITE (nuspecif, yiformat) 'maxsgo  ',maxsgo  ,maxsgo_d  ,'I'
  WRITE (nuspecif, yiformat) 'maxuso  ',maxuso  ,maxuso_d  ,'I'
  WRITE (nuspecif, yiformat) 'maxgpo  ',maxgpo  ,maxgpo_d  ,'I'
  WRITE (nuspecif, yiformat) 'maxtvo  ',maxtvo  ,maxtvo_d  ,'I'
  WRITE (nuspecif, yiformat) 'maxmlv  ',maxmlv  ,maxmlv_d  ,'I'
  WRITE (nuspecif, yiformat) 'mxfrep  ',mxfrep  ,mxfrep_d  ,'I'
  WRITE (nuspecif, yiformat) 'mxfobs  ',mxfobs  ,mxfobs_d  ,'I'
  WRITE (nuspecif, yiformat) 'nolbc   ',nolbc   ,nolbc_d   ,'I'
  WRITE (nuspecif, ylformat) 'lsynop  ',lsynop  ,lsynop_d  ,'L'
  WRITE (nuspecif, ylformat) 'laircf  ',laircf  ,laircf_d  ,'L'
  WRITE (nuspecif, ylformat) 'lsatob  ',lsatob  ,lsatob_d  ,'L'
  WRITE (nuspecif, ylformat) 'ldribu  ',ldribu  ,ldribu_d  ,'L'
  WRITE (nuspecif, ylformat) 'ltemp   ',ltemp   ,ltemp_d   ,'L'
  WRITE (nuspecif, ylformat) 'lpilot  ',lpilot  ,lpilot_d  ,'L'
  WRITE (nuspecif, ylformat) 'lsatem  ',lsatem  ,lsatem_d  ,'L'
  WRITE (nuspecif, ylformat) 'lgps    ',lgps    ,lgps_d    ,'L'
  WRITE (nuspecif, ylformat) 'lscatt  ',lscatt  ,lscatt_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd011  ',lcd011  ,lcd011_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd014  ',lcd014  ,lcd014_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd021  ',lcd021  ,lcd021_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd022  ',lcd022  ,lcd022_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd023  ',lcd023  ,lcd023_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd024  ',lcd024  ,lcd024_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd140  ',lcd140  ,lcd140_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd041  ',lcd041  ,lcd041_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd141  ',lcd141  ,lcd141_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd241  ',lcd241  ,lcd241_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd144  ',lcd144  ,lcd144_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd244  ',lcd244  ,lcd244_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd088  ',lcd088  ,lcd088_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd188  ',lcd188  ,lcd188_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd063  ',lcd063  ,lcd063_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd064  ',lcd064  ,lcd064_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd165  ',lcd165  ,lcd165_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd035  ',lcd035  ,lcd035_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd036  ',lcd036  ,lcd036_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd037  ',lcd037  ,lcd037_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd135  ',lcd135  ,lcd135_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd039  ',lcd039  ,lcd039_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd040  ',lcd040  ,lcd040_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd032  ',lcd032  ,lcd032_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd033  ',lcd033  ,lcd033_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd038  ',lcd038  ,lcd038_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd132  ',lcd132  ,lcd132_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd133  ',lcd133  ,lcd133_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd136  ',lcd136  ,lcd136_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd137  ',lcd137  ,lcd137_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd086  ',lcd086  ,lcd086_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd186  ',lcd186  ,lcd186_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd122  ',lcd122  ,lcd122_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd123  ',lcd123  ,lcd123_d  ,'L'
  WRITE (nuspecif, ylformat) 'lcd096  ',lcd096  ,lcd096_d  ,'L'
  WRITE (nuspecif, yiformat) 'mcdmsg1 ',mcdmsg1 ,mcdmsg1_d ,'I'
  WRITE (nuspecif, yiformat) 'mcdmsg2 ',mcdmsg2 ,mcdmsg2_d ,'I'
  WRITE (nuspecif, yiformat) 'mcdno15 ',mcdno15 ,mcdno15_d ,'I'
  WRITE (nuspecif, yiformat) 'mcdno16 ',mcdno16 ,mcdno16_d ,'I'
  WRITE (nuspecif, yiformat) 'mcdno17 ',mcdno17 ,mcdno17_d ,'I'
  WRITE (nuspecif, yiformat) 'mcdno18 ',mcdno18 ,mcdno18_d ,'I'

  WRITE (nuspecif, '(T11,A)') '-->  2-D analyses'
  WRITE (nuspecif, ylformat) 'lsurfa  ',lsurfa  ,lsurfa_d  ,'L'
  WRITE (nuspecif, ylformat) 'lt2m    ',lt2m    ,lt2m_d    ,'L'
  WRITE (nuspecif, ylformat) 'lrh2m   ',lrh2m   ,lrh2m_d   ,'L'
  WRITE (nuspecif, ylformat) 'lprecp  ',lprecp  ,lprecp_d  ,'L'
  WRITE (nuspecif, ylformat) 'lff10m  ',lff10m  ,lff10m_d  ,'L'
  WRITE (nuspecif, yrformat) 'ht2a    ',ht2a    ,ht2a_d    ,'R'
  WRITE (nuspecif, yrformat) 'ht2i    ',ht2i    ,ht2i_d    ,'R'
  WRITE (nuspecif, yrformat) 'hh2a    ',hh2a    ,hh2a_d    ,'R'
  WRITE (nuspecif, yrformat) 'hh2i    ',hh2i    ,hh2i_d    ,'R'
  WRITE (nuspecif, yrformat) 'hffa    ',hffa    ,hffa_d    ,'R'
  WRITE (nuspecif, yrformat) 'hffi    ',hffi    ,hffi_d    ,'R'
  WRITE (nuspecif, yrformat) 'hprc    ',hprc    ,hprc_d    ,'R'
  WRITE (nuspecif, yrformat) 'raintp  ',raintp  ,raintp_d  ,'R'
  WRITE (nuspecif, '(T8,A,T31,A,T70,A)') 'ydir_lansfc ', TRIM(ydir_lansfc), 'C*250'
  WRITE (nuspecif, '(T8,A,T16,A,T49,A)') 'yform_lansfc', yform_lansfc_d, yform_lansfc

  WRITE (nuspecif, '(T11,A)') '-->  Diagnostic output'
  WRITE (nuspecif, ylformat) 'lpraof  ',lpraof  ,lpraof_d  ,'L'
  WRITE (nuspecif, ylformat) 'lprodr  ',lprodr  ,lprodr_d  ,'L'
  WRITE (nuspecif, ylformat) 'ldiasa  ',ldiasa  ,ldiasa_d  ,'L'
  WRITE (nuspecif, yiformat) 'ionl    ',ionl    ,ionl_d    ,'I'
  WRITE (nuspecif, yiformat) 'jonl    ',jonl    ,jonl_d    ,'I'
  WRITE (nuspecif, yiformat) 'ionl2   ',ionl2   ,ionl2_d   ,'I'
  WRITE (nuspecif, yiformat) 'jonl2   ',jonl2   ,jonl2_d   ,'I'
  WRITE (nuspecif, yiformat) 'noctrq  ',noctrq  ,noctrq_d  ,'I'
  WRITE (nuspecif, yrformat) 'dinlat  ',dinlat  ,dinlat_d  ,'R'
  WRITE (nuspecif, yrformat) 'dislat  ',dislat  ,dislat_d  ,'R'
  WRITE (nuspecif, yrformat) 'diwlon  ',diwlon  ,diwlon_d  ,'R'
  WRITE (nuspecif, yrformat) 'dielon  ',dielon  ,dielon_d  ,'R'
!LHN
  WRITE (nuspecif, '(T11,A)') '-->  Latent heat nudging'
  WRITE (nuspecif, ylformat) 'llhn           ',llhn,llhn_d                  ,'L'
  WRITE (nuspecif, ylformat) 'llhnverif      ',llhnverif,llhnverif_d        ,'L'
  WRITE (nuspecif, ylformat) 'lhn_diag       ',lhn_diag,lhn_diag_d          ,'L'
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                    'nlhn_start ',nlhn_start ,nlhn_start_d ,'I'                &
                                   ,' : hlhn_start ',hlhn_start
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                          'nlhn_end ',nlhn_end ,nlhn_end_d ,'I'                &
                                       ,' : hlhn_end ',hlhn_end
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                    'nlhnverif_start ',nlhnverif_start ,nlhnverif_start_d ,'I' &
                                   ,' : hlhnverif_start ',hlhnverif_start
  WRITE (nuspecif, '(T8,A,T21,I12,T40,I12,T60,A , A ,F8.4)')                   &
                          'nlhnverif_end ',nlhnverif_end ,nlhnverif_end_d ,'I' &
                                       ,' : hlhnverif_end ',hlhnverif_end
  WRITE (nuspecif, yrformat) 'lhn_dt_obs     ',lhn_dt_obs,lhn_dt_obs_d      ,'R'
  WRITE (nuspecif, yrformat) 'lhn_coef       ',lhn_coef,lhn_coef_d          ,'R'
  WRITE (nuspecif, yiformat) 'nradar         ',nradar,nradar_d              ,'I'
  WRITE (nuspecif, yiformat) 'ktop_lhn       ',ktop_lhn,ktop_lhn_d          ,'I'
  WRITE (nuspecif, yiformat) 'kbot_lhn       ',kbot_lhn,kbot_lhn_d          ,'I'
  WRITE (nuspecif, ylformat) 'lhn_search     ',lhn_search,lhn_search_d      ,'L'
  WRITE (nuspecif, ylformat) 'lhn_hum_adj    ',lhn_hum_adj,lhn_hum_adj_d    ,'L'
  WRITE (nuspecif, yiformat) 'rlhn_search     ',rlhn_search,rlhn_search_d   ,'I'
  WRITE (nuspecif, yr1format) 'fac_lhn_search',fac_lhn_search_in,fac_lhn_search_d,'R'
  IF (lhn_logscale) &
     WRITE (nuspecif, yr1format) ' -> logscaled: ',fac_lhn_search
  WRITE (nuspecif, yrformat) 'rad_wobs_lhn   ',rad_wobs_lhn,rad_wobs_lhn_d  ,'R'
  WRITE (nuspecif, yr1format) 'thres_lhn     ',thres_lhn,thres_lhn_d        ,'R'
  WRITE (nuspecif, yr1format) 'rqrsgmax      ',rqrsgmax,rqrsgmax_d          ,'R'
  WRITE (nuspecif, yr1format) 'ktop_temp     ',ktop_temp,ktop_temp_d        ,'R'
  WRITE (nuspecif, yr1format) 'fac_lhn_up    ',fac_lhn_up_in,fac_lhn_up_d      ,'R'
  IF (lhn_logscale) &
     WRITE (nuspecif, yr1format) '  -> logscaled: ',fac_lhn_up
  WRITE (nuspecif, yr1format) 'fac_lhn_down  ',fac_lhn_down_in,fac_lhn_down_d  ,'R'
  IF (lhn_logscale) &
     WRITE (nuspecif, yr1format) '  -> logscaled: ',fac_lhn_down
  WRITE (nuspecif, ylformat) 'lhn_filt       ',lhn_filt,lhn_filt_d          ,'L'
  WRITE (nuspecif, ylformat) 'lhn_relax      ',lhn_relax,lhn_relax_d        ,'L'
  WRITE (nuspecif, yiformat) 'nlhn_relax     ',nlhn_relax,nlhn_relax_d      ,'I'
  WRITE (nuspecif, ylformat) 'lhn_limit      ',lhn_limit,lhn_limit_d        ,'L'
  WRITE (nuspecif, yr1format) 'abs_lhn_lim   ',abs_lhn_lim,abs_lhn_lim_d    ,'R'
  WRITE (nuspecif, ylformat) 'lhn_spqual     ',lhn_spqual,lhn_spqual_d      ,'L'
  WRITE (nuspecif, ylformat) 'lhn_incloud    ',lhn_incloud,lhn_incloud_d    ,'L'
  WRITE (nuspecif, ylformat) 'lhn_black      ',lhn_black   ,lhn_black_d     ,'L'
  WRITE (nuspecif, ylformat) 'lhn_qrs       ',lhn_qrs     ,lhn_qrs_d        ,'L'
  WRITE (nuspecif, ylformat) 'lhn_logscale  ',lhn_logscale,lhn_logscale_d   ,'L'
  WRITE (nuspecif, ylformat) 'lhn_wweight   ',lhn_wweight ,lhn_wweight_d    ,'L'
  WRITE (nuspecif, ylformat) 'lhn_height   ',lhn_height ,lhn_height_d       ,'L'
  WRITE (nuspecif, ylformat) 'lhn_bright   ',lhn_bright ,lhn_bright_d       ,'L'
  WRITE (nuspecif, '(T8,A,T31,A,T70,A)') 'radar_in '     ,radar_in      ,'C*100'
  WRITE (nuspecif, '(T8,A,T31,A,T70,A)') 'blacklist_file',blacklist_file,'C*100'
  WRITE (nuspecif, '(T8,A,T31,A,T70,A)') 'height_file',height_file_d,'C*100'
  
  DO nn=1,n_noobs
     IF (noobs_date(nn) /= '') &
         WRITE (nuspecif, '(T8,A,T20,I4,T31,A,T70,A)') 'noobs_date:',nn,noobs_date(nn),'C*12'
  ENDDO

  WRITE (nuspecif, '(A2)')  '  '

ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_nudging
!==============================================================================

!------------------------------------------------------------------------------
! End of module procedure organize_assimilation
!------------------------------------------------------------------------------

END SUBROUTINE organize_assimilation
