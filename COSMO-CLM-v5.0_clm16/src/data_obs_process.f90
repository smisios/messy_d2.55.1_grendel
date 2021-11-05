!+ Data module for all data, that are used by the data assimilation
!-------------------------------------------------------------------------------

MODULE data_obs_process

!-------------------------------------------------------------------------------
! Description:
!  This module declares and initializes all parametric data and data arrays
!  that are required to read and pre-process observation reports from AOF files
!  (except for those variables that reside in modules data_obs_cdfin or
!   data_obs_lib_cosmo).
!  It contains the following sections: 
!    
!    1. Analysis Observation File (AOF) formats
!    2. Numerical and scaling constants related to AOF 
!    3. Variable extraction and output numbers related to AOF
!    4. Diagnostic arrays and pointers for event counters (for AOF input only)
!    5. Character descriptions of events, flags, observation & code types
!    6. Temporary information buffers and report counters 
!  Moved to module 'data_obs_cdfin.f90':
!    -  Data event counter arrays and diagnostics arrays and formats
!    -  Observation errors
!    -  Limit values
!    -  Further parameters controlling the observation processing
!    -  Variables used for the production of aircraft multi-level reports
!    -  Output buffer and formats
!    -  Global fields
!
! Current Code Owner:  DWD, Christoph Schraff
!  phone:  +49  69  8062 2725
!  fax:    +49  69  8062 3721
!  email:  Christoph.Schraff@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.10       1998/09/29 Michael Buchhold
!  Initial release
! 1.12       1998/10/19 Ulrich Schaettler
!  Changes to introduce nudging in the next release.
! 1.15       1998/11/02 Christoph Schraff
!  All global allocatable arrays moved from module 'src_obs_processing'.
! 1.19       1998/12/11 Christoph Schraff
!  Report and data event counters and descriptions revised. 'nodrln' removed.
!  Character literals used as a Hollerith constant replaced by integers.
! 1.27       1999/03/29 Christoph Schraff
!  Processing of optional groups activated, extreme temperature group
!  introduced, surface analysis limits moved to 'data_obs_record'.
! 1.31       1999/07/01 Christoph Schraff
!  Reduction and redefinition of limit 'rprlim'. 'nmlob', 'nsgob' moved from
!  'src_obs_processing' to section 10.
! 1.36       2000/02/24 Michael Buchhold
!     1. Introduction of ACAR aircraft reports
!     2. Adaptation to 32 bit aof
! 2.5        2001/06/01 Christoph Schraff
!  Thresholds and events for wind shear, lapse rate, and flight track checks.
! 2.6        2001/06/12 Christoph Schraff
!  Variable for message on thinning sequences of aircraft reports.
! 2.13       2002/01/18 Christoph Schraff
!  To allow for the 'hot' compiler option on IBM: Data statements replaced by
!  direct assignment, and parameter attribute added to corresponding variables.
! 2.19       2002/09/27 Michael Buchhold
!  Introduction of new code type for Radar VAD wind profiles
! 3.3        2003/04/22 Maria Tomassini + Christoph Schraff
!  Index for GPS reports introduced.
! 3.6        2003/12/11 Christoph Schraff
!  SATOB only as single-level reports. Correction at data events descriptions.
! 3.12       2004/09/15 Christoph Schraff
!  Bug correction to render RESHAPE (for noct??) consistent with array size.
!  Extension of code type descriptions to include satellite retrievals.
! 3.18       2006/03/03 Christoph Schraff
!  Code type descriptors for statistics updated by MSG and NOAA entries.
! V3_23        2007/03/30 Michael Buchhold
!  Introduction of amdar humidity data
! V4_5         2008/09/10 Christoph Schraff
!  Added section 1 for reading observations from NetCDF files and section 2 for
!  evaluating a black-/whitelist. Revised limits to use significant-level data.
! V4_7         2008/12/12 Christoph Schraff
!  Bug correction for control output (nfmt24, nfmt25 introduced).
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Ulrich Schaettler
!  Changed the code owner
! V4_22        2012/01/31 Christoph Schraff
!  All information that is used for reading obs reports from NetCDF files is
!  moved to new module 'data_obs_cdfin.f90', e.g. event counters, observation
!  errors, some limits values and parameters contolling the obs processing,
!  variables used for the production of aircraft multi-level reports, etc.
!  'nbufhed' increased to 18 to accommodate obs date 'iobdat'.
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!===============================================================================
!
! Modules used: 

!-------------------------------------------------------------------------------

USE netcdf
 
USE data_parameters, ONLY :   &
    ireals,    & ! KIND-type parameter for real variables
    iintegers    ! KIND-type parameter for standard integer variables

!-------------------------------------------------------------------------------

IMPLICIT NONE

!===============================================================================

! Local arrays and scalars:
! -------------------------

! INTEGER (KIND=iintegers), PRIVATE  ::  &
!   i,j              ! loop indices
!-------------------------------------------------------------------------------

! Global (i.e. public) Declarations:
! ----------------------------------

!-------------------------------------------------------------------------------
! Section 1 : Analysis Observation File formats
!             (this section is used ONLY IF obs data are read from AOF file)
!-------------------------------------------------------------------------------
!
!         1.1   file description record format
!               ------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nfddrl = 14  ,& ! fdr record length
    nlrecf =  1  ,& ! fdr record length
    nlprlf =  2  ,& ! length of preliminary array
    nnextf =  3  ,& ! length of next fdr
    nlmaxf =  4  ,& ! maximum record size
    nddr1  =  5  ,& ! record number of first ddr
    ncdatf =  6  ,& ! creation date
    nctimf =  7  ,& ! creation time
    nfltyp =  8  ,& ! file type
    ntdatf =  9  ,& ! number of files
    nsecf  = 10  ,& ! common time unit length (secs)
    nlchrf = 11  ,& ! character description
    nlusrf = 12  ,& ! user area
    nrecf  = 13  ,& ! record number
    nchkwf = 14     ! check word


!         1.1.1  observation file parameters
!                ---------------------------
  INTEGER (KIND=iintegers)  ::    &
    nmxlob =  0  ,& ! maximum length of AOF data record
    nmxbln          ! maxi. length of buffer containing all reports for one PE

  LOGICAL                   ::    &
    lastob =.FALSE. ! .TRUE. if all observations have been processed


!         1.2   data description record format
!               ------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nddrrl = 17  ,& ! ddr record length
    nlrecd =  1  ,& ! length of ddr (data description record)
    nlprld =  2  ,& ! length of preliminary array
    nnextd =  3  ,& ! length of next record
    nlmaxd =  4  ,& ! max data record length
    ndr1   =  5  ,& ! record number of first ddr
    ncdatd =  6  ,& ! creation date
    nctimd =  7  ,& ! creation time
    ndttyp =  8  ,& ! data type
    norec  =  9  ,& ! number of data records
    nsecd  = 10  ,& ! common time unit length (secs)
    ntiop  = 11  ,& ! initial time of observation period
    ntlbx  = 12  ,& ! length of AOF time boxes
    nteop  = 13  ,& ! end of observation period     
    nlchrd = 14  ,& ! character area
    nlusrd = 15  ,& ! user area
    nrecd  = 16  ,& ! record number
    nchkwd = 17     ! check word

 
!         1.3   AOF data record header format
!               -----------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nhdob  = 19  ,& ! header len.
    nrecln =  1  ,& ! record len.
    nprarr =  2  ,& ! len. of prel. array
    nlnnxt =  3  ,& ! len. of next record
    nrecnm =  4  ,& ! record number
    nboxnm =  5  ,& ! anal. box number
    nobtp  =  6  ,& ! observation type
    ncdtp  =  7  ,& ! code type
    nlatit =  8  ,& ! observation latitude
    nlongt =  9  ,& ! observation longitude
    ndate  = 10  ,& ! date of observation
    nsyntm = 11  ,& ! synoptic time of observation
    nextim = 12  ,& ! exact time of observation
    nstid1 = 13  ,& ! part 1 of station id (char 1 to 4)
    nstid2 = 14  ,& ! part 2 of station id (char 5 to 8)
    ndbkey = 15  ,& ! data-base key
    naltit = 16  ,& ! station altitude
    nstcar = 17  ,& ! station characteristics
    ninstr = 18  ,& ! instrument specifications
    nqllta = 19  ,& ! flags on lat./lon./dat./tim./alt.
    nbufhed= 18     ! report header length in send buffer


!         1.3.1  Bit patterns for packed information
!                -----------------------------------
!                station characteristics (nstcar)

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                             word no.
!   --------          -------                             --------
    ntscbp =  5  ,& ! bit pos. for co-loc. satem/tovs     nstcar 
    ntscoc =  1  ,& ! no. of bits occ. by co-loc. ind.      "
    nstcbp =  4  ,& ! bit pos. for st. correct. indic.      "
    nstcoc =  1  ,& ! no. of bits occ. by st. corr. ind.    "
    nstibp =  3  ,& ! bit pos. for important st. indic.     "
    nstioc =  1  ,& ! no. of bits occ. by imp. st. ind.     "
    nscfbp =  2  ,& ! bit pos. for st. confidence indic.    "
    nscfoc =  1  ,& ! no. of bits occ. by st. conf. ind.    "
    ndbfbp =  0  ,& ! bit pos. for data base form. ind.     "
    ndbfoc =  1  ,& ! no. of bits occ. by db form. ind      "
    narabp = 14  ,& ! bit pos. for aircraft roll angle      "
    naraoc =  4  ,& ! no. of bits occ. by airc. roll angle  "
    npafbp = 18  ,& ! bit position for phase of airc. flight"
    npafoc =  4     ! no. of bits occ. by phase of flight   "


!         1.3.2   instrument specifications (ninstr)
!                 ----------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                             word no.
!   --------          -------                             --------
    ninsbp =  0  ,& ! bit pos. for instr. spec. indic.    ninstr
    ninsoc =  7  ,& ! no. of bits occ. by ins. spec. in.    "
    ni1bp  =  7  ,& ! bit position for I1 parameter, name   "
                    !   of country which operates satellite
    ni1oc  =  4  ,& ! no. of bits occ. by I1                "
    ni2bp  = 11  ,& ! bit position for I2I2 parameter       "
                    !   indicator figure for satelite name
    ni2oc  =  7  ,& ! no. of bits occ. by I2                "
    ni3bp  = 22  ,& ! bit position for I3 parameter,        "
                    !   instrument data used in processing
    ni3oc  =  4  ,& ! no. of bits occ. by I1                "
    ni4bp  = 18  ,& ! bit position for I4 parameter,        "
                    !   type of measurement equipment used
                    !   or data processing technique
    ni4oc  =  4     ! no. of bits occ. by I4                "


!         1.3.3   Bit patterns for flag information
!                 ---------------------------------
 
!                 Bit position for each flag

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable             meaning                                      word no.
!   --------             -------                                      --------
                       ! bit pos. for flags on                        nqllta
                       !    lat.  lon.  date  time  alt.                 "
    nf3bps(5) = (/            0,    6,   12,   18,   24 /)                    ,&
    nf3boc    =  6  ,& ! no. of bits occ. by each flag                   "
                       ! inner bit structure (pos./no.) for each flag    "
                       ! hMsubst QCsubs ovrrd  flag QC/hM                "
    nf3ibp(5) = (/            0,    1,    2,    3,    5 /)                    ,&
    nf3ioc(5) = (/            1,    1,    1,    2,    1 /)                    ,&
    nf3rbp    = 30  ,& ! bit pos. for redundancy flag                    "
    nf3roc    =  1     ! no. of bits occ. by redunda. flag               "


!         1.3.4   Constants added on certain parameters (lat/lon/alt)
!                 ---------------------------------------------------
!                     to avoid having negative value

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable             meaning                               word no.
!   --------             -------                               --------
    nltmcn =  9000  ,& ! constant added to lat.                nlatit
    nlgmcn = 18000  ,& ! constant added to lon.                nlongt
    nalmcn =  1000     ! constant added to altit.              naltit
 

!         1.4      AOF data record body format
!                  ---------------------------

!         1.4.1   Level len. arrays and multi/single level report switch
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    ! level 0 length array / level n length array of obs. type:
!                   ! SYNOP   AIREP   SATOB   DRIBU    TEMP   PILOT   SATEM 
    nlev0 (7) = (/        0,      0,      3,      0,      0,      0,      7/) ,&
    nlevn (7) = (/       32,      8,      8,      7,      8,      7,      8/)

  LOGICAL                  , PARAMETER  :: &
                    ! multi/single level report switch of obs. type:
!                   ! SYNOP   AIREP   SATOB   DRIBU    TEMP   PILOT   SATEM
    lmulti(7) = (/   .false.,.false.,.false.,.false., .true., .true., .true./)


!         1.4.2  AOF data record body format

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    !              p d f t td pt z dz w c rl q st
                    !        synop x x x x  x  x x    x x      x
                    !        airep x x x                     x
                    !        satob x x x                       x
                    !        dribu x x x                       x
                    !        temp  x x x x  x    x
                    !        pilot x x x x       x
                    !        satem x     x          x      x x x
    nwrvr (13,7)= & ! body format matrix of: (variable, obs type)
                    !    PP  DD  FF  TT  TD  PT  ZZ  DZ  WW  CL  RL  QQ  ST
             RESHAPE( (/  1,  2,  3,  4,  5,  6,  0,  0,  8,  9,  0,  0,  7,   &
                          1,  2,  3,  4,  6,  0,  5,  0,  0,  0,  0,  0,  0,   &
                          1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  1,   &
                          1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  5,   &
                          1,  2,  3,  4,  5,  0,  6,  0,  0,  0,  0,  0,  0,   &
                          1,  2,  3,  4,  0,  0,  5,  0,  0,  0,  0,  0,  0,   &
                          2,  0,  0,  3,  0,  0,  0,  5,  0,  0,  1,  4,  1/)  &
                    , (/13,7/) )


!         1.4.2   Bit patterns for packed information (weather/cloud group)
!                 --------------------------------------------------------
  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                            word no.
!   --------          -------                            --------
    nwbp  =  0  , & ! bit position  for w                nwrvr( 9,1)
    nwwbp =  7  , & !         "         ww                  "
    nvibp = 14  , & !         "         v                   "
    nwoc  =  7  , & ! no. of bits occupied by w             "
    nwwoc =  7  , & !           "             ww            "
    nvioc = 17  , & !           "             v             "
    nchbp =  0  , & ! bit position for ch                nwrvr(10,1)
    ncmbp =  4  , & !         "        cm                   "
    nclbp = 19  , & !         "        cl                   "
    nnbp  = 23  , & !         "        nn                   "
    nhbp  =  8  , & !         "        h                    "
    nbp   = 27  , & !         "        n                    "
    nchoc =  4  , & ! no. of bits occupied by ch            "
    ncmoc =  4  , & !           "             cm            "
    ncloc =  4  , & !           "             cl            "
    nnoc  =  4  , & !           "             nn            "
    nhoc  = 11  , & !           "             h             "
    noc   =  4      !           "             n             "
 

!         1.4.3   Additional group indic. and bit pattern (only SYNOP)
!                 ---------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                            word no.
!   --------          -------                            --------
    naddgr = 13  ,& ! additional groups indicator          13
    naclbp =  0  ,& ! bit pos. for additional clouds     naddgr
    nacloc =  4  ,& ! no. of bits occ. by above indic.      "
    nhsbp  =  0  ,& ! bit position for hshs              calculated
    ncbp   = 15  ,& !          "       c                    "
    nnsbp  = 19  ,& !          "       ns                   "
    nhsoc  = 15  ,& ! no. of bits occupied h                "
    ncoc   =  4  ,& !              "       c                "
    nnsoc  =  4  ,& !              "       ns               "
    ngrnbp =  4  ,& ! bit pos. for ground group          naddgr
    ngrnoc =  1  ,& ! no. of bits occ. by above indic.      "
    nsbp   =  0  ,& ! bit position for s (snow depth)    calculated
    nebp   =  8  ,& !          "       e                    "      
    ntgbp  = 12  ,& !          "       tgtg                 "      
    nsoc   =  8  ,& ! no. of bits occupied s                "      
    neoc   =  4  ,& !              "       e                "
    ntgoc  = 12  ,& !              "       tgtg             "
    nsppbp =  5  ,& ! bit pos. for special phenomena     naddgr
    nsppoc =  2  ,& ! no. of bits occ. by above indic.      "
    nlspbp =  0  ,& ! bit position for spsp              calculated
    nuspbp =  7  ,& !          "       spsp                 "      
    nlspoc =  7  ,& ! no. of bits occupied spsp             "      
    nuspoc =  7  ,& !              "       spsp             "      
    nicebp =  7  ,& ! bit pos. for ice group             naddgr
    niceoc =  1  ,& ! no. of bits occ. by above indic.      "
    nrsbp  =  0  ,& ! bit position for rs              calculated
    nesbp  =  4  ,& !          "       eses                 "    
    nisbp  = 11  ,& !          "       is                   "    
    nrsoc  =  4  ,& ! no. of bits occupied rs               "    
    nesoc  =  7  ,& !              "       eses             "    
    nisoc  = 11  ,& !              "       is               "    
    nrinbp =  8  ,& ! bit pos. for rain group            naddgr
    nrinoc =  1  ,& ! no. of bits occ. by above indic.      "
    ntrbp  =  0  ,& ! bit position for trtr            calculated
    nrrbp  =  7  ,& !          "       rr                   "
    ntroc  =  7  ,& ! no. of bits occupied trtr             "    
    nrroc  = 14     !              "       rr               "    

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nshpbp =  9  ,& ! bit pos. for ship group            naddgr
    nshpoc =  1  ,& ! no. of bits occ. by above indic.      "
    nvsbp  =  0  ,& ! bit position for vs              calculated
    ndsbp  =  8  ,& !          "       ds                   "     
    nvsoc  =  8  ,& ! no. of bits occupied vs               "     
    ndsoc  = 10  ,& !              "       ds               "     
    nwavbp = 10  ,& ! bit pos. for waves group           naddgr
    nwavoc =  3  ,& ! no. of bits occ. by above indic.      "
    nhwbp  =  0  ,& ! bit position for hwhw            calculated
    npwbp  = 10  ,& !          "       pwpw                 "    
    ndwbp  = 17  ,& !          "       dwdw                 "    
    nhwoc  = 10  ,& ! no. of bits occupied hwhw             "    
    npwoc  =  7  ,& !              "       pwpw             "    
    ndwoc  = 10  ,& !              "       dwdw             "    
    nextbp = 13  ,& ! bit pos. for extreme temp. group   naddgr
    nextoc =  1  ,& ! no. of bits occ. by above indic.      "
    ntnbp  =  0  ,& ! bit position for tntn            calculated
    ntxbp  = 12  ,& !          "       txtx                 "
    ntnoc  = 12  ,& ! no. of bits occupied tntn             "    
    ntxoc  = 12  ,& !              "       txtx             "    
    nradbp = 14  ,& ! bit pos. for radiation group       naddgr
    nradoc =  1  ,& ! no. of bits occ. by above indic.      "
    nffbp  =  0  ,& ! bit position for ffff            calculated
    nffoc  = 12     ! no. of bits occupied ffff             "
 

!         1.4.4   Flag formats (set as 2-d arr. of varb.,obs. typ.)
!                 -------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    !              p d f t td pt z dz w c rl q st
                    !        synop x x x x  x  x x    x x      x
                    !        airep x x x x
                    !        satob x x x                       x
                    !        dribu x x x x                     x
                    !        temp  x x x x  x    x
                    !        pilot x x x x       x
                    !        satem x     x          x      x x x
    nwrvrf(13,7)= & ! flag format matrix of: (variable, obs type)
                    !    PP  DD  FF  TT  TD  PT  ZZ  DZ  WW  CL  RL  QQ  ST
             RESHAPE( (/ 10, 11, 11, 11, 11, 12,  0,  0, 12, 12,  0,  0, 12,   &
                          7,  8,  8,  8,  8,  0,  8,  0,  0,  0,  0,  0,  0,   &
                          7,  7,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0,  3,   &
                          6,  7,  7,  7,  0,  0,  0,  0,  0,  0,  0,  0,  7,   &
                          7,  7,  8,  8,  8,  0,  8,  0,  0,  0,  0,  0,  0,   &
                          6,  6,  7,  7,  0,  0,  7,  0,  0,  0,  0,  0,  0,   &
                          7,  0,  0,  8,  0,  0,  0,  8,  0,  0,  7,  8,  6/)  &
                    , (/13,7/) )

!         1.4.5     6 bit struct. within flag words (2-d arr.)
!                   -----------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    !              p d f t td pt z dz w c rl q st
                    !        synop x x x x  x  x x    x x      x
                    !        airep x x x x       x           x
                    !        satob x x x                       x
                    !        dribu x x x x                     x
                    !        temp  x x x x  x    x
                    !        pilot x x x x       x
                    !        satem x     x          x      x x x
    n6vfb (13,7)= & ! 6 bit flag patterns matrix of: (variable, obs type)
                    !       def. by 'nwrvrf' matrix 
                    !    PP  DD  FF  TT  TD  PT  ZZ  DZ  WW  CL  RL  QQ  ST
             RESHAPE( (/  1,  1,  2,  3,  4,  1,  0,  0,  3,  4,  0,  0,  2,   &
                          1,  1,  2,  3,  5,  0,  4,  0,  0,  0,  0,  0,  0,   &
                          1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  1,   &
                          1,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  4,   &
                          1,  2,  1,  2,  3,  0,  4,  0,  0,  0,  0,  0,  0,   &
                          1,  2,  1,  2,  0,  0,  3,  0,  0,  0,  0,  0,  0,   &
                          2,  0,  0,  1,  0,  0,  0,  3,  0,  0,  1,  2,  1/)  &
                    , (/13,7/) )
 

!         1.4.6    Bit pattern within 6 bit structures
!                  -----------------------------------
  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable             meaning
!   --------             -------
    nf1aoc    =  6  ,& ! no. bits in active part of 6bit pattern
                       ! bit pos. / no. bits within 6bit pattern:     matrix
                       !    HMSUB QCSUB OVRID FLAG  QCHMF
    nf1bps(5) = (/            0,    1,    2,    3,    5 /)                    ,&
    nf1boc(5) = (/            1,    1,    1,    2,    1 /)


!         1.4.7    Other bit patterns for flags
!                  ----------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable          meaning                              word no.
!   --------          -------                              --------
    nwfbp  =  0  ,& ! bit position for flag on w           defined by
    nwwfbp =  2  ,& !             "            ww          nwrvrf matrix
    nvfbp  =  4  ,& !             "            v                "
    nchfbp =  0  ,& !             "            ch   (cloud)     "
    ncmfbp =  2  ,& !             "            cm               "
    nnhfbp =  4  ,& !             "            nh               "
    nclfbp =  6  ,& !             "            cl               "
    nnnfbp =  8  ,& !             "            nn               "
    nnfbp  = 10  ,& !             "            n    (add cld)   "
    nhsfbp =  0  ,& !             "            hshs             "
    ncfbp  =  2  ,& !             "            c                "
    nnsfbp =  4  ,& !             "            ns               "
    nsfbp  =  0  ,& !             "            s    (ground)    "
    nefbp  =  2  ,& !             "            e                "
    ntgfbp =  4  ,& !             "            tgtg             "
    nlspfb =  0  ,& !             "            spsp (special)   "
    nuspfb =  2  ,& !             "            spsp             "
    nrsfbp =  0  ,& !             "            rs   (ice)       "
    nesfbp =  2  ,& !             "            eses             "
    nisfbp =  4  ,& !             "            ises             "
    ntrfbp =  0  ,& !             "            trtr (rain)      "
    nrrfbp =  2  ,& !             "            rr               "
    nvsfbp =  0  ,& !             "            vs   (ship)      "
    ndsfbp =  2  ,& !             "            ds               "
    nhwfbp =  0  ,& !             "            hwhw (wave)      "
    npwfbp =  2  ,& !             "            pwpw             "
    ndwfbp =  4  ,& !             "            dwpw             "
    nttnbp =  0  ,& !             "            tntn (ext.T)     "
    nttxbp =  2  ,& !             "            txtx             "
    nfffbp =  0  ,& !             "            ffff (rad.)      "
    nf2boc =  2     ! no. of bits occupied by each of           "
                    ! above flag (annex f2 of AOF documentation)


!         1.4.8   Bit pattern for level indicator and pressure code
!                 -------------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable             meaning                              word no.
!   --------             -------                              --------
    nlinoc    =  9  ,& ! no. of bits occupied by level id     pressure
                       ! bit position for level id            flag word
!                      ! SYNOP   AIREP   SATOB   DRIBU   TEMP    PILOT   SATEM 
    nlinbp(7) =    (/     -1,     -1,     -1,     -1,     12,     12,     12/),&
    npcdbp    =  6  ,& ! bit position for pressure code           "
    npcdoc    =  4  ,& ! no. of bits occupied by press. code      "
                       ! bit pattern for level id:                " 
                       !   nlidbp(1) = maximum wind
                       !   nlidbp(2) = tropopause
                       !   nlidbp(3) = d part
                       !   nlidbp(4) = c part
                       !   nlidbp(5) = b part
                       !   nlidbp(6) = a part
                       !   nlidbp(7) = surface level
                       !   nlidbp(8) = significant wind
                       !   nlidbp(9) = significant temperature
    nlidbp(9) =    (/  0,  1,  2,  3,  4,  5,  6,  7,  8/)                    ,&
    nlidoc    =  1     ! no. of bits occ. by each of above info


!         1.4.9   Constants added to avoid having negative values
!                 -----------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
!   variable             meaning                                   word no.
!   --------             -------                                   --------
    nptc24  = 5000  ,& ! const. added on 24hr press. tendency      nwrvr(6,1)
    nptc03  =  500     ! const. added on  3hr press. tendency


!         1.5     Extracted parameters and Variables
!                 ----------------------------------

  LOGICAL                   ::    &
    llship     ,& ! ship switch
    lsurfob       ! if true then derive surface report from temp

  INTEGER (KIND=iintegers)  ::    &
!   variable        meaning              
!   --------        -------             
    nobtyp     ,& ! observation type
    ncdtyp     ,& ! code type
    ntime      ,& ! actual time of observation
    nvarib     ,& ! current variable
    nwwgr      ,& ! weather group
    ngclgr     ,& ! general cloud group
    naclgr(4)  ,& ! additonal cloud groups
    ngrngr     ,& ! ground group
    nsppgr(2)  ,& ! special phenomena group
    nicegr     ,& ! ice group
    nraing     ,& ! rain group
    nshpgr     ,& ! ship group
    nwavgr(3)  ,& ! waves group
    nextgr     ,& ! extreme temperature group
    nradgr     ,& ! radiation group flag word   "
 
! flags to be extracted
    nvrfwr     ,& ! current variable flag (bit pattern)
    nflag      ,& ! actual flag for current variable (0-3)
    nlevin     ,& ! level id. (bit pattern)
    nwwfw      ,& ! weather flag word           "
    ngclfw     ,& ! cloud flag word
    naclfw(4)  ,& ! additonal clouds flag word  "
    ngrnfw     ,& ! ground group flag word      "
    nsppfw(2)  ,& ! special phenomena flag word "
    nicefw     ,& ! ice group flag word         "
    nranfw     ,& ! rain group flag word        "
    nshpfw     ,& ! ship movement flag word     "
    nwavfw(3)  ,& ! waves group flag word       "
    nextfw     ,& ! extr. temp. group flag word "
    nradfw        ! radiation group flag word   "

  REAL (KIND=ireals)        ::    &
    taof       ,& ! observation time in forecast hour units
    rppp       ,& ! reported pressure level
    roblat     ,& ! station lat.
    roblon        ! station lon.


!         1.6     AOF time boxes
!                 --------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    naofbox = 1   ! length of AOF time periods to be read
                  ! (to specify shorter periods, replacement of INTEGER
                  !  by REAL only needs adjustment of WRITE statements)

  REAL    (KIND=ireals)     ::    &
    naoflbx    ,& ! initial time of last  AOF time box to be read currently
    naoffbx       ! initial time of first AOF time box to be read currently


 
!-------------------------------------------------------------------------------
! Section 2 : Numerical, scaling and other constants related to AOF
!-------------------------------------------------------------------------------

!         2.1     Numerical constants and scaling constants
!                 -----------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nmdi   =  2147483647     ! missing data indicator for AOF data (2^31 -1)
 
  REAL    (KIND=ireals)    , PARAMETER  :: &
    atol   =  0.0001_ireals  ! area tolerance (because of rounding off problem)
 
  REAL    (KIND=ireals)    , PARAMETER  :: &
                             ! scaling factors
    rfact (16) = (/    10._ireals,  1.0_ireals,  1.0_ireals,  0.1_ireals,      &
!                          PP           U            V            TT    
                       0.1_ireals,  10._ireals,  1.0_ireals,  1.0_ireals,      &
!                          TDT          RPT          ZZ           DZ  
                       1.0_ireals,  1.0_ireals,  10._ireals,  1.0_ireals,      &
!                          WW           CL           RLP          QQ
                       0.1_ireals,  0.01_ireals, 0.1_ireals,  1.0_ireals/)
!                          STT          SNW          RR           TR
 
!         2.2     Pressure and geopotential constants
!                 -----------------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
                  ! pressure multiple (depending on pressure code)
    npfact (12) = (/      1 ,      1 ,      0 ,      0 ,      1 ,      1,      &
                          1 ,      1 ,      1 ,      0 ,      0 ,      0 /)   ,&
                  ! geopot. multiple (depending on pressure code)
    nzfact (12) = (/      0 ,      0 ,      1 ,      1 ,      0 ,      0,      &
                          0 ,      0 ,      0 ,      1 ,      1 ,      1 /)   ,&
                  ! pressure flag multiple
    npflfc (12) = (/      1 ,      1 ,      0 ,      0 ,      1 ,      1,      &
                          1 ,      1 ,      1 ,      0 ,      0 ,      0 /)   ,&
                  ! geopot. flag multiple
    nzflfc (12) = (/      0 ,      0 ,      1 ,      1 ,      0 ,      0,      &
                          0 ,      0 ,      0 ,      1 ,      1 ,      1 /)

  REAL    (KIND=ireals)    , PARAMETER  :: &
                  ! values of pressure for different pressure codes
    rpset  (12) = (/    0._ireals,    0._ireals, 8500._ireals, 7000._ireals,   &
                        0._ireals,    0._ireals,    0._ireals,    0._ireals,   &
                        0._ireals, 9000._ireals,10000._ireals, 5000._ireals/) ,&
                  ! values of geopotential for different pressure codes
    rzset  (12) = (/    0._ireals,    0._ireals,    0._ireals,    0._ireals,   &
                      500._ireals, 1000._ireals, 2000._ireals, 3000._ireals,   &
                     4000._ireals,    0._ireals,    0._ireals,    0._ireals/)

 
!-------------------------------------------------------------------------------
! Section 3 : Variable extraction and output numbers related to AOF
!-------------------------------------------------------------------------------
 
!         3.0.1 Obs. type numbers               --->  moved to 'data_nudge_all'
!               -----------------

!         3.0.2 Observation code type numbers   --->  moved to 'data_nudge_all'
!               -----------------------------
 
!         3.1   Variables extraction numbering
!               ------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
                    ! variables extraction numbering:
                    !   nvar  ( 1) = 1 (pressure)
                    !     "   ( 2) = 2 (wind direction)
                    !     "   ( 3) = 3 (wind speed)
                    !     "   ( 4) = 4 (temperature or 2m temp.)
                    !     "   ( 5) = 5 (dew-point or 2m dew-point)
                    !     "   ( 6) = 6 (pressure tendency)
                    !     "   ( 7) = 7 (geopotential)
                    !     "   ( 8) = 8 (thickness)
                    !     "   ( 9) = 9 (weather group)
                    !     "   (10) =10 (general cloud group)
                    !     "   (11) =11 (precitable water content)
                    !     "   (12) =12 (surface temperature)
    nvar  (13) = (/  1,  2,  3,  4,  5,  6 , 7,  8,  9, 10, 11, 12, 13/)
 

!         3.2    Variables' extraction inventory table
!                -------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nzex  (7) = (/ 1, 1, 0, 1, 1, 1, 0/) ,& ! geopotential extraction table
    nuex  (7) = (/ 1, 1, 1, 1, 1, 1, 0/) ,& ! u-comp. extraction table
    nvex  (7) = (/ 1, 1, 1, 1, 1, 1, 0/) ,& ! v-comp. extraction table
    ntex  (7) = (/ 0, 1, 0, 0, 1, 1, 1/) ,& ! temperature extraction table
    ntdex (7) = (/ 0, 1, 0, 0, 1, 0, 0/) ,& ! dew-point extraction table
    nqex  (7) = (/ 0, 0, 0, 0, 0, 0, 1/) ,& ! water content extraction table
    nt2ex (7) = (/ 1, 0, 0, 1, 0, 0, 0/) ,& ! 2m temperature extraction table
    ntd2ex(7) = (/ 1, 0, 0, 0, 0, 0, 0/)    ! 2m dew-point extraction table
 

!         3.2.1  Variables' extraction inventory table:
!                --------------------------------------
!                SYNOP additional groups
!                -----------------------  

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nwwex  (3) = (/       & ! weather group extraction table
                    0    ,& !   nwwex ( 1)  =  past weather
                    1    ,& !   nwwex ( 2)  =  ww
                    1 /) ,& !   nwwex ( 3)  =  v
    ngclex (6) = (/       & ! general cloud group extraction table
                    0    ,& !   ngclex( 1)  =  ch
                    0    ,& !   ngclex( 2)  =  cm
                    1    ,& !   ngclex( 3)  =  h
                    0    ,& !   ngclex( 4)  =  cl
                    1    ,& !   ngclex( 5)  =  nh
                    1 /) ,& !   ngclex( 6)  =  n
    naclex (3) = (/       & ! additional cloud groups extraction table
                    0    ,& !   naclex( 1)  =  hshs
                    0    ,& !   naclex( 2)  =  c
                    0 /) ,& !   naclex( 3)  =  ns
    ngrnex (3) = (/       & ! ground group extraction table
                    0    ,& !   ngrnex( 1)  =  s
                    0    ,& !   ngrnex( 2)  =  e
                    1 /) ,& !   ngrnex( 3)  =  tgtg
    nsppex (2) = (/       & ! special phenomena group table
                    1    ,& !   nsppex( 1)  =  spsp
                    1 /)    !   nsppex( 2)  =  spsp

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    niceex (3) = (/       & ! ice group extraction table
                    0    ,& !   niceex( 1)  =  rs
                    0    ,& !   niceex( 2)  =  eses
                    0 /) ,& !   niceex( 3)  =  is
    nraine (2) = (/       & ! rain group extraction table
                    1    ,& !   nraine( 1)  =  trtr
                    1 /) ,& !   nraine( 2)  =  rr
    nshipe (2) = (/       & ! ship group extraction table
                    0    ,& !   nshipe( 1)  =  vs
                    0 /) ,& !   nshipe( 2)  =  ds
    nwavee (3) = (/       & ! wave groups extraction table
                    0    ,& !   nwavee( 1)  =  hwhw
                    0    ,& !   nwavee( 2)  =  pwpw
                    0 /) ,& !   nwavee( 3)  =  dwdw
    nextex (2) = (/       & ! extreme temperature group extraction table
                    1    ,& !   nextex( 1)  =  tntn
                    1 /) ,& !   nextex( 2)  =  txtx
    nradex (1) = (/       & ! radiation group extraction table
                    1 /)    !   nradex( 1)  =  ffff
 

!         3.3    Variables' output inventory table
!                ---------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nzot  (7) = (/ 1, 1, 0, 1, 1, 1, 0/) ,& ! geopotential output table
    nuot  (7) = (/ 1, 1, 1, 1, 1, 1, 0/) ,& ! u-comp. output table
    nvot  (7) = (/ 1, 1, 1, 1, 1, 1, 0/) ,& ! v-comp. output table
    ntot  (7) = (/ 0, 1, 0, 0, 1, 1, 1/) ,& ! temperature output table
    ntdot (7) = (/ 0, 1, 0, 0, 1, 0, 0/) ,& ! dew-point output table
    nqot  (7) = (/ 0, 0, 0, 0, 0, 0, 1/) ,& ! water content output table
    nt2ot (7) = (/ 1, 0, 0, 1, 0, 0, 0/) ,& ! 2m temperature output table
    ntd2ot(7) = (/ 1, 0, 0, 0, 0, 0, 0/)    ! dew-point output table
!   nstot (7) = (/ 0, 0, 1, 0, 0, 0, 0/) ,& ! surface temperature output table
!   nptot (7) = (/ 0, 0, 0, 0, 0, 0, 0/) ,& ! pressure tendency output table
!   nrhot (7) = (/ 1, 0, 0, 0, 1, 0, 0/)    ! relative humidity output table

 
!         3.4    Variables' numbering for obs. error extraction
!                ----------------------------------------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nvrpoe (6) = (/ 1, 2, 3, 4, 5, 6/)      ! variable numbering for obs. error
                                            !                      calculation:
                                            ! nvrpoe(1) = 1 (u component)
                                            !   "   (2) = 2 (v component)
                                            !   "   (3) = 3 (geopotential)
                                            !   "   (4) = 4 (thickness)
                                            !   "   (5) = 5 (relative humidity)
                                            !   "   (6) = 6 (temperature)
 
 
!-------------------------------------------------------------------------------
! Section 4 : Diagnostic arrays and pointers for event counters (AOF input only)
!-------------------------------------------------------------------------------
 
  INTEGER (KIND=iintegers) , PARAMETER  :: &
    nmxcdt =  6  ,& ! max. no. of code type for any obs. type
    ntotob =  8     ! total number of obs. types

  INTEGER (KIND=iintegers)  ::    &
                                                   ! 2-d diagnostic arrays of:
    noctpr(8,6) = RESHAPE( (/0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !   processed obs.
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0/)                 , (/8,6/) )        ,&
    noctac(8,6) = RESHAPE( (/0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !   active    obs.
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0/)                 , (/8,6/) )        ,&
    noctps(8,6) = RESHAPE( (/0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !   passive   obs.
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0/)                 , (/8,6/) )        ,&
    noctrj(8,6) = RESHAPE( (/0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !   rejected  obs.
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0,0,0,0,0,0,0,0,0   & !
                            ,0,0,0,0,0,0/)                 , (/8,6/) )        ,&
    nobtpp      =            0     ,& ! array pointer for obs. type
    ncdtpp      =            0        ! array pointer for obs. code type

 
 
!-------------------------------------------------------------------------------
! Section 5 : Character descriptions of events, flags, observation & code types
!-------------------------------------------------------------------------------
 
!         5.1    Character descriptions of events and flags
!                ------------------------------------------

!         5.1.1  Report events
!                -------------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxcrev  = 18    ! length of data event decriptor

!         5.1.2  Data events
!                -----------

  INTEGER (KIND=iintegers) , PARAMETER  :: &
    mxcdev  = 37    ! length of data event decriptor

!         5.1.3  Report flags
!                ------------

  CHARACTER (LEN=50)       , PARAMETER  :: &
                    !  Description of report flags
    crepfl (5)  = (/ 'HUMAN MONITOR SUBSTITUTION                        '     ,&
                     'QUALITY CONTROL SUBSTITUTION                      '     ,&
                     'OVERRIDE BY HUMAN MONITOR                         '     ,&
                     'FLAG ON PARAMETER                                 '     ,&
                     'FLAG SET BY HUMAN MONITOR (1)/BY Q.CONTROL (0)    '/)


!         5.1.4  Data flags
!                ----------

  CHARACTER (LEN=50)       , PARAMETER  :: &
                    !  Description of data flags
    cdatfl (9)  = (/ 'HUMAN MONITOR SUBSTITUTION                        '     ,&
                     'QUALITY CONTROL SUBSTITUTION                      '     ,&
                     'OVERRIDE BY HUMAN MONITOR                         '     ,&
                     'FLAG ON PARAMETER                                 '     ,&
                     'FLAG SET BY HUMAN MONITOR (1)/BY Q.CONTROL (0)    '     ,&
                     'DEGREE OF CONFIDENCE GIVEN BY A PREVIOUS ANALYSIS '     ,&
                     'VARIABLE USED (1)/NOT USED (0) BY PREV. ANALYSIS  '     ,&
                     '                                                  '     ,&
                     '                                                  '/)


!         5.1.5  Upper air level identifier
!                -------------------------

  CHARACTER (LEN=30)       , PARAMETER  :: &
                    !  Description of upper air levels
    clevua (9)  = (/ 'MAX WIND LEVEL                '                         ,&
                     'TROPOPAUSE                    '                         ,&
                     'D PART                        '                         ,&
                     'C PART                        '                         ,&
                     'B PART                        '                         ,&
                     'A PART                        '                         ,&
                     'SURFACE LEVEL                 '                         ,&
                     'SIGNIFICANT WIND LEVEL        '                         ,&
                     'SIGNIFICANT TEMPERATURE LEVEL '/)


!         5.1.6  SYNOP level identifier
!                ----------------------

  CHARACTER (LEN=30)       , PARAMETER  :: &
                    !  Description of synop level
    clevsy (12) = (/ '0  - SEA LEVEL                '                         ,&
                     '1  - STATION LEVEL PRESSURE   '                         ,&
                     '2  - 850MB GEOPOTENTIAL       '                         ,&
                     '3  - 700MB GEOPOTENTIAL       '                         ,&
                     '4  - 500GPM PRESSURE          '                         ,&
                     '5  - 1000GPM PRESSURE         '                         ,&
                     '6  - 2000GPM PRESSURE         '                         ,&
                     '7  - 3000GPM PRESSURE         '                         ,&
                     '8  - 4000GPM PRESSURE         '                         ,&
                     '9  - 900MB GEOPOTENTIAL       '                         ,&
                     '10 - 1000MB GEOPOTENTIAL      '                         ,&
                     '11 - 500MB GEOPOTENTIAL       '/)


!         5.2   Character decriptions of observation types and code types
!               ---------------------------------------------------------

!         5.2.1  Observation types
!                -----------------

  CHARACTER (LEN=8)        , PARAMETER  :: &
                    ! Description of observation types
    chobtp (8)  = (/ 'SYNOP   ','AIREP   ','SATOB   ','DRIBU   ',              &
                     'TEMP    ','PILOT   ','SATEM   ','GPS     ' /)
 

!         5.2.2  Observation code types
!                ----------------------

  CHARACTER (LEN=20)       , PARAMETER  :: &
                    ! Description of observation code types
    chobcd (8,6) =  RESHAPE(                                                   &
                    !
  (/'MANUAL LAND SYNOP   ',    'CODAR AIRCRAFT      ',  'SATOB SATELLITE     ',&
         'BATHY SPHERE        ',    'LAND RADIO SONDE    ',                    &
              'LAND PILOT          ',    'MSG-1 SEVIRI        ',               &
                   'GPZ GPS             ',                                     &
    'AUTOM. LAND SYNOP   ',    'AIREP AIRCRAFT      ',  'SST ANALYSES        ',&
         'TESAC               ',    'SHIP RADIO SONDE    ',                    &
              'SHIP PILOT          ',    'NOAA-15 ATOVS       ',               &
                   '                    ',                                     &
    'MANUAL SHIP SYNOP   ',    'CONST LEV BALLOON   ',  '                    ',&
         'DRIFTING BUOY       ',    'DROP SONDE          ',                    &
              'WIND PROFILER EUR   ',    'NOAA-16 ATOVS       ',               &
                   '                    ',                                     &
    'ABBR. SHIP SYNOP    ',    'AMDAR AIRCRAFT      ',  '                    ',&
         '                    ',    'LAND ROCKET SONDE   ',                    &
              'SODAR/RASS EUR      ',    'NOAA-17 ATOVS       ',               &
                   '                    ',                                     &
    'REDUCED SHIP SYNOP  ',    'ACAR  AIRCRAFT      ',  '                    ',&
         '                    ',    'SHIP ROCKET SONDE   ',                    &
              'PROFILER/RASS USA   ',    'NOAA-18 ATOVS       ',               &
                   '                    ',                                     &
    'AUTOM. SHIP SYNOP   ',    '                    ',  '                    ',&
         '                    ',    '                    ',                    &
              'RADAR VAD           ',    '                    ',               &
                   '                    '                      /)  ,  (/8,6/)  )


!-------------------------------------------------------------------------------
! Section 6 : Temporary information buffers and report counters
!-------------------------------------------------------------------------------

  INTEGER (KIND = iintegers), ALLOCATABLE :: &
    aofbuf    (:)   ,& ! observation report in aof format
    nodbuf    (:)      ! buffer containing all reports for a single node

  INTEGER (KIND = iintegers)  :: &
    nmlob           ,& ! current number of  multi-level report
    nsgob           ,& ! current number of single-level report
    ngpob              ! current number of GPS report


!-------------------------------------------------------------------------------
  END MODULE data_obs_process
