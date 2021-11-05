MODULE mo_memory_g3b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add => add_stream_element, &
                            default_stream_setting,                   &
                            ABOVESUR2, ABOVESUR10, BELOWSUR, HYBRID_H,&
                            HYBRID

  IMPLICIT NONE
  SAVE   ! mz_jb_20030805
  PRIVATE

  PUBLIC :: construct_g3b ! routine to construct the g3b table
  PUBLIC :: destruct_g3b  ! routine to destruct  the g3b table
  PUBLIC :: g3b           ! the g3b table

  ! declaration of predefined fields within this module

  REAL(dp), POINTER, PUBLIC :: geosp(:,:)   ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: tsl(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: ws(:,:)      ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: wl(:,:)      ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: sn(:,:)      ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: slm(:,:)     ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: az0(:,:)     ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: alb(:,:)     ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: forest(:,:)  ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: vgrat(:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: vlt(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: wsmx(:,:)    ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: fao(:,:)     ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: aps(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aprl(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aprc(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aprs(:,:)    ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: ustrgw(:,:)  ! -> OROGW
  REAL(dp), POINTER, PUBLIC :: vstrgw(:,:)  ! -> OROGW
  REAL(dp), POINTER, PUBLIC :: vdisgw(:,:)  ! -> OROGW
#endif
  REAL(dp), POINTER, PUBLIC :: aclcov(:,:) ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: temp2(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: dew2(:,:)    ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: wind10(:,:)  ! -> E5VDIFF
#endif
  REAL(dp), POINTER, PUBLIC :: u10(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: v10(:,:)    ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: srads(:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trads(:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: srad0(:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trad0(:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: vdis(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ustr(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: vstr(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ahfs(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: evap(:,:)   ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ahfl(:,:)   ! -> E5VDIFF
#endif
  REAL(dp), POINTER, PUBLIC :: tslm(:,:)      ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: tslm1(:,:)     ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: emter(:,:,:)   ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trsol(:,:,:)   ! -> RAD
#endif
  REAL(dp), POINTER, PUBLIC :: runoff(:,:)    ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: srad0u(:,:)    ! -> RAD
  REAL(dp), POINTER, PUBLIC :: sradsu(:,:)    ! -> RAD
  REAL(dp), POINTER, PUBLIC :: tradsu(:,:)    ! -> RAD
  REAL(dp), POINTER, PUBLIC :: albedo(:,:)    ! -> RAD
  REAL(dp), POINTER, PUBLIC :: tsurf(:,:)     ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: seaice(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: siced(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: relhum(:,:,:)  ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: wind10w(:,:)   ! -> E5VDIFF
#endif
  REAL(dp), POINTER, PUBLIC :: glac(:,:)      ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: gld(:,:)       ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aclc(:,:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aclcac(:,:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: snmel(:,:)     ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: runtoc(:,:)    ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: apmegl(:,:)    ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: t2max(:,:)     ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: t2min(:,:)     ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: wimax(:,:)     ! -> E5VDIFF (OBSOLETE)
#endif
  REAL(dp), POINTER, PUBLIC :: topmax(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: aclcv(:,:)
  REAL(dp), POINTER, PUBLIC :: qvi(:,:)       ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: xlvi(:,:)      ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: xivi(:,:)      ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: runlnd(:,:)    ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: rgcgn(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: sodif(:,:)     ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: srafs(:,:)     ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trafs(:,:)     ! -> RAD
  REAL(dp), POINTER, PUBLIC :: sraf0(:,:)     ! -> RAD
  REAL(dp), POINTER, PUBLIC :: traf0(:,:)     ! -> RAD
  REAL(dp), POINTER, PUBLIC :: emtef(:,:,:)   ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trsof(:,:,:)   ! -> RAD
  REAL(dp), POINTER, PUBLIC :: emtef0(:,:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: trsof0(:,:,:)  ! -> RAD
  REAL(dp), POINTER, PUBLIC :: tke(:,:,:)     ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: tkem(:,:,:)    ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: tkem1(:,:,:)   ! -> E5VDIFF
#endif
  REAL(dp), POINTER, PUBLIC :: xvar(:,:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: xskew(:,:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: drain(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: grndcapc(:,:) ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: grndhflx(:,:) ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: grndflux(:,:) ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: tsoil(:,:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: grndc(:,:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: grndd(:,:,:)  ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: srad0d(:,:)   ! -> RAD
#endif
  REAL(dp), POINTER, PUBLIC :: acdnc(:,:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: snacl(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: rogl(:,:)     ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: alake(:,:)    ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: aprflux(:,:)  ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: acvtype(:,:)  ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: xtec(:,:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: slf(:,:)      ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: snc(:,:)      ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: rtype(:,:)    ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: rintop(:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: apmeb(:,:)
  REAL(dp), POINTER, PUBLIC :: apmebco(:,:)
  REAL(dp), POINTER, PUBLIC :: qtnew(:,:)
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: abso4(:,:)    ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: so4nat(:,:,:) ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: so4all(:,:,:) ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: ao3(:,:,:)    ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: tropo(:,:)    ! OBSOLETE
#endif
  !
  !  variables for fractional surface coverage
  !
  REAL(dp), POINTER, PUBLIC :: tsi(:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: tsw(:,:)   ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: sni(:,:)   ! op_pj_20160617: d
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: ustri(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: vstri(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ustrw(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: vstrw(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ustrl(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: vstrl(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ahfli(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ahflw(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: ahfll(:,:) ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: az0i(:,:)  ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: az0w(:,:)  ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: az0l(:,:)  ! -> E5VDIFF
  REAL(dp), POINTER, PUBLIC :: alsoi(:,:) ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: alsow(:,:) ! OBSOLETE
  REAL(dp), POINTER, PUBLIC :: alsol(:,:) ! OBSOLETE
#endif
  REAL(dp), POINTER, PUBLIC :: ahfice(:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: qres(:,:)    ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: ahfcon(:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: ahfres(:,:)  ! op_pj_20160617: d
  REAL(dp), POINTER, PUBLIC :: fluxres(:,:) ! op_pj_20160617: d
  !
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: ahfsiac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: ahfswac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: ahfslac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: ahfliac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: ahflwac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: ahfllac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: evapiac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: evapwac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: evaplac(:,:) ! -> E5VDIFF (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: trfllac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: trflwac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: trfliac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: sofllac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: soflwac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: sofliac(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: friac(:,:)   ! -> (OBSOLETE)
#endif
  !
#ifdef OBSOLETE
  !  variables for ocean coupling only
  !
  REAL(dp), POINTER, PUBLIC :: awhea(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: awsol(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: awfre(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: awust(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: awvst(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: awsta(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: aicon(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: aiqre(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: aifre(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: aiust(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: aivst(:,:)   ! -> (OBSOLETE)
#endif
  REAL(dp), POINTER, PUBLIC :: ocu(:,:)
  REAL(dp), POINTER, PUBLIC :: ocv(:,:)
  !
#ifdef OBSOLETE
  !  variables for coupling with HD-model and calving model only
  !
  REAL(dp), POINTER, PUBLIC :: aros(:,:)    ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: adrain(:,:)  ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: disch(:,:)   ! -> (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: apmecal(:,:) ! -> (OBSOLETE)
#endif
  !
  !  variables for sso parametrization
  !
  REAL(dp), POINTER, PUBLIC :: oromea(:,:) ! op_pj_20160617: ioinitial.f90 
  REAL(dp), POINTER, PUBLIC :: orostd(:,:) ! op_pj_20160617: ioinitial.f90 d
  REAL(dp), POINTER, PUBLIC :: orosig(:,:) ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: orogam(:,:) ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: orothe(:,:) ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: oropic(:,:) ! op_pj_20160617: ioinitial.f90
  REAL(dp), POINTER, PUBLIC :: oroval(:,:) ! op_pj_20160617: ioinitial.f90
  !
#ifndef MESSY
  !  variables for mixed layer ocean only
  !
  REAL(dp), POINTER, PUBLIC :: amlcorr(:,:)    ! -> MLOCEAN
  REAL(dp), POINTER, PUBLIC :: amlcorac(:,:)   ! -> MLOCEAN (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: amlheatac(:,:)  ! -> MLOCEAN (OBSOLETE)
  !
  !  variables for 200mb radiation
  !
  REAL(dp), POINTER, PUBLIC :: tradl(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: sradl(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: trafl(:,:) ! -> RAD (OBSOLETE)
  REAL(dp), POINTER, PUBLIC :: srafl(:,:) ! -> RAD (OBSOLETE)
#endif

  ! declaration of table with 2d- and 3d-field entries

  TYPE (t_stream), POINTER     :: g3b

CONTAINS

  SUBROUTINE construct_g3b

#ifdef OBSOLETE
    USE mo_control,    ONLY: lmidatm, lcouple, lhd
#else
    USE mo_control,    ONLY: lmidatm
#endif

    ! set default attributes for the g3b stream

    CALL default_stream_setting (g3b               &
                                ,lrerun=.TRUE.     &
                                ,lpost=.TRUE.      &
                                ,table=128 ,bits=16)

    ! Add fields to the g3b stream.
    ! despite some 3-d fields (lpost=.FALSE.) these fields are written out by default

    CALL add (g3b,'qtnew',    qtnew    ,lpost=.FALSE.,contnorest=.true.)
#ifndef MESSY
    CALL add (g3b,'abso4',    abso4    ,code=235,laccu=.TRUE. ,                  &
         longname='antropogenic sulfur burden'  ,units='kg/m**2' ,contnorest=.true.)
    CALL add (g3b,'so4nat',   so4nat   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='natural sulfate'  ,units='kg/kg' ,contnorest=.true.)
    CALL add (g3b,'so4all',   so4all   ,lpost=.FALSE.,leveltype=HYBRID,          &
         longname='total sulfur'  ,units='kg/m**2' ,contnorest=.true.)
    CALL add (g3b,'ao3',      ao3      ,code=236,leveltype=HYBRID,               &
         longname='ipcc ozone'  ,units='kg/m**2' ,contnorest=.true.)
    CALL add (g3b,'tropo',    tropo    ,code=237 ,contnorest=.true.,             &
         longname='WMO defined tropopause height'          ,units='Pa'       )
    CALL add (g3b,'trfliac',  trfliac  ,code= 91,laccu=.TRUE. ,longname='LW flux over ice'                       ,units='W/m**2'   )
    CALL add (g3b,'trflwac',  trflwac  ,code= 92,laccu=.TRUE. ,longname='LW flux over water'                     ,units='W/m**2'   )
    CALL add (g3b,'trfllac',  trfllac  ,code= 93,laccu=.TRUE. ,longname='LW flux over land'                      ,units='W/m**2'   )
    CALL add (g3b,'sofliac',  sofliac  ,code= 94,laccu=.TRUE. ,longname='SW flux over ice'                       ,units='W/m**2'   )
    CALL add (g3b,'soflwac',  soflwac  ,code= 95,laccu=.TRUE. ,longname='SW flux over water'                     ,units='W/m**2'   )
    CALL add (g3b,'sofllac',  sofllac  ,code= 96,laccu=.TRUE. ,longname='SW flux over land'                      ,units='W/m**2'   )
    CALL add (g3b,'friac',    friac    ,code= 97,laccu=.TRUE. ,longname='ice cover (fraction of grid box)'                         )
#endif
    CALL add (g3b,'tsi',      tsi      ,code=102              ,longname='surface temperature of ice'             ,units='K'        )
    CALL add (g3b,'tsw',      tsw      ,code=103              ,longname='surface temperature of water'           ,units='K'        )
#ifndef MESSY
    CALL add (g3b,'ustri',    ustri    ,code=104              ,longname='zonal      wind stress over ice'        ,units='Pa'       )
    CALL add (g3b,'vstri',    vstri    ,code=105              ,longname='meridional wind stress over ice'        ,units='Pa'       )
    CALL add (g3b,'ustrw',    ustrw    ,code=106              ,longname='zonal      wind stress over water'      ,units='Pa'       )
    CALL add (g3b,'vstrw',    vstrw    ,code=107              ,longname='meridional wind stress over water'      ,units='Pa'       )
    CALL add (g3b,'ustrl',    ustrl    ,code=108              ,longname='zonal      wind stress over land'       ,units='Pa'       )
    CALL add (g3b,'vstrl',    vstrl    ,code=109              ,longname='meridional wind stress over land'       ,units='Pa'       )
    CALL add (g3b,'ahfliac',  ahfliac  ,code=110,laccu=.TRUE. ,longname='latent heat flux over ice'              ,units='W/m**2'   )
    CALL add (g3b,'ahflwac',  ahflwac  ,code=111,laccu=.TRUE. ,longname='latent heat flux over water'            ,units='W/m**2'   )
    CALL add (g3b,'ahfllac',  ahfllac  ,code=112,laccu=.TRUE. ,longname='latent heat flux over land'             ,units='W/m**2'   )
    CALL add (g3b,'evapiac',  evapiac  ,code=113,laccu=.TRUE. ,longname='evaporation over ice'                   ,units='kg/m**2s' )
    CALL add (g3b,'evapwac',  evapwac  ,code=114,laccu=.TRUE. ,longname='evaporation over water'                 ,units='kg/m**2s' )
    CALL add (g3b,'evaplac',  evaplac  ,code=115,laccu=.TRUE. ,longname='evaporation over land'                  ,units='kg/m**2s' )
    CALL add (g3b,'az0i',     az0i     ,code=116              ,longname='roughness length over ice'              ,units='m'        )
    CALL add (g3b,'az0w',     az0w     ,code=117              ,longname='roughness length over water'            ,units='m'        )
    CALL add (g3b,'az0l',     az0l     ,code=118              ,longname='roughness length over land'             ,units='m'        )
    CALL add (g3b,'ahfsiac',  ahfsiac  ,code=119,laccu=.TRUE. ,longname='sensible heat flux over ice'            ,units='W/m**2'   )
    CALL add (g3b,'ahfswac',  ahfswac  ,code=120,laccu=.TRUE. ,longname='sensible heat flux over water'          ,units='W/m**2'   )
    CALL add (g3b,'ahfslac',  ahfslac  ,code=121,laccu=.TRUE. ,longname='sensible heat flux over land'           ,units='W/m**2'   )
    CALL add (g3b,'alsoi',    alsoi    ,code=122              ,longname='albedo of ice'                                            )
    CALL add (g3b,'alsow',    alsow    ,code=123              ,longname='albedo of water'                                          )
    CALL add (g3b,'alsol',    alsol    ,code=124              ,longname='albedo of land'                                           )
#endif
    CALL add (g3b,'ahfice',   ahfice   ,code=125              ,longname='conductive heat flux'                   ,units='W/m**2'   )
    CALL add (g3b,'qres',     qres     ,code=126              ,longname='residual heat flux for melting sea ice' ,units='W/m**2'   )
    CALL add (g3b,'alake',    alake    ,code=127,lpost=.FALSE.,longname='lake fraction of grid box'                                )
    CALL add (g3b,'rintop',   rintop   ,code=128,lpost=.FALSE.,longname='low level inversion      '                                )
    CALL add (g3b,'geosp',    geosp    ,code=129              ,longname='surface geopotential (orography)'       ,units='m**2/s**2')
    !               stp                      130                         temperature                                     K
    !                                        131                         u-velocity                                      m/s
    !                                        132                         v-velocity                                      m/s
    !                                        133                         specific humidity                               kg/kg
    CALL add (g3b,'aps',      aps      ,code=134              ,longname='surface pressure'                       ,units='Pa'       )
    !                                        135                         vertical velocity                               Pa/s
    CALL add (g3b,'acdnc',    acdnc    ,code=136,lpost=.FALSE.,longname='cloud droplet number concentration'     ,units='1/m**3'   )
    CALL add (g3b,'apmeb',    apmeb    ,code=137,laccu=.TRUE. ,longname='vert.integr.tendencies of water',bits=24,units='kg/m**2s' )
    !                         svo            138                         vorticity                                       1/s
    CALL add (g3b,'tslm1',    tslm1    ,code=139              ,longname='surface temperature of land'            ,units='K'        )
    CALL add (g3b,'ws',       ws       ,code=140              ,longname='soil wetness'                           ,units='m'        )
    CALL add (g3b,'sn',       sn       ,code=141              ,longname='snow depth'                             ,units='m'        )
    CALL add (g3b,'aprl',     aprl     ,code=142,laccu=.TRUE. ,longname='large scale precipitation'      ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aprc',     aprc     ,code=143,laccu=.TRUE. ,longname='convective  precipitation'      ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aprs',     aprs     ,code=144,laccu=.TRUE. ,longname='snow fall'                      ,bits=24,units='kg/m**2s' )
#ifndef MESSY
    CALL add (g3b,'vdis',     vdis     ,code=145,laccu=.TRUE. ,longname='boundary layer dissipation'             ,units='W/m**2'   )
    CALL add (g3b,'ahfs',     ahfs     ,code=146,laccu=.TRUE. ,longname='sensible heat flux'                     ,units='W/m**2'   )
    CALL add (g3b,'ahfl',     ahfl     ,code=147,laccu=.TRUE. ,longname='latent heat flux'                       ,units='W/m**2'   )
#endif
    !                                        148                         streamfunction                                  m**2/s
    !                                        149                         velocity potential                              m**2/s
    CALL add (g3b,'xivi',     xivi     ,code=150,laccu=.TRUE. ,longname='vertically integrated cloud ice'        ,units='kg/m**2'  )
    !                                        151                         mean sea level pressure                         Pa
    !                         stp(20)        152                         log surface pressure
    !                         xl             153                         cloud water                                     kg/kg
    !                         xi             154                         cloud ice                                       kg/kg
    !                         sd             155                         divergence                                      1/s
    !                                        156                         geopotential height                             gpm
    CALL add (g3b,'relhum',   relhum   ,code=157              ,longname='relative humidity'                   )
    !                                        158                         tendency of surface pressure                    Pa/s
#ifndef MESSY
    CALL add (g3b,'wind10w', wind10w   ,code=159,lpost=.FALSE.,longname='10m windspeed over water'               ,units='m/s'      )
#endif
    CALL add (g3b,'runoff',  runoff    ,code=160,laccu=.TRUE. ,longname='surface runoff and drainage'    ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'drain',   drain     ,code=161,laccu=.TRUE. ,longname='drainage'                       ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'aclc',    aclc      ,code=162,lpost=.FALSE.,longname='cloud cover'                                              )
    CALL add (g3b,'aclcv',   aclcv     ,code=163,lpost=.FALSE.,longname='total cloud cover'                                        )
    CALL add (g3b,'aclcov',  aclcov    ,code=164,laccu=.TRUE. ,longname='total cloud cover'                                        )
    CALL add (g3b,'u10',     u10       ,code=165              ,longname='10m u-velocity'    ,units='m/s'      )
    CALL add (g3b,'v10',     v10       ,code=166              ,longname='10m v-velocity'    ,units='m/s'      )
#ifndef MESSY
    CALL add (g3b,'temp2',   temp2     ,code=167              ,longname='2m temperature'    ,units='K'        )
    CALL add (g3b,'dew2',    dew2      ,code=168       ,longname='2m dew point temperature' ,units='K'        )
    CALL add (g3b,'tsurf',   tsurf     ,code=169,laccu=.TRUE. ,longname='surface temperature'                    ,units='K'        )
#endif
! mz_jd_20180605+
!!$ CALL add (g3b,'xvar',    xvar      ,code=170,lpost=.FALSE.,longname='variance of total water amount qv+qi+ql',units='kg/kg'    )
    CALL add (g3b,'xvar',    xvar      ,code=170,lpost=.FALSE.,longname='proxy for standard deviation of total water amount qv+qi+ql',units='kg/kg'    )
! mz_jd_20180605-
#ifndef MESSY
    CALL add (g3b,'wind10',  wind10    ,code=171,laccu=.TRUE. ,longname='10m windspeed'     ,units='m/s'      )
#endif
    CALL add (g3b,'slm',     slm       ,code=172              ,longname='land sea mask (1. = land, 0. = sea/lakes)'                )
    CALL add (g3b,'az0',     az0       ,code=173,lpost=.FALSE.,longname='roughness length'                       ,units='m'        )
    CALL add (g3b,'alb',     alb       ,code=174,lpost=.FALSE.,longname='surface background albedo'                                )
#ifndef MESSY
    CALL add (g3b,'albedo',  albedo    ,code=175              ,longname='surface albedo'                                           )
    CALL add (g3b,'srads',   srads     ,code=176,laccu=.TRUE. ,longname='net surface solar radiation'            ,units='W/m**2'   )
    CALL add (g3b,'trads',   trads     ,code=177,laccu=.TRUE. ,longname='net surface thermal radiation'          ,units='W/m**2'   )
    CALL add (g3b,'srad0',   srad0     ,code=178,laccu=.TRUE. ,longname='net top solar radiation'                ,units='W/m**2'   )
    CALL add (g3b,'trad0',   trad0     ,code=179,laccu=.TRUE. ,longname='top thermal radiation (OLR)'            ,units='W/m**2'   )
    CALL add (g3b,'ustr',    ustr      ,code=180,laccu=.TRUE. ,longname='u-stress'                               ,units='Pa'       )
    CALL add (g3b,'vstr',    vstr      ,code=181,laccu=.TRUE. ,longname='v-stress'                               ,units='Pa'       )
    CALL add (g3b,'evap',    evap      ,code=182,laccu=.TRUE. ,longname='evaporation'                    ,bits=24,units='kg/m**2s' )
#endif
! mz_jd_20180605+
!!$ CALL add (g3b,'xskew',   xskew     ,code=183,lpost=.FALSE.,longname='skewness of total water amount qv+qi+ql'                  )
    CALL add (g3b,'xskew',   xskew     ,code=183,lpost=.FALSE.,longname='beta distribution parameter q'                            )
! mz_jd_20180605-
#ifndef MESSY
    CALL add (g3b,'srad0d',  srad0d    ,code=184,laccu=.TRUE. ,longname='top incoming solar radiation'           ,units='W/m**2'   )
    CALL add (g3b,'srafs',   srafs     ,code=185,laccu=.TRUE. ,longname='net surf. solar radiation   (clear sky)',units='W/m**2'   )
    CALL add (g3b,'trafs',   trafs     ,code=186,laccu=.TRUE. ,longname='net surf. thermal radiation (clear sky)',units='W/m**2'   )
    CALL add (g3b,'sraf0',   sraf0     ,code=187,laccu=.TRUE. ,longname='net top solar radiation     (clear sky)',units='W/m**2'   )
    CALL add (g3b,'traf0',   traf0     ,code=188,laccu=.TRUE. ,longname='net top thermal radiation   (clear sky)',units='W/m**2'   )
    !                                        189                         surface solar cloud forcing                     W/m**2
    !                                        190                         surface thermal cloud forcing                   W/m**2
    !                                        191                         SW top cloud forcing (178-187)                  W/m**2
    !                                        192                         LW top cloud forcing (179-188)                  W/m**2
#endif
    CALL add (g3b,'wl',      wl        ,code=193              ,longname='skin reservoir content'                 ,units='m'        )
    CALL add (g3b,'slf',     slf       ,code=194,lpost=.FALSE.,longname='sea land fraction'                                        )
#ifndef MESSY
    CALL add (g3b,'ustrgw',  ustrgw    ,code=195,lpost=.FALSE.,longname='u-gravity wave stress'                  ,units='Pa'       )
    CALL add (g3b,'vstrgw',  vstrgw    ,code=196,lpost=.FALSE.,longname='v-gravity wave stress'                  ,units='Pa'       )
    CALL add (g3b,'vdisgw',  vdisgw    ,code=197              ,longname='gravity wave dissipation'               ,units='W/m**2'   )
#endif
    CALL add (g3b,'vgrat',   vgrat     ,code=198,lpost=.FALSE.,longname='vegetation ratio'                                         )
    CALL add (g3b,'orostd',  orostd    ,code=199,lpost=.FALSE.,longname='orographic standard deviation'          ,units='m'        )
    CALL add (g3b,'vlt',     vlt       ,code=200,lpost=.FALSE.,longname='leaf area index'                                          )
#ifndef MESSY
    CALL add (g3b,'t2max',   t2max     ,code=201,reset=-99.0_dp   ,longname='maximum 2m temperature',units='K'     )
    CALL add (g3b,'t2min',   t2min     ,code=202,reset=999.0_dp   ,longname='minimum 2m temperature',units='K'     )
    CALL add (g3b,'srad0u',  srad0u    ,code=203,laccu=.TRUE. ,longname='top solar radiation upward'             ,units='W/m**2'   )
    CALL add (g3b,'sradsu',  sradsu    ,code=204,laccu=.TRUE. ,longname='surface solar radiation upward'         ,units='W/m**2'   )
    CALL add (g3b,'tradsu',  tradsu    ,code=205,laccu=.TRUE. ,longname='surface thermal radiation upward'       ,units='W/m**2'   )
#endif
    CALL add (g3b,'grndflux',grndflux  ,code=206,laccu=.TRUE. ,longname='surface ground heat flux'               ,units='W/m**2'   )
    CALL add (g3b,'tsoil',   tsoil     ,code=207              ,longname='deep soil temperatures',leveltype=BELOWSUR,units='K'      )
    CALL add (g3b,'ahfcon',  ahfcon    ,code=208,laccu=.TRUE. ,longname='conductive heat flux through ice'       ,units='W/m**2'   )
    CALL add (g3b,'ahfres',  ahfres    ,code=209,laccu=.TRUE. ,longname='melting of ice'                         ,units='W/m**2'   )
    CALL add (g3b,'seaice',  seaice    ,code=210              ,longname='ice cover (fraction of 1-SLM)'                            )
    CALL add (g3b,'siced',   siced     ,code=211              ,longname='ice depth'                              ,units='m'        )
    CALL add (g3b,'forest',  forest    ,code=212,lpost=.FALSE.,longname='forest fraction'                                          )
    CALL add (g3b,'gld',     gld       ,code=213              ,longname='glacier depth'                          ,units='m'        )
    CALL add (g3b,'sni',     sni       ,code=214              ,longname='water equivalent of snow on ice'        ,units='m'        )
    CALL add (g3b,'rogl',    rogl      ,code=215,laccu=.TRUE.,lpost=.FALSE. ,longname='glacier runoff'           ,units='kg/m**2s' )
#ifndef MESSY
    CALL add (g3b,'wimax',   wimax     ,code=216,reset=-99.0_dp   ,longname='maximum 10m-wind speed',units='m/s'  )
#endif
    CALL add (g3b,'topmax',  topmax    ,code=217,reset=99999.0_dp ,longname='maximum height of convective cloud tops',units='Pa'   )
    CALL add (g3b,'snmel',   snmel     ,code=218,laccu=.TRUE. ,longname='snow melt'                              ,units='kg/m**2s' )
#ifndef MESSY
    CALL add (g3b,'runtoc',  runtoc    ,code=219,laccu=.TRUE.,lpost=.FALSE. ,                &
                 longname='surface runoff into ocean',bits=24,units='kg/m**2s' )
    CALL add (g3b,'runlnd',  runlnd    ,code=220,lpost=.FALSE.,longname='surface runnof not running into ocean'  ,units='kg/m**2s' )
#endif
    CALL add (g3b,'apmegl',  apmegl    ,code=221,laccu=.TRUE. ,longname='P-E over land ice'              ,bits=24,units='kg/m**2s' )
    CALL add (g3b,'snacl',   snacl     ,code=222,laccu=.TRUE. ,longname='snow accumulation over land'            ,units='kg/m**2s' )
    CALL add (g3b,'aclcac',  aclcac    ,code=223,laccu=.TRUE. ,longname='cloud cover'                                              )
#ifndef MESSY
    CALL add (g3b,'tke',     tke       ,code=224,lpost=.false.,longname='turbulent kinetic energy'               ,units='m**2/s**2')
    CALL add (g3b,'tkem1',   tkem1     ,code=225,lpost=.false.,longname='turbulent kinetic energy (t-1)'         ,units='m**2/s**2')
#endif
    CALL add (g3b,'fao',     fao       ,code=226,lpost=.false.,longname='FAO data set (soil data flags 0...5.)'                    )
    CALL add (g3b,'rgcgn',   rgcgn     ,code=227,lpost=.false.,longname='volumetric heat capacity of soil and land ice'            &
                                                                                                                ,units='J/(m**3*K)')
    CALL add (g3b,'sodif',   sodif     ,code=228,lpost=.false.,longname='diffusivity  of soil and land ice'      ,units='m**2/s'   )
    CALL add (g3b,'wsmx',    wsmx      ,code=229              ,longname='field capacity of soil'                 ,units='m'        )
    CALL add (g3b,'qvi',     qvi       ,code=230,laccu=.TRUE. ,longname='vertically integrated water vapor'      ,units='kg/m**2'  )
    CALL add (g3b,'xlvi',    xlvi      ,code=231,laccu=.TRUE. ,longname='vertically integrated cloud water'      ,units='kg/m**2'  )
    CALL add (g3b,'glac',    glac      ,code=232              ,longname='fraction of land covered by glaciers'                     )
    CALL add (g3b,'snc',     snc       ,code=233              ,longname='snow depth at the canopy'               ,units='m'        )
#ifndef MESSY
    CALL add (g3b,'rtype',   rtype     ,code=234,lpost=.FALSE.,longname='type of convection 0...3.'                                )
#endif
! op_pj_20160613+
!!$    ! mz_ho_20160412+
!!$#ifdef MESSY
!!$    CALL add (g3b,'az0hi',    az0hi    ,code=238              ,longname='roughness length for heat over ice'     ,units='m'        )
!!$    CALL add (g3b,'az0hw',    az0hw    ,code=239              ,longname='roughness length for heat over water'   ,units='m'        )
!!$    CALL add (g3b,'az0hl',    az0hl    ,code=240              ,longname='roughness length for heat over land'    ,units='m'        )
!!$#endif
!!$    ! mz_ho_20160412-
! op_pj_20160613-
    !                                        259                         windspeed (sqrt(u**2+v**2))
    !                                        260                         total precipitation (142+143)
    !                                        261                         total top radiation (178+179)
    !                                        262                         total surface radiation (176+177)
    !                                        263                         net surface heat flux 
    !                                                                    (146+147+176+177-C*218-208*fice-209); C=3.345E5*fland
    !                                        264                         total surface water (142+143+182-160-221)

    ! Add fields not written to the output stream

    CALL add (g3b,'tsl',     tsl     ,lpost=.FALSE.)
    CALL add (g3b,'tslm',    tslm    ,lpost=.FALSE.)
#ifndef MESSY
    CALL add (g3b,'emter',   emter   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'trsol',   trsol   ,lpost=.FALSE. ,leveltype=HYBRID_H)
    CALL add (g3b,'emtef0',  emtef0  ,lpost=.FALSE. ,leveltype=HYBRID_H,contnorest=.true.)
    CALL add (g3b,'trsof0',  trsof0  ,lpost=.FALSE. ,leveltype=HYBRID_H,contnorest=.true.)
    CALL add (g3b,'emtef',   emtef   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'trsof',   trsof   ,lpost=.FALSE. ,klev=2)
    CALL add (g3b,'tkem',    tkem    ,lpost=.FALSE.)
#endif
    CALL add (g3b,'grndcapc',grndcapc,lpost=.FALSE.)
    CALL add (g3b,'grndhflx',grndhflx,lpost=.FALSE.)
    CALL add (g3b,'grndc',   grndc   ,lpost=.FALSE. ,leveltype=BELOWSUR)
    CALL add (g3b,'grndd',   grndd   ,lpost=.FALSE. ,leveltype=BELOWSUR)
#ifndef MESSY
    CALL add (g3b,'acvtype', acvtype ,lpost=.FALSE.)
#endif
    CALL add (g3b,'xtec',    xtec    ,lpost=.FALSE.)
    !
    ! variables for middle atmosphere only
    !
    IF (lmidatm) &
      CALL add (g3b,'aprflux',  aprflux ,lpost=.FALSE.)
    !
    !  variables for fractional surface coverage
    !
#ifndef MESSY
    CALL add (g3b,'ahfli',    ahfli    ,lpost=.FALSE.)
    CALL add (g3b,'ahflw',    ahflw    ,lpost=.FALSE.)
    CALL add (g3b,'ahfll',    ahfll    ,lpost=.FALSE.)
#endif
    CALL add (g3b,'fluxres',  fluxres  ,lpost=.FALSE.)
    !
    !  variables for mixed layer ocean only
    !
#ifndef MESSY
    CALL add (g3b,'amlcorr',  amlcorr  ,lpost=.FALSE.)
    CALL add (g3b,'amlcorac', amlcorac ,code= 89,laccu=.TRUE.)
    CALL add (g3b,'amlheatac',amlheatac,code= 90,laccu=.TRUE.,lpost=.FALSE.)
#endif
    !
    !  variables for ocean coupling only
    !
    CALL add (g3b,'apmebco',  apmebco  ,lpost=.FALSE.,contnorest=.true.)
    !
#ifndef MESSY
    IF (lcouple) THEN
      CALL add (g3b,'awhea',    awhea   ,lpost=.FALSE.)
      CALL add (g3b,'awsol',    awsol   ,lpost=.FALSE.)
      CALL add (g3b,'awfre',    awfre   ,lpost=.FALSE.)
      CALL add (g3b,'awust',    awust   ,lpost=.FALSE.)
      CALL add (g3b,'awvst',    awvst   ,lpost=.FALSE.)
      CALL add (g3b,'awsta',    awsta   ,lpost=.FALSE.)
      CALL add (g3b,'aicon',    aicon   ,lpost=.FALSE.)
      CALL add (g3b,'aiqre',    aiqre   ,lpost=.FALSE.)
      CALL add (g3b,'aifre',    aifre   ,lpost=.FALSE.)
      CALL add (g3b,'aiust',    aiust   ,lpost=.FALSE.)
      CALL add (g3b,'aivst',    aivst   ,lpost=.FALSE.)
    END IF
#endif
    CALL add (g3b,'ocu',  ocu    ,code= 83,                                             & ! 0 if MPIOM is not used
         longname='ocean eastw. velocity',bits=24,units='m/s',contnorest=.true. )
    CALL add (g3b,'ocv',  ocv    ,code= 84,                                             & ! 0 if MPIOM is not used
         longname='ocean northw. velocity',bits=24,units='m/s',contnorest=.true. )
    !
    !  variables for 200mb radiation
    !
#ifndef MESSY
    CALL add (g3b,'tradl',   tradl     ,code= 85,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'sradl',   sradl     ,code= 86,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'trafl',   trafl     ,code= 87,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
    CALL add (g3b,'srafl',   srafl     ,code= 88,laccu=.TRUE.,units='W/m**2',contnorest=.true.   )
#endif
    !
    !  variables for coupling with HD-model only
    !
#ifndef MESSY    
!mz_ap_20070920+
! defined in MESSY
    IF (lhd) THEN
      CALL add (g3b,'aros',     aros    ,lpost=.FALSE.)
      CALL add (g3b,'adrain',   adrain  ,lpost=.FALSE.)
      CALL add (g3b,'disch',    disch   ,lpost=.FALSE.)
      CALL add (g3b,'apmecal',  apmecal ,lpost=.FALSE.)
    END IF
!mz_ap_20070920-
#endif
    !
    !  variables for sso parametrization
    !
    CALL add (g3b,'oromea', oromea,lpost=.FALSE.)
    CALL add (g3b,'orosig', orosig,lpost=.FALSE.)
    CALL add (g3b,'orogam', orogam,lpost=.FALSE.)
    CALL add (g3b,'orothe', orothe,lpost=.FALSE.)
    CALL add (g3b,'oropic', oropic,lpost=.FALSE.)
    CALL add (g3b,'oroval', oroval,lpost=.FALSE.)

  END SUBROUTINE construct_g3b

  SUBROUTINE destruct_g3b

    CALL delete_stream (g3b)

  END SUBROUTINE destruct_g3b

END MODULE mo_memory_g3b
