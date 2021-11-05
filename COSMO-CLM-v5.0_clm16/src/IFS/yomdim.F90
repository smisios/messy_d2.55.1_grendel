MODULE YOMDIM

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Dimensions of model working arrays

! === COLLOCATION GRID OF THE DYNAMICS ===

! NDGLG  : number of rows of latitudes
! NDGLL  : number of rows of latitudes for which this process is
!          performing Fourier Space calculations
! NDGNH  : number of rows in the northern hemisphere
! NDGSUR : number of additional rows at each pole for horizontal
!          interpolations.
! NDGSAG = -NDGSUR+1
! NDGSAL = Local version of NDGSAG.
! NDGSAH = 1-NSLWIDE in DM version.
! NDGSAFPH=1-NFPWIDE in DM version.
! NDGENG = NDGLG+NDGSUR
! NDGENL = Number of latitude rows for which this process has grid
!          point calculations to perform.
! NDGENH = NDGENL+NSLWIDE in DM version.
! NDGENFPH=NDGENL+NFPWIDE in DM version.
! NDGUNG : first row of the area of interest in Aladin
!        = NDGSAG in the global model
! NDGUXG : last  row of the area of interest in Aladin
!        = NDGENG in the global model
! NDGUNL : local first row in C+I zone in distributed memory Aladin
! NDGUXL : local last row in C+I zone in distributed memory Aladin

! NDLON  : length of a row of latitude near equator
! NDSUR1 : over dimensioning of NDLON for technical reasons (at least 2)
! NDLSUR = NDLON+NDSUR1
! NDLSM  = NDLSUR-1
! NDLUNG : first meridian of the area of interest in Aladin
!        = 1 in the global model
! NDLUXG : last  meridian of the area of interest in Aladin
!        = NDLON in the global model
! NDLUNL : local first meridian in C+I zone in distributed memory Aladin
! NDLUXL : local last meridian in C+I zone in distributed memory Aladin
! NPROMA : working dimension for grid-point computations
! NPROMB : working dimension for a 2nd call to DMN physics
! NPROMC : working dimension for a 2nd call to DMN physics
! NPROMP : working dimension for physics computations
! NPROME : working dimension for ECMWF physics computations
! NPROMM : working dimension for DMN   physics computations
! NPROMV : working dimension for variational gridpoint fields
! NPROMNH: working dimension for non hydrostatic
! NPROMVC: working dimension for some additional grid-point arrays
!          used only when "lvercor=.T.".
! NGPBLKS: number of grid point blocks, i.e. the number of calls
!          to cpg and cpglag per scan2mdm call.It is only used in
!          DM version and defined as (NGPTOT-1)/NPROMA+1
! LOPTPROMA : .TRUE. NPROMA will be optimised
!           : .FALSE. NPROMA will not be optimised (forced by
!           : negative NPROMA in namelist)

! === NUMBER OF FIELDS ===

! NFLEVG : number of levels in grid point space
! NFLEVL : number of levels in Fourier and Legendre space
! NFLEVLMX : maximum NFLEVL among all PEs
! NFTHER : number of spectral thermodynamic variables
! NFPASS : number of spectral passive scalar variables
! NFAUX  : number of auxillary variables in t+dt array
! NF3D   = number of 3D fields in the state of the model
! NFD2D  : number of 2d fields in the dynamics
! NFC2D  : number of 2d fields in the boundaries
! NFLSUR : over dimensioning of NFLEVL for technical reasons, always odd
! NFLSUL : number of additional levels for semi-lagrangian
! NFLSA  = 1    -NFLSUL
! NFLEN  = NFLEVG+NFLSUL
! NFPP3M : maximum number of extra 3-D fields in post-processing
! NPPM   : Number of interpolation methods in post-processing
! NFPPYX : Maximum number of modern dyn.met. post-processed fields
! NFPPYE : as long as SUPP is called after SUDIM/SUALLO,
!            NFPPYE = 1, if LMDYPP = .FALSE. (after SUALLO)
!            NFPPYE = NFPPYX if LMDYPP = .TRUE., but it could become more flex
! NFGPUA : number of (3D) fields in upper air gridpoint
! NFTHERG: number of thermodynamics variables in upper air gridpoint
! NFPASSG: number of passive scalar variables in upper air gridpoint
! NFGPNH : number of (3D) fields in non hydrostatic gridpoint
! NFTC1  : number of 'NFLSUR' transmission coefficients.
! NFTC2  : number of '(0:NFLEVG)*(0:NFLEVG)' transmission coefficients.
! NGTC   : equivalent total number of 2D fields to read or write or
!          transform for transmission coefficients.
! NFLPFLPSUR: odd number overdimension for (NFLEVG+1)*(NFLEVG+1).

! === GRID POINT ARRAYS ===

! NG2D0 :  number of 2D fields in the  t    array
! NG2D1 :  number of 2D fields in the t+dt array
! NG2D9 :  number of 2D fields in the t-dt array
! NG3D0 :  number of 3D fields in the  t   array
! NG3D1 :  number of 3D fields for the t+dt array dimensioning
! NG3D9 :  number of 3D fields in the t-dt array
! NG2C  :  number of 2D fields in the boundary array
! NG2D5 :  number of 2D fields in the trajectory array
! NG3D5 :  number of 3D fields in the  trajectory array
! NGT0   : maximum number of fields in t array
! NGT1   : maximum number of fields in t+dt array
! LVOR  : controls the allocation of vorticity
! LGRAQ : controls the allocation of thermodynamic variables other than T
! LADER : controls the allocation of vor div and derivatives
! LUVDER: controls the allocation of derivatives for u and v

! LSPT  : .TRUE. if temperature variable as spectral field
! LSPQ  : .TRUE. if q as spectral field
! LSPO3 : .TRUE. for ozone in data assimilation (need spectral ozone there)
! LSPL  : .TRUE. if liquid water spectral field
! LSPI  : .TRUE. if ice water as spectral field
! LSPA  : .TRUE. if cloud fraction as spectral field

! LGPUA : .TRUE. if any non-spectral upper air fields
! LGPQ  : .TRUE. if q as upper air grid-point field
! LGPO3 : .TRUE. for ozone as prognostic variable (only grid-point field) (ECMWF)
! LGPL  : .TRUE. if liquid water as upper air grid-point field
! LGPI  : .TRUE. if ice water as upper air grid-point field
! LGPA  : .TRUE. if cloud fraction as upper air grid-point field
! LGPE  : .TRUE. if TKE as upper air grid-point field

! LGPQIN : .TRUE. if Q grid-point upper air fields to be read on input
! LCLDPIN: .TRUE. if (L,I,A) grid-point upper air fields to be read on input

! === SPECTRAL SPACE (TRANSFORMED SPHERE) ===

! NSMAX  : truncation order
! NMSMAX  : truncation order in longitude
! NVARMAX: truncation order in 3d-var distributed direction
!          this is a priori longitude, so that nvarmax = nsmax in Arp/IFS
!          and nvarmax = nmsmax in Aladin
! NSEFRE : number of degrees of freedom in the spectral space
! NSPECG : number of complex spectral coefficients (global)
! NSPEC2G = 2*NSPECG
! NSPEC  : number of complex spectral coefficients (local, i.e. on
!          this PE)
! NSPEC2 = 2*NSPEC
! NSPEC2MX : maximun NSPEC2 among all PEs
! NS3D   : number of 3D fields in spectral space
! NS2D   : number of 2D fields in spectral space
! NSA    : number of antisymmetric Laplace coefficients
! NSS    :   "           symmetric    "        "
! NSAUX  : dimension of auxillary array ( MAX(NFAUX,NFPP3M) )
! NSMIN  : lower troncature for configurations 911 and 912.
! NSPOLEG: number of Legendre polynomials
! NTCMAX : truncation order for transmission coefficients.
! NMTCMAX: truncation order for transmission coefficients
!          in longitude.

! === SPECTRAL SPACE (REAL SPHERE) ===

! NCMAX  : truncation order
! NCPEC  : number of complex spectral coefficients (local)
! NCPEC2 : 2*NCPEC, where NCPEC (local)

! === SPECTRAL SPACE (EXPLICIT NMI TRUNCATION) ===

! NXMAX  : truncation order
! NXPECG : number of complex spectral coefficients (global)
! NXPEC  : number of complex spectral coefficients (local)

! === SPECTRAL SPACE (TENDENCIES, TRANSFORMED SPHERE) ===

! NTMAX  : truncation order for tendencies (on n, m<= NSMAX)
! NTPEC2 : 2*'number of complex spectral coefficients'

! === LEGENDRE TRANSFORM WORKING SPACE (INVERSE TRANSFORM) ===

! NLEI1  = first dimension of input array PIA
! NLEI2  = second    "          "     "    "   and PxOA1
! NLEI3  = first dimension of Legendre polynomials and PxOA1
! NLEI4  = second dimension of Legendre polynomials
! NLTI6  = offset for PIA and PxOA1
! NLTI7  = depth of work for MXMA

! === LEGENDRE TRANSFORM WORKING SPACE (DIRECT TRANSFORM) ===

! NLED1  = first dimension of arrays PSIA and PAIA, second of PLEPO
! NLED2  = second dimension of arrays PSIA, PAIA and POA1
! NLED3  = first dimension of PLEPO
! NLED4  = first dimension of POA1 and POA2
! NLED5  = second dimension POA2
! NLTD1N = 2*NFLEVL*(2+NFTHER+NFPASS+NFAUX)+2  (FOR COMPUTATIONS)

! === DISTRIBUTED MEMORY DIMENSIONS ===

! NUMP  :  Number of spectral waves handled by this processor
!          If LMESSP=.F. then NUMP=NSMAX+1
! NUMXP :  Same as NUMP, but related to NXMAX
! NUMCP :  Same as NUMP, but related to NCMAX
! NUMTP :  Same as NUMP, but related to NTMAX

! === OTHER QUANTITIES ===

! NRLEVX: USED TO DIMENSION NVAUTF IN YOMGEM.
! NDIMDHU:DIMENSION FOR LOCAL ARRAYS IN SPC (UNIFIED HORIZONTAL DIFFUSION).
! NUNDEFLD: index value for unused/undefined fields (default=-9999)
!         : should be set to 1 when using compiler subscript checking

! === POSITIONING IN THE BOUNDARY ARRAY ===


!  Position    Pointee         Description
!     0        GCS             generic
!     1        RCORI           CORIOLIS parameter
!     2        GEMU            Sine of latitude on the real earth
!     3        GSQM2           Cosine of latitude on the real earth
!     4        GELAM           longitude on the real earth
!     5        GELAT           latitude on the real earth
!     6        GECLO           cosinus of the longitude on the real earth
!     7        GESLO             sinus of the longitude on the real earth
!     8        GM              mapping factor
!     9        GOMVRL          2 omega vect. R E-0 component
!    10        GOMVRM          ................N-S..........
!    11        GNORDL          Real North E-O component
!    12        GNORDM          ...........N-S..........
!    13        OROG            orography
!    14        OROGL           orography  E-O Derivative
!    15        OROGM           orography  N-S Derivative
!    16        OROGLL          orography  E-O Second Order Derivative
!    17        OROGMM          orography  N-S Second Order Derivative
!    18        OROGLM          orography  Mixed Second Order Derivative

! === Positioning in the spectral array ===


! Position     Pointee         Description
!    0         SPA3            generic
!    1         SPVOR           vorticity
!    2         SPDIV           divergence
!    3         SPHV            generic name for thermodynamic variables
!    4         SPT             temperature
!    5         SPQ             vapour phase
!    6         SPO3            ozone
!    7         SPW             liquid phase
!    8         SPS             solid phase
!    9         SPSV            passive variables


! Position     Pointee         Description
!    0         SPA2            generic
!    1         SPD2D           2d dynamical fields
!    2         SPSP            surface pressure
!    3         SPC2D           2d boundaries
!    4         SPOR            orography


! Position     Pointee         Description
!    0         SPAUX           generic
!    1         SPR             energy or post-processed variables


INTEGER(KIND=JPIM),ALLOCATABLE:: NDLUNL(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE:: NDLUXL(:,:)


INTEGER(KIND=JPIM) :: NDGLG
INTEGER(KIND=JPIM) :: NDGLL
INTEGER(KIND=JPIM) :: NDGNH
INTEGER(KIND=JPIM) :: NDGSUR
INTEGER(KIND=JPIM) :: NDGSAG
INTEGER(KIND=JPIM) :: NDGSAL
INTEGER(KIND=JPIM) :: NDGENG
INTEGER(KIND=JPIM) :: NDGENL
INTEGER(KIND=JPIM) :: NDLON
INTEGER(KIND=JPIM) :: NDSUR1
INTEGER(KIND=JPIM) :: NDLSUR
INTEGER(KIND=JPIM) :: NDLSM
INTEGER(KIND=JPIM) :: NFLEVL
INTEGER(KIND=JPIM) :: NFLEVLMX
INTEGER(KIND=JPIM) :: NFLEVG
INTEGER(KIND=JPIM) :: NFTHER
INTEGER(KIND=JPIM) :: NFPASS
INTEGER(KIND=JPIM) :: NFAUX
INTEGER(KIND=JPIM) :: NF3D
INTEGER(KIND=JPIM) :: NDGSAH
INTEGER(KIND=JPIM) :: NDGENH
INTEGER(KIND=JPIM) :: NDGSAFPH
INTEGER(KIND=JPIM) :: NDGENFPH
INTEGER(KIND=JPIM) :: NFD2D
INTEGER(KIND=JPIM) :: NFC2D
INTEGER(KIND=JPIM) :: NPROMA
INTEGER(KIND=JPIM) :: NFGPUA
INTEGER(KIND=JPIM) :: NFTHERG
INTEGER(KIND=JPIM) :: NFPASSG
INTEGER(KIND=JPIM) :: NFLSUR
INTEGER(KIND=JPIM) :: NFGPNH
INTEGER(KIND=JPIM) :: NFLSUL
INTEGER(KIND=JPIM) :: NFLSA
INTEGER(KIND=JPIM) :: NFLEN
INTEGER(KIND=JPIM) :: NFPP3M
INTEGER(KIND=JPIM) :: NFLPFLPSUR
INTEGER(KIND=JPIM) :: NFTC1
INTEGER(KIND=JPIM) :: NFTC2
INTEGER(KIND=JPIM) :: NGTC
INTEGER(KIND=JPIM) :: NPPM
INTEGER(KIND=JPIM) :: NG2D0
INTEGER(KIND=JPIM) :: NG2D1
INTEGER(KIND=JPIM) :: NG2D9
INTEGER(KIND=JPIM) :: NG3D0
INTEGER(KIND=JPIM) :: NG3D1
INTEGER(KIND=JPIM) :: NG3D9
INTEGER(KIND=JPIM) :: NG2C
INTEGER(KIND=JPIM) :: NG2D5
INTEGER(KIND=JPIM) :: NG3D5
INTEGER(KIND=JPIM) :: NG3L0
INTEGER(KIND=JPIM) :: NG3L9
INTEGER(KIND=JPIM) :: NGT0
INTEGER(KIND=JPIM) :: NGT1
INTEGER(KIND=JPIM) :: NPROMB
INTEGER(KIND=JPIM) :: NPROMC
INTEGER(KIND=JPIM) :: NPROMP
INTEGER(KIND=JPIM) :: NPROME
INTEGER(KIND=JPIM) :: NPROMM
INTEGER(KIND=JPIM) :: NPROMV
INTEGER(KIND=JPIM) :: NPROMNH
INTEGER(KIND=JPIM) :: NPROMVC
INTEGER(KIND=JPIM) :: NVARMAX
INTEGER(KIND=JPIM) :: NMSMAX
INTEGER(KIND=JPIM) :: NSMAX
INTEGER(KIND=JPIM) :: NSEFRE
INTEGER(KIND=JPIM) :: NSPECG
INTEGER(KIND=JPIM) :: NSPEC
INTEGER(KIND=JPIM) :: NSPEC2G
INTEGER(KIND=JPIM) :: NSPEC2
INTEGER(KIND=JPIM) :: NSPEC2MX
INTEGER(KIND=JPIM) :: NS3D
INTEGER(KIND=JPIM) :: NS2D
INTEGER(KIND=JPIM) :: NSA
INTEGER(KIND=JPIM) :: NSS
INTEGER(KIND=JPIM) :: NSAUX
INTEGER(KIND=JPIM) :: NSMIN
INTEGER(KIND=JPIM) :: NSPOLEG
INTEGER(KIND=JPIM) :: NPMAX
INTEGER(KIND=JPIM) :: NCMAX
INTEGER(KIND=JPIM) :: NTCMAX
INTEGER(KIND=JPIM) :: NMTCMAX
INTEGER(KIND=JPIM) :: NDGUNL
INTEGER(KIND=JPIM) :: NDGUXL
INTEGER(KIND=JPIM) :: NDGUNG
INTEGER(KIND=JPIM) :: NDGUXG
INTEGER(KIND=JPIM) :: NDLUNG
INTEGER(KIND=JPIM) :: NDLUXG
INTEGER(KIND=JPIM) :: NCPEC
INTEGER(KIND=JPIM) :: NCPEC2
INTEGER(KIND=JPIM) :: NXMAX
INTEGER(KIND=JPIM) :: NXPECG
INTEGER(KIND=JPIM) :: NXPEC
INTEGER(KIND=JPIM) :: NTMAX
INTEGER(KIND=JPIM) :: NTPEC2
INTEGER(KIND=JPIM) :: NLEI1
INTEGER(KIND=JPIM) :: NLEI2
INTEGER(KIND=JPIM) :: NLEI3
INTEGER(KIND=JPIM) :: NLEI4
INTEGER(KIND=JPIM) :: NLEI6
INTEGER(KIND=JPIM) :: NLEI7
INTEGER(KIND=JPIM) :: NLED1
INTEGER(KIND=JPIM) :: NLED2
INTEGER(KIND=JPIM) :: NLED3
INTEGER(KIND=JPIM) :: NLED4
INTEGER(KIND=JPIM) :: NLED5
INTEGER(KIND=JPIM) :: NLTD1N
INTEGER(KIND=JPIM) :: NFPPYX
INTEGER(KIND=JPIM) :: NFPPYE
INTEGER(KIND=JPIM) :: NRLEVX
INTEGER(KIND=JPIM) :: NDIMDHU
INTEGER(KIND=JPIM) :: NUMP
INTEGER(KIND=JPIM) :: NUMXP
INTEGER(KIND=JPIM) :: NUMCP
INTEGER(KIND=JPIM) :: NUMTP
INTEGER(KIND=JPIM) :: NGPBLKS
INTEGER(KIND=JPIM) :: NUNDEFLD
LOGICAL LVOR
LOGICAL LGRAQ
LOGICAL LADER
LOGICAL LUVDER
LOGICAL LGPUA
LOGICAL LGPQ
LOGICAL LGPO3
LOGICAL LGPL
LOGICAL LGPI
LOGICAL LGPA
LOGICAL LGPE
LOGICAL LSPT
LOGICAL LSPQ
LOGICAL LSPO3
LOGICAL LSPL
LOGICAL LSPI
LOGICAL LSPA
LOGICAL LGPQIN
LOGICAL LCLDPIN
LOGICAL LOPTPROMA

!     ------------------------------------------------------------------
END MODULE YOMDIM
