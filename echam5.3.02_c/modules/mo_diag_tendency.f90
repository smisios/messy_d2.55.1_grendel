#ifndef MESSY


MODULE mo_diag_tendency
!BOP
  ! !MODULE: mo_diag_tendency

  ! !DESCRIPTION:
  !\begin{verbatim}
  !-------------------------------------------------------------------
  ! overview about the separation of terms in diagnostic equations
  !
  !  LTDIAG      additional switch in &RUNCTL
  !              .true. ==> addional diagnostics
  !
  ! ***************** Interface to ECHAM *************************
  !
  ! *ECHAM*       *DIAGNOSTICS*
  !
  ! CONTROL -+
  !          +--> DIAG_Init(1)                general initialization
  !                  |
  !                  +---> DIAG_Rerun('R')    read rerun data
  !
  ! STEPON --+
  !          +--> DIAG_SumAll                 count total tendency
  !          |
  !          +--> DIAG_Write                  store data in model output
  !          |       |
  !          |       +---> DIAG_Check
  !          |       |
  !          |       +---> DiagWriteSP
  !          |
  !          +--> DIAG_Init(-1)               clear diagnostic memory
  !                  |
  !                  +---> DIAG_Rerun('W')    store rerun data
  !
  ! SCAN1SL -+
  !          +--> DIAG_fftd
  !          |
  !          SI2
  !          |
  !          +--> DIAG_SpecTrans
  !                  |
  !                  +--> DIAG_sym1
  !                  |
  !                  +--> DIAG_ltd
  !
  !---------------------------------------------------------------------
  ! count terms in additional arrays
  !   *spectral*    *Gaussian*
  !
  !                 DYN            adiabatic terms
  !                  |
  !                 TF2            time filter
  !                  |
  !                 TF1            time filter
  !                  |
  !                 GPC
  !                  +--> PHYSC    diabatic terms
  !                  |      |
  !                  |      +---> M_RADHEAT long/short wave
  !                  |
  !                  +--> SI1      semi-implicit terms
  !                  |
  !                 DIAG_fftd
  !                 |
  !               SI2                semi-implicit terms
  !               |
  !             DIAG_SpecTrans
  !                     |
  !                   DIAG_sym1
  !                     |
  !                   DIAG_ltd
  ! SCCD                           semi-implicit terms
  !   |
  ! SCCTP                          semi-implicit terms
  !   |
  ! HDIFF                          horizontal diffusion
  !   |
  ! DIAG_SumAll
  !   |
  ! DIAG_Write
  !
  !**************************************************************
  !\end{verbatim}

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI, October 1998, basic version
  ! I. Kirchner, MPI, May-2000, patch E4v3.22p1
  ! I. Kirchner, MPI, May-2002, revision E5R1.04

  ! !USES:
  USE mo_exception,     ONLY: finish, message
  USE mo_linked_list,   ONLY: t_stream
  USE mo_kind,          ONLY: dp
!BOX
  IMPLICIT NONE

  PRIVATE
!EOX
  PUBLIC  :: DIAG_Init      ! allocate/deallocate memory, initialize diagnostic arrays
  PUBLIC  :: DIAG_SpecTrans ! second part of spectral transform
  PUBLIC  :: DIAG_fftd      ! Fourier transformation of diagnostics arrays
  PRIVATE :: DIAG_sym1      ! composition of symetric/antisymetric part
  PRIVATE :: DIAG_ltd       ! legendre transform of diagnostic arrays
  PUBLIC  :: DIAG_Write     ! store diagnostic fields
  PRIVATE :: DIAG_Check     ! global check of diagnostic arrays
!BOX
  ! special diagnostic fields

  REAL(dp), ALLOCATABLE, TARGET, PUBLIC :: &    ! ARRAYS in spectral space
       pdvor (:,:,:,:), &  !   terms for vorticity equation
       pddiv (:,:,:,:), &  !   terms for divergence equation
       pdtem (:,:,:,:)     !   terms for temperature equation
  REAL(dp), POINTER, PUBLIC :: &
       pdprs (:,:,:), &    !   terms for surface pressure equation
       pdprl (:,:,:), &    !   surface pressure term vertical integration
       p4prs (:,:,:,:)

  REAL(dp), ALLOCATABLE, TARGET :: &    ! old value for accumulation of tendency
       pdovor (:,:,:,:), pdodiv (:,:,:,:), pdotep (:,:,:,:)


  INTEGER, PARAMETER :: NO_PTFH = 3
  REAL(dp), POINTER, PUBLIC ::   &  ! ARRAYS in grid point space
       ptfh1 (:,:,:,:), &  !   time integration parts
       ptfh2 (:,:)
  LOGICAL, SAVE, PUBLIC :: lset_fh1    = .FALSE.

  INTEGER, PARAMETER :: NO_PDIGA=25  ! no. of 3-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDSGA=2   ! no. of 2-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDIGB=10  ! no. of 3-dim fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS=4   ! no. of 2-dim mixed grid point and fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS_FRAC=3   ! no. of 2-dim grid point space fields fraction

  REAL(dp), POINTER, PUBLIC :: pdiga(:,:,:,:), pdsga (:,:,:)
  REAL(dp), ALLOCATABLE, PUBLIC, TARGET :: &    ! ARRAYS for accumulation
                        pdigaa(:,:,:,:), pdigas(:,:,:,:), &
                        pdsgaa(:,:,:),   pdsgas(:,:,:),   &
       pdigb (:,:,:,:), pdigba(:,:,:,:), pdigbs(:,:,:,:), &
       pdigs (:,:,:,:), pdigsa(:,:,:,:), pdigss(:,:,:,:), &
       pdsgs (:,:),     pdsgsa(:,:),     pdsgss(:,:)

!EOX
  TYPE(t_stream), POINTER :: tdiag
!BOX
  INTEGER, PUBLIC         :: dio_index ! event index of the TDIAG stream

  INTEGER, PARAMETER, PUBLIC :: &  ! number of terms in different equations
       NDVOR = 9, &! vorticity
       NDDIV = 9, &! divergence
       NDTEM =14, &! temperature
       NDPRS = 4   ! surface pressure

  INTEGER, PARAMETER, PUBLIC :: &  ! function separation
       IDIAG_INIT     = 1, &
       IDIAG_INI_GBUF = 2, &
       IDIAG_INI_PREV = 3, &
       IDIAG_FREE     = 4

  LOGICAL, SAVE :: &
       laccu_ok    = .FALSE., &! mean of total tendency is fine
       ldiag_start = .TRUE.

  CHARACTER(len=256) :: diag_mess = ''

!EOX
!EOP
CONTAINS
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Init
  ! !INTERFACE:

  SUBROUTINE DIAG_Init(itype)

    ! !DESCRIPTION: 
    ! initialisation of tendency diagnostics

    ! !USES:
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, &
                              add_stream_element,                 &
                              SPECTRAL, GAUSSIAN, GRIB,     &
                              SURFACE
    USE mo_mpi,         ONLY: p_parallel
    USE mo_grib,        ONLY: nudging_table
    USE mo_memory_sp,   ONLY: sd, svo, stp
    USE mo_control,     ONLY: nlev, nlevp1, nsp, ngl, nmp1, nlon

    INTEGER, INTENT(in) :: itype       ! function separator
!EOP
!BOC
!BOX
    REAL(dp), POINTER  :: pdum(:,:,:), p4(:,:,:,:)
    LOGICAL            :: lpdo1 = .FALSE., lpdo2 = .FALSE.
!EOX
    IF (p_parallel) &
         CALL finish('mo_diag_tendency:DIAG_Init','not prepared for parallel mode')
    SELECT CASE(itype)

    CASE(IDIAG_INIT)             ! startup initialization
!BOX
      CALL message('mo_diag_tendency:DIAG_Init','  ------- start initialization -----')
!EOX
      CALL message('','Tendency Diagnostics (E5R1.04) 22-Aug-2002 (kirchner@dkrz.de)')
!BOX
      ! define tendency diagnostic output stream

      CALL new_stream ( tdiag, 'tdiag', lpost=.TRUE., lrerun=.TRUE., filetype=GRIB)
      dio_index = tdiag%post_idx

      CALL default_stream_setting ( tdiag, &
           table=nudging_table, bits=24, units='(X)/s', contnorest=.TRUE., &
           lpost  = .TRUE., lrerun = .TRUE., laccu  = .TRUE.,              &
           repr=SPECTRAL, ldims=(/nlev, 2, nsp/), gdims=(/nlev, 2, nsp/)   )
!EOX
      ! -----------------------------------------------------------
      ! allocate spectral space arrays, these fields are stored

      ALLOCATE (pdvor (nlev,2,nsp,NDVOR)); pdvor(:,:,:,:) = 0.0 !  vorticity equation

      p4 => pdvor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'VOR1',  pdum,p4=p4, code=41, &
           longname = 'VorEQ horizontal advec + press.grad. + cori.term (DYN)')
      p4 => pdvor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'VOR2',  pdum,p4=p4, code=42, &
           longname = 'VorEQ vertical advection (DYN)')
      p4 => pdvor(:,:,:,3:3)
      CALL add_stream_element(tdiag,'VOR3',  pdum,p4=p4, code=43, &
           longname = 'VorEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pdvor(:,:,:,4:4)
      CALL add_stream_element(tdiag,'VOR4',  pdum,p4=p4, code=44, &
           longname = 'VorEQ gravity wave drag (GWDRAG)')
      p4 => pdvor(:,:,:,5:5)
      CALL add_stream_element(tdiag,'VOR5',  pdum,p4=p4, code=45, &
           longname = 'VorEQ moisture mass flux (CUCALL)')
      p4 => pdvor(:,:,:,6:6)
      CALL add_stream_element(tdiag,'VOR6',  pdum,p4=p4, code=46, &
           longname = 'VorEQ timefilter')
      p4 => pdvor(:,:,:,7:7)
      CALL add_stream_element(tdiag,'VOR7',  pdum,p4=p4, code=47, &
           longname = 'VorEQ semi-implicit part of time integration')
      p4 => pdvor(:,:,:,8:8)
      CALL add_stream_element(tdiag,'VOR8',  pdum,p4=p4, code=48, &
           longname = 'VorEQ horizontal diffusion')
      p4 => pdvor(:,:,:,9:9)
      CALL add_stream_element(tdiag,'VORSUM',pdum,p4=p4, code=49, &
           longname = 'VorEQ total tendency')
      
      ALLOCATE (pddiv (nlev,2,nsp,NDDIV)); pddiv(:,:,:,:) = 0.0 !  divergence equation

      p4 => pddiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'DIV1',  pdum,p4=p4, code=61, &
           longname = 'DivEQ horizontal advec + cori. + press.grad. + G-term (DYN)')
      p4 => pddiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'DIV2',  pdum,p4=p4, code=62, &
           longname = 'DivEQ vertical advection')
      p4 => pddiv(:,:,:,3:3)
      CALL add_stream_element(tdiag,'DIV3',  pdum,p4=p4, code=63, &
           longname = 'DivEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pddiv(:,:,:,4:4)
      CALL add_stream_element(tdiag,'DIV4',  pdum,p4=p4, code=64, &
           longname = 'DivEQ gravity wave drag (GWDRAG)')
      p4 => pddiv(:,:,:,5:5)
      CALL add_stream_element(tdiag,'DIV5',  pdum,p4=p4, code=65, &
           longname = 'DivEQ moisture mass flux (CUCALL)')
      p4 => pddiv(:,:,:,6:6)
      CALL add_stream_element(tdiag,'DIV6',  pdum,p4=p4, code=66, &
           longname = 'DivEQ timefilter')
      p4 => pddiv(:,:,:,7:7)
      CALL add_stream_element(tdiag,'DIV7',  pdum,p4=p4, code=67, &
           longname = 'DivEQ semi-implicit part of time integration')
      p4 => pddiv(:,:,:,8:8)
      CALL add_stream_element(tdiag,'DIV8',  pdum,p4=p4, code=68, &
           longname = 'DivEQ horizontal diffusion')
      p4 => pddiv(:,:,:,9:9)
      CALL add_stream_element(tdiag,'DIVSUM',pdum,p4=p4, code=69, &
           longname = 'DivEQ total tendency')

      ALLOCATE (pdtem (nlev,2,nsp,NDTEM)); pdtem(:,:,:,:) = 0.0_dp ! temperature

      p4 => pdtem(:,:,:,1:1)
      CALL add_stream_element(tdiag,'TEM01', pdum,p4=p4,  code=81, &
           longname = 'TempEQ horizontal advection (DYN)')
      p4 => pdtem(:,:,:,2:2)
      CALL add_stream_element(tdiag,'TEM02', pdum,p4=p4,  code=82, &
           longname = 'TempEQ vertical advection (DYN)')
      p4 => pdtem(:,:,:,3:3)
      CALL add_stream_element(tdiag,'TEM03', pdum,p4=p4,  code=83, &
           longname = 'TempEQ energy conversion (DYN)')
      p4 => pdtem(:,:,:,4:4)
      CALL add_stream_element(tdiag,'TEM04', pdum,p4=p4,  code=84, &
           longname = 'TempEQ radiation (RADHEAT)')
      p4 => pdtem(:,:,:,5:5)
      CALL add_stream_element(tdiag,'TEM05', pdum,p4=p4,  code=85, &
           longname = 'TempEQ vertical diffusion due to turbulence (VDIFF)')
      p4 => pdtem(:,:,:,6:6)
      CALL add_stream_element(tdiag,'TEM06', pdum,p4=p4,  code=86, &
           longname = 'TempEQ gravity wave drag (GWDRAG)')
      p4 => pdtem(:,:,:,7:7)
      CALL add_stream_element(tdiag,'TEM07', pdum,p4=p4,  code=87, &
           longname = 'TempEQ convection (CUCALL)')
      p4 => pdtem(:,:,:,8:8)
      CALL add_stream_element(tdiag,'TEM08', pdum,p4=p4,  code=88, &
           longname = 'TempEQ large scale cloud processes (COND)')
      p4 => pdtem(:,:,:,9:9)
      CALL add_stream_element(tdiag,'TEM09', pdum,p4=p4,  code=89, &
           longname = 'TempEQ timefilter')
      p4 => pdtem(:,:,:,10:10)
      CALL add_stream_element(tdiag,'TEM10', pdum,p4=p4, code=90, &
           longname = 'TempEQ semi-implicit part of time integration')
      p4 => pdtem(:,:,:,11:11)
      CALL add_stream_element(tdiag,'TEM11', pdum,p4=p4, code=91, &
           longname = 'TempEQ horizontal diffusion')
      p4 => pdtem(:,:,:,12:12)
      CALL add_stream_element(tdiag,'TEM12', pdum,p4=p4, code=92, &
           longname = 'TempEQ longwave radiation')
      p4 => pdtem(:,:,:,13:13)
      CALL add_stream_element(tdiag,'TEM13', pdum,p4=p4, code=93, &
           longname = 'TempEQ shortwave radiation')
      p4 => pdtem(:,:,:,14:14)
      CALL add_stream_element(tdiag,'TEMSUM',pdum,p4=p4, code=94, &
           longname = 'TempEQ total tendency')
      
      ! divergence for each layer
      CALL add_stream_element(tdiag,'PRL',   pdprl, code=100, &
           longname = 'PresEQ convergence in each layer')

      ALLOCATE (p4prs (1,2,nsp,NDPRS)); p4prs(:,:,:,:) = 0.0_dp !  pressure equation
!BOX
      pdprs => p4prs(1,:,:,:)
      CALL default_stream_setting ( tdiag, &
            table=nudging_table, bits=24, &
            lpost  = .TRUE., lrerun = .TRUE., laccu  = .TRUE., &
            repr=SPECTRAL, leveltype=SURFACE, &
            ldims=(/1, 2, nsp/), gdims=(/1, 2, nsp/),&
            units='(X)/s', no_default=.TRUE.)
!EOX
      p4 => p4prs(1:1,:,:,1:1)
      CALL add_stream_element(tdiag,'PRS1',  pdum, code=101, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ vertical integrated convergence')
      p4 => p4prs(1:1,:,:,2:2)
      CALL add_stream_element(tdiag,'PRS2',  pdum, code=102, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ timefilter')
      p4 => p4prs(1:1,:,:,3:3)
      CALL add_stream_element(tdiag,'PRS3',  pdum, code=103, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ semi-implicit part of time integration')
      p4 => p4prs(1:1,:,:,4:4)
      CALL add_stream_element(tdiag,'PRSSUM',pdum, code=104, klev=1, p4=p4, &
           ldims=(/1, 2, nsp/) ,gdims=(/1, 2, nsp/), &
           longname = 'PresEQ total step to step')
!BOX
      !-------------------------------------------------------------------------
      ! fields used for restart

      CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., lpost = .FALSE., laccu = .FALSE., &
            repr=SPECTRAL, &
            ldims=(/nlev, 2, nsp/), gdims=(/nlev, 2, nsp/), &
            units='(X)', no_default=.TRUE.)

      ALLOCATE (pdovor(nlev,2,nsp,2)); pdovor(:,:,:,:) = 0.0_dp
      p4 => pdovor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOV1',  pdum,p4=p4, &
           longname = 'PDO Vor time level 1')
      p4 => pdovor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOV2',  pdum,p4=p4, &
           longname = 'PDO Vor time level 2')

      ALLOCATE (pdodiv(nlev,2,nsp,2)); pdodiv(:,:,:,:) = 0.0_dp
      p4 => pdodiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOD1',  pdum,p4=p4, &
           longname = 'PDO Div time level 1')
      p4 => pdodiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOD2',  pdum,p4=p4, &
           longname = 'PDO Div time level 2')

      ALLOCATE (pdotep(nlevp1,2,nsp,2)); pdotep(:,:,:,:) = 0.0_dp
      p4 => pdotep(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOTP1', pdum,p4=p4, &
           longname = 'PDO STP time level 1',klev=nlevp1)
      p4 => pdotep(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOTP2', pdum,p4=p4, &
           longname = 'PDO STP time level 2',klev=nlevp1)

!EOX       
      ! -----------------------------------------------------------
      ! allocate gridpoint space arrays
!BOX
       CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., repr=GAUSSIAN, &
            ldims=(/nlon, nlev, ngl/), gdims=(/nlon, nlev, ngl/), &
            no_default=.TRUE.)
!EOX
      !
      !   PTFH1   memory of timefilter
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH1', ptfh1, &
           longname = 'time filter buffer 1 ', &
           ldims=(/nlon, nlev,NO_PTFH,ngl/) ,gdims=(/nlon, nlev,NO_PTFH,ngl/) )
!EOX
      !           1 ... vorticity
      !           2 ... divergence
      !           3 ... temperature
      !   PTFH2   pressure
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH2', ptfh2, &
           longname = 'time filter buffer 2', &
           ldims=(/nlon,ngl/) ,gdims=(/nlon, ngl/))
!EOX
      !
      ! -----------------------------------------------------------
      !
      !     g .... accumulated in grid point space
      !     gs ... parts from gridpoint space, accumulated in spectral space
      !     s .... accumulated in spectral space
      !
      ! workspace for accumulation of terms in grid point space
      !
      ! **** PDIGA
      !
      !     Index
      !     ............. vorticity and divergenc equation
      !     VOM = Fu, VOL = Fv
      !
      !     Fu (VOM parts in GPC)
      ! DYN    1   coriolisterm and pressure tendency of Fu
      ! DYN    2   vertical advection of Fu
      ! PHYSC  3   diffusion due to impuls Pu VDIFF
      ! PHYSC  4   diffusion due to gravity wave drag Pu GWDRAG
      ! PHYSC  5   diffusion due to mass flux Pu CUCALL
      !
      !     Fv (VOL parts in GPC)
      ! DYN    6   coriolisterm and pressure tendency of Fv
      ! DYN    7   vertical advection of Fv
      ! PHYSC  8   diffusion due to impuls Pv VDIFF
      ! PHYSC  9   diffusion due to gravity wave drag Pv GWDRAG
      ! PHYSC  10  diffusion due to mass flux Pv CUCALL
      !
      ! DYN    11  potential and kinetic energy term G
      !
      !     ............. temperature tendency equation
      ! DYN    12  horizontal advection term
      ! DYN    13  vertical advection term
      ! DYN    14  energy conversion
      ! PHYSC  15  RADHEAT radiation tendency
      ! PHYSC  16  VDIFF turbulence
      ! PHYSC  17  GWDRAG gravity wave drag
      ! PHYSC  18  CUCALL convective condensation
      ! PHYSC  19  COND large+scale condensation
      !
      !     ............. pressure tendency equation
      ! DYN    20 level dependend divergence part, 
      !
      ! TF2    21  timefilter of vorticity
      ! TF2    22  timefilter of divergence
      ! TF2    23  timefilter of temperature
      !
      ! RADHEAT 24 longwave radiation
      ! RADHEAT 25 shortwave radiation
      !
!BOX
      CALL add_stream_element(tdiag,'PDIGA', pdiga, &
           longname = 'diagnostic term DIGA', &
           ldims=(/nlon, nlev,NO_PDIGA,ngl/) ,gdims=(/nlon, nlev,NO_PDIGA,ngl/))
      ALLOCATE (pdigaa (nlev,2,nmp1,NO_PDIGA))
      ALLOCATE (pdigas (nlev,2,nmp1,NO_PDIGA))
!EOX
      !
      ! **** PDSGA
      !
      ! one level arrays of surface pressure
      !      1     last level total integral
      !      2     timefilter for pressure
      !
!BOX
      CALL add_stream_element(tdiag,'PDSGA', pdsga, &
           longname = 'diagnostic term DSGA', &
           ldims=(/nlon,NO_PDSGA,ngl/) ,gdims=(/nlon,NO_PDSGA,ngl/))
      ALLOCATE (pdsgaa (2,nmp1,NO_PDSGA))
      ALLOCATE (pdsgas (2,nmp1,NO_PDSGA))
!EOX
      !
      ! ***** PDIGB
      !
      !     array for d/dlambda derivatives calculated in fourierspace
      !     corresponding to pdiga(*,1...10,*)
      !
!BOX
      ALLOCATE (pdigb    (nlon,nlev,NO_PDIGB,ngl))
      ALLOCATE (pdigba (nlev,2,nmp1,NO_PDIGB))
      ALLOCATE (pdigbs (nlev,2,nmp1,NO_PDIGB))
!EOX
      !
      ! ***** PDIGS
      !
      !     local memory buffer for accumulation of semi-implicit parts
      !     fraction solved in grid point space
      !
      !    1 vorticity and implicit and explicit part (L)
      !    2 divergence
      !    3 temperatur and pressure
      !    4 implicit part of vorticity (M) used in Fourierspace
!BOX
      ALLOCATE (pdigs    (nlon,nlev,NO_PDIGS,ngl))
      ALLOCATE (pdigsa (nlev,2,nmp1,NO_PDIGS))
      ALLOCATE (pdigss (nlev,2,nmp1,NO_PDIGS))
!EOX
      !
      ! ****** PDSGS
      !
      !     semi-implicit part of log surface pressure
!BOX
      ALLOCATE (pdsgs    (nlon,ngl))
      ALLOCATE (pdsgsa (2,nmp1))
      ALLOCATE (pdsgss (2,nmp1))
      
      CALL message('mo_diag_tendency:DIAG_Init','  ------- end of initialization -----')
!EOX
    CASE (IDIAG_FREE)        ! get memory free
!BOX
      DEALLOCATE (pdvor, pddiv, pdtem, p4prs)
      DEALLOCATE (pdovor, pdodiv, pdotep)
      DEALLOCATE        (pdigaa, pdigas)
      DEALLOCATE        (pdsgaa, pdsgas)
      DEALLOCATE (pdigb, pdigba, pdigbs)
      DEALLOCATE (pdigs, pdigsa, pdigss)
      DEALLOCATE (pdsgs, pdsgsa, pdsgss)
!EOX
    CASE (IDIAG_INI_GBUF)    ! reset local accumulation buffer
!BOX
      pdigb (:,:,:,:) = 0.0_dp
      pdigs (:,:,:,:) = 0.0_dp
      pdsgs (:,:)     = 0.0_dp
!EOX
    CASE (IDIAG_INI_PREV)    ! prepare reference values
!BOX
      ! detect restart fields using the global mean of temperature
      lpdo1 = pdotep(1,1,1,1) > 2._dp*EPSILON(1.0_dp)
      lpdo2 = pdotep(1,1,1,2) > 2._dp*EPSILON(1.0_dp)

      ! calculate tendencies
      IF (lpdo1) THEN
        ! mean is correct, also for the first call
        IF (ldiag_start) laccu_ok = .TRUE.

        ! memory is filled, the tendency can be calculated as [X(t+1) - X(t-1)]
        pdvor(:,:,:,NDVOR) = pdvor(:,:,:,NDVOR)+(svo(:,:,:)      - pdovor(:,:,:,1))
        pddiv(:,:,:,NDDIV) = pddiv(:,:,:,NDDIV)+(sd (:,:,:)      - pdodiv(:,:,:,1))
        pdtem(:,:,:,NDTEM) = pdtem(:,:,:,NDTEM)+(stp(1:nlev,:,:) - pdotep(1:nlev,:,:,1))
        pdprs(  :,:,NDPRS) = pdprs(  :,:,NDPRS)+(stp(nlevp1,:,:) - pdotep(nlevp1,:,:,1))
      END IF

      IF (lpdo2) THEN
        ! move (t) into (t-1) memory
        pdovor(:,:,:,1) = pdovor(:,:,:,2)
        pdodiv(:,:,:,1) = pdodiv(:,:,:,2)
        pdotep(:,:,:,1) = pdotep(:,:,:,2)
      END IF

      ! store (t+1) for next integration loop
      ! it is the unfiltered valu ein the spectral space
      ! in TF1 the same field is available at grid points
      pdovor(:,:,:,2) = svo(:,:,:)
      pdodiv(:,:,:,2) = sd (:,:,:)
      pdotep(:,:,:,2) = stp(:,:,:)

      ldiag_start = .FALSE.  ! switch after the first pass

    CASE default
      WRITE (diag_mess,*) 'type not implemented, ITYPE= ',itype
      CALL finish('mo_diag_tendency:DIAG_Init',diag_mess)

    END SELECT
    
  END SUBROUTINE DIAG_Init
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_SpecTrans
  ! !INTERFACE:

  SUBROUTINE DIAG_SpecTrans

    ! !DESCRIPTION: 
    ! second part of spectral transform
    !
    ! insert in SCAN1SL after CALL SI2

    ! !USES:
    USE mo_legendre,      ONLY: legmod
    USE mo_control,       ONLY: nhgl
!EOP
!BOC
!BOX
    INTEGER :: i

    north_south_loop: DO i  = 1, nhgl

      CALL DIAG_sym1(i)

      CALL legmod   (i)
      CALL DIAG_ltd

    END DO north_south_loop

  END SUBROUTINE DIAG_SpecTrans
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_fftd 
  ! !INTERFACE:
  SUBROUTINE DIAG_fftd

    ! !DESCRIPTION: 
    ! transform the diagnostic arrays into fourierspace
    !
    ! insert in SCAN1SL after CALL FFTD

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE mo_fft992,        ONLY: fft992
    USE mo_control,       ONLY: nlp2, nlon, nlev, ngl
    USE mo_spectral,      ONLY: compress, expand
!EOP
!BOC
!BOX
    INTEGER, PARAMETER :: inc = 1, isign = -1
    INTEGER            :: iamount
    REAL(dp), TARGET :: zinp(nlp2*(nlev*NO_PDIGA+NO_PDSGA)*ngl) ! work space for transformation
    REAL(dp) , POINTER :: f1d(:), f4d(:,:,:,:)

    ! transform 3-dim fields semi-implicit parts counted in grid space
    iamount = nlev*NO_PDIGS_FRAC*ngl
    f1d => zinp( (iamount*nlp2+1): )
    f4d => pdigs(:,:,1:NO_PDIGS_FRAC,:)
    CALL expand  (zinp, f4d  ,2)
    CALL expand  (f1d,  pdsgs,2)
    CALL fft992(zinp,inc, nlp2, nlon, iamount+ngl, isign)
    CALL compress(zinp, f4d  ,2)
    CALL compress(f1d,  pdsgs,2)

    ! transform other terms during output time step
    IF (l_putdata(dio_index)) THEN
      iamount = nlev*NO_PDIGA*ngl
      f1d => zinp( (iamount*nlp2+1): )
      CALL expand  (zinp, pdiga ,2)
      CALL expand  (f1d,  pdsga ,2)
      CALL fft992(zinp,inc,nlp2,nlon, iamount+NO_PDSGA*ngl,isign)
      CALL compress(zinp, pdiga ,2)
      CALL compress(f1d,  pdsga ,2)
    ENDIF

  END SUBROUTINE DIAG_fftd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_sym1
  ! !INTERFACE:

  SUBROUTINE DIAG_sym1(ihrow)

    ! !DESCRIPTION: 
    ! separate the diagnostics arrays into symetric and asymetric part
    !
    ! insert in SCAN1SL after CALL SYM1

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE mo_control,       ONLY: nmp1, nlev, ngl

    INTEGER, INTENT(in) :: ihrow   ! latitude index
!EOP
!BOC
!BOX
    INTEGER :: jl, jm, irow_n, irow_s

    !     even and odd components
    irow_n = ihrow
    irow_s = ngl + 1 - ihrow

    spec_loop1: DO jm = 1,nmp1

      level_loop1: DO jl = 1,nlev

        pdigss(jl,:,jm,:) = 0.5_dp*(pdigs(2*jm-1:2*jm,jl,:,irow_n) + pdigs(2*jm-1:2*jm,jl,:,irow_s))
        pdigsa(jl,:,jm,:) = 0.5_dp*(pdigs(2*jm-1:2*jm,jl,:,irow_n) - pdigs(2*jm-1:2*jm,jl,:,irow_s))

      ENDDO level_loop1

      pdsgss(:,jm) = 0.5_dp*(pdsgs(2*jm-1:2*jm,irow_n) + pdsgs(2*jm-1:2*jm,irow_s))
      pdsgsa(:,jm) = 0.5_dp*(pdsgs(2*jm-1:2*jm,irow_n) - pdsgs(2*jm-1:2*jm,irow_s))

    ENDDO spec_loop1

    ! transform only during output step
    IF (l_putdata(dio_index)) THEN

      spec_loop2: DO jm = 1,nmp1

        level_loop2: DO jl = 1,nlev

          pdigas(jl,:,jm,:) = 0.5_dp*(pdiga(2*jm-1:2*jm,jl,:,irow_n) + pdiga(2*jm-1:2*jm,jl,:,irow_s))
          pdigaa(jl,:,jm,:) = 0.5_dp*(pdiga(2*jm-1:2*jm,jl,:,irow_n) - pdiga(2*jm-1:2*jm,jl,:,irow_s))
          pdigbs(jl,:,jm,:) = 0.5_dp*(pdigb(2*jm-1:2*jm,jl,:,irow_n) + pdigb(2*jm-1:2*jm,jl,:,irow_s))
          pdigba(jl,:,jm,:) = 0.5_dp*(pdigb(2*jm-1:2*jm,jl,:,irow_n) - pdigb(2*jm-1:2*jm,jl,:,irow_s))

        ENDDO level_loop2

        pdsgas(:,jm,:) = 0.5_dp*(pdsga(2*jm-1:2*jm,:,irow_n) + pdsga(2*jm-1:2*jm,:,irow_s))
        pdsgaa(:,jm,:) = 0.5_dp*(pdsga(2*jm-1:2*jm,:,irow_n) - pdsga(2*jm-1:2*jm,:,irow_s))

      ENDDO spec_loop2

    ENDIF

  END SUBROUTINE DIAG_sym1
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_ltd 
  ! !INTERFACE:
  SUBROUTINE DIAG_ltd

    ! !DESCRIPTION: 
    ! perform legendre transform for diagnostic arrays
    !
    ! insert in SCAN1SL after CALL LTD

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE mo_truncation,    ONLY: nmp, nnp
    USE mo_legendre,      ONLY: pnmd, anmd, rnmd
    USE mo_control,       ONLY: nmp1
!EOP
!BOC
!BOX
    INTEGER :: jm, ims, ins, is, jn, jh, iu

    north_south_loop: DO jh = 1,2 ! 1: north, 2:south
      iu = 2-jh

      spec_loop: DO jm = 1,nmp1
        ims = nmp(jm)-iu
        ins = nnp(jm)+iu

        sym_loop: DO jn=2,ins,2
          is = ims+jn
          IF (jh == 1) THEN   !     calculations for northern hemisphere
            !
            ! semi-implicit parts of vorticity
            pdvor(:,:,is, 7) = pdvor(:,:,is, 7) + pdigss(:,:,jm,1)*pnmd(is) &
                                                - pdigsa(:,:,jm,4)*anmd(is)

            ! explicit part of divergence, temperature and pressure
            pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigss(:,:,jm,2)*rnmd(is)
            pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigss(:,:,jm,3)*pnmd(is)
            pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgss  (:,jm  )*pnmd(is)
            !
            IF (l_putdata(dio_index)) THEN
              ! dynamic and physical tendencies
              ! remark: ANMD contains the minus sign
              ! vorticity
              pdvor(:,:,is,1:5) = pdvor (:,:,is,1:5) + pdigbs(:,:,jm,6:10)*pnmd(is) &
                                                     + pdigaa(:,:,jm,1:5 )*anmd(is)
              pdvor(:,:,is,  6) = pdvor (:,:,is,  6) + pdigas(:,:,jm,21  )*pnmd(is)
              !
              ! divergence
              pddiv(:,:,is,1:5) = pddiv (:,:,is,1:5) + pdigbs(:,:,jm,1:5 )*pnmd(is) &
                                                     - pdigaa(:,:,jm,6:10)*anmd(is)
              pddiv(:,:,is,  6) = pddiv (:,:,is,  6) + pdigas(:,:,jm,22  )*pnmd(is)
              ! laplacian operation with G-term
              pddiv(:,:,is,  1) = pddiv (:,:,is,  1) + pdigas(:,:,jm,11  )*rnmd(is)
              !
              ! temperature
              pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigas(:,:,jm,12:19)*pnmd(is)
              pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigas(:,:,jm,23   )*pnmd(is)
              pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigas(:,:,jm,24:25)*pnmd(is)
              !
              ! pressure
              pdprl(:,:,is)     = pdprl(:,:,is)      + pdigas(:,:,jm,20)*pnmd(is)
              pdprs  (:,is,1:2) = pdprs  (:,is,1:2)  + pdsgas  (:,jm, :)*pnmd(is)
            ENDIF
            !
          ELSE                ! calculations for southern hemisphere
            !
            ! semi-implicit parts of vorticity
            pdvor(:,:,is, 7) = pdvor (:,:,is, 7) + pdigsa(:,:,jm,1)*pnmd(is) &
                                                 - pdigss(:,:,jm,4)*anmd(is)
            ! explicit part divergence, temperature and pressure
            pddiv(:,:,is, 7) = pddiv(:,:,is, 7) + pdigsa(:,:,jm,2)*rnmd(is)
            pdtem(:,:,is,10) = pdtem(:,:,is,10) + pdigsa(:,:,jm,3)*pnmd(is)
            pdprs  (:,is, 3) = pdprs  (:,is, 3) + pdsgsa  (:,jm  )*pnmd(is)
            !
            IF (l_putdata(dio_index)) THEN
              ! dynamic and physical tendencies
              ! vorticity
              pdvor(:,:,is,1:5) = pdvor(:,:,is,1:5) + pdigba(:,:,jm,6:10)*pnmd(is) &
                                                    + pdigas(:,:,jm,1:5 )*anmd(is)
              pdvor(:,:,is,  6) = pdvor(:,:,is,  6) + pdigaa(:,:,jm,21  )*pnmd(is)
              !
              ! divergence
              pddiv(:,:,is,1:5) = pddiv(:,:,is,1:5) + pdigba(:,:,jm,1:5 )*pnmd(is) &
                                                    - pdigas(:,:,jm,6:10)*anmd(is)
              pddiv(:,:,is,  6) = pddiv(:,:,is,  6) + pdigaa(:,:,jm,22  )*pnmd(is)
              pddiv(:,:,is,  1) = pddiv(:,:,is,  1) + pdigaa(:,:,jm,11  )*rnmd(is)
              !
              ! temperature
              pdtem(:,:,is, 1: 8) = pdtem(:,:,is, 1: 8) + pdigaa(:,:,jm,12:19)*pnmd(is)
              pdtem(:,:,is,    9) = pdtem(:,:,is,    9) + pdigaa(:,:,jm,23   )*pnmd(is)
              pdtem(:,:,is,12:13) = pdtem(:,:,is,12:13) + pdigaa(:,:,jm,24:25)*pnmd(is)
              !
              ! pressure
              pdprl(:,:,is)     = pdprl(:,:,is)     + pdigaa(:,:,jm,20)*pnmd(is)
              pdprs  (:,is,1:2) = pdprs  (:,is,1:2) + pdsgaa  (:,jm, :)*pnmd(is)
            ENDIF
            !
          ENDIF

        ENDDO sym_loop

      ENDDO spec_loop

     ENDDO north_south_loop

  END SUBROUTINE DIAG_ltd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Write
  ! !INTERFACE:

  SUBROUTINE DIAG_Write

    ! !DESCRIPTION: 
    ! correction of output buffer, reset local buffer
!EOP
!BOC
!BOX

    IF (laccu_ok) THEN
      CALL DIAG_Check
    ELSE
      pdvor(:,:,:,NDVOR) = 0.0_dp
      pddiv(:,:,:,NDDIV) = 0.0_dp
      pdtem(:,:,:,NDTEM) = 0.0_dp
      pdprs(  :,:,NDPRS) = 0.0_dp
    END IF

    ! correction of all terms due to the leap frog scheme
    pdvor(:,:,:,:) = 0.5_dp*pdvor(:,:,:,:)
    pddiv(:,:,:,:) = 0.5_dp*pddiv(:,:,:,:)
    pdtem(:,:,:,:) = 0.5_dp*pdtem(:,:,:,:)
    pdprl(:,:,:)   = 0.5_dp*pdprl(:,:,:)
    p4prs(:,:,:,:) = 0.5_dp*p4prs(:,:,:,:)

    ! reset accumulated grid point arrays
    pdiga(:,:,:,:) = 0.0_dp
    pdsga(:,:,:)   = 0.0_dp

    ! for the next postprocessing cycle the accumulation should be fine
    IF ( pdotep(1,1,1,1) > 2._dp*EPSILON(1.0_dp) ) laccu_ok = .TRUE.

  END SUBROUTINE DIAG_Write
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Check
  ! !INTERFACE:

  SUBROUTINE DIAG_Check
    ! !DESCRIPTION: 
    ! The procedure diagnoses the tendency calculation. The sum of all terms
    ! is compared with the total tendency. For the first accumulation interval
    ! the correlation can not be used. But for all following intervals all 
    ! correlations must be 1.00, except for nudging mode. In nudging mode
    ! the nudging term will not be put into the account, therefore the
    ! diagnostics are not complete.
    !
    
    ! !USES:
    USE mo_control,  ONLY: nsp, nlev
    USE mo_spectral, ONLY: corrsp
!EOP
!BOC
!BOX
    REAL(dp)          :: ccv, ccd, cct, ccp
    REAL(dp), POINTER :: p3(:,:,:), p2a(:,:), p2b(:,:)
    REAL(dp), TARGET  :: wrk(nlev,2,nsp)
    INTEGER           :: i

    wrk = 0.0_dp
    DO i=1,NDVOR-1
       wrk(:,:,:) = wrk(:,:,:) + pdvor(:,:,:,i)
    END DO
    p3 => pdvor(:,:,:,NDVOR)
    ccv = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDDIV-1
       wrk(:,:,:) = wrk(:,:,:) + pddiv(:,:,:,i)
    END DO
    p3 => pddiv(:,:,:,NDDIV)
    ccd = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDTEM-3
       wrk(:,:,:) = wrk(:,:,:) + pdtem(:,:,:,i)
    END DO
    p3 => pdtem(:,:,:,NDTEM)
    cct = corrsp (wrk,p3)

    wrk = 0.0_dp
    DO i=1,NDPRS-1
       wrk(1,:,:) = wrk(1,:,:) + pdprs(:,:,i)
    END DO
    p2a => wrk(1,:,:)
    p2b => pdprs(:,:,NDPRS)
    ccp = corrsp (p2a,p2b)

    WRITE(diag_mess,'(4(a,f8.5))') 'VOR ',ccv,' DIV ',ccd,' TEM ',cct,' PRS ',ccp
    CALL message('mo_diag_tendency:DIAG_Check',diag_mess)

  END SUBROUTINE DIAG_Check
!EOX
!EOC
END MODULE mo_diag_tendency


#else


MODULE mo_diag_tendency
!BOP
  ! !MODULE: mo_diag_tendency

  ! !DESCRIPTION:
  !\begin{verbatim}
  !-------------------------------------------------------------------
  ! overview about the separation of terms in diagnostic equations
  !
  !  LTDIAG      additional switch in &RUNCTL
  !              .true. ==> addional diagnostics
  !
  ! ***************** Interface to ECHAM *************************
  !
  ! *ECHAM*       *DIAGNOSTICS*
  !
  ! CONTROL -+
  !          +--> DIAG_Init(1)                general initialization
  !                  |
  !                  +---> DIAG_Rerun('R')    read rerun data
  !
  ! STEPON --+
  !          +--> DIAG_SumAll                 count total tendency
  !          |
  !          +--> DIAG_Write                  store data in model output
  !          |       |
  !          |       +---> DIAG_Check
  !          |       |
  !          |       +---> DiagWriteSP
  !          |
  !          +--> DIAG_Init(-1)               clear diagnostic memory
  !                  |
  !                  +---> DIAG_Rerun('W')    store rerun data
  !
  ! SCAN1SL -+
  !          +--> DIAG_fftd
  !          |
  !          SI2
  !          |
  !          +--> DIAG_SpecTrans
  !                  |
  !                  +--> DIAG_sym1
  !                  |
  !                  +--> DIAG_ltd
  !
  !---------------------------------------------------------------------
  ! count terms in additional arrays
  !   *spectral*    *Gaussian*
  !
  !                 DYN            adiabatic terms
  !                  |
  !                 TF2            time filter
  !                  |
  !                 TF1            time filter
  !                  |
  !                 GPC
  !                  +--> PHYSC    diabatic terms
  !                  |      |
  !                  |      +---> M_RADHEAT long/short wave
  !                  |
  !                  +--> SI1      semi-implicit terms
  !                  |
  !                 DIAG_fftd
  !                 |
  !               SI2                semi-implicit terms
  !               |
  !             DIAG_SpecTrans
  !                     |
  !                   DIAG_sym1
  !                     |
  !                   DIAG_ltd
  ! SCCD                           semi-implicit terms
  !   |
  ! SCCTP                          semi-implicit terms
  !   |
  ! HDIFF                          horizontal diffusion
  !   |
  ! DIAG_SumAll
  !   |
  ! DIAG_Write
  !
  !**************************************************************
  !\end{verbatim}

  ! !REVISION HISTORY: 
  ! I. Kirchner, MPI, October 1998, basic version
  ! I. Kirchner, MPI, May-2000, patch E4v3.22p1
  ! I. Kirchner, MPI, May-2002, revision E5R1.04

  ! LAST CHANGES: 
  ! parallelisation and vectorisation by Joachim Buchholz, Max Planck Institute
  !   for Chemistry, Mainz, Germany, August 2004  


  ! !USES:
  USE mo_exception,     ONLY: finish, message
  USE mo_linked_list,   ONLY: t_stream
  USE mo_kind,          ONLY: dp
  USE mo_control,       ONLY: nsp, nlev

!BOX
  IMPLICIT NONE
  
  INTRINSIC null

  PRIVATE
!EOX
  PUBLIC  :: DIAG_Init      ! allocate/deallocate memory
                            ! and initialize diagnostic arrays
  PUBLIC  :: DIAG_SpecTrans ! second part of spectral transform
  PUBLIC  :: DIAG_fftd      ! Fourier transformation of diagnostics arrays
  PUBLIC  :: DIAG_Write     ! store diagnostic fields
  ! mz_ht_20050310+
  !PRIVATE  :: DIAG_sp2gp   ! coverts spectral tendencies into gridpoint 
                            ! called just before output in stepon.f90
  ! mz_ht_20050310-

!BOX
  

  LOGICAL, SAVE, PUBLIC :: lset_fh1    = .FALSE.

  INTEGER, PARAMETER :: NO_PTFH = 3
  INTEGER, PARAMETER :: NO_PDIGA=25  ! no. of 3-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDSGA=2   ! no. of 2-dim grid point space fields
  INTEGER, PARAMETER :: NO_PDIGB=10  ! no. of 3-dim fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS=4   ! no. of 2-dim mixed grid point 
                                     ! and fourier space fields
  INTEGER, PARAMETER :: NO_PDIGS_FRAC=3   ! no. of 2-dim grid point 
                                          ! space fields fraction


  ! variables in gridpoint space
  REAL(dp), POINTER, PUBLIC, SAVE :: &   ! time integration parts
       ptfh1(:,:,:,:)=>NULL(), ptfh2(:,:)=>NULL()
  REAL(dp), POINTER, PUBLIC, SAVE :: &   ! stream elements
    pdigs(:,:,:,:)=>NULL(), pdiga(:,:,:,:)=>NULL(), pdsga (:,:,:,:)=>NULL(), &
    pdsgs(:,:,:)=>NULL()
  ! mz_pj_20050712+
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: mpdiga(:,:,:,:,:)
  ! mz_pj_20050712-

  ! allocatable arrays in Fourier space
  REAL(dp), ALLOCATABLE, PUBLIC, SAVE :: &
    pdigs_fourier(:,:,:,:), &
    pdsgs_fourier(:,:,:),   &
    pdiga_fourier(:,:,:,:), &
    pdigb_fourier(:,:,:,:), &
    pdsga_fourier(:,:,:,:)
  
  ! allocatable arrays in asymmetric/symmetric decomposition
  REAL(dp), ALLOCATABLE, SAVE :: &
    pdigaa(:,:,:,:,:), pdigas(:,:,:,:,:), &
    pdsgaa(:,:,:,:,:), pdsgas(:,:,:,:,:), &
    pdigba(:,:,:,:,:), pdigbs(:,:,:,:,:), &
    pdigsa(:,:,:,:,:), pdigss(:,:,:,:,:), &
    pdsgsa(:,:,:,:),   pdsgss(:,:,:,:)

  ! variables in spectral space
  REAL(dp), ALLOCATABLE, TARGET, PUBLIC, SAVE :: &
    pdvor(:,:,:,:), &  !   terms for vorticity equation
    pddiv(:,:,:,:), &  !   terms for divergence equation
    pdtem(:,:,:,:)     !   terms for temperature equation
  REAL(dp), POINTER, PUBLIC, SAVE :: &
    pdprs(:,:,:,:)=>NULL(), &  ! terms for surface pressure equation
    pdprl(:,:,:)  =>NULL(), &  ! surface pressure term vertical integration
    p4prs(:,:,:,:)=>NULL()
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: &   ! old values
    pdovor(:,:,:,:), pdodiv(:,:,:,:), pdotep(:,:,:,:)
  ! mz_ht_20050310+
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: diag_vor_ten_gp(:,:,:,:)
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: diag_div_ten_gp(:,:,:,:)
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: diag_tem_ten_gp(:,:,:,:)
  REAL(dp), ALLOCATABLE, TARGET, SAVE :: diag_sfp_ten_gp(:,:,:,:)

  REAL(dp), POINTER, SAVE  :: diag_pressl_ten_gp(:,:,:) => NULL()
  TYPE(t_stream), POINTER :: tdiag_gp=>NULL()

  ! lpost_sp is a switch for the additional grib output (original) of all 
  ! tendencies in spectral space which is less consistent than the online 
  ! calculated grid point output........ See Warning below !!!!

  ! both can be switched of, however the correlation ocuuring in the log file
  ! will still be calculated and reset every default output time event

  LOGICAL, PARAMETER  :: lpost_sp =.FALSE.
  LOGICAL, PARAMETER  :: lpost_gp =.TRUE.
  ! mz_ht_20050310-


!EOX
  TYPE(t_stream), POINTER :: tdiag=>NULL()
!BOX
  INTEGER, PUBLIC         :: dio_index ! event index of the TDIAG stream

  INTEGER, PARAMETER, PUBLIC :: &  ! number of terms in different equations
       NDVOR = 9, &! vorticity
       NDDIV = 9, &! divergence
       NDTEM =14, &! temperature
       NDPRS = 4   ! surface pressure

  INTEGER, PARAMETER, PUBLIC :: &  ! function separation
       IDIAG_INIT     = 1, &
       IDIAG_INI_GBUF = 2, &
       IDIAG_INI_PREV = 3, &
       IDIAG_FREE     = 4

  LOGICAL, SAVE :: &
       laccu_ok    = .FALSE., &! mean of total tendency is fine
       ldiag_start = .TRUE.

  CHARACTER(len=256) :: diag_mess = ''

!EOX
!EOP
CONTAINS
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Init
  ! !INTERFACE:

  SUBROUTINE DIAG_Init(itype)

    ! !DESCRIPTION: 
    ! initialisation of tendency diagnostics

    ! !USES:
    USE mo_memory_base, ONLY: new_stream, default_stream_setting, &
                              add_stream_element,                 &
                              SPECTRAL, GAUSSIAN, GRIB, NETCDF,   &
                              SURFACE
    USE mo_grib,        ONLY: nudging_table
    USE mo_memory_sp,   ONLY: sd, svo, stp
    USE mo_decomposition, &
      ONLY: dcl=>local_decomposition

    INTEGER, INTENT(in) :: itype       ! function separator
!EOP
!BOC
!BOX
    REAL(dp), POINTER  :: pdum(:,:,:)=>NULL(), p4(:,:,:,:)=>NULL()
    LOGICAL            :: lpdo1 = .FALSE., lpdo2 = .FALSE.
    INTEGER            :: i
    REAL(dp), POINTER  :: p3(:,:,:) => NULL()
    CHARACTER(LEN=2)   :: istr
!EOX
    INTRINSIC epsilon, maxval, min

    SELECT CASE(itype)

    CASE(IDIAG_INIT)             ! startup initialization
!BOX
      CALL message('mo_diag_tendency:DIAG_Init', &
        '  ------- start initialization -----')
!EOX
      CALL message('mo_diag_tendency', &
        'Tendency Diagnostics (E5R1.04) 22-Aug-2002 (kirchner@dkrz.de)')
!BOX
      ! define tendency diagnostic output stream

      CALL new_stream ( tdiag, 'tdiag', lpost=lpost_sp, lrerun=.TRUE., &
        filetype=GRIB)
      dio_index = tdiag%post_idx

      CALL default_stream_setting ( tdiag, &
           table=nudging_table, bits=24, units='(X)/s', contnorest=.TRUE., &
           lpost  = lpost_sp, lrerun = .TRUE., laccu  = .TRUE.,              &
           repr=SPECTRAL)
!EOX
      ! -----------------------------------------------------------
      ! allocate spectral space arrays, these fields are stored

      ! vorticity equation: 
      ALLOCATE (pdvor (dcl%nlev,2,dcl%snsp,NDVOR)); pdvor(:,:,:,:) = 0.0_dp 

      p4 => pdvor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'VOR1',  pdum,p4=p4, code=41, &
           longname = 'VorEQ horizontal advec + press.grad. + cori.term (DYN)')
      p4 => pdvor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'VOR2',  pdum,p4=p4, code=42, &
           longname = 'VorEQ vertical advection (DYN)')
      p4 => pdvor(:,:,:,3:3)
      CALL add_stream_element(tdiag,'VOR3',  pdum,p4=p4, code=43, &
           longname = 'VorEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pdvor(:,:,:,4:4)
      CALL add_stream_element(tdiag,'VOR4',  pdum,p4=p4, code=44, &
           longname = 'VorEQ gravity wave drag (GWDRAG)')
      p4 => pdvor(:,:,:,5:5)
      CALL add_stream_element(tdiag,'VOR5',  pdum,p4=p4, code=45, &
           longname = 'VorEQ moisture mass flux (CUCALL)')
      p4 => pdvor(:,:,:,6:6)
      CALL add_stream_element(tdiag,'VOR6',  pdum,p4=p4, code=46, &
           longname = 'VorEQ timefilter')
      p4 => pdvor(:,:,:,7:7)
      CALL add_stream_element(tdiag,'VOR7',  pdum,p4=p4, code=47, &
           longname = 'VorEQ semi-implicit part of time integration')
      p4 => pdvor(:,:,:,8:8)
      CALL add_stream_element(tdiag,'VOR8',  pdum,p4=p4, code=48, &
           longname = 'VorEQ horizontal diffusion')
      p4 => pdvor(:,:,:,9:9)
      CALL add_stream_element(tdiag,'VORSUM',pdum,p4=p4, code=49, &
           longname = 'VorEQ total tendency')
           
      ! divergence equation
      ALLOCATE (pddiv (dcl%nlev,2,dcl%snsp,NDDIV)); pddiv(:,:,:,:) = 0.0_dp 

      p4 => pddiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'DIV1',  pdum,p4=p4, code=61, &
       longname = 'DivEQ horizontal advec + cori. + press.grad. + G-term (DYN)')
      p4 => pddiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'DIV2',  pdum,p4=p4, code=62, &
           longname = 'DivEQ vertical advection')
      p4 => pddiv(:,:,:,3:3)
      CALL add_stream_element(tdiag,'DIV3',  pdum,p4=p4, code=63, &
           longname = 'DivEQ vertical diffusion due to impuls (VDIFF)')
      p4 => pddiv(:,:,:,4:4)
      CALL add_stream_element(tdiag,'DIV4',  pdum,p4=p4, code=64, &
           longname = 'DivEQ gravity wave drag (GWDRAG)')
      p4 => pddiv(:,:,:,5:5)
      CALL add_stream_element(tdiag,'DIV5',  pdum,p4=p4, code=65, &
           longname = 'DivEQ moisture mass flux (CUCALL)')
      p4 => pddiv(:,:,:,6:6)
      CALL add_stream_element(tdiag,'DIV6',  pdum,p4=p4, code=66, &
           longname = 'DivEQ timefilter')
      p4 => pddiv(:,:,:,7:7)
      CALL add_stream_element(tdiag,'DIV7',  pdum,p4=p4, code=67, &
           longname = 'DivEQ semi-implicit part of time integration')
      p4 => pddiv(:,:,:,8:8)
      CALL add_stream_element(tdiag,'DIV8',  pdum,p4=p4, code=68, &
           longname = 'DivEQ horizontal diffusion')
      p4 => pddiv(:,:,:,9:9)
      CALL add_stream_element(tdiag,'DIVSUM',pdum,p4=p4, code=69, &
           longname = 'DivEQ total tendency')

      ALLOCATE (pdtem (dcl%nlev,2,dcl%snsp,NDTEM))
      pdtem(:,:,:,:) = 0.0_dp ! temperature

      p4 => pdtem(:,:,:,1:1)
      CALL add_stream_element(tdiag,'TEM01', pdum,p4=p4,  code=81, &
           longname = 'TempEQ horizontal advection (DYN)')
      p4 => pdtem(:,:,:,2:2)
      CALL add_stream_element(tdiag,'TEM02', pdum,p4=p4,  code=82, &
           longname = 'TempEQ vertical advection (DYN)')
      p4 => pdtem(:,:,:,3:3)
      CALL add_stream_element(tdiag,'TEM03', pdum,p4=p4,  code=83, &
           longname = 'TempEQ energy conversion (DYN)')
      p4 => pdtem(:,:,:,4:4)
      CALL add_stream_element(tdiag,'TEM04', pdum,p4=p4,  code=84, &
           longname = 'TempEQ radiation (RADHEAT)')
      p4 => pdtem(:,:,:,5:5)
      CALL add_stream_element(tdiag,'TEM05', pdum,p4=p4,  code=85, &
           longname = 'TempEQ vertical diffusion due to turbulence (VDIFF)')
      p4 => pdtem(:,:,:,6:6)
      CALL add_stream_element(tdiag,'TEM06', pdum,p4=p4,  code=86, &
           longname = 'TempEQ gravity wave drag (GWDRAG)')
      p4 => pdtem(:,:,:,7:7)
      CALL add_stream_element(tdiag,'TEM07', pdum,p4=p4,  code=87, &
           longname = 'TempEQ convection (CUCALL)')
      p4 => pdtem(:,:,:,8:8)
      CALL add_stream_element(tdiag,'TEM08', pdum,p4=p4,  code=88, &
           longname = 'TempEQ large scale cloud processes (COND)')
      p4 => pdtem(:,:,:,9:9)
      CALL add_stream_element(tdiag,'TEM09', pdum,p4=p4,  code=89, &
           longname = 'TempEQ timefilter')
      p4 => pdtem(:,:,:,10:10)
      CALL add_stream_element(tdiag,'TEM10', pdum,p4=p4, code=90, &
           longname = 'TempEQ semi-implicit part of time integration')
      p4 => pdtem(:,:,:,11:11)
      CALL add_stream_element(tdiag,'TEM11', pdum,p4=p4, code=91, &
           longname = 'TempEQ horizontal diffusion')
      p4 => pdtem(:,:,:,12:12)
      CALL add_stream_element(tdiag,'TEM12', pdum,p4=p4, code=92, &
           longname = 'TempEQ longwave radiation')
      p4 => pdtem(:,:,:,13:13)
      CALL add_stream_element(tdiag,'TEM13', pdum,p4=p4, code=93, &
           longname = 'TempEQ shortwave radiation')
      p4 => pdtem(:,:,:,14:14)
      CALL add_stream_element(tdiag,'TEMSUM',pdum,p4=p4, code=94, &
           longname = 'TempEQ total tendency')
      
      ! divergence for each layer
      CALL add_stream_element(tdiag,'PRL',   pdprl, code=100, &
           ldims=(/dcl%nlev, 2, dcl%snsp/), & ! mz_pj_20050712
           gdims=(/nlev, 2, nsp/),          & ! mz_pj_20050712
           longname = 'PresEQ convergence in each layer')

      ! pressure equation: 
      ALLOCATE (p4prs (1,2,dcl%snsp,NDPRS)); p4prs(:,:,:,:) = 0.0_dp 
!BOX
      pdprs => p4prs(:,:,:,:)
      CALL default_stream_setting ( tdiag, &
            table=nudging_table, bits=24, &
            lpost  = lpost_sp, lrerun = .TRUE., laccu  = .TRUE., &
            repr=SPECTRAL, leveltype=SURFACE, &
            units='(X)/s', no_default=.TRUE.)
!EOX
      p4 => p4prs(1:1,:,:,1:1)
      CALL add_stream_element(tdiag,'PRS1',  pdum, code=101, klev=1, p4=p4, &
           longname = 'PresEQ vertical integrated convergence')
      p4 => p4prs(1:1,:,:,2:2)
      CALL add_stream_element(tdiag,'PRS2',  pdum, code=102, klev=1, p4=p4, &
           longname = 'PresEQ timefilter')
      p4 => p4prs(1:1,:,:,3:3)
      CALL add_stream_element(tdiag,'PRS3',  pdum, code=103, klev=1, p4=p4, &
           longname = 'PresEQ semi-implicit part of time integration')
      p4 => p4prs(1:1,:,:,4:4)
      CALL add_stream_element(tdiag,'PRSSUM',pdum, code=104, klev=1, p4=p4, &
           longname = 'PresEQ total step to step')
!BOX
      !-------------------------------------------------------------------------
      ! fields used for restart

      CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., lpost = .FALSE., laccu = .FALSE., &
            repr=SPECTRAL, &
            units='(X)', no_default=.TRUE.)

      ALLOCATE (pdovor(dcl%nlev,2,dcl%snsp,2)); pdovor(:,:,:,:) = 0.0_dp
      p4 => pdovor(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOV1',  pdum,p4=p4, &
           longname = 'PDO Vor time level 1')
      p4 => pdovor(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOV2',  pdum,p4=p4, &
           longname = 'PDO Vor time level 2')

      ALLOCATE (pdodiv(dcl%nlev,2,dcl%snsp,2)); pdodiv(:,:,:,:) = 0.0_dp
      p4 => pdodiv(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOD1',  pdum,p4=p4, &
           longname = 'PDO Div time level 1')
      p4 => pdodiv(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOD2',  pdum,p4=p4, &
           longname = 'PDO Div time level 2')

      ALLOCATE (pdotep(dcl%nlev+1,2,dcl%snsp,2)); pdotep(:,:,:,:) = 0.0_dp
      p4 => pdotep(:,:,:,1:1)
      CALL add_stream_element(tdiag,'PDOTP1', pdum,p4=p4, &
           klev=dcl%nlev+1, longname = 'PDO STP time level 1')
      p4 => pdotep(:,:,:,2:2)
      CALL add_stream_element(tdiag,'PDOTP2', pdum,p4=p4, &
           klev=dcl%nlev+1, longname = 'PDO STP time level 2')

!EOX       
      ! -----------------------------------------------------------
      ! allocate gridpoint space arrays
!BOX
       CALL default_stream_setting ( tdiag, &
            lrerun=.TRUE., repr=GAUSSIAN, &
            no_default=.TRUE.)
!EOX
      !
      !   PTFH1   memory of timefilter
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH1', ptfh1, &
           longname = 'time filter buffer 1 ', &
           ktrac = NO_PTFH)
!EOX
      !           1 ... vorticity
      !           2 ... divergence
      !           3 ... temperature
      !   PTFH2   pressure
      ! updated in TF1
!BOX
      CALL add_stream_element(tdiag,'PTFH2', ptfh2, &
           longname = 'time filter buffer 2')
!EOX
      !
      ! -----------------------------------------------------------
      !
      !     g .... accumulated in grid point space
      !     gs ... parts from gridpoint space, accumulated in spectral space
      !     s .... accumulated in spectral space
      !
      ! workspace for accumulation of terms in grid point space
      !
      ! **** PDIGA
      !
      !     Index
      !     ............. vorticity and divergenc equation
      !     VOM = Fu, VOL = Fv
      !
      !     Fu (VOM parts in GPC)
      ! DYN    1   coriolisterm and pressure tendency of Fu
      ! DYN    2   vertical advection of Fu
      ! PHYSC  3   diffusion due to impuls Pu VDIFF
      ! PHYSC  4   diffusion due to gravity wave drag Pu GWDRAG
      ! PHYSC  5   diffusion due to mass flux Pu CUCALL
      !
      !     Fv (VOL parts in GPC)
      ! DYN    6   coriolisterm and pressure tendency of Fv
      ! DYN    7   vertical advection of Fv
      ! PHYSC  8   diffusion due to impuls Pv VDIFF
      ! PHYSC  9   diffusion due to gravity wave drag Pv GWDRAG
      ! PHYSC  10  diffusion due to mass flux Pv CUCALL
      !
      ! DYN    11  potential and kinetic energy term G
      !
      !     ............. temperature tendency equation
      ! DYN    12  horizontal advection term
      ! DYN    13  vertical advection term
      ! DYN    14  energy conversion
      ! PHYSC  15  RADHEAT radiation tendency
      ! PHYSC  16  VDIFF turbulence
      ! PHYSC  17  GWDRAG gravity wave drag
      ! PHYSC  18  CUCALL convective condensation
      ! PHYSC  19  COND large+scale condensation
      !
      !     ............. pressure tendency equation
      ! DYN    20 level dependend divergence part, 
      !
      ! TF2    21  timefilter of vorticity
      ! TF2    22  timefilter of divergence
      ! TF2    23  timefilter of temperature
      !
      ! RADHEAT 24 longwave radiation
      ! RADHEAT 25 shortwave radiation
      !
!BOX
      ! mz_pj_20050712+
      !CALL add_stream_element(tdiag,'PDIGA', pdiga, &
      !     longname = 'diagnostic term DIGA', &
      !     ktrac=NO_PDIGA)
      ALLOCATE(mpdiga(dcl%nproma,dcl%nlev,NO_PDIGA,dcl%ngpblks,1))
      mpdiga(:,:,:,:,:) = 0.0_DP
      pdiga => mpdiga(:,:,:,:,1)
      ! mz_pj_20050712-

      DO i=1, NO_PDIGA
         istr = ''
         IF (i<10) THEN
            WRITE(istr, '(i1)') i
         ELSE
            WRITE(istr, '(i2)') i
         END IF
         !p4 => pdiga(:,:,i:i,:)  ! mz_pj_20050712
         p4 => mpdiga(:,:,i,:,:)  ! mz_pj_20050712
         CALL add_stream_element(tdiag,'PDIGA'//TRIM(istr), p3, &
              p4=p4, longname = 'diagnostic term DIGA')
      END DO
      ALLOCATE (pdiga_fourier(dcl%nlon+2, dcl%nflev, NO_PDIGA, dcl%nflat))
      ALLOCATE (pdigaa(dcl%nllev,2,dcl%nlm,NO_PDIGA, dcl%nlat/2))
      ALLOCATE (pdigas(dcl%nllev,2,dcl%nlm,NO_PDIGA, dcl%nlat/2))
!EOX
      !
      ! **** PDSGA
      !
      ! one level arrays of surface pressure
      !      1     last level total integral
      !      2     timefilter for pressure
      !
!BOX
      CALL add_stream_element(tdiag,'PDSGA', pdsga, &
           longname = 'diagnostic term DSGA', &
           ktrac=NO_PDSGA, klev=1)
      ALLOCATE (pdsga_fourier(dcl%nlon+2,min(1,dcl%nflev),NO_PDSGA,dcl%nflat))
      ALLOCATE (pdsgaa(min(1,dcl%nllev),2,dcl%nlm,NO_PDSGA,dcl%nlat/2))
      ALLOCATE (pdsgas(min(1,dcl%nllev),2,dcl%nlm,NO_PDSGA,dcl%nlat/2))
      pdigaa = 0._dp
      pdigas = 0._dp
!EOX
      !
      ! ***** PDIGB
      !
      !     array for d/dlambda derivatives calculated in fourierspace
      !     corresponding to pdiga(*,1...10,*)
      !
!BOX
      ALLOCATE (pdigb_fourier(dcl%nlon+2,dcl%nflev,NO_PDIGB,dcl%nflat))
      ALLOCATE (pdigba(dcl%nllev,2,dcl%nlm,NO_PDIGB,dcl%nlat/2))
      ALLOCATE (pdigbs(dcl%nllev,2,dcl%nlm,NO_PDIGB,dcl%nlat/2))
!EOX
      !
      ! ***** PDIGS
      !
      !     local memory buffer for accumulation of semi-implicit parts
      !     fraction solved in grid point space
      !
      !    1 vorticity and implicit and explicit part (L)
      !    2 divergence
      !    3 temperatur and pressure
      !    4 implicit part of vorticity (M) used in Fourierspace
!BOX
      ALLOCATE (pdigs(dcl%nproma,dcl%nlev,NO_PDIGS_FRAC,dcl%ngpblks))
      ALLOCATE (pdigs_fourier(dcl%nlon+2,dcl%nflev,NO_PDIGS,dcl%nflat))
      ALLOCATE (pdigsa(dcl%nllev,2,dcl%nlm,NO_PDIGS,dcl%nlat/2))
      ALLOCATE (pdigss(dcl%nllev,2,dcl%nlm,NO_PDIGS,dcl%nlat/2))
!EOX
      !
      ! ****** PDSGS
      !
      !     semi-implicit part of log surface pressure
!BOX
      CALL add_stream_element(tdiag, 'PDSGS', pdsgs, &
           longname = 'diagnostic term DSGS', &
           klev=1, lpost=.false.)
      ALLOCATE (pdsgs(dcl%nproma,1,dcl%ngpblks))
      ALLOCATE (pdsgs_fourier(dcl%nlon+2,min(1,dcl%nflev),dcl%nflat))
      ALLOCATE (pdsgsa(min(1,dcl%nllev),2,dcl%nlm,dcl%nlat/2))
      ALLOCATE (pdsgss(min(1,dcl%nllev),2,dcl%nlm,dcl%nlat/2))
      
      ! mz_ht_20050310+
      ! setting up of new stream for diagnostic output of diag_tendency 
      ! in gridpoint space
      
!=============================================================================

! *********************** WARNING ********************************************

      ! The output in gridpoint stream does not give identical results compared
      ! to offline transformed spectral fields (afterburner)
      ! because transformation algorithms of ECHAM5 and AFTERBURN differ.
      ! Check difference between sum and respective contributions
      ! of vorticity, divergence and temperature.
      ! These differences are smaller in this gridpoint stream than in the
      ! offline transformed fields.
      ! in BRIEF:   USE THIS GRIDPOINT STREAM FOR DIAGNOSTICS !!!!!!!!!!!!!

!=============================================================================

      CALL message('mo_diag_tendency',' Building GP-Stream')

      CALL NEW_STREAM(tdiag_gp, 'tdiag_gp', lpost=lpost_gp, lrerun=.FALSE., &
                                            filetype=NETCDF)

      CALL default_stream_setting ( tdiag_gp, contnorest=.TRUE.,        &
                                    lpost  = lpost_gp, lrerun=.FALSE.,     &
                                    laccu  = .TRUE., repr = GAUSSIAN)

      ! the laccu switch is in this case only used to trigger division 
      ! by time interval before output 
      ! Since the gridpoint field is only a transformed instantaneous
      ! "copy" of the accumulated spectral fields, explicit accumulation
      ! is not required.

      ! vorticity      
      ALLOCATE(diag_vor_ten_gp(dcl%nproma,dcl%nlev,dcl%ngpblks, NDVOR))
      diag_vor_ten_gp(:,:,:,:) = 0._dp
      p4 => diag_vor_ten_gp(:,:,:,1:1)
      CALL add_stream_element(tdiag_gp, 'VOR1', pdum, p4=p4,                 &
        longname = 'VorEQ horizontal advec + press.grad. + cori.term (DYN)', &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,2:2)
      CALL add_stream_element(tdiag_gp, 'VOR2', pdum, p4=p4,                 &
        longname = 'VorEQ vertical advection (DYN)',                         &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,3:3)
      CALL add_stream_element(tdiag_gp, 'VOR3', pdum, p4=p4,                 &
        longname = 'VorEQ vertical diffusion due to impuls (VDIFF)',         &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,4:4)
      CALL add_stream_element(tdiag_gp, 'VOR4', pdum ,p4=p4,                 &
        longname = 'VorEQ gravity wave drag (GWDRAG)',                       &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,5:5)
      CALL add_stream_element(tdiag_gp, 'VOR5', pdum, p4=p4,                 &
        longname = 'VorEQ moisture mass flux (CUCALL)',                      &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,6:6)
      CALL add_stream_element(tdiag_gp, 'VOR6', pdum, p4=p4,                 &
        longname = 'VorEQ timefilter',                                       &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,7:7)
      CALL add_stream_element(tdiag_gp, 'VOR7', pdum, p4=p4,                 &
        longname = 'VorEQ semi-implicit part of time integration',           &
        units = '[VOR]/s')
      p4 => diag_vor_ten_gp(:,:,:,8:8)
      CALL add_stream_element(tdiag_gp, 'VOR8', pdum, p4=p4,                 &
        longname = 'VorEQ horizontal diffusion',                             &
        units = '[VOR]/s') 
      p4 => diag_vor_ten_gp(:,:,:,9:9)
      CALL add_stream_element(tdiag_gp,'VORSUM', pdum, p4=p4,                &
        longname = 'VorEQ total tendency',                                   &
        units = '[VOR]/s') 

      !      divergence
      ALLOCATE(diag_div_ten_gp(dcl%nproma,dcl%nlev,dcl%ngpblks, NDDIV))
      diag_div_ten_gp(:,:,:,:) = 0._dp
      p4 => diag_div_ten_gp(:,:,:,1:1)
      CALL add_stream_element(tdiag_gp, 'DIV1', pdum, p4=p4,                 &
        longname = 'DIVEQ horizontal advec + press.grad. + cori.term (DYN)', &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,2:2)
      CALL add_stream_element(tdiag_gp, 'DIV2', pdum, p4=p4,                 &
        longname = 'DIVEQ vertical advection (DYN)',                         &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,3:3)
      CALL add_stream_element(tdiag_gp, 'DIV3', pdum, p4=p4,                 &
        longname = 'DIVEQ vertical diffusion due to impuls (VDIFF)',         &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,4:4)
      CALL add_stream_element(tdiag_gp, 'DIV4', pdum ,p4=p4,                 &
        longname = 'DIVEQ gravity wave drag (GWDRAG)',                       &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,5:5)
      CALL add_stream_element(tdiag_gp, 'DIV5', pdum, p4=p4,                 &
        longname = 'DIVEQ moisture mass flux (CUCALL)',                      &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,6:6)
      CALL add_stream_element(tdiag_gp, 'DIV6', pdum, p4=p4,                 &
        longname = 'DIVEQ timefilter',                                       &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,7:7)
      CALL add_stream_element(tdiag_gp, 'DIV7', pdum, p4=p4,                 &
        longname = 'DIVEQ semi-implicit part of time integration',           &
        units = '[DIV]/s')
      p4 => diag_div_ten_gp(:,:,:,8:8)
      CALL add_stream_element(tdiag_gp, 'DIV8', pdum, p4=p4,                 &
        longname = 'DIVEQ horizontal diffusion',                             &
        units = '[DIV]/s') 
      p4 => diag_div_ten_gp(:,:,:,9:9)
      CALL add_stream_element(tdiag_gp,'DIVSUM', pdum, p4=p4,                &
        longname = 'DIVEQ total tendency',                                   &
        units = '[DIV]/s') 

      !     temperature
      ALLOCATE (diag_tem_ten_gp (dcl%nproma,dcl%nlev,dcl%ngpblks, NDTEM))
      diag_tem_ten_gp(:,:,:,:) = 0.0_dp ! temperature

      p4 => diag_tem_ten_gp(:,:,:,1:1)
      CALL add_stream_element(tdiag_gp,'TEM01', pdum, p4=p4,                 &
           longname = 'TempEQ horizontal advection (DYN)',                   &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,2:2)
      CALL add_stream_element(tdiag_gp,'TEM02', pdum, p4=p4,                 &
           longname = 'TempEQ vertical advection (DYN)',                     &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,3:3)
      CALL add_stream_element(tdiag_gp,'TEM03', pdum, p4=p4,                 &
           longname = 'TempEQ energy conversion (DYN)',                      &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,4:4)
      CALL add_stream_element(tdiag_gp,'TEM04', pdum, p4=p4,                 &
           longname = 'TempEQ radiation (RADHEAT)',                          &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,5:5)
      CALL add_stream_element(tdiag_gp,'TEM05', pdum, p4=p4,                 &
           longname = 'TempEQ vertical diffusion due to turbulence (VDIFF)', &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,6:6)
      CALL add_stream_element(tdiag_gp,'TEM06', pdum, p4=p4,                 &
           longname = 'TempEQ gravity wave drag (GWDRAG)',                   &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,7:7)
      CALL add_stream_element(tdiag_gp,'TEM07', pdum, p4=p4,                 &
           longname = 'TempEQ convection (CUCALL)',                          &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,8:8)
      CALL add_stream_element(tdiag_gp,'TEM08', pdum, p4=p4,                 &
           longname = 'TempEQ large scale cloud processes (COND)',           &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,9:9)
      CALL add_stream_element(tdiag_gp,'TEM09', pdum, p4=p4,                 &
           longname = 'TempEQ timefilter',                                   &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,10:10)
      CALL add_stream_element(tdiag_gp,'TEM10', pdum, p4=p4,                 &
           longname = 'TempEQ semi-implicit part of time integration',       &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,11:11)
      CALL add_stream_element(tdiag_gp,'TEM11', pdum, p4=p4,                 &
           longname = 'TempEQ horizontal diffusion',                         &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,12:12)
      CALL add_stream_element(tdiag_gp,'TEM12', pdum, p4=p4,                 &
           longname = 'TempEQ longwave radiation',                           &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,13:13)
      CALL add_stream_element(tdiag_gp,'TEM13', pdum, p4=p4,                 &
           longname = 'TempEQ shortwave radiation',                          &
           units = 'K/s')
      p4 => diag_tem_ten_gp(:,:,:,14:14)
      CALL add_stream_element(tdiag_gp,'TEMSUM', pdum, p4=p4,                &
           longname = 'TempEQ total tendency',                               &
           units = 'K/s')

      ! divergence for each layer
      CALL add_stream_element(tdiag_gp,'PRL',   diag_pressl_ten_gp,          &
           longname = 'PresEQ convergence in each layer', units= 'Pa/s')

      ! pressure equation: 
      ALLOCATE (diag_sfp_ten_gp (dcl%nproma,1,dcl%ngpblks, NDPRS))
      diag_sfp_ten_gp(:,:,:,:) = 0.0_dp 

      p4 => diag_sfp_ten_gp(:,:,:,1:1)
      CALL add_stream_element(tdiag_gp,'PRS1', pdum, klev=1, p4=p4,     &
        longname = 'PresEQ vertical integrated convergence',            &
        units = 'Pa/s')
      p4 => diag_sfp_ten_gp(:,:,:,2:2)
      CALL add_stream_element(tdiag_gp,'PRS2', pdum, klev=1, p4=p4,     &
           longname = 'PresEQ timefilter',                              &
           units = 'Pa/s')
      p4 =>  diag_sfp_ten_gp(:,:,:,3:3)
      CALL add_stream_element(tdiag_gp,'PRS3', pdum, klev=1, p4=p4,     &
           longname = 'PresEQ semi-implicit part of time integration',  &
           units = 'Pa/s')
      p4 =>  diag_sfp_ten_gp(:,:,:,4:4)
      CALL add_stream_element(tdiag_gp,'PRSSUM', pdum, klev=1, p4=p4,   &
           longname = 'PresEQ total step to step',                      &
           units = 'Pa/s')
      CALL message('mo_diag_tendency',' Finished GP-Stream')
    ! mz_ht_20050310-

      CALL message('mo_diag_tendency:DIAG_Init', &
        '  ------- end of initialization -----')
!EOX
    CASE (IDIAG_FREE)        ! get memory free
!BOX
      DEALLOCATE (pdvor, pddiv, pdtem, p4prs)
      DEALLOCATE (pdovor, pdodiv, pdotep)
      DEALLOCATE        (pdiga_fourier, pdigaa, pdigas)
      DEALLOCATE        (pdsga_fourier, pdsgaa, pdsgas)
      DEALLOCATE        (pdigb_fourier, pdigba, pdigbs)
      DEALLOCATE (pdigs, pdigs_fourier, pdigsa, pdigss)
      DEALLOCATE (pdsgs, pdsgs_fourier, pdsgsa, pdsgss)
      ! mz_ht_20050310+
      DEALLOCATE (diag_vor_ten_gp, diag_div_ten_gp, &
           diag_tem_ten_gp, diag_sfp_ten_gp)
      ! mz_ht_20050310-
      ! mz_pj_20050712+
      DEALLOCATE(mpdiga)
      ! mz_pj_20050712-

!EOX
    CASE (IDIAG_INI_GBUF)    ! reset local accumulation buffer
!BOX
      pdigb_fourier(:,:,:,:) = 0.0_dp
      pdigs(:,:,:,:)         = 0.0_dp
      pdigs_fourier(:,:,:,:) = 0.0_dp
      pdsgs(:,:,:)           = 0.0_dp
      pdsgs_fourier(:,:,:)   = 0.0_dp
!EOX
    CASE (IDIAG_INI_PREV)    ! prepare reference values
!BOX
      ! detect restart fields using the global mean of temperature
      lpdo1 = maxval(ABS(pdotep(:,:,:,1))) > 2.0_dp*EPSILON(1.0_dp)
      lpdo2 = maxval(ABS(pdotep(:,:,:,2))) > 2.0_dp*EPSILON(1.0_dp)

      ! calculate tendencies
      IF (lpdo1) THEN
        ! mean is correct, also for the first call
        IF (ldiag_start) laccu_ok = .TRUE.

        ! memory is filled, the tendency can be calculated as [X(t+1) - X(t-1)]
        pdvor(:,:,:,NDVOR) = pdvor(:,:,:,NDVOR)+(svo(:,:,:) &
                            -pdovor(:,:,:,1))
        pddiv(:,:,:,NDDIV) = pddiv(:,:,:,NDDIV)+(sd (:,:,:) &
                            -pdodiv(:,:,:,1))
        pdtem(:,:,:,NDTEM) = pdtem(:,:,:,NDTEM)+(stp(1:dcl%nlev,:,:) &
                            -pdotep(1:dcl%nlev,:,:,1))
        pdprs(1,:,:,NDPRS) = pdprs(1,:,:,NDPRS)+(stp(dcl%nlev+1,:,:) &
                            -pdotep(dcl%nlev+1,:,:,1))
      END IF

      IF (lpdo2) THEN
        ! move (t) into (t-1) memory
        pdovor(:,:,:,1) = pdovor(:,:,:,2)
        pdodiv(:,:,:,1) = pdodiv(:,:,:,2)
        pdotep(:,:,:,1) = pdotep(:,:,:,2)
      END IF

      ! store (t+1) for next integration loop
      ! it is the unfiltered valu ein the spectral space
      ! in TF1 the same field is available at grid points
      pdovor(:,:,:,2) = svo(:,:,:)
      pdodiv(:,:,:,2) = sd (:,:,:)
      pdotep(:,:,:,2) = stp(:,:,:)

      ldiag_start = .FALSE.  ! switch after the first pass

    CASE default
      WRITE (diag_mess,*) 'type not implemented, ITYPE= ',itype
      CALL finish('mo_diag_tendency:DIAG_Init',diag_mess)

    END SELECT
    
  END SUBROUTINE DIAG_Init
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_SpecTrans
  ! !INTERFACE:

  SUBROUTINE DIAG_SpecTrans

    ! !DESCRIPTION: 
    ! second part of spectral transform
    !
    ! insert in SCAN1SL after CALL SI2

!EOP
!BOC
!BOX

      CALL DIAG_sym1
      CALL DIAG_ltd

  END SUBROUTINE DIAG_SpecTrans
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_fftd 
  ! !INTERFACE:
  SUBROUTINE DIAG_fftd

    ! !DESCRIPTION: 
    ! transform the diagnostic arrays into fourierspace
    !
    ! Step 1: 
    !   transpose pdigs, pdsgs, pdiga, pdsga from gridpoint space into a
    !   gridpoint/Fourier space intermediate decomposition with the subroutine
    !   trp_gp_fs_3d from messy_main_transform_bi.f90
    ! Step 2: 
    !   transform gridpoint/Fourier space intermediate decomposition variables
    !   pdigs_gf, pdsgs_gf, pdiga_gf, pdsga_gf into Fourier space
    ! Step 3: 
    !   copy Fourier space values from pdigs_gf, pdsgs_gf, pdiga_gf, pdsga_gf
    !   into Fourier space variables pdigs_fourier, pdsgs_fourier, 
    !   pdiga_fourier, pdsga_fourier
    !
    ! insert in SCAN1SL after CALL FFTD

    ! !USES:
    USE messy_main_transform_bi, &
      ONLY: trp_gp_fs_3d, trf_fftd
    USE mo_decomposition, &
      ONLY: dcl=>local_decomposition
    USE mo_time_control,  ONLY: l_putdata
    
!EOP
!BOC
!BOX

    INTRINSIC min

    REAL(dp) :: pdigs_gf(dcl%nlon+2,dcl%nflev,dcl%nflat), &
                pdsgs_gf(dcl%nlon+2,min(1,dcl%nflev),dcl%nflat),           &
                pdiga_gf(dcl%nlon+2,dcl%nflev,dcl%nflat), &
                pdsga_gf(dcl%nlon+2,min(1,dcl%nflev),dcl%nflat)
    INTEGER :: count

    DO count = 1, NO_PDIGS_FRAC
      pdigs_gf(:,:,:) = 0.0_dp
      CALL trp_gp_fs_3d(1,pdigs(:,:,count,:),pdigs_gf(:,:,:))
      CALL trf_fftd(pdigs_gf(:,:,:))
      pdigs_fourier(:,:,count,:) = pdigs_gf(:,:,:)
    END DO

    pdsgs_gf(:,:,:) = 0.0_dp
    CALL trp_gp_fs_3d(1,pdsgs(:,:,:),pdsgs_gf(:,:,:))
    CALL trf_fftd(pdsgs_gf(:,:,:))
    pdsgs_fourier(:,:,:) = pdsgs_gf(:,:,:)

    ! transform other terms during output time step
    IF (l_putdata(dio_index)) THEN
      DO count = 1, NO_PDIGA
        pdiga_gf(:,:,:) = 0.0_dp
        CALL trp_gp_fs_3d(1,pdiga(:,:,count,:),pdiga_gf(:,:,:))
        CALL trf_fftd(pdiga_gf(:,:,:))
        pdiga_fourier(:,:,count,:) = pdiga_gf(:,:,:)
      END DO

      DO count = 1, NO_PDSGA
        pdsga_gf(:,:,:) = 0.0_dp
        CALL trp_gp_fs_3d(1,pdsga(:,:,count,:),pdsga_gf(:,:,:))
        CALL trf_fftd(pdsga_gf(:,:,:))
        pdsga_fourier(:,:,count,:) = pdsga_gf(:,:,:)
      END DO
    ENDIF

  END SUBROUTINE DIAG_fftd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_sym1
  ! !INTERFACE:

  SUBROUTINE DIAG_sym1

    ! !DESCRIPTION: 
    ! separate the diagnostics arrays into symetric and asymetric part
    !
    ! Step 1: 
    !   transpose pdigs_fourier, pdsgs_fourier, pdiga_fourier, pdigb_fourier,
    !   pdsga_fourier from Fourier space into a Fourier/asymmetric-symmetric
    !   intermediate decomposition with the subroutine trp_fs_ls_3d from
    !   messy_main_transform_bi.f90
    ! Step 2: 
    !   transpose Fourier/asymmetric-symmetric intermediate decomposition
    !   variables pdigs_fas, pdsgs_fas, pdiga_fas, pdigb_fas, pdsga_fas into 
    !   asymmetric/symmetric decomposition
    !
    ! insert in SCAN1SL after CALL SYM1

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE messy_main_transform_bi, &
      ONLY: trp_fs_ls_3d, trp_fs_as
    USE mo_decomposition, &
      ONLY: dcl=>local_decomposition

!EOP
!BOC
!BOX
      
    INTRINSIC min

    REAL(dp) :: pdigs_fas(dcl%nlm*2,dcl%nflev,dcl%nlat), &
                pdsgs_fas(dcl%nlm*2,min(1,dcl%nflev),dcl%nlat),           &
                pdiga_fas(dcl%nlm*2,dcl%nflev,dcl%nlat), &
                pdigb_fas(dcl%nlm*2,dcl%nflev,dcl%nlat), &
                pdsga_fas(dcl%nlm*2,min(1,dcl%nflev),dcl%nlat)    
    INTEGER :: count

    DO count=1,NO_PDIGS
      pdigs_fas(:,:,:) = 0.0_dp
      CALL trp_fs_ls_3d(1,pdigs_fourier(:,:,count,:),pdigs_fas(:,:,:))
      CALL trp_fs_as(1,pdigs_fas(:,:,:), &
                     pdigsa(:,:,:,count,:),pdigss(:,:,:,count,:))
    END DO

    pdsgs_fas(:,:,:) = 0.0_dp
    CALL trp_fs_ls_3d(1,pdsgs_fourier(:,:,:),pdsgs_fas(:,:,:))
    CALL trp_fs_as(1,pdsgs_fas(:,:,:), &
                   pdsgsa(:,:,:,:),pdsgss(:,:,:,:))

    ! transform only during output step
    IF (l_putdata(dio_index)) THEN

      DO count=1,NO_PDIGA
        pdiga_fas(:,:,:) = 0.0_dp
        CALL trp_fs_ls_3d(1,pdiga_fourier(:,:,count,:),pdiga_fas(:,:,:))
        CALL trp_fs_as(1,pdiga_fas(:,:,:), &
                       pdigaa(:,:,:,count,:),pdigas(:,:,:,count,:))
      END DO
      
      DO count=1,NO_PDIGB
        pdigb_fas(:,:,:) = 0.0_dp
        CALL trp_fs_ls_3d(1,pdigb_fourier(:,:,count,:),pdigb_fas(:,:,:))
        CALL trp_fs_as(1,pdigb_fas(:,:,:), &
                       pdigba(:,:,:,count,:),pdigbs(:,:,:,count,:))
      END DO
      
      DO count=1,NO_PDSGA
        pdsga_fas(:,:,:) = 0.0_dp
        CALL trp_fs_ls_3d(1,pdsga_fourier(:,:,count,:),pdsga_fas(:,:,:))
        CALL trp_fs_as(1,pdsga_fas(:,:,:), &
                       pdsgaa(:,:,:,count,:),pdsgas(:,:,:,count,:))
      END DO

    ENDIF

  END SUBROUTINE DIAG_sym1
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_ltd 
  ! !INTERFACE:
  SUBROUTINE DIAG_ltd

    ! !DESCRIPTION: 
    ! perform legendre transform for diagnostic arrays
    !
    ! Step 1: 
    !   transform the asymmetric/symmetric variable couples pdigsa/pdigss, 
    !   pdsgsa/pdsgss, pdigaa/pdigas, pdigba/pdigbs, pdsga/pdsgs into Legendre
    !   space with the subroutine trf_ltd from messy_main_transform_bi.f90
    ! Step 2: 
    !   transpose the Legendre space variables pdigs_legendre, pdsgs_legendre, 
    !   pdiga_legendre, pdigb_legendre, pdsga_legendre into spectral space
    !   with the subroutine trp_ls_sp_3d from messy_main_transform_bi.f90
    ! Step 3: 
    !   add values from the spectral space variables pdigs_spectral, 
    !   pdsgs_spectral, pdiga_spectral, pdigb_spectral, pdsga_spectral to 
    !   the spectral space variables pdvor, pddiv, pdtem, pdprl, pdprs
    !
    ! insert in SCAN1SL after CALL LTD

    ! !USES:
    USE mo_time_control,  ONLY: l_putdata
    USE messy_main_transform_bi, &
                          ONLY: trf_ltd, trp_ls_sp_3d, trf_ltd_a, trf_ltd_r
    USE mo_decomposition, &
      ONLY: dcl=>local_decomposition

!EOP
!BOC
!BOX

    INTRINSIC min

    REAL(dp) :: pdigs_legendre(dcl%nllev,2,dcl%lnsp), &
                pdsgs_legendre(min(1,dcl%nllev),2,dcl%lnsp),         &
                pdiga_legendre(dcl%nllev,2,dcl%lnsp), &
                pdigb_legendre(dcl%nllev,2,dcl%lnsp), &
                pdsga_legendre(min(1,dcl%nllev),2,dcl%lnsp)
    REAL(dp) :: pdigs_spectral(dcl%nlev,2,dcl%snsp,NO_PDIGS), &
                pdsgs_spectral(1,2,dcl%snsp),        &
                pdiga_spectral(dcl%nlev,2,dcl%snsp,NO_PDIGA), &
                pdigb_spectral(dcl%nlev,2,dcl%snsp,NO_PDIGB), &
                pdsga_spectral(1,2,dcl%snsp,NO_PDSGA)
    INTEGER :: count

    pdigs_spectral(:,:,:,:) = 0.0_dp
    CALL trf_ltd(pdigs_legendre(:,:,:), &
      pdigsa(:,:,:,1,:),pdigss(:,:,:,1,:))
    CALL trp_ls_sp_3d(1,pdigs_legendre(:,:,:), &
      pdigs_spectral(:,:,:,1))
    CALL trf_ltd_r(pdigs_legendre(:,:,:), &
      pdigsa(:,:,:,2,:),pdigss(:,:,:,2,:))
    CALL trp_ls_sp_3d(1,pdigs_legendre(:,:,:), &
      pdigs_spectral(:,:,:,2))
    CALL trf_ltd(pdigs_legendre(:,:,:), &
      pdigsa(:,:,:,3,:),pdigss(:,:,:,3,:))
    CALL trp_ls_sp_3d(1,pdigs_legendre(:,:,:), &
      pdigs_spectral(:,:,:,3))
    CALL trf_ltd_a(pdigs_legendre(:,:,:), &
      pdigss(:,:,:,4,:),pdigsa(:,:,:,4,:))        ! CAUTION !
    CALL trp_ls_sp_3d(1,pdigs_legendre(:,:,:), &
      pdigs_spectral(:,:,:,4))

    pdsgs_spectral(:,:,:) = 0.0_dp
    CALL trf_ltd(pdsgs_legendre(:,:,:), &
                 pdsgsa(:,:,:,:),pdsgss(:,:,:,:))
    CALL trp_ls_sp_3d(1,pdsgs_legendre(:,:,:), &
                      pdsgs_spectral(:,:,:))

    ! transform only during output step
    IF (l_putdata(dio_index)) THEN
      pdiga_spectral(:,:,:,:) = 0.0_dp
      DO count=1,10
        CALL trf_ltd_a(pdiga_legendre(:,:,:), &
                     pdigas(:,:,:,count,:),pdigaa(:,:,:,count,:)) ! CAUTION !
        CALL trp_ls_sp_3d(1,pdiga_legendre(:,:,:), &
                          pdiga_spectral(:,:,:,count))
      END DO
        CALL trf_ltd_r(pdiga_legendre(:,:,:), &
                     pdigaa(:,:,:,11,:),pdigas(:,:,:,11,:))
        CALL trp_ls_sp_3d(1,pdiga_legendre(:,:,:), &
                          pdiga_spectral(:,:,:,11))
      DO count=12,NO_PDIGA
        CALL trf_ltd(pdiga_legendre(:,:,:), &
                     pdigaa(:,:,:,count,:),pdigas(:,:,:,count,:))
        CALL trp_ls_sp_3d(1,pdiga_legendre(:,:,:), &
                          pdiga_spectral(:,:,:,count))
      END DO

      pdigb_spectral(:,:,:,:) = 0.0_dp
      DO count=1,NO_PDIGB
        CALL trf_ltd(pdigb_legendre(:,:,:), &
                     pdigba(:,:,:,count,:),pdigbs(:,:,:,count,:))
        CALL trp_ls_sp_3d(1,pdigb_legendre(:,:,:), &
                          pdigb_spectral(:,:,:,count))
      END DO

      pdsga_spectral(:,:,:,:) = 0.0_dp
      DO count=1,NO_PDSGA
        CALL trf_ltd(pdsga_legendre(:,:,:), &
                     pdsgaa(:,:,:,count,:),pdsgas(:,:,:,count,:))
        CALL trp_ls_sp_3d(1,pdsga_legendre(:,:,:), &
                          pdsga_spectral(:,:,:,count))
      END DO
    END IF

    pdvor(:,:,:,7) = pdvor(:,:,:,7) + pdigs_spectral(:,:,:,1) &
                                    - pdigs_spectral(:,:,:,4)
    pddiv(:,:,:,7) = pddiv(:,:,:,7) + pdigs_spectral(:,:,:,2)
    pdtem(:,:,:,10)= pdtem(:,:,:,10)+ pdigs_spectral(:,:,:,3)
    pdprs(:,:,:,3) = pdprs(:,:,:,3) + pdsgs_spectral(:,:,:)

    IF (l_putdata(dio_index)) THEN

      pdvor(:,:,:,1:5)   =  pdvor(:,:,:,1:5)           &
                           +pdigb_spectral(:,:,:,6:10) &
                           +pdiga_spectral(:,:,:,1:5)
      pdvor(:,:,:,6)     =  pdvor(:,:,:,6)             &
                           +pdiga_spectral(:,:,:,21)
      pddiv(:,:,:,1:5)   =  pddiv(:,:,:,1:5)           &
                           +pdigb_spectral(:,:,:,1:5)  &
                           -pdiga_spectral(:,:,:,6:10)
      pddiv(:,:,:,6)     =  pddiv(:,:,:,6)             &
                           +pdiga_spectral(:,:,:,22) 
      pddiv(:,:,:,1)     =  pddiv(:,:,:,1)             &
                           +pdiga_spectral(:,:,:,11)
      pdtem(:,:,:,1:8)   =  pdtem(:,:,:,1:8)           &
                           +pdiga_spectral(:,:,:,12:19)
      pdtem(:,:,:,9)     =  pdtem(:,:,:,9)             &
                           +pdiga_spectral(:,:,:,23)
      pdtem(:,:,:,12:13) =  pdtem(:,:,:,12:13)         &
                           +pdiga_spectral(:,:,:,24:25)
      pdprl(:,:,:)       =  pdprl(:,:,:)               &
                           +pdiga_spectral(:,:,:,20)
      pdprs(:,:,:,1:2)   =  pdprs(:,:,:,1:2)           &
                           +pdsga_spectral(:,:,:,:)

    END IF

  END SUBROUTINE DIAG_ltd
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Write
  ! !INTERFACE:

  SUBROUTINE DIAG_Write

    ! !DESCRIPTION: 
    ! correction of output buffer, reset local buffer
    
    INTRINSIC epsilon, maxval
        
!EOP
!BOC
!BOX

    IF (laccu_ok) THEN
      CALL DIAG_Check
    ELSE
      pdvor(:,:,:,NDVOR) = 0.0_dp
      pddiv(:,:,:,NDDIV) = 0.0_dp
      pdtem(:,:,:,NDTEM) = 0.0_dp
      pdprs(:,:,:,NDPRS) = 0.0_dp
    END IF

    ! correction of all terms due to the leap frog scheme
    pdvor(:,:,:,:) = 0.5_dp*pdvor(:,:,:,:)
    pddiv(:,:,:,:) = 0.5_dp*pddiv(:,:,:,:)
    pdtem(:,:,:,:) = 0.5_dp*pdtem(:,:,:,:)
    pdprl(:,:,:)   = 0.5_dp*pdprl(:,:,:)
    p4prs(:,:,:,:) = 0.5_dp*p4prs(:,:,:,:)

    ! reset accumulated grid point arrays
    pdiga(:,:,:,:) = 0.0_dp
    pdsga(:,:,:,:) = 0.0_dp

    ! for the next postprocessing cycle the accumulation should be fine
    IF ( maxval(ABS(pdotep(:,:,:,1))) > 2.0_dp*EPSILON(1.0_dp) ) &
         laccu_ok = .TRUE.

    ! mz_ht_20050310+
    CALL DIAG_sp2gp

    if (.not. lpost_sp) then
      pdvor(:,:,:,:) = 0._dp
      pddiv(:,:,:,:) = 0._dp
      pdtem(:,:,:,:) = 0._dp
      pdprl(:,:,:)   = 0._dp
      p4prs(:,:,:,:) = 0._dp
    endif
    ! mz_ht_20050310-

  END SUBROUTINE DIAG_Write
!EOX
!EOC
!======================================================================
!BOP
  ! !IROUTINE:  DIAG_Check
  ! !INTERFACE:

  SUBROUTINE DIAG_Check
    ! !DESCRIPTION: 
    ! The procedure diagnoses the tendency calculation. The sum of all terms
    ! is compared with the total tendency. For the first accumulation interval
    ! the correlation can not be used. But for all following intervals all 
    ! correlations must be 1.00, except for nudging mode. In nudging mode
    ! the nudging term will not be put into the account, therefore the
    ! diagnostics are not complete.
    !
    
    ! !USES:
    USE mo_spectral,      ONLY: corrsp
    USE mo_decomposition, ONLY: dcl=>local_decomposition &
                              , dcg=>global_decomposition
    USE mo_mpi,           ONLY: p_parallel_io
    USE mo_control,       ONLY: nlev, nsp
    USE mo_transpose,     ONLY: gather_sp
!EOP
!BOC
!BOX
    REAL(dp)          :: ccv, ccd, cct, ccp
    REAL(dp), POINTER :: p3(:,:,:)=>NULL(), p2a(:,:)=>NULL(), p2b(:,:)=>NULL()
    REAL(dp), TARGET  :: wrk(dcl%nlev,2,dcl%snsp)
    INTEGER           :: i

      ! mz_ht_20050419+
    REAL(dp), POINTER :: g_wrk(:,:,:)
    REAL(dp), POINTER :: g_wrk_2(:,:,:)

!    REAL(dp), POINTER :: wrk(:,:,:)
!    ALLOCATE  (wrk(nlev,2,nsp))

!   correlation for global field for vorticity

    if (p_parallel_io) then
      ALLOCATE(g_wrk(nlev,2,nsp))
      ALLOCATE(g_wrk_2(nlev,2,nsp))
    else
      NULLIFY(g_wrk)
      NULLIFY(g_wrk_2)
    endif

    wrk = 0.0_dp
    DO i=1,NDVOR-1
       wrk(:,:,:) = wrk(:,:,:) + pdvor(:,:,:,i)
    END DO

    CALL gather_sp(g_wrk, wrk, dcg)
    CALL gather_sp(g_wrk_2, pdvor(:,:,:,NDVOR), dcg)
    
    if (p_parallel_io) then
      !    p3 => g_wrk(:,:,:,NDVOR)
      !    ccv = corrsp(wrk,p3)
      ccv = corrsp (g_wrk,g_wrk_2)
    ENDIF

    if (associated(g_wrk)) deallocate(g_wrk)
    if (associated(g_wrk_2)) deallocate(g_wrk_2)

!   correlation for global field for divergence

    if (p_parallel_io) then
      ALLOCATE(g_wrk(nlev,2,nsp))
      ALLOCATE(g_wrk_2(nlev,2,nsp))
    else
      NULLIFY(g_wrk)
      NULLIFY(g_wrk_2)
    endif

    wrk = 0.0_dp
    DO i=1,NDDIV-1
      wrk(:,:,:) = wrk(:,:,:) + pddiv(:,:,:,i)
    END DO

    CALL gather_sp(g_wrk, wrk, dcg)
    CALL gather_sp(g_wrk_2, pddiv(:,:,:,NDDIV), dcg)
    
    if (p_parallel_io) then
      !    p3 => g_wrk(:,:,:,NDDIV)
      !    ccv = corrsp(wrk,p3)
      ccd = corrsp (g_wrk,g_wrk_2)
    ENDIF

    if (associated(g_wrk)) deallocate(g_wrk)
    if (associated(g_wrk_2)) deallocate(g_wrk_2)


!   correlation for global field for temperature

    if (p_parallel_io) then
      ALLOCATE(g_wrk(nlev,2,nsp))
      ALLOCATE(g_wrk_2(nlev,2,nsp))
    else
      NULLIFY(g_wrk)
      NULLIFY(g_wrk_2)
    endif

    wrk = 0.0_dp
    DO i=1,NDTEM-3
      wrk(:,:,:) = wrk(:,:,:) + pdtem(:,:,:,i)
    END DO

    CALL gather_sp(g_wrk, wrk, dcg)
    CALL gather_sp(g_wrk_2, pdtem(:,:,:,NDTEM), dcg)
    
    if (p_parallel_io) then
      !    p3 => g_wrk(:,:,:,NDTEM)
      !    ccv = corrsp(wrk,p3)
      cct = corrsp (g_wrk,g_wrk_2)
    ENDIF

    if (associated(g_wrk)) deallocate(g_wrk)
    if (associated(g_wrk_2)) deallocate(g_wrk_2)


!   correlation for global field for pressure

    if (p_parallel_io) then
      ALLOCATE(g_wrk(nlev,2,nsp))
      ALLOCATE(g_wrk_2(1,2,nsp))
    else
      NULLIFY(g_wrk)
      NULLIFY(g_wrk_2)
    endif

    wrk = 0.0_dp
    DO i=1,NDPRS-1
      wrk(1,:,:) = wrk(1,:,:) + pdprs(1,:,:,i)
    END DO

    CALL gather_sp(g_wrk, wrk, dcg)
    CALL gather_sp(g_wrk_2, pdprs(:,:,:,NDPRS), dcg)
    
    if (p_parallel_io) then
      p2a => g_wrk(1,:,:)
      p2b => g_wrk_2(1,:,:)
      ccp = corrsp (p2a,p2b)  
    ENDIF

    if (associated(g_wrk)) deallocate(g_wrk)
    if (associated(g_wrk_2)) deallocate(g_wrk_2)

    if (p_parallel_io) WRITE(*,*)  'mo_diag_tendency: ', &
      'VOR ',ccv,' DIV ',ccd,  ' TEM ',cct,' PRS ',ccp

  END SUBROUTINE DIAG_Check
!EOX
!EOC

  ! mz_ht_20050310+
  SUBROUTINE DIAG_SP2GP

    ! this subroutine should transform the spectral diagnostic
    ! tendency fields to the gridpoint tendency fields

    USE MESSY_MAIN_TRANSFORM_BI,   ONLY: sp2gp

    IMPLICIT NONE

    INTEGER :: kvor, kdiv, ktem, kpres

    do kvor = 1, ndvor
      call sp2gp(pdvor(:,:,:,kvor), diag_vor_ten_gp(:,:,:,kvor), .false. )
    enddo
    do kdiv = 1, nddiv
      call sp2gp(pddiv(:,:,:,kdiv), diag_div_ten_gp(:,:,:,kdiv), .false. )
    enddo
    do ktem = 1, ndtem
      call sp2gp(pdtem(:,:,:,ktem), diag_tem_ten_gp(:,:,:,ktem), .false. )
    enddo
    call sp2gp(pdprl(:,:,:), diag_pressl_ten_gp(:,:,:), .false. )
    do kpres = 1, ndprs
      call sp2gp(pdprs(:,:,:,kpres), diag_sfp_ten_gp(:,:,:,kpres), .false. )
    enddo

  END SUBROUTINE DIAG_SP2GP
  ! mz_ht_20050310-

END MODULE mo_diag_tendency
#endif
