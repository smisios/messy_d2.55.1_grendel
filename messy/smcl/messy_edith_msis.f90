MODULE messy_edith_msis

  !=============================================================================
  !
  !- Description:
  !
  !  This module allows to get MSIS model profiles interpolated to a
  !  specified prssure level grid.
  !
  !  This module contains:
  !  B) Data of MSIS model for global mean vertical structure of
  !     composition of middle atmosphere and lower thermosphere.
  !     Profiles are given at 45 levels from surface to top for:
  !     - atmospheric molar weight
  !     - volume mixing ratio of O, N2, O2
  !     - specific heat at constant pressure
  !     - gravity
  !     - gas 'constant' R
  !  C) Arrays to store interpolated MSIS profiles for specified pressure grid
  !     and related profiles
  !
  !  D) Subroutines to allocate and deallocate arrays
  !
  !  E) Subroutine to interpolate MSIS data to specified pressure grid
  !
  !  F) Interpolation function
  !
  !-----------------------------------------------------------------------------
  !
  !- Author:
  !
  !  V. Fomichev,    November, 1997: original source
  !  M.A. Giorgetta, MPI, June 2001: rewrite for ECHAM5
  !
  !=============================================================================

  !USE mo_kind     , ONLY: dp
  USE messy_main_constants_mem, ONLY: DP
  !USE mo_constants, ONLY: avo
  USE messy_main_constants_mem, ONLY: avo => N_A

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xo,                         & ! coordinate x=log(1000hPa/p)
       &    muo, oo, n2o, o2o, cpo, go, & ! MSIS climatology
       &    ro, aro,                    & ! MSIS climatology
       &    l_msis_alloc,               & ! status of allocation
       &    l_msis_def,                 & ! status of initialization
       &    msis_allocate,              & ! subroutine to allocate arrays
       &    msis_deallocate,            & ! subroutine to deallocate arrays
       &    msis_clim,                  & ! subroutine for MSIS climatology
       &    pvf,xvf,amuvf,ovf,cn2vf,    & ! MSIS climatology interpolated
       &    o2vf,cpvf,gvf,tox,tn2,to2,  & ! to specified pressure levels
       &    rvf,arvf
  
  !===========================================================================

  ! A) status indicators for allocation and initialization
  ! ======================================================

  LOGICAL               :: l_msis_alloc=.FALSE.
  LOGICAL               :: l_msis_def  =.FALSE.

  !===========================================================================

  ! B) MSIS model
  ! =============

  ! number of levels of MSIS model
  INTEGER , PARAMETER :: nmsis=45

  ! vertical coordinate xo=log(1000/po)
  REAL(dp), PARAMETER   :: xo(nmsis) =                                       &
       & (/  0.0,     0.5,     1.0,     1.5,     2.0,     2.30000, 2.55259,  &
       &     3.05259, 3.55259, 4.05259, 4.55259, 5.05259, 5.55259, 6.05259,  &
       &     6.55259, 7.05259, 7.55259, 8.05259, 8.55259, 9.05259, 9.55259,  &
       &    10.05260,10.55260,11.05260,11.55260,12.05260,12.55260,13.05260,  &
       &    13.55260,14.05260,14.55260,15.05260,15.55260,16.05260,16.55260,  &
       &    17.05260,17.55260,18.05260,18.55260,19.05260,19.55260,20.05260,  &
       &    20.55260,21.05260,21.30000/)

  ! atmospheric molar mass (g/mol)
  REAL(dp), PARAMETER   :: muo(nmsis) =                                      &
       & (/ 28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,     &
       &    28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,28.95,     &
       &    28.95,28.95,28.94,28.94,28.94,28.93,28.92,28.89,28.83,28.73,     &
       &    28.58,28.38,28.13,27.80,27.39,26.89,26.30,25.65,24.97,24.27,     &
       &    23.53,22.73,21.89,21.03,20.60/)

  ! O volume mixing ratio (m3/m3)
  REAL(dp), PARAMETER   :: oo(nmsis) =                                       &
       & (/ 0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,     &
       &    0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,     &
       &    0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,0.000e+00,     &
       &    0.000e+00,0.000e+00,0.000e+00,0.240e-08,0.127e-06,0.176e-05,     &
       &    0.177e-04,0.125e-03,0.615e-03,0.217e-02,0.575e-02,0.122e-01,     &
       &    0.219e-01,0.350e-01,0.520e-01,0.736e-01,0.101e+00,0.135e+00,     &
       &    0.177e+00,0.224e+00,0.275e+00,0.329e+00,0.388e+00,0.452e+00,     &
       &    0.519e+00,0.589e+00,0.623e+00/)

  ! N2 volume mixing ratio (m3/m3)
  REAL(dp), PARAMETER   :: n2o(nmsis) =                                      &
       & (/ 0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,     &
       &    0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,0.781,     &
       &    0.782,0.782,0.782,0.783,0.783,0.784,0.785,0.786,0.786,0.785,     &
       &    0.782,0.778,0.774,0.767,0.758,0.745,0.726,0.698,0.662,0.618,     &
       &    0.568,0.512,0.452,0.389,0.357/)

  ! O2 volume mixing ratio (m3/m3)
  REAL(dp), PARAMETER   :: o2o(nmsis) =                                      &
       & (/ 0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.210,     &
       &    0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.210,0.209,     &
       &    0.209,0.209,0.209,0.208,0.208,0.207,0.205,0.203,0.200,0.195,     &
       &    0.188,0.179,0.168,0.153,0.136,0.116,0.094,0.076,0.061,0.051,     &
       &    0.043,0.036,0.029,0.022,0.020/)

  ! Ar volume mixing ratio (m3/m3)
  REAL(dp), PARAMETER   :: aro(nmsis) =                                      &
       & (/ 9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,          &
       &    9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,          &
       &    9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,          &
       &    9.3432E-03,9.3432E-03,9.3432E-03,9.3432E-03,9.3082E-03,          &
       &    9.2436E-03,9.1793E-03,9.1279E-03,9.0820E-03,9.0250E-03,          &
       &    8.9516E-03,8.8531E-03,8.7139E-03,8.5128E-03,8.2260E-03,          &
       &    7.8307E-03,7.3092E-03,6.6546E-03,5.8799E-03,5.0228E-03,          &
       &    4.1379E-03,3.2804E-03,2.4976E-03,1.8317E-03,1.3070E-03,          &
       &    9.1593E-04,6.3086E-04,4.2447E-04,2.7725E-04,2.2180E-04/)

  ! cp for dry air (J/K/kg)
  REAL(dp), PARAMETER   :: cpo(nmsis) =                                      &
       & (/ 1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,  &
       &    1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,1002.5,  &
       &    1002.5,1002.6,1002.6,1002.7,1002.8,1002.9,1003.0,1003.2,1003.5,  &
       &    1004.1,1005.2,1006.9,1009.4,1012.8,1017.3,1023.0,1030.3,1039.2,  &
       &    1049.7,1061.3,1073.3,1085.8,1099.4,1114.8,1132.1,1151.0,1160.9/)

  ! gravity g (m/s2)
  REAL(dp), PARAMETER   :: go(nmsis) =                                       &
       & (/ 9.8,9.79,9.78,9.77,9.765,9.76,9.75,9.74,9.73,9.72,9.71,9.70,     &
       &    9.69,9.68,9.67,9.66,9.64,9.63,9.62,9.61,9.60,9.59,9.58,9.57,     &
       &    9.56,9.55,9.55,9.54,9.53,9.52,9.51,9.51,9.50,9.49,9.48,9.46,     &
       &    9.45,9.43,9.40,9.37,9.33,9.29,9.23,9.18,9.15/)

  ! gas 'constant' (J/kg/K)
  REAL(dp), PARAMETER   :: ro(nmsis) =                                       &
       & (/ 288.33,288.33,288.33,288.33,288.33,288.33,288.33,288.33,         &
       &    288.33,288.33,288.33,288.33,288.33,288.33,288.33,288.33,         &
       &    288.33,288.33,288.33,288.34,288.35,288.36,288.38,288.40,         &
       &    288.42,288.48,288.65,289.08,290.02,291.68,294.16,297.51,         &
       &    301.85,307.37,314.39,323.11,333.49,345.06,357.24,370.00,         &
       &    383.76,398.63,414.36,430.39,438.23/)

  !===========================================================================

  ! C) storage arrays for interpolated MSIS profiles and derived profiles
  ! =====================================================================

  REAL(dp), ALLOCATABLE :: pvf(:)   ! pressure grid (hPa, surface to top)
  REAL(dp), ALLOCATABLE :: xvf(:)   ! x=log(1000hPa/pvf)
  REAL(dp), ALLOCATABLE :: amuvf(:) ! muo interpolated to xvf
  REAL(dp), ALLOCATABLE :: ovf(:)   ! oo  interpolated to xvf
  REAL(dp), ALLOCATABLE :: cn2vf(:) ! n2o interpolated to xvf
  REAL(dp), ALLOCATABLE :: o2vf(:)  ! o2o interpolated to xvf
  REAL(dp), ALLOCATABLE :: cpvf(:)  ! cpo interpolated to xvf
  REAL(dp), ALLOCATABLE :: rvf(:)   ! ro  interpolated to xvf
  REAL(dp), ALLOCATABLE :: gvf(:)   ! go  interpolated to xvf
  REAL(dp), ALLOCATABLE :: arvf(:)  ! aro interpolated to xvf
  REAL(dp), ALLOCATABLE :: tox(:)   ! O  number density integrated from top
  REAL(dp), ALLOCATABLE :: tn2(:)   ! N2 number density integrated from top
  REAL(dp), ALLOCATABLE :: to2(:)   ! O2 number density integrated from top

CONTAINS

  !=============================================================================

  ! D) subroutines to allocate and deallocate arrays
  ! ================================================

  SUBROUTINE msis_allocate(klev)

    INTEGER, INTENT(in)  :: klev

    if (.not. allocated(pvf)) ALLOCATE(pvf(klev))
    if (.not. allocated(xvf)) ALLOCATE(xvf(klev))
    if (.not. allocated(amuvf)) ALLOCATE(amuvf(klev))
    if (.not. allocated(ovf)) ALLOCATE(ovf(klev))
    if (.not. allocated(cn2vf)) ALLOCATE(cn2vf(klev))
    if (.not. allocated(o2vf)) ALLOCATE(o2vf(klev))
    if (.not. allocated(cpvf)) ALLOCATE(cpvf(klev))
    if (.not. allocated(rvf)) ALLOCATE(rvf(klev))
    if (.not. allocated(gvf)) ALLOCATE(gvf(klev))
    if (.not. allocated(arvf)) ALLOCATE(arvf(klev))
    if (.not. allocated(tox)) ALLOCATE(tox(klev))
    if (.not. allocated(tn2)) ALLOCATE(tn2(klev))
    if (.not. allocated(to2)) ALLOCATE(to2(klev))

    l_msis_alloc=.TRUE.

  END SUBROUTINE msis_allocate


  SUBROUTINE msis_deallocate

    DEALLOCATE(pvf)
    DEALLOCATE(xvf)
    DEALLOCATE(amuvf)
    DEALLOCATE(ovf)
    DEALLOCATE(cn2vf)
    DEALLOCATE(o2vf)
    DEALLOCATE(cpvf)
    DEALLOCATE(rvf)
    DEALLOCATE(gvf)
    DEALLOCATE(arvf)
    DEALLOCATE(tox)
    DEALLOCATE(tn2)
    DEALLOCATE(to2)

    l_msis_alloc=.FALSE.

  END SUBROUTINE msis_deallocate

  !=============================================================================

  ! E) compute MSIS climatology for specified grid
  ! ==============================================

  SUBROUTINE msis_clim(ppf)

    ! full level pressure in Pa from top to surface
    REAL(dp), INTENT(in) :: ppf(:)

    ! log of MSIS profiles
    REAL(dp) :: zlogmuo(nmsis)
    REAL(dp) :: zlogoo(nmsis)
    REAL(dp) :: zlogn2o(nmsis)
    REAL(dp) :: zlogo2o(nmsis)
    REAL(dp) :: zlogcpo(nmsis)
    REAL(dp) :: zlogro(nmsis)
    REAL(dp) :: zloggo(nmsis)
    REAL(dp) :: zlogaro(nmsis)

    INTEGER  :: klev
    INTEGER  :: jk

    REAL(dp) :: ana, g1, g2

    klev=SIZE(ppf)

    ! p coordinate in hPa, from surface to top
    DO jk=1,klev
       pvf(jk)=ppf(klev+1-jk)/100._dp
    END DO

    ! log-p coordinate
    xvf(:) = LOG(1000._dp/pvf(:))

    ! log(MSIS profiles)
    zlogmuo(:)= LOG(muo(:))
    zlogoo(:) = LOG(oo(:)+1._dp)
    zlogn2o(:)= LOG(n2o(:))
    zlogo2o(:)= LOG(o2o(:))
    zlogcpo(:)= LOG(cpo(:))
    zlogro(:) = LOG(ro(:))
    zloggo(:) = LOG(go(:))
    zlogaro(:)= LOG(aro(:))

    ! interpolate MSIS profiles to specified pressure levels
    DO jk=1,SIZE(ppf)
       amuvf(jk)= EXP(a18lin(xvf(jk),xo(:),zlogmuo(:),1,nmsis))
       ovf(jk)  = EXP(a18lin(xvf(jk),xo(:),zlogoo(:) ,1,nmsis))-1._dp
       cn2vf(jk)= EXP(a18lin(xvf(jk),xo(:),zlogn2o(:),1,nmsis))
       o2vf(jk) = EXP(a18lin(xvf(jk),xo(:),zlogo2o(:),1,nmsis))
       cpvf(jk) = EXP(a18lin(xvf(jk),xo(:),zlogcpo(:),1,nmsis))
       rvf(jk)  = EXP(a18lin(xvf(jk),xo(:),zlogro(:) ,1,nmsis))
       gvf(jk)  = EXP(a18lin(xvf(jk),xo(:),zloggo(:) ,1,nmsis))
       arvf(jk) = EXP(a18lin(xvf(jk),xo(:),zlogaro(:),1,nmsis))
    END DO

    ! compute column number densities (molecules/cm2)

    ana=avo*10._dp

    ! to calculate tO2(cm-2)

    g1=o2vf(klev)*ana*pvf(klev)/(amuvf(klev)*gvf(klev))
    to2(klev)=g1*amuvf(klev)/32._dp

    DO jk=klev-1,1,-1
       g2=o2vf(jk)*ana*pvf(jk)/(amuvf(jk)*gvf(jk))
       to2(jk)=to2(jk+1)+(g1+g2)*0.5_dp*(xvf(jk+1)-xvf(jk))
       g1=g2
    END DO

    ! to calculate tN2(cm-2)

    g1=cn2vf(klev)*ana*pvf(klev)/(amuvf(klev)*gvf(klev))
    tn2(klev)=g1*amuvf(klev)/28._dp

    DO jk=klev-1,1,-1
       g2=cn2vf(jk)*ana*pvf(jk)/(amuvf(jk)*gvf(jk))
       tn2(jk)=tn2(jk+1)+(g1+g2)*0.5_dp*(xvf(jk+1)-xvf(jk))
       g1=g2
    END DO

    ! to calculate tOX(cm-2)

    g1=ovf(klev)*ana*pvf(klev)/(amuvf(klev)*gvf(klev))
    tox(klev)=g1*amuvf(klev)/16._dp

    DO jk=klev-1,1,-1
       g2=ovf(jk)*ana*pvf(jk)/(amuvf(jk)*gvf(jk))
       tox(jk)=tox(jk+1)+(g1+g2)*0.5_dp*(xvf(jk+1)-xvf(jk))
       g1=g2
    END DO

    l_msis_def=.TRUE.

  END SUBROUTINE msis_clim


  ! F) linear interpolation function
  ! ================================

  FUNCTION a18lin(x,xn,yn,m,n)

    ! linear interpolation
    ! input:
    !  X - argument for which a value of function should be found
    !  XN(N),YN(N) - values of function YN(N) at XN(N) grid. X(N) should be
    !                ordered so that X(I-1) < X(I).
    ! output:
    !  A18LIN - value of function for X
    !
    ! V. Fomichev, original source

    ! modification:
    !  The interpolation y(x)is done only if xn(m)<x<xn(n),
    !  otherwise y is set to yn(m) if x<xn(m) or yn(n) if x>xn(n)
    !
    ! M.A. Giorgetta


    INTEGER , INTENT(in) :: m, n   ! first and last index of array range
    ! to be used for interpolation

    REAL(dp), INTENT(in) :: xn(n)  ! coordinate array xn
    REAL(dp), INTENT(in) :: yn(n)  ! value array yn=y(xn)

    REAL(dp), INTENT(in) :: x      ! new coordinate x
    REAL(dp)             :: y      ! interpolated value y(x)
    REAL(dp)             :: a18lin ! function result

    INTEGER :: i,k

    IF(x<=xn(m))      THEN         ! y=yn(m) if x<=xn(m)
       y=yn(m)
    ELSE IF(x>=xn(n)) THEN         ! y=yn(n) if x>=xn(n)
       y=yn(n)
    ELSE                           ! y is interpolated if xn(m)<x<xn(n)
       k=m-1
       DO i=m,n
          k=k+1
          IF(x<=xn(i)) EXIT
       END DO
       IF(k == 1) k=2

       ! k has been found so that xn(k).le.x.lt.xn(k+1)

       y=(yn(k)-yn(k-1))/(xn(k)-xn(k-1))*(x-xn(k))+yn(k)
    END IF

    a18lin=y

  END FUNCTION a18lin

  !=============================================================================

END MODULE messy_edith_msis
