! **********************************************************************
MODULE messy_cloudopt
! **********************************************************************

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_tools,         ONLY: t_reset_par

  PRIVATE
  SAVE

  PUBLIC :: DP

  INTEGER, PARAMETER, PUBLIC :: JPBAND = 16 ! number of bands of LW scheme
  INTEGER, PARAMETER, PUBLIC :: NSW    =  4 ! number of bands of SW scheme

  ! correction for asymmetry factor of ice clouds
  REAL(dp), PUBLIC :: asic  
  ! cloud inhom. factor   
  REAL(DP), PUBLIC :: zinhomi
  REAL(DP) :: zinpar
  ! Switch to reset the corresponding parameters
  TYPE(t_reset_par), PUBLIC :: rset_zinhoml = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par), PUBLIC :: rset_zinhomi = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par), PUBLIC :: rset_zinpar  = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par), PUBLIC :: rset_asic    = t_reset_par(.FALSE.,0.0_dp)

  LOGICAL          , PUBLIC :: loice = .TRUE.

  REAL(dp),PARAMETER :: diff  =1.66_dp   ! diffusivity factor

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'cloudopt'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1b'

  PUBLIC :: calc_clcv
  PUBLIC :: cloud_set_param
  PUBLIC :: cloudopt_read_nml_ctrl
  PUBLIC :: rad_std

  PUBLIC :: cloud_opt_prelwsw ! routine to calculate input for 
                              ! final lw/sw computation
  PUBLIC :: cloud_opt_lw1     ! only final lw computation
  PUBLIC :: cloud_opt_sw1     ! only final sw computation

CONTAINS

  ! =========================================================================
  SUBROUTINE calc_clcv(kproma, kbdim, klev, pclfr, clcv)

    IMPLICIT NONE
    INTRINSIC :: MAX, MIN, EPSILON

    ! I/O
    INTEGER, INTENT(IN) :: kproma, kbdim, klev
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: pclfr 
    REAL(DP), DIMENSION(kbdim),       INTENT(OUT) :: clcv

    ! LOCAL
    INTEGER :: jk

    ! total cloud cover
    ! EPSILON(1.) avoids 0/0 in the diagnostic of total cloud cover.
    clcv(:) = 1._dp - pclfr(:,1)
    DO jk = 2, klev
       clcv(1:kproma) = clcv(1:kproma)                              &
            *(1._dp-MAX(pclfr(1:kproma,jk), pclfr(1:kproma,jk-1)))  &
            /(1._dp-MIN(pclfr(1:kproma,jk-1),1._dp-EPSILON(1._dp)))
    END DO
    clcv(1:kproma) = 1._dp-clcv(1:kproma)

  END SUBROUTINE calc_clcv
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE cloud_opt_prelwsw(kproma,kbdim,klev, &  !IN
       pclfr,totclfr, pxl,pxi,                &  !IN 
       zdp,                                   &  !IN
       lcouple,                               &  !IN 
       icldlyr,                               &  !OUT
       zlwp,ziwp,                             &  !OUT  !RENAME????
       zinhoml)                                  !OUT
     

    USE messy_main_constants_mem, ONLY:  g

    IMPLICIT NONE
    INTRINSIC :: MAX, INT, EPSILON

    ! I/O
    INTEGER,                          INTENT(IN)  :: kproma,kbdim,klev
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: pclfr
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: totclfr
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: pxl
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: pxi
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zdp ! use directly as it is calculated in si
    LOGICAL,                          INTENT(IN)  :: lcouple
    ! OUTPUT variables
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(OUT) :: zlwp    ! needed for zlwpt/ lw and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(OUT) :: ziwp    ! needed for lw and sw
    REAL(DP), DIMENSION(kbdim),       INTENT(OUT) :: zinhoml ! needed in lw and sw comp
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(OUT) :: icldlyr
    
    ! LOCAL
    ! cloud optical properties
    INTEGER :: jk,jl
    REAL(DP), DIMENSION(kbdim,klev)  :: zclwa   ! for zlwgkg
    REAL(DP), DIMENSION(kbdim,klev)  :: zclic   ! pre; for ziwgkg
    REAL(DP), DIMENSION(kbdim,klev)  :: zclfr   ! needed pre ziwgkg and zlwgkg

    REAL(DP), DIMENSION(kbdim,klev)  :: zlwgkg  ! needed for zlwp
    REAL(DP), DIMENSION(kbdim,klev)  :: ziwgkg  ! pre; needed for ziwp
    REAL(DP), DIMENSION(kbdim)       :: zlwpt   ! needed to calc zinhoml

! Calculate diverse cloud properties liquid water path/ice water path....etc (for zinhoml)

    ! cloud liquid water
    zclwa(1:kproma,:)=MAX(pxl(1:kproma,:),0._dp)

    ! cloud ice
    zclic(1:kproma,:)=MAX(pxi(1:kproma,:),0._dp)

    zclfr(1:kproma,:)=pclfr(1:kproma,:)

    ! Clear/cloudy index
    WHERE (zclfr(1:kproma,:)>EPSILON(1.0_dp))
       icldlyr(1:kproma,:)=1._dp
    ELSEWHERE
       icldlyr(1:kproma,:)=0._dp
    END WHERE

    ! Secure cloud fraction
    zclfr(1:kproma,:)=MAX(pclfr(1:kproma,:),EPSILON(1.0_dp))

    WHERE (INT(icldlyr(1:kproma,:))==1)
       ! Specific ice water content, g/kg
       ziwgkg(1:kproma,:)=zclic(1:kproma,:)*1000._dp/zclfr(1:kproma,:)
       ! Ice water path, g/m2
       ziwp(1:kproma,:)=(ziwgkg(1:kproma,:)*zdp(1:kproma,:)/g)* &
            (pclfr(1:kproma,:)/totclfr(1:kproma,:))
       ! Specific liquid water content, g/kg
       zlwgkg(1:kproma,:)=zclwa(1:kproma,:)*1000._dp/zclfr(1:kproma,:)
       ! Liquid water path, g/m2
       zlwp(1:kproma,:)=(zlwgkg(1:kproma,:)*zdp(1:kproma,:)/g)* &
            (pclfr(1:kproma,:)/totclfr(1:kproma,:))
    ELSEWHERE
       ! Specific ice water content, g/kg
       ziwgkg(1:kproma,:)=0._dp
       ! Ice water path, g/m2
       ziwp(1:kproma,:)=0._dp
       ! Specific liquid water content, g/kg
       zlwgkg(1:kproma,:)=0._dp
       ! Liquid water path, g/m2
       zlwp(1:kproma,:) = 0.0_dp
    END WHERE

    ! Initialisation of asic has moved to messy_rad_e5:rad_set_param
    ! to enable namelist input.
    ! fb_mk_20120917-

    IF(lcouple) THEN
       zinhoml(:) = 0.7_dp
    ELSE
       zlwpt(:) = 0._dp
       DO jk = 1, klev
          DO jl = 1, kproma
             zlwpt(jl) = zlwpt(jl)+zlwp(jl,jk)
          END DO
       END DO
       DO jl = 1, kproma
          IF(zlwpt(jl)>1._dp) THEN
             zinhoml(jl) = zlwpt(jl)**(-zinpar)
          ELSE
             zinhoml(jl) = 1.0_dp
          END IF
       END DO
    END IF

    ! Overwrite the initialization of zinhoml with the
    ! the value provided by the namelist
    IF (rset_zinhoml%l) THEN
       zinhoml = rset_zinhoml%v
    END IF

  END SUBROUTINE cloud_opt_prelwsw
  ! =========================================================================

 ! =========================================================================
  SUBROUTINE cloud_opt_lw1(kproma,kbdim,klev, &  !IN Q? REMOVE kbdim in the end
       zlwp,ziwp,                              &  !IN    
       zradlp, zradip,                         &  !IN
       zinhoml,                                &  !IN
       taucld_lw)                                 !OUT 

    USE messy_main_constants_mem, ONLY: ccwmin

    IMPLICIT NONE
    INTRINSIC :: EXP

    ! I/O
    INTEGER,                          INTENT(IN)  :: kproma,kbdim,klev

    ! cloud optical properties
   
    ! INPUT 
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zlwp    ! needed for zlwpt/ lw and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: ziwp    ! needed for lw and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zradlp  ! needed for zmacl(LW) and SW
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zradip  ! needed for zmaci(LW) and SW
    REAL(DP), DIMENSION(kbdim),       INTENT(IN)  :: zinhoml ! needed in lw and sw comp

    ! start OUT lw variables
    REAL(DP), DIMENSION(kbdim,klev,jpband), INTENT(OUT) :: taucld_lw

    ! start parameters lw computation

    !  Coefficients for the parameterization of the emissivity
    !  of cloud particles. Coefficients are given for cloud liquid droplets
    !  and ice crystals.
    !
    !  Coefficients are given for 16 band resolution of the LW spectrum
    !  as used by the RRTM scheme.
    !
    !  Reference : follows ECMWF-CY23R1,
    !              Smith and Shi (1992), Ebert and Curry (1992)

    ! repcug = C7 in mass absorption coefficient formula
    REAL(dp), PARAMETER, DIMENSION(jpband) :: & 
         rebcug = (/ 0.718_dp, 0.726_dp, 1.136_dp, 1.320_dp, &
         &           1.505_dp, 1.290_dp, 0.911_dp, 0.949_dp, &
         &           1.021_dp, 1.193_dp, 1.279_dp, 0.626_dp, &
         &           0.647_dp, 0.668_dp, 0.690_dp, 0.690_dp /)

    ! repcuh = C8 in mass absorption coefficient formula
    REAL(dp), PARAMETER, DIMENSION(jpband) :: & 
         rebcuh = (/ 0.0069_dp, 0.0060_dp, 0.0024_dp, 0.0004_dp, &
         &          -0.0016_dp, 0.0003_dp, 0.0043_dp, 0.0038_dp, &
         &           0.0030_dp, 0.0013_dp, 0.0005_dp, 0.0054_dp, &
         &           0.0052_dp, 0.0050_dp, 0.0048_dp, 0.0048_dp /)
    
    ! parameters lw computation


    !LW + SW

    INTEGER :: jk,jl,jb

   ! LW-Calculated
    REAL(DP), DIMENSION(kbdim,klev)  ::zmacl  ! LW
    REAL(DP), DIMENSION(kbdim,klev)  ::zmaci  ! LW    
    REAL(DP) :: zmsald,zmsaid   ! LW

    !========================================================================
    !           LW
    !========================================================================

    DO jk=1,klev
       DO jl=1,kproma
          zmacl(jl,jk)=0.025520637_dp   &
               +0.2854650784_dp*EXP(-0.088968393014_dp*zradlp(jl,jk))
          zmaci(jl,jk)=0.020219423_dp   &
               +0.2058619832_dp*EXP(-0.067631070625_dp*zradip(jl,jk))
       END DO
    END DO

    DO jb=1,jpband
       DO jk=1,klev
          DO jl=1,kproma
             IF (zlwp(jl,jk)+ziwp(jl,jk)>ccwmin) THEN
                zmsald=zmacl(jl,jk)
                IF (loice) THEN
                   ! ice cloud emissivity after Ebert and Curry (1992)
                   ! with diffusivity factor diff
                   zmsaid=(rebcuh(jb)+rebcug(jb)/zradip(jl,jk))*diff
                ELSE
                   ! ice cloud emissivity after Rockel et al. (1991)
                   zmsaid=zmaci(jl,jk)
                END IF
                ! combine
                taucld_lw(jl,jk,jb)=zmsald*zlwp(jl,jk)*zinhoml(jl) &
                     +zmsaid*ziwp(jl,jk)*zinhomi
             ELSE
                taucld_lw(jl,jk,jb)=0._dp
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE cloud_opt_lw1
  ! =========================================================================

 ! =========================================================================
  SUBROUTINE cloud_opt_sw1(kproma,kbdim,klev,  &  !IN
       icldlyr,                                &  !IN          
       zlwp,ziwp,                              &  !IN    
       zradlp, zradip,                         &  !IN
       zinhoml,                                &  !IN
       tau_sw,tau_sw_raw,omega_sw,gamma_sw, &    ! OUT 
       valid)                                    ! OUT

    USE messy_main_constants_mem, ONLY:  ccwmin

    IMPLICIT NONE
    INTRINSIC :: INT

       ! cloud optical properties
    INTEGER,                          INTENT(IN)  :: kproma,kbdim,klev   

    ! INPUT 
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: icldlyr ! needed for sw calculation
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zlwp    ! needed for zlwpt/ lw and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: ziwp    ! needed for lw and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zradlp  ! needed for zmacl(lW) and sw
    REAL(DP), DIMENSION(kbdim,klev),  INTENT(IN)  :: zradip  ! needed for zmaci(lW) and sw
    REAL(DP), DIMENSION(kbdim),       INTENT(IN)  :: zinhoml ! needed in lw and sw comp

    ! start OUT sw variables
    REAL(DP), DIMENSION(kbdim,klev,nsw)   , INTENT(OUT) :: tau_sw
    REAL(DP), DIMENSION(kbdim,klev,nsw)   , INTENT(OUT) :: tau_sw_raw
    REAL(DP), DIMENSION(kbdim,klev,nsw)   , INTENT(OUT) :: omega_sw
    REAL(DP), DIMENSION(kbdim,klev,nsw)   , INTENT(OUT) :: gamma_sw
    REAL(DP), DIMENSION(kbdim,klev)       , INTENT(OUT) :: valid
    ! end OUT sw variables

!SW


    ! start parameters lw computation
    ! note that only 4 sw bands are possible

    !  Coefficients for the parameterization of the extinction (tau),
    !  single scattering albedo (omega) and asymmetry factor (gamma)
    !  of cloud particles. Coefficients are given for cloud liquid droplets
    !  and ice crystals.
    !
    !  Coefficients are given for 2 or 4 band resolution of the SW spectrum.
    !  Reference : Rockel et al. (1991).
    !  Note: only 4 bands are used!
    !
    !  2 bands: as in ECHAM4
    !           visible : band 1 : 0.25 - 0.68 micrometer
    !           near IR : band 2 : 0.68 - 4.0  micrometer
    !
    !  4 bands: as in ECMWF model
    !           visible : band 1 : 0.25 - 0.69 micrometer
    !           near IR : band 2 : 0.69 - 1.19 micrometer
    !                     band 3 : 1.19 - 2.38 micrometer
    !                     band 4 : 2.38 - 4.00 micrometer
    !
    !  Coefficients: xNBPI, N = 2 or 4 for SW resolution
    !                       B = 1,..N band index
    !                       P = 1 for liquid, 2 for ice phase
    !                       I = coefficient index
    !
    !  extinction  : x=a , single scattering : x=b , asymmetry factor  : x=c

    !  4 bands
    !   extinction
    !    band 1
    !     liquid
    REAL(dp), PARAMETER :: a41l0= 1.8362_dp
    REAL(dp), PARAMETER :: a41l1=-1.0665_dp
    !     ice
    REAL(dp), PARAMETER :: a41i0= 1.9787_dp
    REAL(dp), PARAMETER :: a41i1=-1.0365_dp
    !    band 2
    !     liquid
    REAL(dp), PARAMETER :: a42l0= 2.0731_dp
    REAL(dp), PARAMETER :: a42l1=-1.1079_dp
    !     ice
    REAL(dp), PARAMETER :: a42i0= 2.1818_dp
    REAL(dp), PARAMETER :: a42i1=-1.0611_dp
    !    band 3
    !     liquid
    REAL(dp), PARAMETER :: a43l0= 1.8672_dp
    REAL(dp), PARAMETER :: a43l1=-1.0420_dp
    !     ice
    REAL(dp), PARAMETER :: a43i0= 1.9608_dp
    REAL(dp), PARAMETER :: a43i1=-1.0212_dp
    !    band 4
    !     liquid
    REAL(dp), PARAMETER :: a44l0= 1.0787_dp
    REAL(dp), PARAMETER :: a44l1=-0.79772_dp
    !     ice
    REAL(dp), PARAMETER :: a44i0= 1.2558_dp
    REAL(dp), PARAMETER :: a44i1=-0.88622_dp
    !   single scattering
    !    band 1
    !     liquid
    REAL(dp), PARAMETER :: b41l0= 1._dp
    REAL(dp), PARAMETER :: b41l1=-2.2217e-7_dp
    !     ice
    REAL(dp), PARAMETER :: b41i0= 1._dp
    REAL(dp), PARAMETER :: b41i1=-1.143e-7_dp
    !    band 2
    !     liquid
    REAL(dp), PARAMETER :: b42l0= 1._dp
    REAL(dp), PARAMETER :: b42l1=-1.6712e-5_dp
    !     ice
    REAL(dp), PARAMETER :: b42i0= 0.99999_dp
    REAL(dp), PARAMETER :: b42i1=-7.9238e-6_dp
    !    band 3
    !     liquid
    REAL(dp), PARAMETER :: b43l0= 0.99936_dp
    REAL(dp), PARAMETER :: b43l1=-0.0013632_dp
    !     ice
    REAL(dp), PARAMETER :: b43i0= 0.99975_dp
    REAL(dp), PARAMETER :: b43i1=-0.001662_dp
    REAL(dp), PARAMETER :: b43i2= 6.9726e-6_dp
    !    band 4
    !     liquid
    REAL(dp), PARAMETER :: b44l0= 0.90032_dp
    REAL(dp), PARAMETER :: b44l1=-0.091955_dp
    !     ice
    REAL(dp), PARAMETER :: b44i0= 0.89779_dp
    REAL(dp), PARAMETER :: b44i1=-0.0802_dp
    !   asymmetry factor
    !    band 1
    !     liquid
    REAL(dp), PARAMETER :: c41l0= 0.78063_dp
    REAL(dp), PARAMETER :: c41l1= 0.12600_dp
    REAL(dp), PARAMETER :: c41l2=-0.042412_dp
    !     ice
    REAL(dp), PARAMETER :: c41i0= 0.79602_dp
    REAL(dp), PARAMETER :: c41i1= 0.10183_dp
    REAL(dp), PARAMETER :: c41i2=-0.028648_dp
    !    band 2
    !     liquid
    REAL(dp), PARAMETER :: c42l0= 0.74102_dp
    REAL(dp), PARAMETER :: c42l1= 0.16315_dp
    REAL(dp), PARAMETER :: c42l2=-0.050268_dp
    !     ice
    REAL(dp), PARAMETER :: c42i0= 0.77176_dp
    REAL(dp), PARAMETER :: c42i1= 0.11995_dp
    REAL(dp), PARAMETER :: c42i2=-0.030557_dp
    !    band 3
    !     liquid
    REAL(dp), PARAMETER :: c43l0= 0.70730_dp
    REAL(dp), PARAMETER :: c43l1= 0.18299_dp
    REAL(dp), PARAMETER :: c43l2=-0.045693_dp
    !     ice
    REAL(dp), PARAMETER :: c43i0= 0.74691_dp
    REAL(dp), PARAMETER :: c43i1= 0.13514_dp
    REAL(dp), PARAMETER :: c43i2=-0.027140_dp
    !    band 4
    !     liquid
    REAL(dp), PARAMETER :: c44l0= 0.70554_dp
    REAL(dp), PARAMETER :: c44l1= 0.88798_dp
    REAL(dp), PARAMETER :: c44l2=-1.8214_dp
    REAL(dp), PARAMETER :: c44l3= 1.5775_dp
    REAL(dp), PARAMETER :: c44l4=-0.46293_dp
    !     ice
    REAL(dp), PARAMETER :: c44i0= 0.77614_dp
    REAL(dp), PARAMETER :: c44i1= 0.15186_dp
    REAL(dp), PARAMETER :: c44i2=-0.031452_dp
    ! end parameters lw computation


    !LW + SW

    INTEGER :: jk,jl,jb
  
! SW   
    REAL(DP) :: zasic           ! SW
    REAL(DP) :: zrl,zri,zlp,zip ! SW
    REAL(DP), DIMENSION(nsw):: ztol,ztoi,ztaumx ! SW
    REAL(DP), DIMENSION(nsw):: zol,zoi,zomgmx   ! SW
    REAL(DP), DIMENSION(nsw):: zgl,zgi,zasymx   ! SW



    ! Initialisation of asic has moved to messy_rad_e5:rad_set_param
    ! to enable namelist input.
    zasic = asic

    !========================================================================
    !           SW
    !========================================================================
    ! cloud optical properties: extinction tau,
    !                           single scattering albedo omega,
    !                           asymmetry factor gamma

    ! a) set values for cloud free conditions
    tau_sw(1:kbdim,1:klev,1:nsw)     = 0._dp
    tau_sw_raw(1:kbdim,1:klev,1:nsw) = 0._dp
    omega_sw(1:kbdim,1:klev,1:nsw)   = 0._dp
    gamma_sw(1:kbdim,1:klev,1:nsw)   = 0._dp
    valid(1:kbdim,1:klev)            = 0._dp

    ! b) compute values for cloudy grid boxes
    DO jk=1,klev
       DO jl=1,kproma
          IF (INT(icldlyr(jl,jk))==1 .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN
             ! effective radii
             zrl=zradlp(jl,jk)
             zri=zradip(jl,jk)
             ! log10 of effective radii
             zlp=LOG10(zrl)
             zip=LOG10(zri)
             ! optical depth tau of cloud layer
             ! extinction is approximated by exponential of reff
             ! tau liquid
             ztol(1)=zlwp(jl,jk)*a41l0*zrl**a41l1
             ztol(2)=zlwp(jl,jk)*a42l0*zrl**a42l1
             ztol(3)=zlwp(jl,jk)*a43l0*zrl**a43l1
             ztol(4)=zlwp(jl,jk)*a44l0*zrl**a44l1
             ! tau ice
             ztoi(1)=ziwp(jl,jk)*a41i0*zri**a41i1
             ztoi(2)=ziwp(jl,jk)*a42i0*zri**a42i1
             ztoi(3)=ziwp(jl,jk)*a43i0*zri**a43i1
             ztoi(4)=ziwp(jl,jk)*a44i0*zri**a44i1
             ! single scattering albedo omega
             ! approximated by polynomials of reff for bands 1 to 3
             ! approximated by exponential of reff for band 4
             ! omega liquid
             zol(1)=b41l0+zrl*b41l1
             zol(2)=b42l0+zrl*b42l1
             zol(3)=b43l0+zrl*b43l1
             zol(4)=b44l0*zrl**b44l1
             !       ! omega ice
             zoi(1)=b41i0+zri*b41i1
             zoi(2)=b42i0+zri*b42i1
             zoi(3)=b43i0+zri*(b43i1+zri*b43i2)
             zoi(4)=b44i0*zri**b44i1
             ! asymmetry factor gamma
             ! approximated  by polynomials of log10(reff)
             ! gamma liquid: approximated by polynomials of log10(reff)
             zgl(1)=c41l0+zlp*(c41l1+zlp*c41l2)
             zgl(2)=c42l0+zlp*(c42l1+zlp*c42l2)
             zgl(3)=c43l0+zlp*(c43l1+zlp*c43l2)
             zgl(4)=c44l0+zlp*(c44l1+zlp*(c44l2+zlp*(c44l3+zlp*c44l4)))
             ! gamma ice with empirical factor zasic
             zgi(1)=zasic*(c41i0+zip*(c41i1+zip*c41i2))
             zgi(2)=zasic*(c42i0+zip*(c42i1+zip*c42i2))
             zgi(3)=zasic*(c43i0+zip*(c43i1+zip*c43i2))
             zgi(4)=zasic*(c44i0+zip*(c44i1+zip*c44i2))
             DO jb=1,nsw
                ztaumx(jb)=ztol(jb)+ztoi(jb)   ! tau liquid and ice
                zomgmx(jb)=ztol(jb)*zol(jb)+ztoi(jb)*zoi(jb) 
                zasymx(jb)=ztol(jb)*zol(jb)*zgl(jb)+ztoi(jb)*zoi(jb)*zgi(jb) 
! Note: this is uncommented to allow aggregaton of perturbation
!       contributions outside in SMIL:
!                zasymx(jb)=zasymx(jb)/zomgmx(jb)    ! gamma liquid and ice
!                zomgmx(jb)=zomgmx(jb)/ztaumx(jb)    ! omega liquid and ice
                ! ------ input arrays for SW ----------------------
                tau_sw(jl,jk,jb)=ztol(jb)*zinhoml(jl)+ztoi(jb)*zinhomi
                tau_sw_raw(jl,jk,jb)=ztaumx(jb)
                omega_sw(jl,jk,jb)=zomgmx(jb)
                gamma_sw(jl,jk,jb)=zasymx(jb)
             END DO
             valid(jl,jk)=1.0_dp
          END IF
       END DO
    END DO

  END SUBROUTINE cloud_opt_sw1
  ! =========================================================================

  ! =========================================================================
  ELEMENTAL SUBROUTINE rad_std(                &
       ppf, pdp, ptf,                          &   ! IN
       pxi, pxl, pclfr, pcdnc,                 &   ! IN
       laland, laglac,                         &   ! IN
       radip, radlp                            &   ! OUT 
       )

    USE messy_main_constants_mem, ONLY: rd, ceffmin, ceffmax, pi, rho_H2O, g

    IMPLICIT NONE
    INTRINSIC :: MAX, MIN, EPSILON, INT

    ! I/O
    REAL(DP), INTENT(IN)  :: ppf
    REAL(DP), INTENT(IN)  :: pdp
    REAL(DP), INTENT(IN)  :: ptf
    REAL(DP), INTENT(IN)  :: pxi
    REAL(DP), INTENT(IN)  :: pxl
    REAL(DP), INTENT(IN)  :: pclfr 
    REAL(DP), INTENT(IN)  :: pcdnc
    LOGICAL,  INTENT(IN)  :: laland
    LOGICAL,  INTENT(IN)  :: laglac
    !
    REAL(DP), INTENT(OUT) :: radip
    REAL(DP), INTENT(OUT) :: radlp

    ! LOCAL
    REAL(DP) :: zclic, zclfr, icldlyr, ziwgkg, ziwc
    REAL(DP) :: zkap, zrex, zref, zcdnc, ziwp, zlwgkg, zclwa, zlwc, zlwp


    ! cloud liquid water
    zclwa = MAX(pxl, 0._dp)

    ! cloud ice
    zclic = MAX(pxi,0._dp)

    zclfr = pclfr

    ! Clear/cloudy index
    IF (zclfr > EPSILON(1.0_dp)) THEN
       icldlyr = 1._dp
    ELSE
       icldlyr = 0._dp
    END IF

    ! Secure cloud fraction
    zclfr = MAX(pclfr, EPSILON(1.0_dp))

    ! Specific ice water content, g/kg
    IF (INT(icldlyr) == 1) THEN
       ziwgkg = zclic * 1000._dp/zclfr
    ELSE
       ziwgkg = 0._dp
    END IF

    ! Ice water content per volume g/m3
    ziwc = ziwgkg*ppf/ptf/rd

    ! Ice water path, g/m2
    ziwp = ziwgkg * pdp/g

    ! Specific liquid water content, g/kg
    IF (INT(icldlyr) == 1) THEN
       zlwgkg = zclwa*1000._dp/zclfr
    ELSE
       zlwgkg = 0._dp
    END IF

    ! Liquid water content per volume, g/m3
    zlwc = zlwgkg*ppf/ptf/rd

    ! Liquid water path, g/m2
    zlwp = zlwgkg*pdp/g

    ! Effective radii for ice and liquid particles, micrometer
    ! Boucher/Lohmann (1995) and Moss et al. (1995)
    !
    ! - ice
    radip = MAX(ceffmin, MIN(ceffmax,83.8_dp*ziwc**0.216_dp))

    ! - liquid
    IF (laland .AND. (.NOT.laglac)) THEN
       zkap = 1.143_dp
    ELSE
       zkap = 1.077_dp
    END IF
    zrex = 1._dp/3._dp
    zref = 1.e6_dp * (3.e-9_dp / (4._dp * pi * rho_H2O))**zrex
    zcdnc = pcdnc * 1.e-6_dp

    if (zcdnc > 0._dp) then
       radlp = zref * zkap * (zlwc/zcdnc)**zrex
    else
       radlp = 0._dp
    endif
    radlp = MAX(4._dp, MIN(24._dp,radlp))

  END SUBROUTINE rad_std
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE cloud_set_param(status,nn,nlev,lmidatm,lcouple)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: status                !error status
    INTEGER, INTENT(IN)  :: nn,nlev
    LOGICAL, INTENT(IN)  :: lmidatm,lcouple     

    !LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloud_set_param'

    status=1

    !
    !                11 Level, no middle atmosphere
    !
    IF(nlev == 11 .AND. .NOT. lmidatm) THEN
       !    IF(nn == 21) THEN
       IF((nn == 21).OR.(nn == 10)) THEN
          zinpar=0.1_dp
          zinhomi=0.7_dp
       ELSE IF(nn == 31) THEN
          zinpar=0.1_dp
          zinhomi=0.7_dp
       ELSE
          RETURN
       ENDIF
       !
       !                19 Level, no middle atmosphere
       !
    ELSE IF(nlev == 19 .AND. .NOT. lmidatm) THEN
       !    IF(nn == 21) THEN
       IF((nn == 21).OR.(nn == 10)) THEN
          zinpar=0.05_dp
          zinhomi=0.8_dp
       ELSE IF(nn == 31) THEN
          zinpar=0.06_dp
          zinhomi=0.8_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.70_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE
          RETURN
       ENDIF
       !
       !                31 Level, no middle atmosphere
       !
    ELSE IF(nlev == 31  .AND. .NOT. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.1_dp
          zinhomi=0.85_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.80_dp  ! mz_ap_20090519
       ELSE IF(nn == 85) THEN
          zinpar=0.1_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=1.0_dp
       ELSE
          RETURN
       ENDIF
    ELSE IF(nlev == 41  .AND. .NOT. lmidatm) THEN
       IF (nn == 10) THEN
          zinpar = 0.05_dp
          zinhomi = 0.80_dp
       ELSE IF(nn == 21) THEN
          zinpar = 0.05_dp
          zinhomi = 0.80_dp
       ELSE IF(nn == 31) THEN
          zinpar = 0.08_dp
          zinhomi = 0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar = 0.08_dp
          zinhomi = 0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar = 0.10_dp
          zinhomi = 0.85_dp
       ELSE IF(nn == 85) THEN
          zinpar = 0.10_dp
          zinhomi = 0.90_dp
       ELSE IF(nn == 106) THEN
          zinpar = 0.12_dp
          zinhomi = 0.90_dp
       ELSE IF(nn == 159) THEN
          zinpar = 0.12_dp
          zinhomi = 0.90_dp
       ELSE IF(nn == 255) THEN
          zinpar = 0.12_dp
          zinhomi = 0.90_dp
       ELSE
          RETURN
       ENDIF
       !
       !                60 Level
       !
    ELSE IF(nlev == 60) THEN
       IF (nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE
          RETURN
       ENDIF
       !
       !                39 Level, middle atmosphere
       !
    ELSE IF(nlev == 39 .AND. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.06_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE
          RETURN
       ENDIF
       !
       !                47 Level, middle atmosphere
       !
    ELSE IF(nlev == 47  .AND. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.80_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.1_dp
          zinhomi=0.85_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.80_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.1_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=1.0_dp
       ELSE
          RETURN
       ENDIF
       !
       !                49 Level, middle atmosphere
       !
    ELSE IF(nlev == 49  .AND. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.1_dp
          zinhomi=0.85_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.80_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.1_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=1.0_dp
       ELSE
          RETURN
       ENDIF
       !
       !                90 Level, middle atmosphere
       !
    ELSE IF(nlev == 90  .AND. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.1_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.1_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE
          RETURN
       ENDIF
       !
       !                95 Level, middle atmosphere
       !
    ELSE IF(nlev == 95  .AND. lmidatm) THEN
       IF(nn == 31) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 42) THEN
          zinpar=0.08_dp
          zinhomi=0.85_dp
       ELSE IF(nn == 63) THEN
          zinpar=0.1_dp
          zinhomi=0.85_dp
          ! coupled model:
          IF(lcouple) zinhomi=0.80_dp
       ELSE IF(nn == 85) THEN
          zinpar=0.1_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=1.0_dp
       ELSE
          RETURN
       ENDIF
       !
       !                191 Level, middle atmosphere
       !
    ELSE IF(nlev == 191  .AND. lmidatm) THEN
       IF (nn == 106) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 159) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 255) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE IF(nn == 319) THEN
          zinpar=0.12_dp
          zinhomi=0.9_dp
       ELSE
          RETURN
       ENDIF
       !
       !       74 Level, middle atmosphere and extended above
       !
    ELSE IF(nlev == 74 .AND. lmidatm) THEN
       IF(nn == 42) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE
          RETURN
       ENDIF
       !       84 Level, middle atmosphere and extended above
       !
    ELSE IF(nlev == 84 .AND. lmidatm) THEN
       IF(nn == 42) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE
          RETURN
       ENDIF
       !       94 Level, middle atmosphere and extended above
       !
    ELSE IF(nlev == 94 .AND. lmidatm) THEN
       IF(nn == 42) THEN
          zinpar=0.07_dp
          zinhomi=0.85_dp
       ELSE
          RETURN
       ENDIF
       ! CESM1:
    ELSE IF(nlev == 26) THEN
       zinpar=0.1_dp
       zinhomi=0.85_dp
    ELSE IF(nlev == 30) THEN
       zinpar=0.1_dp
       zinhomi=0.85_dp
    ELSE IF(nlev == 51 .AND. lmidatm) THEN
       zinpar=0.1_dp
       zinhomi=0.85_dp
    ELSE
       RETURN
    ENDIF

    !
    ! Overwrite the initialization of zinhomi and zinpar with the
    ! the value provided by the namelist.
    IF (rset_zinhomi%l) zinhomi = rset_zinhomi%v
    IF (rset_zinpar%l)  zinpar  = rset_zinpar%v

    !
    ! Initialisation of asic has moved from SUBROUTINE rad_rad_smcl to here,
    ! to enable the change of asic via namelist.
    !
    IF(lcouple) THEN
       asic = 0.89_dp
    ELSE
       asic = 0.85_dp
    ENDIF
    ! Overwrite the initialization of asic with the
    ! the value provided by the namelist.
    IF (rset_asic%l) asic = rset_asic%v

    status=0 !NO ERROR

  END SUBROUTINE cloud_set_param
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE cloudopt_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ rset_asic, rset_zinhomi, rset_zinhoml, rset_zinpar, loice

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='cloudopt_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    !IF (p_parallel_io) THEN
    !   WRITE (*,*) 'SUBROUTINE cloud_set_param:'
    !   WRITE (*,*) ' cloud inhomogeneity factor (ice): zinhomi = ',zinhomi
    !   WRITE (*,*) '                                   zinpar  = ',zinpar
    !   WRITE (*,*) ' parameter to correct the asymmetry factor of ice clouds:'
    !   WRITE (*,*) '                                   asic    = ',asic
    !END IF
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    !
  END SUBROUTINE cloudopt_read_nml_ctrl
  ! =========================================================================

  ! **********************************************************************
END MODULE messy_cloudopt
! **********************************************************************

