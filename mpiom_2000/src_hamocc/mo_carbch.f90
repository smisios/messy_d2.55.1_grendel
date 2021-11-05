      MODULE mo_carbch
!***********************************************************************
!
!**** *MODULE mo_carbch* - Variables for inorganic carbon cycle.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!
!     Patrick Wetzel    *MPI-Met, HH*    16.04.02
!     - new: atm, atdifv, suppco2
!     - changed: chemc(:,:,:) to chemcm(:,:,:,:)
!     - new: bgcmean(:,:,:,:)
!
!     Purpose
!     -------
!     - declaration and memory allocation
!
!     *ocetra*       *REAL*  - .
!     *hi*           *REAL*  - .
!     *co3*          *REAL*  - .
!     *chemcm*       *REAL*  - .
!     *co2flu*       *REAL*  - .
!     *co2fluacc     *REAL*  - .
!     *akw3*         *REAL*  - .
!     *akb3*         *REAL*  - .
!     *ak13*         *REAL*  - .
!     *ak23*         *REAL*  - .
!     *aksp*         *REAL*  - .
!
!**********************************************************************

      USE mo_kind, ONLY: wp
      implicit none

      REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: ocetra
      REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE :: chemcm

      REAL(wp), DIMENSION (:,:,:),   ALLOCATABLE, TARGET :: atm

      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: co3
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: hi
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: akw3
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: akb3
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: ak13
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: ak23
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE,TARGET :: aksp
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: satoxy
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: satn2o
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: atdifv
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: suppco2
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: sedfluxi                ! never used
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: sedfluxo

      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: co2flux     ! sea-air C-flux, ingetrated over whole simulation period
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: o2flux      ! sea-air O2-flux, ingetrated over whole simulation period
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: n2flux      ! sea-air N2-flux, ingetrated over whole simulation period
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: n2oflux      ! sea-air N2-flux, ingetrated over whole simulation period

#ifdef __c_isotopes
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: c13flux     ! sea-air C13-flux, ingetrated over whole simulation period
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: o14flux     ! sea-air C14-flux, ingetrated over whole simulation period
#endif

#ifdef __cpl_co2
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: co2trans    ! transfer coefficient for CO2 atmosph/ocean
                                                          ! exchange, still to be multiplied by v**2
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: co2conc
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: co2flux_cpl
#endif
#ifdef __cpl_dust
      ! for input of dust fields (on-line from ECHAM)
      ! reintroduced js 19062006 */
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE :: dustdep                 ! + updated allocation
#endif
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: dusty
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: c14pool

      REAL(wp) :: dmspar(6)

! decay coefficient for sco214 plus more for 13C/14C
      REAL(wp) :: c14dec, D14Catm, Ratm, D14Cocn, Roc14
      REAL(wp) :: roc13, atcoa
      REAL(wp) :: ozkoa, c14prod, c14inva, eins, c14ret, c14prosta

      REAL(wp) :: atm_co2, atm_o2, atm_n2
!                              .35e-3 * 5.1e14*12. .35e-3 * 5.1e14
      REAL(wp) :: globalmean_co2, globalmean_o2, globalmean_n2
      REAL(wp) :: atm_c13, atm_c14, atmacon,           atmacmol
      REAL(wp) :: ppm2con,contppm
      REAL(wp) :: ems_per_step

      REAL(wp) :: totalcarbon_old

      REAL(wp) :: calcinp, orginp, silinp
      REAL(wp),TARGET :: calcinpglint, orginpglint, silinpglint

      REAL(wp) :: totarea

!hh#ifdef PCFC
      REAL(wp) ::           cfc11_atm_1s,cfc11_atm_2s,cfc11_atm_3s        &
     &                 ,cfc11_atm_1n,cfc11_atm_2n,cfc11_atm_3n        &
     &                 ,cfc12_atm_1s,cfc12_atm_2s,cfc12_atm_3s        &
     &                 ,cfc12_atm_1n,cfc12_atm_2n,cfc12_atm_3n

      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: cfc_int
!hh#endif
!hh#ifdef ANTC14
      REAL(wp) :: D14C_north, D14C_south, D14C_equator
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: Rbomb
!hh#endif

! n2budget closes the mass balance for the alkalinity for biogenic induced changes in N2
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: n2budget
! h2obudget closes the mass balance for the oxygen for biogenic induced changes in N2
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: h2obudget

      CONTAINS

      SUBROUTINE ALLOC_MEM_CARBCH(kpie,kpje,kpke)

      use mo_control_bgc
      use mo_param1_bgc

      INTEGER :: kpie,kpje,kpke

        WRITE(io_stdo_bgc,*)'Memory allocation for variable ocetra ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nocetra

        ALLOCATE (ocetra(kpie,kpje,kpke,nocetra))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable hi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (hi(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (co3(kpie,kpje,kpke))
        ALLOCATE (c14pool(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable chemcm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',8
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',12

        ALLOCATE (chemcm(kpie,kpje,8,12))

!        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxi ...'
!        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
!        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
!        WRITE(io_stdo_bgc,*)'Third dimension    : ',3
!
!        ALLOCATE (sedfluxi(kpie,kpje,3))


        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedfluxo ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',npowtra

        ALLOCATE (sedfluxo(kpie,kpje,npowtra))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable satn2o ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (satn2o(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable aksp ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (aksp(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable satoxy ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (satoxy(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable ak23 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (ak23(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable ak13 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (ak13(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable akb3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (akb3(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable akw3 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (akw3(kpie,kpje,kpke))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',natm

        ALLOCATE (atm(kpie,kpje,natm))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable atdifv ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (atdifv(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable suppco2 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (suppco2(kpie,kpje))
        suppco2(:,:) = 0.0_wp

! from esm version start
!       WRITE(io_stdo_bgc,*)'Memory allocation for variable dmsflux ...'
!       WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
!       WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

!       ALLOCATE (dmsflux(kpie,kpje))
!       dmsflux(:,:)=0.0

!       WRITE(io_stdo_bgc,*)'Memory allocation for variable dms ...'
!       WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
!       WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

!       ALLOCATE (dms(kpie,kpje))
!       dms(:,:)=0.0

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2flux(kpie,kpje))
        co2flux(:,:) = 0.0_wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable o2flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (o2flux(kpie,kpje))
        o2flux(:,:) = 0.0_wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable n2flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (n2flux(kpie,kpje))
        n2flux(:,:) = 0.0_wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable n2oflux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (n2oflux(kpie,kpje))
        n2oflux(:,:) = 0.0_wp

#ifdef __c_isotopes
        WRITE(io_stdo_bgc,*)'Memory allocation for variable c13flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (c13flux(kpie,kpje))
        c13flux(:,:) = 0.0_wp
        WRITE(io_stdo_bgc,*)'Memory allocation for variable c14flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (c14flux(kpie,kpje))
        c14flux(:,:) = 0.0_wp
#endif

#ifdef __cpl_co2

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2trans ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2trans(kpie,kpje))
        co2trans(:,:) = 0.0_wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2conc ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2conc(kpie,kpje))
        co2conc(:,:) = 0.0_wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2flux_cpl ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2flux_cpl(kpie,kpje))
        co2flux_cpl(:,:) = 0.0_wp
#endif
#ifdef __cpl_dust

        WRITE(io_stdo_bgc,*)'Memory allocation for variable dustdep ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (dustdep(kpie,kpje))
        dustdep(:,:) = 0.0_wp

! js changed from else to endif for supervolcano exp.
#endif

        WRITE(io_stdo_bgc,*)'Memory allocation for variable dusty ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',12

        ALLOCATE (dusty(kpie,kpje,12))
! end insert from esm version

#ifdef PCFC
        WRITE(io_stdo_bgc,*)'Memory allocation for variable cfc_int ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (cfc_int(kpie,kpje))
#endif
#ifdef ANTC14
        WRITE(io_stdo_bgc,*)'Memory allocation for variable Rbomb ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (Rbomb(kpie,kpje))
#endif

        WRITE(io_stdo_bgc,*)'Memory allocation for variable n2budget ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (n2budget(kpie,kpje,kpke))
        n2budget = 0._wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable h2obudget ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

        ALLOCATE (h2obudget(kpie,kpje,kpke))
        h2obudget(:,:,:) = 0._wp


      END SUBROUTINE ALLOC_MEM_CARBCH
      END MODULE mo_carbch
