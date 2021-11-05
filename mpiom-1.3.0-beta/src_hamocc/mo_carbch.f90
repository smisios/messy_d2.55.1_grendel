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
      
      implicit none
      
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ocetra
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: chemcm
!hh#ifdef DIFFAT      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: atm      
!hh#endif
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: co3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: hi
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: akw3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: akb3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ak13
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ak23
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: aksp
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: satoxy
      REAL, DIMENSION (:,:),   ALLOCATABLE :: satn2o
      REAL, DIMENSION (:,:),   ALLOCATABLE :: atdifv
      REAL, DIMENSION (:,:),   ALLOCATABLE :: suppco2
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedfluxi                ! never used
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedfluxo
#ifdef __cpl_co2
      REAL, DIMENSION (:,:),   ALLOCATABLE :: co2flux
      REAL, DIMENSION (:,:),   ALLOCATABLE :: co2conc
#endif
#ifdef __cpl_dust      ! for input of dust fields (on-line from ECHAM)! reintroduced js 19062006
      REAL, DIMENSION (:,:),   ALLOCATABLE :: dustdep                 ! + updated allocation
#endif
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: dusty
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: c14pool
      
      REAL :: dmspar(6)
            
! decay coefficient for sco214 plus more for 13C/14C
      REAL :: c14dec, D14Catm, Ratm, D14Cocn, Roc14
      REAL :: roc13, atcoa
      REAL :: ozkoa, c14prod, c14inva, eins, c14ret, c14prosta

!hh#ifndef DIFFAT            
      REAL :: atm_co2, atm_o2, atm_n2 
!                              .35e-3 * 5.1e14*12. .35e-3 * 5.1e14
      REAL :: atm_c13, atm_c14, atmacon,           atmacmol
!hh#endif  

!hh#ifdef DIFFAT
      REAL :: emission, ems_per_step
!hh#else
      REAL :: co2_atm_1,co2_atm_2,co2_atm_3
!hh#endif      


!hh#ifdef PCFC
      REAL ::           cfc11_atm_1s,cfc11_atm_2s,cfc11_atm_3s        &
     &                 ,cfc11_atm_1n,cfc11_atm_2n,cfc11_atm_3n        &
     &                 ,cfc12_atm_1s,cfc12_atm_2s,cfc12_atm_3s        &
     &                 ,cfc12_atm_1n,cfc12_atm_2n,cfc12_atm_3n        

      REAL, DIMENSION (:,:), ALLOCATABLE :: cfc_int
!hh#endif
!hh#ifdef ANTC14
      REAL :: D14C_north, D14C_south, D14C_equator
      REAL, DIMENSION (:,:), ALLOCATABLE :: Rbomb
!hh#endif


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
!       WRITE(io_stdo_bgc,*)'Third dimension    : ',3
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
#ifdef DIFFAT      
      
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atm ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',natm

        ALLOCATE (atm(kpie,kpje,natm))
#endif      
        WRITE(io_stdo_bgc,*)'Memory allocation for variable atdifv ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (atdifv(kpie,kpje))
      
      WRITE(io_stdo_bgc,*)'Memory allocation for variable suppco2 ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (suppco2(kpie,kpje))

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

#ifdef __cpl_co2

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2flux ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2flux(kpie,kpje))
        co2flux(:,:)=0.0

        WRITE(io_stdo_bgc,*)'Memory allocation for variable co2conc...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (co2conc(kpie,kpje))
        co2conc(:,:)=0.0   
#endif
#ifdef __cpl_dust

        WRITE(io_stdo_bgc,*)'Memory allocation for variable dustdep ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (dustdep(kpie,kpje))
        dustdep(:,:)=0.0

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

      END SUBROUTINE ALLOC_MEM_CARBCH

      END MODULE mo_carbch
