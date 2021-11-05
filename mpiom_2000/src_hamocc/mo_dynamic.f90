      MODULE mo_dynamic

!***********************************************************************
!
!**** *MODULE mo_dynamic* - Variables for dynamic MLD sampling.
!
!     Patrick Wetzel    *MPI-Met, HH*    28.05.04
!
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      USE mo_kind, ONLY: wp
      implicit none

        INTEGER, PARAMETER ::                                  &
     &          kdphyto  =1,                                   &
     &          kdgrazer =2,                                   &
     &          kddoc    =3,                                   &
     &          kddic    =4,                                   &
     &          kdphosph =5,                                   &
     &          kdoxygen =6,                                   &
     &          kdiron   =7,                                   &
     &          kdano3   =8,                                   &
     &          kdalkali =9,                                   &
     &          kdsilica =10,                                  &
     &          kdtemp   =11,                                  &
     &          kdsal    =12,                                  &
#ifdef PANTHROPOCO2
     &          kddic_ant=13,                                  &
     &          kdalk_ant=14,                                  &
     &          nbgcdyn  =14
#else  /*  PANTHROPOCO2 */
     &          nbgcdyn  =12
#endif /*  PANTHROPOCO2 */


        INTEGER, PARAMETER ::                                &
     &          kdadv  =1,                                   &
     &          kddif  =2,                                   &
     &          kdpre  =3,                                   &
     &          kdgmp  =4,                                   &
     &          kdbio  =5,                                   &
     &          kdtot  =5

        INTEGER, PARAMETER ::                                &
     &          kdyndiff  =1,                                &
     &          kdynsave  =2

      INTEGER :: nc_dyn_id

      INTEGER, DIMENSION (:,:),   ALLOCATABLE :: nbgc_mld

      REAL(wp), DIMENSION (:,:),   ALLOCATABLE ::     bgc_zmld
      REAL(wp), DIMENSION (:,:),   ALLOCATABLE ::     bgc_nmld

      REAL(wp), DIMENSION (:,:,:,:),   ALLOCATABLE :: bgcdyn
      REAL(wp), DIMENSION (:,:,:),   ALLOCATABLE :: bgcdyntmp

      CONTAINS

      SUBROUTINE ALLOC_MEM_DYNAMIC(kpie,kpje,kpke)

      use mo_control_bgc
      use mo_bgcmean
      use mo_param1_bgc

      INTEGER :: kpie,kpje,kpke

        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcdyn ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcdyn
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',kdtot

        ALLOCATE (bgcdyn(kpie,kpje,nbgcdyn,kdtot))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcdyntmp ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcdyn

        ALLOCATE (bgcdyntmp(kpie,kpje,nbgcdyn))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_zmld ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (bgc_zmld(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_nmld ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (bgc_nmld(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable nbgc_mld ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (nbgc_mld(kpie,kpje))

      END SUBROUTINE ALLOC_MEM_DYNAMIC
      END MODULE mo_dynamic
