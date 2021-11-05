      MODULE mo_biomod
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_biomod.f90,v $\\
!$Revision: 1.2.22.1 $\\
!$Date: 2006/04/03 11:27:49 $\\

!***********************************************************************
!
!**** *MODULE mo_biomod* - Variables for marine biology.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     note: kbo,bolay shall be moved to mo_control_bgc
!
!     Purpose
!     -------
!     - declaration and memory allocation.
!
!     *kbo*         *INTEGER*  - k-index of bottom layer (2d)
!     *bolay*          *REAL*  - height of bottom cell.
!
!**********************************************************************
      USE mo_kind, ONLY: wp
      implicit none

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: kbo
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: bolay

!     REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: sutrao
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: expoor
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: expoca
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: exposi

      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: strahl

      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: alar1max
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: TSFmax
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: TMFmax

      REAL(wp) :: phytomi,grami,grazra,rrrcl
      REAL(wp) :: remido,dyphy,zinges,epsher,spemor,gammap,gammaz,ecan
      REAL(wp) :: ro2ut,rcar,rnit,rnoi,rnit23,rnit13,rcalc,ropal,bluefix
      REAL(wp) :: bkphy,bkzoo,bkopal,bifr13,bifr14,plafr13,plafr14
      REAL(wp) :: wpoc,wcal,wopal,drempoc,dremdoc,dremcalc,dremn2o,dremsul
      REAL(wp) :: dphymor,dzoomor,dremopal,calmax, gutc
      REAL(wp) :: psedi,csedi,ssedi
      REAL(wp) :: perc_diron, riron, fesoly, relaxfe, wdust,bolaymin
      REAL(wp) :: pi_alpha
      REAL(wp) :: nitdem,n2prod,ro2nitri

#ifdef AGG
      REAL(wp) :: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
      REAL(wp) :: alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
      REAL(wp) :: vsmall,safe,pupper,plower,zdis
#endif
      REAL(wp) :: dustd1,dustd2,dustd3,dustsink

      CONTAINS

      SUBROUTINE ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

      use mo_control_bgc
      use mo_param1_bgc

      INTEGER :: kpie,kpje,kpke

         WRITE(io_stdo_bgc,*)'Memory allocation for variable expoor ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (expoor(kpie,kpje))
!        ALLOCATE (sutrao(kpie,kpje,nocetra))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable expoca ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (expoca(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable exposi ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (exposi(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable kbo ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (kbo(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable strahl ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (strahl(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable bolay ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (bolay(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable alar1max'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (alar1max(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable TSFmax'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (TSFmax(kpie,kpje))

         WRITE(io_stdo_bgc,*)'Memory allocation for variable TMFmax'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

         ALLOCATE (TMFmax(kpie,kpje))

      END SUBROUTINE ALLOC_MEM_BIOMOD
      END MODULE mo_biomod
