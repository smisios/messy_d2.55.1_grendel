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
      implicit none

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: kbo
      REAL, DIMENSION (:,:), ALLOCATABLE :: bolay

!     REAL, DIMENSION (:,:,:), ALLOCATABLE :: sutrao
      REAL, DIMENSION (:,:), ALLOCATABLE :: expoor
      REAL, DIMENSION (:,:), ALLOCATABLE :: expoca
      REAL, DIMENSION (:,:), ALLOCATABLE :: exposi

      REAL, DIMENSION (:,:), ALLOCATABLE :: strahl

      REAL, DIMENSION (:,:), ALLOCATABLE :: alar1max
      REAL, DIMENSION (:,:), ALLOCATABLE :: TSFmax
      REAL, DIMENSION (:,:), ALLOCATABLE :: TMFmax

#ifdef FB_BGC_OCE     
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: abs_oce
#endif 

      REAL :: phytomi,grami,grazra,rrrcl
      REAL :: remido,dyphy,zinges,epsher,spemor,gammap,gammaz,ecan
      REAL :: ro2ut,rcar,rnit,rnoi,rnit23,rnit13,rcalc,ropal,bluefix
      REAL :: bkphy,bkzoo,bkopal,bifr13,bifr14,plafr13,plafr14
      REAL :: wpoc,wcal,wopal,drempoc,dremdoc,dremcalc,dremn2o
      REAL :: dphymor,dzoomor,dremopal,calmax, gutc
      REAL :: psedi,csedi,ssedi
      REAL :: perc_diron, riron, fesoly, relaxfe, wdust,bolaymin  
      REAL :: ctochl, atten_w, atten_c, atten_f
      REAL :: pi_alpha


#ifdef AGG      
      REAL :: SinkExp, FractDim, Stick, cellmass, cellsink, fsh, fse
      REAL :: alow1, alow2,alow3,alar1,alar2,alar3,TSFac,TMFac
      REAL :: vsmall,safe,pupper,plower,zdis
      REAL :: dustd1,dustd2,dustd3,dustsink
#endif

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

#ifdef FB_BGC_OCE 
         WRITE(io_stdo_bgc,*)'Memory allocation for variable abs_oce ...'
         WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
         WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
         WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke

         ALLOCATE (abs_oce(kpie,kpje,kpke))
#endif 

      END SUBROUTINE ALLOC_MEM_BIOMOD

      END MODULE mo_biomod
