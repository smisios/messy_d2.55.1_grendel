       MODULE mo_sedmnt

!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_sedmnt.f90,v $\\
!$Revision: 1.2.2.1.20.1 $\\
!$Date: 2006/04/03 11:27:49 $\\

!***********************************************************************
!
!**** *MODULE mo_sedmnt* - Variables for sediment modules.
!
!     S.Legutke,        *MPI-MaD, HH*    31.10.01
!
!     Modified
!     --------
!     js add z_sed (depth of sediment layers)
!     
!     Purpose
!     -------
!     - declaration and memory allocation
!
!     *sedlay*         *REAL*  - .
!     *sedla1*         *REAL*  - .
!     *sedtot*         *REAL*  - .
!     *sedtoa*         *REAL*  - .
!     *seffel*         *REAL*  - .
!     *sedhpl*         *REAL*  - .   H+
!     *powtra*         *REAL*  - .
!     *prorca*         *REAL*  - .
!     *prcaca*         *REAL*  - .
!     *silpro*         *REAL*  - .
!     *porwat*         *REAL*  - .
!     *porsol*         *REAL*  - .
!     *seddzi*         *REAL*  - .
!     *dzs*            *REAL*  - .
!     *porwah*         *REAL*  - .
!     *seddw*          *REAL*  - .
!     *sedict*         *REAL*  - .
!     *rno3*           *REAL*  - .
!     *calcon*         *REAL*  - .
!     *ansed*          *REAL*  - .
!     *o2ut*           *REAL*  - .
!
!**********************************************************************
      implicit none

      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: sedlay
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: powtra

      REAL, DIMENSION (:,:,:), ALLOCATABLE :: sedhpl

      REAL, DIMENSION (:), ALLOCATABLE :: seddw
      REAL, DIMENSION (:), ALLOCATABLE :: porsol
      REAL, DIMENSION (:), ALLOCATABLE :: porwah
      REAL, DIMENSION (:), ALLOCATABLE :: porwat

      REAL, DIMENSION (:), ALLOCATABLE :: dzs
      REAL, DIMENSION (:), ALLOCATABLE :: seddzi
      REAL, DIMENSION (:), ALLOCATABLE :: z_sed

      REAL, DIMENSION (:,:), ALLOCATABLE :: silpro
      REAL, DIMENSION (:,:), ALLOCATABLE :: prorca
      REAL, DIMENSION (:,:), ALLOCATABLE :: prcaca
      REAL, DIMENSION (:,:), ALLOCATABLE :: pror13
      REAL, DIMENSION (:,:), ALLOCATABLE :: prca13
      REAL, DIMENSION (:,:), ALLOCATABLE :: pror14
      REAL, DIMENSION (:,:), ALLOCATABLE :: prca14
      REAL, DIMENSION (:,:), ALLOCATABLE :: produs
      
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: burial


      REAL :: sedict,calcon,rno3,o2ut,ansed,sedac,sedifti
      REAL :: calcwei, opalwei, orgwei
      REAL :: calcdens, opaldens, orgdens, claydens
      REAL :: calfa, oplfa, orgfa, clafa, solfu

      CONTAINS

      SUBROUTINE ALLOC_MEM_SEDMNT(kpie,kpje)

      use mo_control_bgc
      use mo_param1_bgc 

      INTEGER :: kpie,kpje
      


        WRITE(io_stdo_bgc,*)' '
        WRITE(io_stdo_bgc,*)'******************************************'
        WRITE(io_stdo_bgc,*)' '
        WRITE(io_stdo_bgc,*)'Memory allocation for sediment modules :'

        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedlay ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nsedtra   !js should be nsedtra (number of sediment tracers)

        ALLOCATE (sedlay(kpie,kpje,ks,nsedtra))


        WRITE(io_stdo_bgc,*)'Memory allocation for variable sedhpl ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks

        ALLOCATE (sedhpl(kpie,kpje,ks))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable burial ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nsedtra

        ALLOCATE (burial(kpie,kpje,nsedtra))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable powtra ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',npowtra

        ALLOCATE (powtra(kpie,kpje,ks,npowtra))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable silpro ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (silpro(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable prorca ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (prorca(kpie,kpje))
        ALLOCATE (pror13(kpie,kpje))
        ALLOCATE (pror14(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable prcaca ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (prcaca(kpie,kpje))
        ALLOCATE (prca13(kpie,kpje))
        ALLOCATE (prca14(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable produs ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje

        ALLOCATE (produs(kpie,kpje))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable dzs ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ksp

        ALLOCATE (dzs(ksp))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable seddzi ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ksp

        ALLOCATE (seddzi(ksp))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable z_sed ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks

        ALLOCATE (z_sed(ks))


        WRITE(io_stdo_bgc,*)'Memory allocation for variable seddw ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks

        ALLOCATE (seddw(ks))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable porsol ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks

        ALLOCATE (porsol(ks))


        WRITE(io_stdo_bgc,*)'Memory allocation for variable porwah ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks

        ALLOCATE (porwah(ks))

        WRITE(io_stdo_bgc,*)'Memory allocation for variable porwat ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',ks

        ALLOCATE (porwat(ks))


      END SUBROUTINE ALLOC_MEM_SEDMNT

    END MODULE mo_sedmnt
