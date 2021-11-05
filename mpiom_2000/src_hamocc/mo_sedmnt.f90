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
      USE mo_kind, ONLY: wp
      implicit none

      REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE,TARGET :: sedlay
      REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE,TARGET :: powtra

      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE,TARGET :: sedhpl

      REAL(wp), DIMENSION (:), ALLOCATABLE :: seddw
      REAL(wp), DIMENSION (:), ALLOCATABLE :: porsol
      REAL(wp), DIMENSION (:), ALLOCATABLE :: porwah
      REAL(wp), DIMENSION (:), ALLOCATABLE :: porwat

      REAL(wp), DIMENSION (:), ALLOCATABLE :: dzs
      REAL(wp), DIMENSION (:), ALLOCATABLE :: seddzi
      REAL(wp), DIMENSION (:), ALLOCATABLE :: z_sed

      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: silpro
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: prorca
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: prcaca
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: pror13
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: prca13
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: pror14
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: prca14
      REAL(wp), DIMENSION (:,:), ALLOCATABLE :: produs

      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: burial

! pown2bud closes the mass balance for the alkalinity for biogenic induced changes in N2 in sediment
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: pown2bud
! powh2obud closes the mass balance for oxygen
      REAL(wp), DIMENSION (:,:,:), ALLOCATABLE :: powh2obud


      REAL(wp) :: sedict,calcon,rno3,o2ut,ansed,sedac,sedifti
      REAL(wp) :: calcwei, opalwei, orgwei
      REAL(wp) :: calcdens, opaldens, orgdens, claydens
      REAL(wp) :: calfa, oplfa, orgfa, clafa, solfu

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

        WRITE(io_stdo_bgc,*)'Memory allocation for variable pown2bud ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks

        ALLOCATE (pown2bud(kpie,kpje,ks))
        pown2bud = 0._wp

        WRITE(io_stdo_bgc,*)'Memory allocation for variable powh2obud ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks

        ALLOCATE (powh2obud(kpie,kpje,ks))
        powh2obud(:,:,:) = 0._wp

      END SUBROUTINE ALLOC_MEM_SEDMNT
    END MODULE mo_sedmnt
