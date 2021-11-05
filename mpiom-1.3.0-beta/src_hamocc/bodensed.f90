      SUBROUTINE BODENSED(kpie,kpje,kpke,pddpo)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/bodensed.f90,v $\\
!$Revision: 1.2.2.1.4.1.2.2.4.1.2.2.2.4.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!**********************************************************************
!
!**** *BODENSED* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     j. segschneider  04.08.2005  introduced z_sed (depth of sediment layers)
!                                  for output in bgcmean_sed
!
!     Purpose
!     -------
!     set up of sediment layers. + kbo + several constants as in beleg_bgc.f90
!
!     Method:
!     ------
!     -
!
!     *CALL*       *BODENSED
!     called by ini_bgc.f90
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *dtbgc*     - ocean model time step [sec].
!
!     Externals
!     ---------
!     none.
!**********************************************************************
! evaluates the min depth of all layers, to be checked vs. sinking rate
!   in routine BELEG_BGC, to assure that flux rate is < 1 per time step

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod

      USE mo_control_bgc
      use mo_param1_bgc 

      use mo_parallel

implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: sumsed

! define sediment layer thickness [m]
      dzs(1) = 0.001
      dzs(2) = 0.003
      dzs(3) = 0.005
      dzs(4) = 0.007
      dzs(5) = 0.009
      dzs(6) = 0.011
      dzs(7) = 0.013
      dzs(8) = 0.015
      dzs(9) = 0.017
      dzs(10) = 0.019
      dzs(11) = 0.021
      dzs(12) = 0.023
      dzs(13) = 0.025

      WRITE(io_stdo_bgc,*)  ' '
      WRITE(io_stdo_bgc,*)  'Sediment layer thickness [m] : '
      WRITE(io_stdo_bgc,'(5F9.3)') dzs
      WRITE(io_stdo_bgc,*)  ' '
    
      porwat(1) = 0.85
      porwat(2) = 0.83
      porwat(3) = 0.8
      porwat(4) = 0.79
      porwat(5) = 0.77
      porwat(6) = 0.75
      porwat(7) = 0.73
      porwat(8) = 0.7
      porwat(9) = 0.68
      porwat(10) = 0.66
      porwat(11) = 0.64
      porwat(12) = 0.62

      WRITE(io_stdo_bgc,*)  'Pore water in sediment: ',porwat

      
      seddzi(1)=500.
      DO 1131 k=1,ks
         porsol(k)=1.-porwat(k)
         IF(k.EQ.1) porwah(k)=0.5*(1.+porwat(1))
         IF(k.GE.2) porwah(k)=0.5*(porwat(k)+porwat(k-1))
         seddzi(k+1)=1./dzs(k+1)                           ! inverse depth of sediment layer
         seddw(k)=0.5*(dzs(k)+dzs(k+1))
1131  CONTINUE

! ******************************************************************
! the following section is to include the SEDIMENT ACCELERATION
! mechanism to accelerate the sediment:

      sedac = 1./float(isac)

! determine total solid sediment thickness sumsed
! and reduced sediment volume
      sumsed = 0.
      do k = 1, ks
        porwat(k) = porwat(k) * sedac
        porsol(k) = porsol(k) * sedac
        sumsed = sumsed + seddw(k)
      enddo

! depth of sediment layers for output in bgcmean
      z_sed(1) = dzs(1)
      do k = 2, ks
        z_sed(k) = z_sed(k-1) + dzs(k)
      enddo

! determine reduced diffusion rate sedict,
! and scaling for bottom layer ventilation, sedifti

      sedict=1.e-9*dtbgc                           ! units? m/s2 ?
      sedifti = sedict / (sumsed**2)               ! not used
      sedict=sedict*sedac

      WRITE(io_stdo_bgc,*)  'sediment acceleration factor: ',isac
      WRITE(io_stdo_bgc,*)  'sediment area reduction: ',sedac
      WRITE(io_stdo_bgc,*)  'new diffusion is: ',sedict
      WRITE(io_stdo_bgc,*)  'total sediment thickness [m]: ',sumsed
      WRITE(io_stdo_bgc,*)  'sedict / sumsed**2: ',sedifti
      WRITE(io_stdo_bgc,*)  'new pore water fraction: ',porwat


! ******************************************************************

! ******************************************************************
! densities etc. for SEDIMENT SHIFTING

! define (wei)ght of calcium carbonate, opal, and poc [kg/kmol]
      calcwei = 100. ! 40+12+3*16 kg/kmol C
      opalwei = 60.  ! 28 + 2*16  kg/kmol Si
      orgwei  = 30.  ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                     ! after Alldredge, 1998: POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

! define densities of caco3, opal, caco3 [kg/m3]
      calcdens = 2600. 
      opaldens = 2200. 
      orgdens  = 1000. 
      claydens = 2600. !quartz
      
! define volumes occupied by solid constituents [m3/kmol]
      calfa = calcwei/calcdens         ! 3.85e-2
      oplfa = opalwei/opaldens         ! 2.73e-2
      orgfa = orgwei/orgdens           ! 3.0e-2
      clafa = 1./claydens              ! 3.85e-4 (clay is calculated in kg/m3)

! determine total solid column integrated sediment volume  (1-d)
      solfu = 0.
      DO k=1,ks
         solfu = solfu + seddw(k)*porsol(k)
      ENDDO

! ******************************************************************

      k=kpke
      DO 1321 j=1,kpje
      DO 1321 i=1,kpie
         kbo(i,j)=1
         bolay(i,j)=0.
         IF(pddpo(i,j,k).GT.0.5) THEN
            bolay(i,j)=pddpo(i,j,k)                  ! thickness of bottom layer
            kbo(i,j)=k
         ENDIF
1321  CONTINUE

! evaluate min depth of last layer 
      bolaymin=8000.
      DO 1322 k=kpke-1,1,-1
      DO 1322 j=1,kpje
      DO 1322 i=1,kpie
         IF(pddpo(i,j,k).GT.0.5 .AND. pddpo(i,j,k+1).LE.0.5) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
            bolaymin = min(bolaymin,bolay(i,j))
         ENDIF
1322  CONTINUE
      CALL global_min(bolaymin)
      WRITE(io_stdo_bgc,*)  'bolaymin=', bolaymin

      RETURN
      END
