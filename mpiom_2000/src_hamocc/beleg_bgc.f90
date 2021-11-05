      SUBROUTINE BELEG_BGC(kpie,kpje,kpke,pddpo)
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/beleg_bgc.f90,v $\\
!$Revision: 1.2.2.1.16.1.2.1.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\

!****************************************************************
!
!**** *BELEG_BGC* - initialize bgc variables.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     - set start values for  bgc variables.
!
!     Method
!     -------
!     - bgc variables are initialized. They might be overwritten
!       if read from restart by call of AUFR_BGC.
!     - physical constants are defined
!     - fields used for budget calculations (should be in extra SBR!)
!       and as boundary conditions are set to zero.
!
!
!**   Interface.
!     ----------
!
!     *CALL*       *BELEG_BGC*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1, ONLY: je1
      use mo_param1_bgc
      use mo_diagnosis , only : ibek
      USE mo_param1, ONLY: je1
      use mo_commo1
      USE mo_planetary_constants, ONLY: g,rhoref_water
      use mo_parallel
      USE mo_boundsexch, ONLY: bounds_exch

#ifdef PDYNAMIC_BGC
      USE mo_dynamic
#endif /* PDYNAMIC_BGC */


      IMPLICIT NONE

      INTEGER :: i,j,k,l,m,kpie,kpje,kpke,kmon

      INTEGER :: jj
#if defined(PCFC) || defined(ANTC14)
      INTEGER :: ii
#endif
      REAL(wp) :: north
      REAL(wp) :: south

      REAL(wp) :: pddpo(kpie,kpje,kpke)

      REAL(wp), ALLOCATABLE  ::ocean_area(:,:)   ! help field to sum up global ocean area ; result stored in totarea


      REAL(wp) :: phosat,phosmed,phospac,oxyat,oxymed,oxypac,silat,silmed,silpac

#ifdef AGG
      REAL(wp) :: shear,zmini,talar1,snow, checksink
#endif
!
! Initialize overall time step counter.
!
      ldtbgc = 0
!
!  Initialize 2D and 3D time step counter.
!
      meancnt_bgc_2D = 0
      meancnt_bgc_3D = 0
!
      if(diffat) then
! all concentrations will be calculated in carchm
      else
         atm_co2 = 278._wp
#ifdef __c_isotopes
         atm_c13 = 276.2_wp ! preindustrial atmospheric d13c=-6.5 permil --> 276.2ppm? test js 15082006 ok.
                        ! piston velocity ~8m/yr -> equilibration time atmosphere/ocean ~500 years
         atm_c14 = 274._wp
#endif
         atcoa=atm_co2                                    ! emr for use in maschk

         atm_o2  = 196800._wp
         atm_n2  = 802000._wp
      endif
         ozkoa=3.e15_wp     ! ozean kohlenstoff alt, start with dummy value in maschk
         atmacon =0.35e-3_wp * 5.1e14_wp*12._wp                    ! factor for atmospheric concentration of CO2 in g C -> ppm
         atmacmol=0.35e-3_wp * 5.1e14_wp                        ! factor for atmospheric concentration of CO2 in mol C -> ppm
                                                          ! 5.1e14 = total surface of Earth in [m^2]
! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
         ppm2con=0.35e-3_wp
         contppm=1._wp/0.35e-3_wp                               ! inverse of conversion factor molC/m**2 to ppm (in AT)

!        write(io_stdo_bgc,*)'atmacon=',atmacon,atmacmol  ! 2.14e+12, 1.78e+11


!
! Biology
!
!ik note that the unit is kmol/m3!
      phytomi=1.e-11_wp            !i.e. 1e-5 mmol P/m3 minimum concentration of phytoplankton
      grami=1.e-11_wp              !i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton | test e-11 ->e-8 very slow decay

!ik addded parameter definition; taken from OCPROD.F (dtb= 1./ steps/day)
!      remido=0.004*dtb      !1/d -remineralization rate (of DOM)
            remido = 0.008_wp * dtb      !KS, JS, EMR 12/07/2007
      dyphy = 0.008_wp * dtb       !1/d -mortality rate of phytoplankton
      zinges = 0.6_wp            !dimensionless fraction -assimilation efficiency of zooplankton
      epsher = 0.8_wp            !dimensionless fraction - (1-epsher)=fraction of grazing egested
      grazra = 1.0_wp * dtb        !1/d -grazing rate [emr: 0.6-0.9]
      spemor = 3.e6_wp * dtb       !1/d -mortality rate of zooplankton
      gammap = 0.03_wp * dtb       !1/d -exudation rate
      gammaz = 0.06_wp * dtb       !1/d -excretion rate
      ecan = 0.95_wp             ! fraction of mortality as PO_4
      pi_alpha = 0.02_wp * 0.4_wp     ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
!new emr
#ifdef __c_isotopes
! fractionation for photosynthesis plafr13 and plafr14 valid for particles
! for avoiding too many tracers, surface gross rates work with reduced
! values bifr13 and bifr14
       plafr13 = 1._wp                 ! planktonic fractionatio 13C   (never used) (HAMOCC3.1: 0.98)
       plafr14 = 1._wp
       bifr13 = 0.98_wp                ! biogenic fractionation ?
       bifr14 = 0.98_wp
#endif

! half sat. constants, note that the units are kmol/m3 ! (conc0 in hamocc3.1)
      bkphy  = 1.e-8_wp !i.e. 0.04 mmol P/m3 |js: 0.01 vs. 0.04? check 0.4 16.9.
      bkzoo  = 4.e-8_wp !i.e. 0.04 mmol P/m3
      bkopal = 1.e-6_wp !i.e. 1.0  mmol Si/m3
!js: no bkiron? (fitzwater et al 1996) et al. 200 120 nMol/m3, moore et al 80nMol/m3)

!sinking speed
      wpoc  = 5._wp * dtb  !m/d  iris : 5. (for dtb=.1 (10 time steps / day) -> wpoc=1) pw:10
      wcal  = 30._wp * dtb !m/d
      wopal = 30._wp * dtb !m/d  iris : 60, pw: 30, emr: 3 (20070506), emr: 10 (20070701)

! water column remineralisation constants

      drempoc  = 0.025_wp * dtb !1/d      ! 0.75/month. H3.1: 0.2/month k=1, 0.1/month k>1
      dremdoc  = 0.004_wp * dtb !1/d
      dphymor  = 0.1_wp * dtb   !1/d
      dzoomor  = 0.02_wp * dtb  !1/d
      dremopal = 0.0075_wp * dtb !1/d      ! 0.01 -> 0.001 js10072006 : slightly overdone -->0.0075
      dremopal = 0.01_wp * dtb

      dremn2o  = 0.01_wp*dtb       !1/d
      dremsul  = 0.005_wp*dtb
      dremcalc = 0.075_wp *dtb              ! 0.2 -> 0.02 js10072006 : slightly overdone --> 0.075

#ifdef AGG
      drempoc  = 0.1_wp *dtb      !1/d       re-introduced 09062006
      dremopal = 3.3e-3_wp *dtb   ! js 4.7.2006 0.0033 from .01/3. (60/20 m/d)
      dphymor  = 0.2_wp *dtb      !1/d
#endif

! nitrogen fixation by blue green algae (cyano.f90)

      bluefix = 0.005_wp * dtb     !1/d

! extended redfield ratio declaration
! Note: stoichiometric ratios are based on Takahashi etal. (1985)
! P:N:C:-O2 + 1:16:122:172

      ro2ut = 172._wp
      rcar = 122._wp
      rnit = 16._wp
      rnoi = 1._wp/rnit
! N consumption of denitrification corrected after Paulmier etal, 2009)
      nitdem = 121.6_wp      ! nitrate demand to remin. 1 mol P in suboxic water
      n2prod = 68.8_wp       ! N2 production for 1 mol P remineralized in suboxic water

      rcalc = 35._wp ! iris 40 !calcium carbonate to organic phosphorous production ratio
#ifdef __cpl_co2
      rcalc = 20._wp ! emr !calcium carbonate to organic phosphorous production ratio
#endif
      ropal = 25._wp ! iris 25 !opal to organic phosphorous production ratio
      calmax=0.20_wp ! max fraction of CaCO3 production (AGG) (max rain ratio)
      gutc = 0._wp   ! fraction of caco3 that dissolves during passage through zooplankton guts (not used yet)

!ik for interaction with sediment module
      o2ut=172._wp
      rno3=16._wp

!ik weight percent iron in dust deposition (0.035) times Fe solubility (0.01) /55.85 g--> mol
      perc_diron = 0.035_wp * 0.01_wp / 55.85_wp
! the three values below are from Johnson et al., 1997 Mar. Chemistry 57, 137-161
!     riron   = 5.*rcar*1.e-6       ! 'Q'=5.*122.*1.e-6 = 6.1e-4   (iron to phosphate stoichiometric ratio * umol->mol)
!     riron   = 2.*rcar*1.e-6       ! 14.2.06: 5.->2. (2.5umol P/nmol Fe) emr: 3.*rcar*1.e-6
      riron   = 3._wp * rcar*1.e-6_wp       ! 06.3.06: 2.->3. coex90=12.2GTC/a before change
      fesoly  = 0.6_wp *1.e-9_wp          ! global mean/max (?) dissolved iron concentration in deep water (for relaxation) [mol/m3]
                                    ! 'fesoly' stands for 'apparent solubility value' of iron
!      relaxfe = 0.05/365.*dtb       ! relaxation time for iron to fesoly ->1.37e-3*dtb(iris' paper states 0.005!?)back15206js
      relaxfe = 0.005_wp/365._wp*dtb       ! changed as suggested by emr
                                    ! (should probably be named scavenge_fe or remove_fe)
                                    ! (relaxfe is only used to remove iron, corresponding to Johnson's 200yr residence time)
                                    ! but 200 yrs might be too long, as iron concentrations in the Atlantic are too high.

! decay parameter for sco214, half life time = 5730years
      c14dec=(LOG(2._wp)/(5730._wp*365._wp))*dtb
      eins=1._wp
      c14ret=eins-c14dec

! Ratm: normalized atmospheric ratio of C-14/C-12, D14Catm: atmospheric Delta C-14
      D14Catm=0._wp             ! D14Catm=0. only for equilibrium runs
      Ratm = 1._wp + D14Catm/1000._wp

!
! Set constants for calculation of DMS ( mo_carbch )
! Parameters are a result from Kettle optimisation (fit to Kettle data, version emr, not tinka)

       dmspar(6)=0.100000000E-07_WP       !0 half saturation microbial
       dmspar(5)=1.25_wp*0.107638502E+00_WP  !2*1.3e-5 production with delsil (diatoms)
       dmspar(4)=1.25_wp*0.109784522E-01_WP  !2*0.02   production with delcar (coccolithoporides)
       dmspar(3)=0.0096_wp  !4.8e-5       !2*1.6e-3 microbial consumption
       dmspar(2)=0.0075_wp  !0.0003       !2*0.005  photolysis (uv-destruction)
       dmspar(1)=10._wp                   !2*5.     temperature dependent release by phytoplankton


      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'* Values of BELEG_BGC variables : '
      IF (.NOT. diffat) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_co2      = ',atm_co2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_o2       = ',atm_o2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_n2       = ',atm_n2
      ENDIF
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              phytomi      = ',phytomi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grami        = ',grami
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              remido       = ',remido
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dyphy        = ',dyphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              zinges       = ',zinges
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              epsher       = ',epsher
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grazra       = ',grazra
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              spemor       = ',spemor
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammap       = ',gammap
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammaz       = ',gammaz
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ecan         = ',ecan
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkphy        = ',bkphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkzoo        = ',bkzoo
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkopal       = ',bkopal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wpoc         = ',wpoc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wcal         = ',wcal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wopal        = ',wopal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              drempoc      = ',drempoc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremdoc      = ',dremdoc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremopal     = ',dremopal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dphymor      = ',dphymor
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dzoomor      = ',dzoomor
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bluefix      = ',bluefix
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ro2ut        = ',ro2ut
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcar         = ',rcar
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnit         = ',rnit
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnoi         = ',rnoi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcalc        = ',rcalc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ropal        = ',ropal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gutc         = ',gutc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              o2ut         = ',o2ut
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rno3         = ',rno3
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              perc_diron   = ',perc_diron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              riron        = ',riron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              fesoly       = ',fesoly
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              relaxfe      = ',relaxfe
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              c14dec       = ',c14dec
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              D14Catm      = ',D14Catm
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              Ratm         = ',Ratm
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(1)    = ',dmspar(1)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(2)    = ',dmspar(2)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(3)    = ',dmspar(3)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(4)    = ',dmspar(4)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(5)    = ',dmspar(5)

#ifndef AGG
      dustd1 = 0.0001_wp !cm = 1 um, boundary between clay and silt
      dustd2 = dustd1*dustd1
      dustsink = (g * 86400._wp / 18._wp                  &  ! g * sec per day / 18. | js: Stoke's law for small particles
     &         * (claydens - rhoref_water) / 1.567_wp * 1000._wp    &  ! excess density / dyn. visc. | -> cm/s to m/day
     &         * dustd2 * 1.e-4_wp)*dtb                      ! *diameter**2 |*1000 *1.e-4?
      wdust = dustsink

      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd1       = ',dustd1
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd2       = ',dustd2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustsink     = ',dustsink
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wdust        = ',wdust
#endif

      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'

#ifdef AGG
! parameters needed for the aggregation module (see Kriest 2002, DSR I vol.49, p. 2133-2162)

      SinkExp = 0.62_wp          ! exponent of the sinking speed vs. diameter relationship
      FractDim = 1.62_wp         ! exponent of the diameter vs. phosphorous content relationship
      Stick = 0.40_wp            ! maximum stickiness
      cellmass = 0.012_wp/rnit   ! [nmol P]   minimum mass of a particle in phosphorous units (rnit=16)
!ik   cellmass = 0.0039_wp/rnit  ! [nmol P] for 10 um diameter
      cellsink = 1.40_wp *dtb    ! [m/d] see Kriest 2002, Table 2 Ref 8 (from Stokes' formula, delta rho 0.052 g/cm3)
!ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
      shear = 86400._wp          ! wind induced shear in upper 100m , 1 d^-1
      fsh = 0.163_wp * shear *dtb !  turbulent shear (used for aggregation)
      fse = 0.125_wp * 3.1415927_wp * cellsink * 100._wp ! differential settling (used for aggregation) (100=10**2 [d**2])
      alow1 = 0.002_wp           !diameter of smallest particle [cm]
!ik      alow1 = 0.001_wp        !diameter of smallest particle [cm]
      alow2 = alow1 * alow1
      alow3 = alow2 * alow1
      alar1 = 1.0_wp             !diameter of largest particle for size dependend aggregation and sinking [cm]
      vsmall = 1.e-9_wp
      safe = 1.e-6_wp
      pupper = safe/((FractDim+safe)*cellmass)         ! upper boundary for cells per aggregate (?)
      plower = 1._wp/(1.1_wp*cellmass)                       ! lower   --------------"------------------
      zdis = 0.01_wp / ((FractDim + 0.01_wp)*cellmass)

!ik check max possible sinking speed in relation to min.
!ik layer thinkness and time step for all standard layers, except
!ik the bottom layer.
!ik if max possible sinking speed (per time step) is greater
!ik than min layer thickness, decrease max. length for sinking and
!ik aggregation

      zmini = 8000._wp
      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,kbo(i,j)-1
        if(pddpo(i,j,k).gt.0.5_wp) then
         zmini=min(pddpo(i,j,k),zmini)
        endif
      ENDDO
      ENDDO
      ENDDO

      CALL global_min(zmini)

      checksink =(zmini/cellsink)**(1._wp/SinkExp)*alow1

      if(checksink.lt.alar1) then

      write(io_stdo_bgc,*) 'Allowed max. length for sinking'              &
     & ,' with min. depth of '                                            &
     & , zmini, ' m for layers 1-(kbo-1) and time step of ',dtb           &
     & ,' days is' , checksink                                            &
     & ,'cm, which is smaller than prescribed value of', alar1, ' cm'

        talar1 = alar1
        alar1 = checksink
        write(io_stdo_bgc,*) 'Set max. length for sinking and aggregation &
     &  from ',talar1,' to ', alar1

      endif

      alar2 = alar1 * alar1
      alar3 = alar2 * alar1
      TSFac = (alar1/alow1)**SinkExp
      TMFac = (alar1/alow1)**FractDim

!ik check the maximum possible sinking speed for the bottom layer (which
!ik may be smaller than zmini, and write to array alar1max, tsfmax, tmfmax

      DO j=1,kpje
      DO i=1,kpie
         alar1max(i,j) = alar1
         TSFmax(i,j) = TSFac
         TMFmax(i,j) = TMFac
         if(pddpo(i,j,kbo(i,j)).gt.0.5_wp) then

!ik evaluate safe length scale for size dependent sinking and
!ik aggregation, and the resulting sinking rate and aggregation rate.

          checksink = (pddpo(i,j,kbo(i,j))/cellsink)**(1._wp/SinkExp)        &
     &                    *alow1
          if(checksink.lt.alar1) then
           alar1max(i,j) = checksink
           TSFmax(i,j) = (checksink/alow1)**SinkExp
           TMFmax(i,j) = (checksink/alow1)**FractDim
           write(io_stdo_bgc,*) 'resetting alar1 to',checksink,'at i =', &
     &     i,' j = ',j,' k = ', kbo(i,j), ' with dz = ',                 &
     &     pddpo(i,j,kbo(i,j))
          endif
        ENDIF
      ENDDO
      ENDDO

! for shear aggregation of dust:
      dustd1 = 0.0001_wp ![cm] = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustd3=dustd2*dustd1
      dustsink = (g * 86400._wp / 18._wp                & ! g * sec per day / 18.
     &         * (claydens - rhoref_water) / 1.567_wp * 1000._wp  & !excess density / dyn. visc.
     &         * dustd2 * 1.e-4_wp)*dtb                   ! --> 4.73e-2 m/d
      write(io_stdo_bgc,*) 'dust diameter (cm)', dustd1
      write(io_stdo_bgc,*) 'dust sinking speed (m/d)', dustsink / dtb
      if(dustsink.gt.cellsink) then
        write(io_stdo_bgc,*) 'dust sinking speed greater than cellsink'
        dustsink=cellsink
        write(io_stdo_bgc,*) 'set dust sinking speed to cellsink'
      endif

#endif /*AGG*/

      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,8
      DO  kmon=1,12
         chemcm(i,j,k,kmon)=rmasko
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO  i=1,kpie
      DO  j=1,kpje
      DO  k=1,kpke
         aksp(i,j,k)=rmasko
      ENDDO
      ENDDO
      ENDDO

!   initialisation of solid sediment fraction
      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,ks
      DO  l=1,nsedtra
         sedlay(i,j,k,l)=0._wp ! note that initialisation with zero causes 'clay digging from below' in sedshi.f90
         burial(i,j,l)=0._wp   ! deepest sediment layer (burial layer)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!   initialising flux to sediment (also used to sum up export production!)
      DO  j=1,kpje
      DO  i=1,kpie
!     IF(pddpo(i,j,k) .GT. 0.5) THEN     ! wet points (approach does not work due to double use [bottom mask and 90m mask])
         expoor(i,j)=0._wp
         expoca(i,j)=0._wp
         exposi(i,j)=0._wp
!     ELSE
!        expoor(i,j)=rmasko
!        expoca(i,j)=rmasko
!        exposi(i,j)=rmasko
!     ENDIF
      ENDDO
      ENDDO

!
!  Initial values for aquatic (advected) ocean tracers (for new start)
!
      phosat=1.e-6_wp
      phospac=2.5e-6_wp
      phosmed=2.e-7_wp
      silat=2.e-5_wp
      silmed=5.e-6_wp
      silpac=1.e-4_wp
      oxypac=1.5e-4_wp
      oxyat=2.5e-4_wp
      oxymed=2.e-4_wp
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

      IF(pddpo(i,j,k) .GT. 0.5_wp) THEN     ! wet points
         m=ibek(i,j)
         ocetra(i,j,k,isco212)=2.27e-3_wp     ! [kmol/m3]
         ocetra(i,j,k,ialkali)=2.37e-3_wp

         if(m.eq.6.or.m.eq.7.or.m.eq.8)then
         ocetra(i,j,k,iphosph)=phospac
         ocetra(i,j,k,ioxygen)=oxypac
         ocetra(i,j,k,isilica)=silpac
         endif
         if(m.le.5)then
         ocetra(i,j,k,iphosph)=phosat
         ocetra(i,j,k,ioxygen)=oxyat
         ocetra(i,j,k,isilica)=silat
         endif
         if(m.eq.9)then
         ocetra(i,j,k,iphosph)=phosmed
         ocetra(i,j,k,ioxygen)=oxymed
         ocetra(i,j,k,isilica)=silmed
         endif
         ocetra(i,j,k,igasnit)=0.0_wp
         ocetra(i,j,k,iano3)  =rnit*ocetra(i,j,k,iphosph)
         ocetra(i,j,k,idoc)   =1.e-10_wp
         ocetra(i,j,k,iphy)   =1.e-8_wp
         ocetra(i,j,k,izoo)   =1.e-8_wp
         ocetra(i,j,k,idet)   =1.e-8_wp
         ocetra(i,j,k,icalc)  =0._wp
         ocetra(i,j,k,iopal)  =1.e-8_wp
#ifdef __c_isotopes
         ocetra(i,j,k,isco214)=0.75_wp*2.27e-3_wp !Paddy: oldest water: 1600y --> X0.83
#endif
         ocetra(i,j,k,ian2o)  =0._wp
         ocetra(i,j,k,idms)   =0._wp
         ocetra(i,j,k,ifdust) =0._wp
         ocetra(i,j,k,iiron)  =0.6_wp*1.e-9_wp
!        ocetra(i,j,k,ibeten) =0.6*1.e-12
         hi(i,j,k)            =3.e-9_wp
         co3(i,j,k)           =0._wp                          ! this good for initialisation -> 2.e-4?
#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=2.27e-3_wp     ! adjusted to reference ratio 13C/12C=1 (*100)!
         ocetra(i,j,k,isco214)=2.27e-3_wp
         ocetra(i,j,k,idet13) =1.e-8_wp
         ocetra(i,j,k,idet14) =1.e-8_wp
         ocetra(i,j,k,icalc13)=0._wp
         ocetra(i,j,k,icalc14)=0._wp
#endif

#ifdef AGG
! calculate initial numbers from mass, to start with appropriate size distribution
         snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6_wp
         ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1._wp)
         ocetra(i,j,k,iadust) =0._wp
#endif /*AGG*/

#ifdef PANTHROPOCO2
         ocetra(i,j,k,isco2_ant)=ocetra(i,j,k,isco212)
         ocetra(i,j,k,ialk_ant) =ocetra(i,j,k,ialkali)
         ocetra(i,j,k,icalc_ant)=ocetra(i,j,k,icalc)
         hi_ant(i,j,k)          =hi(i,j,k)
         co3_ant(i,j,k)         =co3(i,j,k)
#endif
#ifdef ANTC14
         ocetra(i,j,k,iantc14)=ocetra(i,j,k,isco214)
#endif
#ifdef PCFC
         ocetra(i,j,k,icfc11)=0._wp
         ocetra(i,j,k,icfc12)=0._wp
#endif
      ELSE                              ! dry points
         ocetra(i,j,k,iphosph)=rmasko
         ocetra(i,j,k,isilica)=rmasko
         ocetra(i,j,k,ioxygen)=rmasko
         ocetra(i,j,k,ialkali)=rmasko
         ocetra(i,j,k,isco212)=rmasko
         ocetra(i,j,k,iano3)  =rmasko
         ocetra(i,j,k,igasnit)=rmasko
         ocetra(i,j,k,idoc)   =rmasko
         ocetra(i,j,k,iphy)   =rmasko
         ocetra(i,j,k,izoo)   =rmasko
         ocetra(i,j,k,idet)   =rmasko
         ocetra(i,j,k,icalc)  =rmasko
         ocetra(i,j,k,iopal)  =rmasko
         ocetra(i,j,k,ian2o)  =rmasko
         ocetra(i,j,k,idms)   =rmasko
         ocetra(i,j,k,ifdust) =rmasko
         ocetra(i,j,k,iiron)  =rmasko
!        ocetra(i,j,k,ibeten) =rmasko
         hi(i,j,k)            =rmasko
         co3(i,j,k)           =rmasko
#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=rmasko
         ocetra(i,j,k,isco214)=rmasko
         ocetra(i,j,k,idet13) =rmasko
         ocetra(i,j,k,icalc13)=rmasko
         ocetra(i,j,k,idet14) =rmasko
         ocetra(i,j,k,icalc14)=rmasko
#endif
#ifdef AGG
         ocetra(i,j,k,inos)   =rmasko
         ocetra(i,j,k,iadust) =rmasko
#endif /*AGG*/
#ifdef PANTHROPOCO2
         ocetra(i,j,k,isco2_ant)=rmasko
         ocetra(i,j,k,ialk_ant) =rmasko
         ocetra(i,j,k,icalc_ant)=rmasko
         hi_ant(i,j,k)          =rmasko
         co3_ant(i,j,k)         =rmasko
#endif
#ifdef ANTC14
         ocetra(i,j,k,iantc14)  =rmasko
#endif
#ifdef PCFC
         ocetra(i,j,k,icfc11)   =rmasko
         ocetra(i,j,k,icfc12)   =rmasko
#endif
      ENDIF

      ENDDO
      ENDDO
      ENDDO

#ifdef BGCLEVI

      CALL BGC_LEV(kpie,kpje,kpke,pddpo)     ! js: bgc_lev not in current code ?

#endif

#ifndef __cpl_dust   /* from IPCC_HAM, may be changed for super volcano */
      CALL GET_DUST(kpie,kpje,kpke,pddpo)
#endif

!
!  Initial values for sediment pore water tracers. (solid components?)
!
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie
      IF(bolay(i,j) .GT. 0._wp) THEN
         powtra(i,j,k,ipowaic)=ocetra(i,j,kbo(i,j),isco212)
         powtra(i,j,k,ipowaal)=ocetra(i,j,kbo(i,j),ialkali)
         powtra(i,j,k,ipowaph)=ocetra(i,j,kbo(i,j),iphosph)
         powtra(i,j,k,ipowaox)=ocetra(i,j,kbo(i,j),ioxygen)
         powtra(i,j,k,ipown2) =0._wp
         powtra(i,j,k,ipowno3)=ocetra(i,j,kbo(i,j),iano3)
         powtra(i,j,k,ipowasi)=ocetra(i,j,kbo(i,j),isilica)

         sedlay(i,j,k,issso12)=1.e-8_wp
         sedlay(i,j,k,isssc12)=1.e-8_wp
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=1.e-8_wp
         sedlay(i,j,k,isssc13)=1.e-8_wp
         sedlay(i,j,k,issso14)=1.e-8_wp
         sedlay(i,j,k,isssc14)=1.e-8_wp
#endif

         sedlay(i,j,k,issster)=30._wp
         sedlay(i,j,k,issssil)=3._wp

         sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
      ELSE
         powtra(i,j,k,ipowno3)=rmasks   ! pore water
         powtra(i,j,k,ipown2) =rmasks
         powtra(i,j,k,ipowaic)=rmasks
         powtra(i,j,k,ipowaal)=rmasks
         powtra(i,j,k,ipowaph)=rmasks
         powtra(i,j,k,ipowaox)=rmasks
         powtra(i,j,k,ipowasi)=rmasks
         sedlay(i,j,k,issso12)=rmasks   ! solid sediment
         sedlay(i,j,k,isssc12)=rmasks
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=rmasks
         sedlay(i,j,k,isssc13)=rmasks
         sedlay(i,j,k,issso14)=rmasks
         sedlay(i,j,k,isssc14)=rmasks
#endif
         sedlay(i,j,k,issssil)=rmasks
         sedlay(i,j,k,issster)=rmasks
         sedhpl(i,j,k)        =rmasks
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
! Values for bgcmean
!
      DO  l=1,nbgct2d
      DO  j=1,kpje
      DO  i=1,kpie
         bgct2d(i,j,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO


      DO  l=1,nbgcm2d
      DO  j=1,kpje
      DO  i=1,kpie
         bgcm2d(i,j,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nbgcm3d
      DO k=1,kwrbioz
      DO j=1,kpje
      DO i=1,kpie
         bgcm3d(i,j,k,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO l=1,nbgct3d
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
         bgct3d(i,j,k,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!js: output of sediment data

      DO l=1,nbgct_sed
      DO k=1,ks
      DO j=1,kpje
      DO i=1,kpie
         bgct_sed(i,j,k,l) = 0._wp
      ENDDO
      ENDDO
      ENDDO
      ENDDO

!
! values for sediment fluxes
!

      DO j=1,kpje
      DO i=1,kpie
      DO k=1,npowtra
        sedfluxo(i,j,k)=0._wp
      ENDDO
      ENDDO
      ENDDO
!
! fluxes of organic carbon, caco3, opal, and dust to the sediment
! (js: still old [bad] nomenclature from HAMOCC3....)
!
! change default from 0. to rmasko

      DO j=1,kpje
      DO i=1,kpie
        prorca(i,j)=0._wp
        prcaca(i,j)=0._wp
        silpro(i,j)=0._wp
        produs(i,j)=0._wp
      ENDDO
      ENDDO

!
!#ifdef __cpl_co2
!
!  Rate of continental carbon input per unit area
!
!
! Compute total ocean area
!
      ALLOCATE(ocean_area(kpie,kpje))

      ocean_area(:,:)=dlxp(:,:)*dlyp(:,:)*weto(:,:,1)
      CALL global_sum_2d_pio(ocean_area,totarea)

      calcinp = deltacalc * dtbgc / totarea
      orginp  = deltaorg  * dtbgc / totarea
      silinp  = deltasil  * dtbgc / totarea

      DEALLOCATE(ocean_area)

!
!  integrated continental carbon input (overwritten from restart file)
!
      calcinpglint = 0.0_wp
      orginpglint  = 0.0_wp
      silinpglint  = 0.0_wp

!#endif


!
!  atmospheric concentrations (overwritten from restart file)
!

      DO  j=1,kpje
      DO  i=1,kpie
         atm(i,j,iatmco2) = 278._wp
         atm(i,j,iatmo2)  = 196800._wp
         atm(i,j,iatmn2)  = 802000._wp
#ifdef __c_isotopes
         atm(i,j,iatmc13) = 278._wp - (278._wp * 0.0065_wp)
         atm(i,j,iatmc14) = 278._wp - (278._wp*0.0065_wp)**2
#endif
      ENDDO
      ENDDO

      IF (diffat) THEN
        atdifv(:,:)=1._wp

!  finding the equator   (p_ioff defined in ../src_oce/mo_parallel.f90)
!js this could be within IFDEF DIFFAT? (seems to be used only in atmotr, which is called only if DIFFAT (bgc.f90))
      write(io_stdo_bgc,*)'setting equatorial diffusion'
      DO  i=1,kpie
        !ii=1+(i+p_ioff-1)*2 ! global i-index for giph_g   | giph_g global latitude (parallelization)
        north=1._wp
        south=1._wp
        DO  j=1,kpje
          jj=1+(j+p_joff-1)*2 ! global j-index for giph_g
          IF(jj<=2) CYCLE
          ! was once north = giph_g(ii,jj-2), see doc/giph_and_giph_g.txt
          north=giph(2*i - 1, 2*j - 3)
          ! was once south=giph_g(ii,jj), see doc/giph_and_giph_g.txt
          south=giph(2*i - 1, 2*j - 1)
          if ((north .ge. 0._wp).and.(south.le.0._wp)) then
            atdifv(i,j)=0.01_wp
            if(j<=je-1) atdifv(i,j+1)=0.02_wp
            if(j>=   2) atdifv(i,j-1)=0.02_wp
            if(j<=je-2) atdifv(i,j+2)=0.05_wp
            if(j>=   3) atdifv(i,j-2)=0.05_wp
          endif
          if ((north .ge. 30._wp).and.(south .le. 30._wp)) then
            atdifv(i,j)=0.1_wp
            if(j<=je-1) atdifv(i,j+1)=0.2_wp
            if(j>=   2) atdifv(i,j-1)=0.2_wp
          endif
          if ((north .ge. -30._wp).and.(south .le. -30._wp)) then
            atdifv(i,j)=0.1_wp
            if(j<=je-1) atdifv(i,j+1)=0.2_wp
            if(j>=   2) atdifv(i,j-1)=0.2_wp
          endif
        ENDDO
      ENDDO

      ! tjahns: hope atdifv contains p-points, if not please correct
      CALL bounds_exch(1, 'p', atdifv, 'beleg_bgc diffat') ! for safety only

! put ENDIF for DIFFAT here?
!
!     no diffusion into the poles
!
      DO i=1,kpie
         if(have_g_js) atdifv(i,1)=0._wp
!        if(have_g_js) atdifv(i,2)=0._wp
         if(have_g_je) atdifv(i,kpje)=0._wp
         if(have_g_je) atdifv(i,je1)=0._wp
      ENDDO
      ENDIF

#ifdef PCFC
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)*2 ! global i-index for giph_g
        north=1._wp
        DO  j=1,kpje
          jj=1+(j+p_joff-1)*2 ! global j-index for giph_g
          IF(jj<=2) CYCLE
          north=giph_g(ii,jj)
          if (north .gt. 10._wp) then
            cfc_int(i,j) = 1._wp
          endif
          if (north .lt. -10._wp) then
            cfc_int(i,j) = 0._wp
          endif
          if ((north .le. 10._wp).and.(north .ge. -10._wp)) then
            cfc_int(i,j) = (north +10._wp)/20._wp
          endif
!         WRITE(io_stdo_bgc,*)'cfc_int: ',i,j,north,cfc_int(i,j)
        ENDDO
      ENDDO
#endif
#ifdef ANTC14
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)*2 ! global i-index for giph_g
        north=1._wp
        DO  j=1,kpje
          jj=1+(j+p_joff-1)*2 ! global j-index for giph_g
          north=giph_g(ii,jj)
          if (north .gt. 20._wp) then
            Rbomb(i,j) = D14C_north
          endif
          if (north .lt. -20._wp) then
            Rbomb(i,j) = D14C_south
          endif
          if ((north .le. 20._wp).and.(north .ge. -20._wp)) then
            Rbomb(i,j) = D14C_equator
          endif
!         WRITE(io_stdo_bgc,*)'Rbomb: ',i,j,north,Rbomb(i,j)
        ENDDO
      ENDDO
#endif
#ifdef PDYNAMIC_BGC
        DO k=1,nbgcdyn
        DO j=1,kpje
        DO i=1,kpie
          bgcdyntmp(i,j,k) = 0._wp
        ENDDO
        ENDDO
        ENDDO

        DO l=1,kdtot
        DO k=1,nbgcdyn
        DO j=1,kpje
        DO i=1,kpie
          bgcdyn(i,j,k,l) = 0._wp
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        DO j=1,kpje
        DO i=1,kpie
          bgc_zmld(i,j) = 0._wp
          bgc_nmld(i,j) = 0._wp
          nbgc_mld(i,j) = 0
        ENDDO
        ENDDO
#endif /* PDYNAMIC_BGC */


      END
