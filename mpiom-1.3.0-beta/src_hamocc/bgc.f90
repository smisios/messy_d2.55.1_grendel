      SUBROUTINE BGC(kpie,kpje,kpke,                                   &
     &    pfswr,psicomo,ptho,psao,pddpo,pdlxp,pdlyp,ptiestu,pdpio,     &
     &    pfu10,pwo,ppao,kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday)
     
!**********************************************************************
!
!**** *BGC* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     - time step biogeochemestry.
!
!     Method:
!     ------
!     kchck=1 is used to check max/min of bgc arrays on wet/dry cells.
!     Note: prosil is used only for k=1,2. It is adressed, however, for
!           k=1,4 in loop 100. To save memory, 
!
!     *CALL*  *BGC(kpie,kpje,kpke,pfswr,ptho,pddpo,pdlxp,pdlyp,pdpio)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h*        - .
!     *COMMON*     *COMMO1_BGC.h*        - .
!     *COMMON*     *CONTRBGC.h*   - control variables.
!     *COMMON*     *UNITS_BGC.h*        - .
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pfswr*   - solar radiation [W/m**2].
!     *REAL*    *psicomo* - sea ice concentration
!     *REAL*    *ptho*    - potential temperature [deg C].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ppao*    - sea level pressure [Pascal].
!     *REAL*    *pddpo*   - size of scalar grid cell (depth) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (longitudinal) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (latitudinal) [m].
!     *REAL*    *pdpio*   - inverse size of grid cell (depth)[m].
!     *INTEGER* *kldtmon* - time step within current month (from mpi-om) (1=1st time step of month, 2=2nd etc)
!     *INTEGER* *kldtday* - time step within current day (from mpi-om)
!
!     Externals
!     ---------
!     OCPROD CYANO CARCHM POWACH SEDSHI
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      USE mo_timeser_bgc
      use mo_param1_bgc 
#ifdef PDYNAMIC_BGC
      use mo_dynamic
#endif /* PDYNAMIC_BGC */ 
      
implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l

      REAL,intent(in) :: pfswr  (kpie,kpje)
      REAL,intent(in) :: psicomo(kpie,kpje)
      REAL,intent(in) :: pfu10(kpie,kpje)

      REAL,intent(in) :: ptho (kpie,kpje,kpke)
      REAL,intent(in) :: psao (kpie,kpje,kpke)
      REAL,intent(in) :: pddpo(kpie,kpje,kpke)
      REAL,intent(in) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL,intent(in) :: pdpio(kpie,kpje,kpke)
      REAL,intent(in) :: pwo(kpie,kpje,kpke+1)
      REAL,intent(in) :: ppao(kpie,kpje)            
      REAL,intent(in) :: ptiestu(kpke+1)
      INTEGER ,intent(in) :: kplyear,kplmon,kplday,kmonlen,kldtmon,kldtday


      

      
!      WRITE(io_stdo_bgc,*) 'BGC mit :  ',                                        &
!     &            kplyear,kplmon,kplday,kmonlen,kldtmon,               &
!     &            pfswr(50,50),psicomo(50,50),pfu10(50,50),            &
!     &            psao(50,50,1),pddpo(50,50,1),pdlxp(50,50),           &
!     &            pdlyp(50,50),ptiestu(1),ptho(50,50,1),               &            
!     &            pdpio(50,50,1),ppao(50,50),pwo(50,50,1)

!
! Increment bgc time step counter of run (initialized in INI_BGC).
!
      ldtrunbgc = ldtrunbgc + 1
!
! Increment bgc time step counter of experiment (initialized if IAUFR=0).
!
      ldtbgc = ldtbgc + 1
!
! Increment timeseries-1 sample counter (initialized in INI_TIMSER_BGC).
!

!
!--------------------------------------------------------------------
! Net solar radiation: multiply  with sea ice concentration

      DO  j=1,kpje
      DO  i=1,kpie
        strahl(i,j)=pfswr(i,j)*(1.-psicomo(i,j))
      ENDDO
      ENDDO

#ifdef PBGC_CK_TIMESTEP   
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'before BGC: call INVENTORY'
      call INVENTORY_BGC(kpie,kpje,kpke)
#endif
       call maschk(kpie,kpje,kpke,-1)
!
!--------------------------------------------------------------------
! Get flux or partial pressure of atmospheric CO2.
!
!      IF (MOD(ldtrunbgc-1,ndtdaybgc) .EQ. 0) THEN
!         CALL GET_ATMCO2
!      ENDIF
!_______________________________________________________________
!
! First timestep of the month --> kldtmon eq 1
! js check here annual output of 2d/bioz fields for glacial runs
      
      IF (kldtmon .EQ. 1) THEN 
                   
! Calculate the chemical constants for the current month.
!         WRITE(io_stdo_bgc,*) 'CHEMCON gerufen bei kldtmon: ',          &
!     &                         kplmon,kplday,kmonlen,kldtmon
!         CALL CHEMCON(-26,kpie,kpje,kpke,psao,ptho,                     &
!     &                 pddpo,pdlxp,pdlyp,ptiestu,kplmon)     

!#ifdef PCOMPONENT_ANALYSIS
!         CALL component_analysis                                        &
!     &       (kpie,kpje,kpke,pddpo,pdlxp,pdlyp,psao,ptho,psicomo,pfu10, &
!     &        kplyear,kplmon,kplday,kmonlen)
!#endif /* PCOMPONENT_ANALYSIS */

      ENDIF ! First timestep of the month


!      
!-------------------------------------------------------------------------      
! First timestep of the day --> kldtday eq 1     
      IF (kldtday .EQ. 1) THEN 

         CALL CHEMCON(-26,kpie,kpje,kpke,psao,ptho,                     &
              pddpo,pdlxp,pdlyp,ptiestu,kplmon)  
            
      ENDIF ! First timestep of the day
           
      meancnt_bgc_2D = meancnt_bgc_2D + 1 
      meancnt_bgc_3D = meancnt_bgc_3D + 1 
      
!_______________________________________________________________
!
!
!     Biogeochemistry

#ifdef PDYNAMIC_BGC
      call compute_dyn_diff(kpie,kpje,kpke,pddpo,psao,ptho,kdynsave,kdbio)
#endif /* PDYNAMIC_BGC */  
      
      CALL OCPROD(kpie,kpje,kpke,ptho,pddpo,pdlxp,pdlyp,pdpio,kplmon)
 
      do l=1,nocetra
      do K=1,kpke
      do J=1,kpje
      do I=1,kpie
        OCETRA(I,J,K,L)=MAX(0.,OCETRA(I,J,K,L))
      enddo
      enddo
      enddo
      enddo
       
#ifdef PDYNAMIC_BGC
      call compute_dyn_diff(kpie,kpje,kpke,pddpo,psao,ptho,kdyndiff,kdbio)
#endif /* PDYNAMIC_BGC */       

      call maschk(kpie,kpje,kpke,21)

      CALL CYANO(kpie,kpje,kpke,pddpo)

      call maschk(kpie,kpje,kpke,22)

      CALL CARCHM(kpie,kpje,kpke,pddpo,psao,ptho,psicomo,            &
           pfu10,kplmon,kplday,kmonlen,pdlxp,pdlyp)
 
      call maschk(kpie,kpje,kpke,23)

#ifdef PANTHROPOCO2 
      CALL CARCHM_ANT                                                  &
     &(kpie,kpje,kpke,pddpo,psao,ptho,psicomo,                         &
     & pfu10,kplmon,kplday,kmonlen)             ! check consistency (call of pdlxp,yp in carchm)
    
#if defined DIFFAT 
#ifndef __cpl_co2
      CALL ATMOTR(kpie,kpje,kpke,kplmon,pdlxp,pdlyp)   ! distribution of CO2 emissions to atm(iantco2) 
#endif
#endif      
#endif    

!_______________________________________________________________
!
!     Sediment module

      CALL POWACH(kpie,kpje,kpke,pdlxp,pdlyp,psao,pwo)

      call maschk(kpie,kpje,kpke,24)

      IF(MOD(ldtrunbgc,ndtdaybgc) .EQ. 0 ) THEN
         WRITE(io_stdo_bgc,*)                                          &
     &   'Sediment shifting at runstep ',ldtrunbgc,' ...'

        CALL SEDSHI(kpie,kpje)

      ENDIF

      IF (MOD(ldtrunbgc,nfreqts1).EQ.0) THEN
      
        CALL AVRG_TIMESER_BGC
      ENDIF


#ifdef PBGC_CK_TIMESTEP 
      WRITE(io_stdo_bgc,*)' '
      WRITE(io_stdo_bgc,*)'after BGC: call INVENTORY'
      call INVENTORY_BGC(kpie,kpje,kpke)
#endif       
      call maschk(kpie,kpje,kpke,101)

      RETURN
      END
