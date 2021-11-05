    SUBROUTINE READ_NAMELIST

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/read_namelist.f90,v $\\
!$Revision: 1.2.2.1.4.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!****************************************************************
!
!**** *READ_NAMELIST* - set namelist variables.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - set default values of namelist variables
!     - read from namelist
!     - print values of namelist variables
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!     called by ini_bgc
!
!     *CALL*       *READ_NAMELIST
!
!     *COMMON*     *PARAM1_BGC.h*   - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h*   - ocean/sediment tracer arrays.
!     *COMMON*     *NAMELIST_BGC.h* - bgc modules namelist.
!     *COMMON*     *UNITS_BGC.h*    - std I/O logical units.
!     *MODULE*     *mo_timeser_bgc* - parameter and memory for time series.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************
!ik added trap depths for time series

      USE mo_control_bgc
      USE mo_timeser_bgc
      use mo_param1_bgc 
      use mo_bgcmean 

      use mo_mpi

#ifdef PANTHROPOCO2
      USE mo_carbch
#endif

implicit none      

      INTEGER :: L

      NAMELIST /BGCCTL/ rmasko,rmasks,kchck,isac                      &
     &                 ,mean_2D_freq,mean_3D_freq                     &      
     &                 ,io_stdo_bgc,io_stdi_bgc                       &
     &                 ,io_rsti_bgc,io_rsto_bgc,io_timeser_bgc        &
     &                 ,rlonts1,rlatts1,nfreqts1                      &
#ifndef __cpl_co2                         /*for millenium version*/
#ifdef PANTHROPOCO2
#ifdef DIFFAT
     &                 ,emission                                      &
#else     
     &                 ,co2_atm_1,co2_atm_2,co2_atm_3                 &
#endif /*DIFFAT*/
#endif  /*PANTHROPOCO2*/  
#endif /*! __cpl_co2*/
#ifdef ANTC14
     &                 ,D14C_north,D14C_south,D14C_equator            &
#endif
#ifdef PCFC
     &                 ,cfc11_atm_1s,cfc11_atm_2s,cfc11_atm_3s        &
     &                 ,cfc11_atm_1n,cfc11_atm_2n,cfc11_atm_3n        &
     &                 ,cfc12_atm_1s,cfc12_atm_2s,cfc12_atm_3s        &
     &                 ,cfc12_atm_1n,cfc12_atm_2n,cfc12_atm_3n        &    
#endif
     &                 ,rdep1ts1,rdep2ts1,rdep3ts1,lspinbgc


! Coordinates etc. of time series 1 
      nfreqts1=0
      DO l=1,nts
         rlonts1(l)=-1.0
         rlatts1(l)=-1.0
         rdep1ts1(l) = 100.
         rdep2ts1(l) = 200.
         rdep3ts1(l) = 300.
      ENDDO

!                        
! Read BGCCTL namelist
!
      IF(p_pe==p_io) THEN
      !mz_ap_20070515
        OPEN(io_stdi_bgc,FILE='NAMELIST_BGC.nml',STATUS='UNKNOWN',           &
     &        ACCESS='SEQUENTIAL',FORM='FORMATTED')

        READ(io_stdi_bgc,BGCCTL)

        CLOSE(io_stdi_bgc)
      ENDIF

      CALL p_bcast(lspinbgc,p_io)
      CALL p_bcast(rmasko,p_io)
      CALL p_bcast(rmasks,p_io)
      CALL p_bcast(kchck,p_io)
      CALL p_bcast(isac,p_io)
      CALL p_bcast(mean_2D_freq,p_io)
      CALL p_bcast(mean_3D_freq,p_io)
      CALL p_bcast(io_stdo_bgc,p_io)
      CALL p_bcast(io_stdi_bgc,p_io)
      CALL p_bcast(io_rsti_bgc,p_io)
      CALL p_bcast(io_rsto_bgc,p_io)
      CALL p_bcast(io_timeser_bgc,p_io)
      CALL p_bcast(rlonts1,p_io)
      CALL p_bcast(rlatts1,p_io)
      CALL p_bcast(nfreqts1,p_io)
      CALL p_bcast(rdep1ts1,p_io)
      CALL p_bcast(rdep2ts1,p_io)
      CALL p_bcast(rdep3ts1,p_io)
#ifndef __cpl_co2
#ifdef PANTHROPOCO2
#ifdef DIFFAT
      CALL p_bcast(emission,p_io)
#else
      CALL p_bcast(co2_atm_1,p_io)
      CALL p_bcast(co2_atm_2,p_io)
      CALL p_bcast(co2_atm_3,p_io)
#endif /*DIFFAT*/
#endif  /*PANTHROPOCO2*/
#endif /*__cpl_co2*/
#ifdef ANTC14
      CALL p_bcast(D14C_north,p_io)
      CALL p_bcast(D14C_south,p_io)
      CALL p_bcast(D14C_equator,p_io)
#endif
#ifdef PCFC
      CALL p_bcast(cfc11_atm_1s,p_io)
      CALL p_bcast(cfc11_atm_2s,p_io)
      CALL p_bcast(cfc11_atm_3s,p_io)
      CALL p_bcast(cfc11_atm_1n,p_io)
      CALL p_bcast(cfc11_atm_2n,p_io)
      CALL p_bcast(cfc11_atm_3n,p_io)
      CALL p_bcast(cfc12_atm_1s,p_io)
      CALL p_bcast(cfc12_atm_2s,p_io)
      CALL p_bcast(cfc12_atm_3s,p_io)
      CALL p_bcast(cfc12_atm_1n,p_io)
      CALL p_bcast(cfc12_atm_2n,p_io)
      CALL p_bcast(cfc12_atm_3n,p_io)
#endif


!
!  Open std out file
!

      CALL OPEN_STDOUT(io_stdo_bgc,'bgcout')

      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'* Values of BGCCTL namelist variables : '
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              kchck      = ',kchck
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rmasko     = ',rmasko
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rmasks     = ',rmasks
      WRITE(io_stdo_bgc,*)                                             &
     &'* Time series :'
      WRITE(io_stdo_bgc,*)                                             &
     &'*                                nts      = ',nts,'  elements.'
      WRITE(io_stdo_bgc,*)                                             &
     &'*     sampled at / averaged over nfreqts1 = ',nfreqts1          &
     &                                             ,' model time steps.'
      DO l=1,nts
        WRITE(io_stdo_bgc,*)                                           &
     &'*         Position of element ',l,' at    = '                   &
     &                                          ,rlonts1(l),rlatts1(l)
      ENDDO

      WRITE(io_stdo_bgc,*)                                             &
     &'* Standard output unit        io_stdo_bgc = ',io_stdo_bgc
      WRITE(io_stdo_bgc,*)                                             &
     &'* '

      WRITE(io_stdo_bgc,*) '* sediment acceleration factor isac read'
      WRITE(io_stdo_bgc,*) '* sediment acceleration factor is', isac
      If (isac.eq.0) then
        WRITE(io_stdo_bgc,*)                                           &
     &      '* sediment acceleration factor isac not defined'
        WRITE(io_stdo_bgc,*)                                           &
     &      '* isac is set to 1'
        WRITE(io_stdo_bgc,*)'*' 
      isac = 1
      endif     

      WRITE(io_stdo_bgc,*) '* mean_2D_freq is set to  : ', mean_2D_freq
      WRITE(io_stdo_bgc,*) '* mean_3D_freq is set to: ', mean_3D_freq


#ifndef __cpl_co2
#ifdef PANTHROPOCO2
#ifdef DIFFAT
      emission = emission/1000. ! million metric tons --> GigaTons
      ems_per_step=(1.e12/12.)*emission/ndtrunbgc
      WRITE(io_stdo_bgc,*)                                             &
     &'* '      
      WRITE(io_stdo_bgc,*)                                             &
     &      '* PgC CO2 added to atmosphere in run:    ',emission
      WRITE(io_stdo_bgc,*)                                             &
     &      '* kmol CO2 added to atmosphere per step, ',ems_per_step
#else /*DIFFAT*/
      WRITE(io_stdo_bgc,*)                                             &
     &'* '      
      WRITE(io_stdo_bgc,*)                                             &
     &      '* Atmospheric CO2 start  of year:   ',co2_atm_1
      WRITE(io_stdo_bgc,*)                                             &
     &      '* Atmospheric CO2 middle of year:   ',co2_atm_2
      WRITE(io_stdo_bgc,*)                                             &      
     &      '* Atmospheric CO2 end    of year:   ',co2_atm_3
#endif /*DIFFAT*/
#endif /*ANTHROPOCO2*/
#endif /*__cpl_co2*/
#ifdef ANTC14
      WRITE(io_stdo_bgc,*)                                             &
     &      '* Anthropogenic DC14 north:        ',D14C_north
      WRITE(io_stdo_bgc,*)                                             &
     &      '* Anthropogenic DC14 south:        ',D14C_south
      WRITE(io_stdo_bgc,*)                                             & 
     &      '* Anthropogenic DC14 equator:      ',D14C_equator 
#endif
#ifdef PCFC
      WRITE(io_stdo_bgc,*)                                             &
     &      '* CFC11 south:' ,cfc11_atm_1s,cfc11_atm_2s,cfc11_atm_3s      
      WRITE(io_stdo_bgc,*)                                             &
     &      '* CFC11 north:' ,cfc11_atm_1n,cfc11_atm_2n,cfc11_atm_3n      
      WRITE(io_stdo_bgc,*)                                             &
     &      '* CFC12 south:' ,cfc12_atm_1s,cfc12_atm_2s,cfc12_atm_3s      
      WRITE(io_stdo_bgc,*)                                             &
     &      '* CFC12 north:' ,cfc12_atm_1n,cfc12_atm_2n,cfc12_atm_3n      

#endif


      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'

      RETURN
      END
