      MODULE mo_param1_bgc


!***********************************************************************
!
!**** *MODULE mo_param1_bgc* - bgc tracer parameters.
!
!     Patrick Wetzel    *MPI-Met, HH*    01.09.03
!  
!     Purpose
!     -------
!     - declaration and memory allocation
!
!**********************************************************************
      implicit none
      
      INTEGER, PARAMETER :: ks=12,ksp=ks+1                       ! sediment layers
      INTEGER, PARAMETER :: kwrbioz=8                            ! euphotic layers


! advected tracers
      INTEGER, PARAMETER :: i_base_adv=14,                              &
     &                      isco212  =1,                                &
     &                      ialkali  =2,                                &
     &                      iphosph  =3,                                &
     &                      ioxygen  =4,                                &
     &                      igasnit  =5,                                &
     &                      iano3    =6,                                &
     &                      isilica  =7,                                &
     &                      idoc     =8,                                &
     &                      iphy     =9,                                &
     &                      izoo     =10,                               &
     &                      ian2o    =11,                               &
     &                      idms     =12,                               &
     &                      iiron    =13,                               &
     &                      ifdust   =14  
#ifdef __c_isotopes
      INTEGER, PARAMETER ::                                             &
     &                      isco213  =i_base_adv+1,                     &
     &                      isco214  =i_base_adv+2,                     &
     &                      i_iso_adv=2                    
#else 
      INTEGER, PARAMETER ::                                             &
     &                      i_iso_adv=0
#endif

      INTEGER, PARAMETER ::                                             &
#ifdef PCFC  
     &                      i_cfc_adv= 3,                               &
     &                      icfc11   = i_base_adv+i_iso_adv+1,          &
     &                      icfc12   = i_base_adv+i_iso_adv+2,          & 
                            iantc14  = i_base_adv+i_iso_adv+3
#else 
     &                      i_cfc_adv= 0
#endif
      INTEGER, PARAMETER ::                                             &
#ifdef AGG
     &                      i_agg_adv= 2,                               &
     &                      inos     = i_base_adv+i_iso_adv+i_cfc_adv+1,&
     &                      iadust   = i_base_adv+i_iso_adv+i_cfc_adv+2
#else 
                            i_agg_adv= 0
#endif

! total number of advected tracers
      INTEGER, PARAMETER :: ntraad=i_base_adv+i_iso_adv+i_cfc_adv+i_agg_adv

! non-advected (fast sinking) tracers
      INTEGER, PARAMETER ::                                             &
     &                      idet     =ntraad+1,                         &
     &                      icalc    =ntraad+2,                         &
     &                      iopal    =ntraad+3,                         &
     &                      i_base   =3 
                            
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      idet13   =ntraad+i_base+1,                  &
     &                      icalc13  =ntraad+i_base+2,                  &      
     &                      idet14   =ntraad+i_base+3,                  &
     &                      icalc14  =ntraad+i_base+4,                  &
     &                      i_iso    =4  
#else
     &                      i_iso    =0

#endif

     INTEGER, PARAMETER :: nocetra = ntraad+i_base+i_iso

! 10Be to be done later    &                     ,ibeten=24                              


! atmosphere
      INTEGER, PARAMETER :: iatmco2=1,                                  &
     &                      iatmo2 =2,                                  &
     &                      iatmn2 =3,                                  &    
     &                      i_base_atm=3
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      iatmc13 = i_base_atm+1,                     &
     &                      iatmc14 = i_base_atm+2,                     &
     &                      i_iso_atm = 2
#else
     &                      i_iso_atm = 0
#endif

     INTEGER, PARAMETER ::  natm=i_base_atm+i_iso_atm

! sediment
      INTEGER, PARAMETER :: issso12=1,                                  &
     &                      isssc12=2,                                  &
     &                      issssil=3,                                  &
     &                      issster=4,                                  &
     &                      nsss_base=4
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      issso13=nsss_base+1,                        &
     &                      issso14=nsss_base+2,                        &
     &                      isssc13=nsss_base+3,                        &
     &                      isssc14=nsss_base+4,                        &
     &                      nsss_iso = 4
#else
                            nsss_iso = 0
#endif
      INTEGER, PARAMETER :: nsedtra=nsss_base+nsss_iso
     
! pore water tracers, index must be the same as for ocetra otherwise problems in dipowa.f90!
      INTEGER, PARAMETER :: ipowaic=1,                                  &
     &                      ipowaal=2,                                  &
     &                      ipowaph=3,                                  &
     &                      ipowaox=4,                                  &
     &                      ipown2 =5,                                  &
     &                      ipowno3=6,                                  &
     &                      ipowasi=7,                                  &
     &                      npowa_base = 7
      INTEGER, PARAMETER ::                                             &
#ifdef __c_isotopes
     &                      ipowc13=npowa_base+1,                       &
     &                      ipowc14=npowa_base+2,                       &
     &                      npowa_iso = 2
#else
     &                      npowa_iso = 0
#endif
      INTEGER, PARAMETER :: npowtra=npowa_base+npowa_iso

      END MODULE mo_param1_bgc
