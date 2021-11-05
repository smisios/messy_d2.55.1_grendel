      MODULE mo_bgcmean


!***********************************************************************
!
!**** *MODULE mo_bgcmean* - Variables for bgcmean.
!
!     Patrick Wetzel    *MPI-Met, HH*    09.12.02
!  
!     Purpose
!     -------
!     - declaration and memory allocation for output fields
!
!**********************************************************************
      implicit none

! 2d fields -fluxes-

      INTEGER, PARAMETER ::                                   &
     &          ibase_2d=9,                                   &
     &          jco2flux  =1,                                 &
     &          jco214f   =2,                                 &         
     &          jo2flux   =3,                                 &
     &          jn2flux   =4,                                 &
     &          jn2oflux  =5,                                 &
     &          jprorca   =6,                                 &
     &          jprcaca   =7,                                 &
     &          jsilpro   =8,                                 &
     &          jprodus   =9

      INTEGER, PARAMETER ::                                   &
#ifdef __c_isotopes
     &          jc13flux  =ibase_2d+1,                        &
     &          jc14flux  =ibase_2d+2,                        &
     &          iiso_2d   =2
#else
     &          iiso_2d   =0
#endif


      INTEGER, PARAMETER :: nbgct2d=ibase_2d+iiso_2d    ! used for dimension of field bgct2d


!----------------------------------------------------------------      
! 2d fields for monthly means
      
      INTEGER, PARAMETER ::                                            &
     &         i_bsc_m2d=23,                                           &
     &         jkwco2   =1,                                            &
     &         jpco2    =2,                                            &
     &         jdms     =3,                                            &
     &         jdmsflux =4,                                            &
     &         jdmsprod =5,                                            &
     &         jdms_bac =6,                                            &
     &         jdms_uv  =7,                                            &
     &         jco2fxd  =8,                                            &
     &         jco2fxu  =9,                                            &
     &         joxflux  =10,                                           &
     &         jniflux  =11,                                           &
     &         jcoex90  =12,                                           &
     &         jopex90  =13,                                           &
     &         jcaex90  =14,                                           &
     &         jcoex1000=15,                                           &
     &         jopex1000=16,                                           &
     &         jcaex1000=17,                                           &
     &         jcoex2000=18,                                           &
     &         jopex2000=19,                                           &
     &         jcaex2000=20,                                           &
     &         jexport  =21,                                           &
     &         jexpoca  =22,                                           &
     &         jexposi  =23


      INTEGER, PARAMETER ::                                            &
#ifdef DIFFAT
     &     i_atm_m2d=3,                                                &
     &     jatmco2  =i_bsc_m2d+1,                                      &
     &     jatmo2   =i_bsc_m2d+2,                                      &
     &     jatmn2   =i_bsc_m2d+3
#else
     &     i_atm_m2d=0
#endif


      INTEGER, PARAMETER ::                                            &
#ifdef PCFC
     &     i_cfc_m2d=5,                                                &
     &     jac14fx  =i_bsc_m2d+i_atm_m2d+1,                  &
     &     jcfc11fx =i_bsc_m2d+i_atm_m2d+2,                  &
     &     jcfc12fx =i_bsc_m2d+i_atm_m2d+3,                  & 
     &     jpcfc11  =i_bsc_m2d+i_atm_m2d+4,                  & 
     &     jpcfc12  =i_bsc_m2d+i_atm_m2d+5
#else
     &     i_cfc_m2d=0
#endif         

      INTEGER, PARAMETER :: nbgcm2d=i_bsc_m2d+i_atm_m2d+i_cfc_m2d ! for dimension of bgcm2d

!until here 280906js
!----------------------------------------------------------------
! 2d: jdms,jexport, jexpoca, jexposi
! out : jpoc, jatten, jcalc, jopal (js: out=removed?)

      INTEGER, PARAMETER ::                                   &
     &          jphyto  =1,                                   &
     &          jgrazer =2,                                   &  
     &          jdoc    =3,                                   & 
     &          jphosy  =4,                                   &  
     &          jphosph =5,                                   &  
     &          joxygen =6,                                   &  
     &          jiron   =7,                                   &
     &          jano3   =8,                                   &
     &          jalkali =9,                                   & 
     &          jsilica =10,                                  &   
     &          jdic    =11,                                  &
     &          i_base_m3d=11

! conditional 3d fields
      INTEGER, PARAMETER ::                                   &
#ifdef __c_isotopes
     &          jdic13  =i_base_m3d+1,                        &
     &          jdic14  =i_base_m3d+2,                        & 
     &          i_iso_m3d=2                                  
#else
     &          i_iso_m3d=0
#endif

      INTEGER, PARAMETER ::                                   &
#ifdef AGG
     &          i_agg_m3d=1,                                  &
     &          jnos     =i_base_m3d+i_iso_m3d+1               
!    &          jwmass   =i_base_m3d+i_iso_m3d+2       ! js 22.2.2006 needs update of two above lines
#else 
     &          i_agg_m3d=0
#endif 

      INTEGER, PARAMETER :: nbgcm3d = i_base_m3d + i_iso_m3d + i_agg_m3d

!----------------------------------------------------------------

! js: sediment
!
      INTEGER, PARAMETER ::                                   &
     &          i_bsc_sed = 11,                               &
     &          jpowaic   =  1,                               &
     &          jpowaal   =  2,                               &
     &          jpowaph   =  3,                               &
     &          jpowaox   =  4,                               &
     &          jpown2    =  5,                               &
     &          jpowno3   =  6,                               &
     &          jpowasi   =  7,                               &
     &          jssso12   =  8,                               &
     &          jssssil   =  9,                               &
     &          jsssc12   = 10,                               &
     &          jssster   = 11          


      INTEGER, PARAMETER :: nbgct_sed = i_bsc_sed
     
!----------------------------------------------------------------
!  3d 'total' fields
     
      INTEGER, PARAMETER ::                                   &
     &          jphosph_t =1,                                 &  
     &          jano3_t   =2,                                 &
     &          jsilica_t =3,                                 & 
     &          jiron_t   =4,                                 &  
     &          joxygen_t =5,                                 &  
     &          jalkali_t =6,                                 & 
     &          jdic_t    =7,                                 & 
     &          jdoc_t    =8,                                 &  
     &          jpoc_t    =9,                                 & 
     &          jcalc_t   =10,                                & 
     &          jopal_t   =11
#ifdef __c_isotopes
      INTEGER, PARAMETER ::                                   &
     &          jdic13_t  =12,                                &
     &          jdic14_t  =13,                                &
     &          i_bsc_t3d=13
#else
      INTEGER, PARAMETER ::                                   &
     &          i_bsc_t3d=11
#endif

#ifdef PCFC
      INTEGER, PARAMETER ::                                   &
     &          i_cfc_t3d=3,                                  &         
     &          jcfc11_t=i_bsc_t3d+1,               & 
     &          jcfc12_t=i_bsc_t3d+2,               & 
     &          jac14_t =i_bsc_t3d+3
#else
      INTEGER, PARAMETER :: i_cfc_t3d = 0
#endif
      
      INTEGER, PARAMETER :: nbgct3d =i_bsc_t3d+i_cfc_t3d
     
      
      REAL, DIMENSION (:),       ALLOCATABLE :: stepspm      
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgct2d               ! 2d summed-up field (js, ?[patrick's 'total')
      REAL, DIMENSION (:,:,:),   ALLOCATABLE :: bgcm2d               ! 2d averaged  field
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgcm3d               ! 3d averaged  field
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct3d               ! 3d summed-up field
      REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: bgct_sed             ! 3d sediment  field
      
      INTEGER :: meancnt_bgc_2D,meancnt_bgc_3D

      INTEGER :: mean_2D_freq,mean_3D_freq
      INTEGER :: meantime_2d, nmeantime_2d
      INTEGER :: meantime_3d, nmeantime_3d

      INTEGER :: n90depth, n1000depth, n2000depth
      
      INTEGER :: nc_2d_id, nc_bioz_id, nc_3d_id, nc_sed_id
      
      
      CONTAINS

      SUBROUTINE ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

      use mo_control_bgc
      use mo_param1_bgc 
      
      INTEGER :: kpie,kpje,kpke
      
      WRITE(io_stdo_bgc,*)'Memory allocation for variable stepspm ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',nmeantime_2d
      
        ALLOCATE (stepspm(nmeantime_2d))      

      WRITE(io_stdo_bgc,*)'Memory allocation for variable bgct2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgct2d
      
        ALLOCATE (bgct2d(kpie,kpje,nbgct2d))      

      WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm2d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',nbgcm2d
      
        ALLOCATE (bgcm2d(kpie,kpje,nbgcm2d))      
      
      WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcm3d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kwrbioz
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgcm3d
      
        ALLOCATE (bgcm3d(kpie,kpje,kwrbioz,nbgcm3d))            

      WRITE(io_stdo_bgc,*)'Memory allocation for variable bgct3d ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
      WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgct3d
      
        ALLOCATE (bgct3d(kpie,kpje,kpke,nbgct3d))            

        WRITE(io_stdo_bgc,*)'Memory allocation for variable bgctsed ...'
        WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
        WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
        WRITE(io_stdo_bgc,*)'Third dimension    : ',ks
        WRITE(io_stdo_bgc,*)'Forth dimension    : ',nbgct_sed

        ALLOCATE (bgct_sed(kpie,kpje,ks,nbgct_sed))

      END SUBROUTINE ALLOC_MEM_BGCMEAN

      END MODULE mo_bgcmean
