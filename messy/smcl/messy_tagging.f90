!*************************************************************************
MODULE messy_tagging
!*************************************************************************

  ! MODULE FOR CALCULATIN DIAGNOSTIC OZONE ORIGIN TRACER
  !
  ! MESSy-SMCL
  !
  ! Description: Grewe, 2004, Technical Note, ACP4.
  !              Grewe et al., 2010, GMD.
  !              Grewe, 2013, GMD.
  !  See messy_tagging_si for more detail.
  !
  ! Authors:
  !    Volker Grewe, DLR, August 2008.
  !    Eleni Tsati, DLR, 2013.

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, SP, I4, I8
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG, STRLEN_XLONG
 
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP, SP, I4, I8

  ! ----------- <



  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'tagging'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.1'

  INTEGER,PUBLIC  :: i_integrate=1,      &    ! Integration method=Euler backward
                     i_tracer_init=0,    &    ! Internal initialisation
                     i_species=3,        &    ! full mechanism
                     i_advect_scaling=1       ! adjust possible errors due to 
                     !i_method=2               ! numericl diffusion of advection
  LOGICAL, PUBLIC :: l_adv_err_diag=.false.   ! diagnose the impact of these errors.
  

  INTEGER :: i_tagging_ctrl

  !op_mm_20171610 (Number of categories is not controlled via namelist!)
  !INTEGER, PARAMETER, PUBLIC ::   N_tag_species=10
                                ! Number of categories for tagged species
                                ! In principle 8 Emissions + CH4 
                                !       + O2 photolosis for O3

                                !op_mm_20171610 (Not set in SMIL) 
!  INTEGER, PARAMETER, PUBLIC :: ilig=1, ibio=2, isoi=3, iind=4, itra=5, &    
!                                ishp=6, iair=7, in2o=8, ich4=9, istr=10
 
  REAL(DP), PARAMETER        :: r_small=1.e-20

!op_mm_20171610+
  !set in namelist nr of tagging categories  
  INTEGER, public  ::   N_tag_species=0

  !ids of special catagories (set after reading the namelists)  
  INTEGER,public :: ilig=0, isoi=0, in2o=0, ich4=0, istr=0

  TYPE T_CATEGORY_IN
     CHARACTER(LEN=STRLEN_ULONG)  :: abr  = '' ! abreviation
     CHARACTER(LEN=STRLEN_ULONG)  :: long_name  = '' ! long_name
     CHARACTER(LEN=STRLEN_ULONG)  :: cat_type  = '' ! category_taypes
  END TYPE T_CATEGORY_IN
#ifdef LF
  PUBLIC :: T_CATEGORY_IN
#endif
!op_mm_20171610-
 




  INTEGER, PARAMETER,PUBLIC                           :: N_MAX_CATEGORY = 30
  TYPE(T_CATEGORY_IN), DIMENSION(N_MAX_CATEGORY),PUBLIC, SAVE :: CATEGORY_IN





  ! SUBROUTINES
  Public :: tagging_chemistry                      ! calculate short lived species
  PUBLIC :: tagging_get_tendencies                 ! calculate the rest of the 
                                                   ! tendencies by closing the budget
  PUBLIC :: tagging_tendencies                     ! adapt tendencies to integration scheme
  PUBLIC :: tagging_HOx_tendencies                 ! correction factors for HOx tendencies
                                                   !
  PUBLIC :: tagging_read_nml
  PUBLIC :: parse_catstr
  !op_mm_20171610+
  
  
CONTAINS

! ----------------------------------------------------------------------
  
  SUBROUTINE tagging_tendencies(kproma,nlev,dt,X1,X2,P,D,alpha,beta)

    ! Adapt Production and Loss terms to the chosen tagging integration scheme.
    ! I.e. the MECCA chemistry uses whatever is chosen. 
    !      and produces some Loss and Prod. rates. for a species X
    !      In the tagging these rates are used but not the same integration scheme
    !      hence there are inconsistencies. The same species wouldn't give the
    !      same results. Therefore this routine fixes it
    !      The idea is to calculate new tendencies, which are as close as possible
    !      to the old ones.
    !      Pnew=alpha*P   and Dnew=beta*D, such that
    !          X_2=Integration_scheme(X_1,P,D)
    !             and
    !         f(a,b)=((alpha-1)*P)^2+((beta-1)*D)^2 = Min       a=alpha b=beta
    !
    ! Integration scheme: X2=(X1+P*DT)/(1+D/X2*dt)
    !      Find a and b such that
    !          X2=(X1+a*P*DT)/(1+b*D/X2*dt)
    !       i.e. b(a)=P/D*a + (X1-X2)/D/dt                   (1)
    !         f'(a,b(a))=4P^2 a + 2P ((x1-x2)/dt-D-P)
    !         f''(a,b(a)) = 4P^2 > 0
    !         f'=0  => b=1-P/D (a-1)                         (2)
    !    (1)-(2)       a=((x2-x1)/dt+P+D)/(2P)
    !    =>            b=((x1-x2)/dt+P+D)/(2D)
    ! (a,b) valid?
    !       a,b must be non-negative hence
    !       from b>=0 and (1) we obtain a>=(x2-x1)/(P dt)=a_min (or 0 if term neg.)
    !       from b>=0 and (2) we obtain a<=1+D/P         =a_max
    !       -In the case of a_max<a_min optimisation fails 
    !           it follows that X2-X1> (P+D) dt > 0 ! 
    !                           and b<0 always
    !                =>    b=0 and a=(x2-x1)/(P dt)>0
    !       - In case a<a_min b or a are negative. 
    !                 The nearest solution is a=a_min and b=(1)
    !       - In case a>a_min b is negative 
    !                 The nearest solution is a=a_max and b=(1)
    !                 Note that in this case (2) doesn't hold and 
    !                 b is not necessarily 0
    !                 (Does it imply that there is a better solution?)
    ! Author Volker Grewe, DLR-OP, August 2008, durig NCAR visit, Boulder

  IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(IN)    :: kproma,nlev
    REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: X1, X2        ! species at two timesteps 
                                                              !           in [mol/mol]
    REAL(DP),                  INTENT(IN)    :: dt            ! timestep  [s]
    REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: P,D           ! respective loss and 
                                                              ! prod rate [mol/mol/s]
    REAL(DP), DIMENSION(:,:),  INTENT(OUT)   :: alpha,beta    ! Correction factors 
                                                              ! alpha and beta >=0

    ! Local vars
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_tendencies'

    REAL(DP)                                 :: a,b           ! interims alpha, beta
    REAL(DP)                                 :: a_min,a_max   ! a_min<=a<=a_max 
                                                              !  valid range
     
    INTEGER                                  :: i,j
                                                            
    do j=1,nlev
    do i=1,kproma
!------
! Make sure that P,D <> 0.
!------
        if (P(i,j).lt.r_small .and. D(i,j).lt.r_small) then 
            a=1.
            b=1.
        elseif (P(i,j).lt.r_small) then
            a=1.
            b= ((X2(i,j)-X1(i,j))/dt-P(i,j))/D(i,j)
            b=MAX(0._dp,b)
        elseif (D(i,j).lt.r_small) then
            b=1.
            a=((X2(i,j)-X1(i,j))/dt+D(i,j))/P(i,j)
            a=MAX(0._dp,a)
        else
!------
! Calculate a,b, and valid range for a
!------
        a_min=MAX(0._dp,(X2(i,j)-X1(i,j))/P(i,j)/dt)
        a_max=1+D(i,j)/P(i,j)
        a=((X2(i,j)-X1(i,j))/dt+P(i,j)+D(i,j))/(2*P(i,j))
        b=((X1(i,j)-X2(i,j))/dt+P(i,j)+D(i,j))/(2*D(i,j))
!------
! Test validity of solution and adjust it accordingly
!------
        if (a_max.lt.a_min) then   ! b always < 0 for all a
           a=(X2(i,j)-X1(i,j))/P(i,j)/dt
           b=0._dp
        elseif (a.lt.a_min) then 
           a=a_min
           b=P(i,j)/D(i,j)*a+(x1(i,j)-x2(i,j))/(D(i,j)*dt)
        elseif (a.gt.a_max) then
           a=a_max
           b=P(i,j)/D(i,j)*a+(x1(i,j)-x2(i,j))/(D(i,j)*dt)
        endif
        endif ! PD, small
!------
! Update Prod and Dest rates after call of this subroutine
!------
!       P=a*P
!       D=b*D
        alpha(i,j)=a
        beta(i,j)=b
    enddo
    enddo
     
  END SUBROUTINE tagging_tendencies

!------------------------------------------------------------------------------------------------

  SUBROUTINE tagging_HOx_tendencies(kproma,nlev,OH,HO2,P,D1,D2,alpha,beta,gamma)
                           
    ![OH]t = P2.../D2...  , [HO2]t = P1... /D1...
    ![OH]t = alpha*P2.../beta*D2... (1a) , [HO2]t = delta*P1... /gamma*D1... (1b)
    !alpha*P2=beta*D2 (1a), delta*P1=gamma*D1  (1b)
    !min. the error: (alpha-1)**2+(beta-1)**2+(gamma-1)**2+(delta-1)**2  (2)
    !alpha = (1+P2/D2)/(1+P2**2/D2**2) (2a), delta=(1+P1/D1)/(1+P1**2/D1**2)  (2b)
    !q2=P2/D2  (2a), q1=P1/D1 (2b)
    !
    !!only if: Loss_OH=Prod_HO2 in total
    ! Note:  P,D1,D2 in mol/mol/s
    ! Calculate correction terms so that    
    !            alpha P                     beta D1
    !      OH= ------------            HO2=--------------  
    !          beta (D1/OH)                 gamma D2/HO2
    ! or
    ! (1)  alpha P = beta D1 = gamma D2
    ! by minimizing (alpha-1)**2+(beta-1)**2+(gamma-1)**2  (2)
    !                alpha P          alpha P
    ! (1) gives beta=--------   gamma=--------
    !                   D1              D2
    ! (1) in (2) + derivate=0:
    ! alpha = (1+P/D1+P/D2)/(1+P*P/D1/D1+P*P/D2/D2)
    ! 
    ! Author Volker Grewe, DLR-OP, August 2008, durig NCAR visit, Boulder

! USE??  useless

  IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(IN)    :: kproma,nlev
    REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: OH,HO2
    REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: P,D1,D2       ! respective loss and 
                                                !P1,P2        ! prod rate [mol/mol/s]
    REAL(DP), DIMENSION(:,:),  INTENT(OUT)   :: alpha,beta,&  ! Correction factors >=0
                                                gamma !,delta         ! 

    ! Local vars
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_HOx_tendencies'
    INTEGER                                  :: i,j
    REAL(DP)                                 :: a,b,c    ! interims alpha, beta, gamma
                                                         ! d  delta 
    REAL(DP)                                 :: q1,q2    ! P/D1 P/D2

    do j=1,nlev
    do i=1,kproma
       if (D1(i,j).lt.r_small .or. D2(i,j).lt.r_small & 
                              .or. P(i,j).lt.r_small) then 

     ! if (D1(i,j).lt.r_small .or. D2(i,j).lt.r_small & 
     !  .or. P1(i,j).lt.r_small .or. P2(i,j).lt.r_small) then 
  
           a=1._dp
           b=1._dp
           c=1._dp
          !d=1._dp

       else
           q1=P(i,j)/D1(i,j)
           q2=P(i,j)/D2(i,j)
           a=(1+q1+q2)/(1+q1*q1+q2*q2)
          !q1=P1(i,j)/D1(i,j)
          !q2=P2(i,j)/D2(i,j)
          !a=(1+q2)/(1+q2*q2)
          !d=(1+q1)/(1+q1*q1)
           !b=a*q2
          !c=d*q1
           b=a*q1
           c=a*q2
       endif
       alpha(i,j)=a
       beta(i,j)=b
       gamma(i,j)=c
       !delta(i,j)=d
    enddo
    enddo

! WRITE(*,*) 'HOX a,b,s',MINVAL(alpha(:,:)),MAXVAL(alpha(:,:)),MINVAL(beta(:,:)),MAXVAL(beta(:,:)),&
!       MINVAL(beta(:,:)),MAXVAL(beta(:,:))

  END SUBROUTINE tagging_HOx_tendencies

!----------------------------------------------------------------------------------

  SUBROUTINE tagging_get_tendencies(kproma,nlev,dt,            &
             ch4_1,ch4_2,co_1,co_2,pan_1,pan_2,nmhc_1,nmhc_2,  &
             ch4loss,coprod,coloss)

    ! Re-Calculate the poduction and loss terms of the C-compounds based on 
    ! MECCA output
    ! X_1 values before MECCA
    ! X_2 values after MECCA
    ! In this concept both prod and loss terms are non-negative   
    ! ch4loss=-(CH4_2-CH4_1), assuming that there is no CH4 chemical production
    !                      i.e. neglecting C3H6+O3 -> ... + 0.06 CH4 
    ! coloss=CH4_2+NMHC_2+2*PAN_2+CO_2 - (CH4_1+NMHC_1+2*PAN_1+CO_1),
    !                      i.e. C can only be lost via CO,
    !                      i.e. neglecting PA+CH3O2 -> CO2 and PA + OH-> 2 CO2 
    ! coprod=CO_2-CO_1+coloss
    ! 
    ! For some reasons loss or prod terms may not be strictly positive and 
    !                                                            negative,
    !    respectively. This might have 2 reasons
    !    - Chemistry code is inaccurate
    !    - Reactions might be included, which screw up the whole idea
    !    First one cannot be fixed
    !    Second are assumed to be small. Just cut off
    !    Note: 
    !    CO:     CO2 +  hv  =  CO + O3P   not accounted 
    !            coloss negative at some grid points in boundary layer 
    !            (factor 1000 lower than 'normal coloss')
    !    CH4: Tests showed that it is fine
    ! Author Volker Grewe, DLR-OP, August 2008, durig NCAR visit, Boulder

! USE??    or useless

  IMPLICIT NONE

    ! I/O
    INTEGER,                   INTENT(IN)    :: kproma,nlev
    REAL(DP),                  INTENT(IN)    :: dt            ! timestep  [s]
    REAL(DP), DIMENSION(:,:),  INTENT(IN)    :: CH4_1,CH4_2,& ! species at 
                                                              ! two timesteps
                                                              ! in [mol/mol]
                                                CO_1,CO_2,  &
                                                NMHC_1,NMHC_2, &
                                                PAN_1,PAN_2
    REAL(DP), DIMENSION(:,:),  INTENT(OUT)   :: ch4loss, coloss, coprod

    ! Local vars
!   CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_tendencies'

! Determine ch4loss
! in case ch4 is produced somewhere simply neglect it
    ch4loss(:,:)=-(CH4_2(:,:)-CH4_1(:,:))/dt
    ch4loss(:,:)=MAX(ch4loss(:,:),0._dp)

! Determine co loss rate
    coloss(:,:)=-(  (CH4_2(:,:)+NMHC_2(:,:)+2._dp*PAN_2(:,:)+CO_2(:,:))  &
                   -(CH4_1(:,:)+NMHC_1(:,:)+2._dp*PAN_1(:,:)+CO_1(:,:)) ) / dt
! Occasionally this may be negatve !
    coloss(:,:)=MAX(coloss(:,:),0._dp)
! and occationally the coloss rate may not be sufficient to 
! account for the total chemical CO decrease.
    coloss(:,:)=MAX(coloss(:,:),-((CO_2(:,:)-CO_1(:,:))/dt))
! Now we have a well defined (=large enough to account for co change)
!              and non-negative coloss
    coprod(:,:)= (CO_2(:,:)-CO_1(:,:))/dt+coloss(:,:)
    
  END SUBROUTINE tagging_get_tendencies

! ----------------------------------------------------------------------


SUBROUTINE tagging_chemistry(kproma,nlev,dt,                &
        o3_tagging,noy_tagging,pan_tagging,                 &
        nmhc_tagging,co_tagging,                            &
        co_tagging_sum,nmhc_tagging_sum,noy_tagging_sum,    &
        o3_tagging_sum,pan_tagging_sum,                     &
        oh, ho2, do3,dnoy,dpan,dnmhc,dco,                   &
        oh_tag, ho2_tag,                                    &
        o3prod_ho2,o3prod_ro2,o3prod_o2,                    &
        o3loss_oh,o3loss_ho2,o3loss_ro,                     &
        o3loss_no,o3loss_xo,o3loss,                         &
        panprod,panloss,coprod,                             &
        coloss,ch4loss,ch4loss_tag,                         &
        LossG2100,LossG2103,LossG2104,                      &
        LossG2105,LossG2106,LossG2107,                      &
        LossG2109,LossG2110,LossG2111,                      &
        LossG2112,LossG3200,LossG3201,                      &
        LossG3202,LossG3203,LossG3207,                      &
        LossG4101,LossG4110,LossG6203,                      &
        LossG6204,LossG7201,LossJ3200,                      &
        LossJ3201,LossJ4101b,                               &
        OHlossNMHC,HO2lossNMHC,                             &
        OHlossHO2prodNMHC,HO2prodNMHCNOy,                   &
        HO2prodNMHCphoto,frac_o3,                           &
        po3_ho2, po3_ro2, do3_ro,                           &
        do3_oh, do3_ho2, do3_xo, do3_no,                    &
        status)



! Purpose:  Calculate impact of emissions on short-lived species 
!           like OH, HO2, RO2 and NO.   
! Detailed: 
!   i=tagged category i
! NO 
! - NOi = NO * NOyi/NOy    i.e. NOi/NO = NOyi/NOy  
!         this approach is valid for 'large' (whatever that is ateast T42-T30)
!         grid boxes since a local instanteneous emission enhances the actual 
!         grid box concentration only to some extend.
!         Also the exchange between 2 grid boxes is limited compared to
!         mass of the grid box. See Grewe (2004) box model study
!
! RO
! - ROi = RO *NMHCi/NMHC    as NO reasonable ?
!
!
! HOx tagging
! theory is described in Grewe et al. 2017 GMD and Rieger et al. 2017 GMD
!
! - considers the following reactions:
!   reaction rates:    reactions:
!   LossG2100          H + O2 -> HO2
!   LossG2103          OH + O3P -> H + O2
!   LossG2104          OH + O3 -> HO2 + O2
!   LossG2105          OH + H2 -> H2O + H
!   LossG2106          OH + H2 -> H2O + H
!   LossG2107          HO2 + O3 -> OH + 2 O2
!   LossG2109          HO2 + OH -> H2O + O2
!   LossG2110          HO2 + HO2 -> H2O2 + O2
!   LossG2111          H2O + O1D -> 2 OH
!   LossG2112          H2O2 + OH -> H2O + HO2
!   LossG3200          NO + OH -> HONO
!   LossG3201          NO + HO2 -> NO2 + OH
!   LossG3202          NO2 + OH -> HNO3
!   LossG3203          NO2 + HO2 -> HNO4
!   LossG3207          HNO4 -> NO2 + HO2
!   LossG4101          CH4 + OH -> CH3O2 + H2O
!   LossG6203          ClO + OH -> .94 Cl + .94 HO2 + .06 HCl + .06 O2
!   LossG6204          ClO + HO2 -> HOCl + O2
!   LossG7201          BrO + HO2 -> HOBr + O2
!   LossJ2101          H2O2 + hv -> 2 OH  (neglected here!)
!   LossJ3200          HONO + hv -> NO + OH
!   LossJ3201          HNO3 + hv -> NO2 + OH
!   LossJ4101b         HCHO + hv -> H + CO + HO2
!   LossJ6201          HOCl + hv -> Cl + OH  (neglected here!)
!   LossJ7200          HOBr + hv -> Br + OH  (neglected here!)
!   OHlossNMHC         NMHC + OH -> NMHC
!   HO2lossNMHC        NMHC + HO2 -> NMHC
!   OHlossHO2prodNMHC  NMHC + OH -> NMHC + HO2
!   HO2prodNMHCNOy     NMHC + NOy -> HO2 + NMHC + NOy
!   HO2prodNMHCphoto   NMHC + hv -> NMHC + HO2
!
!   caution: reaction rate considers only how often the species undergo
!            this reaction and does not consider the quantity of produced 
!            or destroyed HOx molecules
!
!
! - Steady-State assumption 
!   See PhD Tsati and Grewe et al., 2014 (REACT4C V1.0)
! - Poh und DCH4 DNMHC in mol/mol/s
!
!
! Integration
! 
!
! Structure
!   1. Set variables global -> local
!   2. HOx Tagging
!      - Calculate steady state OH and HO2
!      - set global OH and HO2 variables for output
!      - choose OH HO2 acording to ispecies settings 
!   3. Advected species
!      - Calculated P and L 
!      - Integrate according to i_integrate
!      - set tendencies
!
! Author: Volker Grewe, DLR-Oberpfaffenhofen, August 2008 during NCAR visit, Boulder
!         Eleni Tsati, DLR-Oberpfaffenhofen, August 2013.
!         Vanessa Rieger, DLR-Oberpfaffenhofen, March 2017


   IMPLICIT NONE
   
   ! Variables I/O

   INTEGER,                               INTENT(IN)  :: kproma, nlev
   REAL(DP),                              INTENT(IN)  :: dt         

   REAL(DP), DIMENSION(:,:,:),            INTENT(IN)  :: o3_tagging,noy_tagging,pan_tagging, &
                                                         nmhc_tagging,co_tagging

   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: nmhc_tagging_sum,co_tagging_sum, &
                                                         noy_tagging_sum,o3_tagging_sum,pan_tagging_sum, &
                                                         oh, ho2 

   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: o3prod_ro2,o3prod_o2,o3prod_ho2 
   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: o3loss_oh,o3loss_ho2,o3loss_ro
   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: o3loss_no,o3loss_xo,o3loss
   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: panprod,panloss

   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: coloss,coprod    

   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: ch4loss  

   REAL(DP), DIMENSION(:,:),              INTENT(IN)  :: LossG2100,LossG2103,LossG2104,     &
                                                         LossG2105,LossG2106,LossG2107,     &
                                                         LossG2109,LossG2110,LossG2111,     &
                                                         LossG2112,LossG3200,LossG3201,     &
                                                         LossG3202,LossG3203,LossG3207,     &
                                                         LossG4101,LossG4110,LossG6203,     &
                                                         LossG6204,LossG7201,LossJ3200,     &
                                                         LossJ3201,LossJ4101b,              &
                                                         OHlossNMHC,HO2lossNMHC,            &
                                                         OHlossHO2prodNMHC,HO2prodNMHCNOy,  &
                                                         HO2prodNMHCphoto


   REAL(DP), DIMENSION(:,:,:),            INTENT(OUT) :: ch4loss_tag    ! tagged CH4 loss (needed for change
                                                                        ! of CH4 lifetime)
  
   REAL(DP), DIMENSION(:,:,:),            INTENT(OUT) :: do3,dnoy,dco,dnmhc,dpan  ! tendencies [mol/mol/s]
   INTEGER,                               INTENT(OUT) :: status ! 0=ok 
                                                                ! 1=i_species error
                                                                ! 2=i_method error
                                                                ! 3=i_species-i?method conflict


   REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)  :: Do3_oh, Do3_ho2, Do3_no, Do3_ro, &
                                                 Do3_xo, PO3_ho2, PO3_ro2
   REAL(DP), DIMENSION(:,:,:), INTENT(INOUT)  :: OH_tag,HO2_tag

   REAL(DP),  DIMENSION(kproma,nlev,N_tag_species), INTENT(OUT) :: frac_o3


  ! Local variables

   CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_chemistry'

   ! Define limits for contribution calculation. Necessary, e.g. during night, when OH is absent.                                          
   REAL(DP), PARAMETER                             :: limit_chem = 1.e-15_dp, limit_pan = 1.e-18_dp
        
   REAL(DP), DIMENSION(kproma,nlev,N_tag_species)  ::  o3_tag2,co_tag2,nmhc_tag2, &
                                                       noy_tag2,pan_tag2
                                                      
   REAL(DP),  DIMENSION(kproma,nlev)               :: o3loss_ua
   REAL(DP),  DIMENSION(kproma,nlev,N_tag_species) :: P_O3,D_O3,P_noy,D_noy,P_pan, & 
                                                      D_pan,P_nmhc,D_nmhc,         &
                                                      P_co,D_co,frac_ho2, frac_oh, &
                                                      Numerator_ho2, Numerator_oh

   REAL(DP),  DIMENSION(kproma,nlev)               :: P_OH, L_OH, P_HO2, L_HO2
   REAL(DP),  DIMENSION(kproma,nlev)               :: redOHprod, redOHloss,  &
                                                      redHO2prod,redHO2loss
   REAL(DP),  DIMENSION(kproma,nlev)               :: resOH, resHO2, resH           ! correction terms

   REAL(DP),  DIMENSION(kproma,nlev)               :: o3_tag_sum, co_tag_sum, pan_tag_sum,   & 
                                                      noy_tag_sum, nmhc_tag_sum, OH_tag_sum, HO2_tag_sum

   REAL(DP),  DIMENSION(kproma,nlev,N_tag_species) :: B_tag, A_tag 

   REAL(DP),  DIMENSION(kproma,nlev,N_tag_species) :: frac_noy, frac_pan, frac_nmhc, frac_co

   REAL(DP),  DIMENSION(kproma,nlev)               ::  Par_fra_nsum

   INTEGER                                         :: i,j,k

   INTRINSIC                                       :: SHAPE, MAXVAL, MINVAL



   !write(*,*) 'vgvg IN CHEM: Bounds o3_tag_in', SHAPE(o3_tag_in)
   !write(*,*) 'vgvg IN CHEM: kproma,nlev,dt',kproma,nlev,dt
   !write(*,*) 'src_kproma,nlev, time step, o3_tagging',kproma, nlev,dt, o3_tagging(1,1,1)

   !PRINT*, "In tagging chem"


   status=0 ! no error

   do3(1:kproma,:,:)=0.0_dp
   dnoy(1:kproma,:,:)=0.0_dp
   dco(1:kproma,:,:)=0.0_dp
   dnmhc(1:kproma,:,:)=0.0_dp
   dpan(1:kproma,:,:)=0.0_dp



!--------------------------------------
! 1. Set variables
!------

    ! Set variables, calculate according to system and feed chem. tendencies  back
    o3_tag2(1:kproma,:,:)    = o3_tagging(1:kproma,:,:)
    noy_tag2(1:kproma,:,:)   = noy_tagging(1:kproma,:,:)
    pan_tag2(1:kproma,:,:)   = pan_tagging(1:kproma,:,:)
    nmhc_tag2(1:kproma,:,:)  = nmhc_tagging(1:kproma,:,:)
    co_tag2(1:kproma,:,:)    = co_tagging(1:kproma,:,:)
    
    o3_tag_sum(1:kproma,:)    = o3_tagging_sum(1:kproma,:)
    noy_tag_sum(1:kproma,:)   = noy_tagging_sum(1:kproma,:)
    pan_tag_sum(1:kproma,:)   = pan_tagging_sum(1:kproma,:)
    nmhc_tag_sum(1:kproma,:)  = nmhc_tagging_sum(1:kproma,:)
    co_tag_sum(1:kproma,:)    = co_tagging_sum(1:kproma,:)
    
    OH_tag_sum(1:kproma,:)  = oh(1:kproma,:)
    HO2_tag_sum(1:kproma,:)  = ho2(1:kproma,:)



   ! unaffected o3 loss
   SELECT CASE (i_species)
   CASE (0)
      co_tag2(1:kproma,:,:)=0._dp
      nmhc_tag2(1:kproma,:,:)=0._dp
      pan_tag2(1:kproma,:,:)=0._dp
      o3loss_ua(1:kproma,:)=o3loss(1:kproma,:)     ! ozone loss unaffected
   CASE (1)
      co_tag2(1:kproma,:,:)=0._dp
      nmhc_tag2(1:kproma,:,:)=0._dp
      pan_tag2(1:kproma,:,:)=0._dp
      o3loss_ua(1:kproma,:)=o3loss(1:kproma,:)-o3loss_no(1:kproma,:)
   CASE (2) 
      o3loss_ua(1:kproma,:)=o3loss_xo(1:kproma,:)   & 
                     +o3loss_oh(1:kproma,:)   & 
                     +o3loss_ho2(1:kproma,:)
   CASE (3) 
      o3loss_ua(1:kproma,:)=o3loss_xo(1:kproma,:)
   CASE DEFAULT
      status=1
   END SELECT

   
!------------------------------------------------------





!------------------------------------------------------
! 2. HOx Tagging
!------

   ! determine total production and loss rates of reduced HOx reaction system
   redOHprod(1:kproma,:) =   LossG2106(1:kproma,:)         &
                           + LossG2107(1:kproma,:)         &
                           + 2._dp * LossG2111(1:kproma,:) &
                           + LossG3201(1:kproma,:)         &
                           + LossJ3200(1:kproma,:)         &
                           + LossJ3201(1:kproma,:)         

   redOHloss(1:kproma,:) =   LossG2103(1:kproma,:)         &
                           + LossG2104(1:kproma,:)         &
                           + LossG2105(1:kproma,:)         &
                           + LossG2109(1:kproma,:)         &
                           + LossG2112(1:kproma,:)         &
                           + LossG4110(1:kproma,:)         &
                           + LossG4101(1:kproma,:)         &
                           + LossG3200(1:kproma,:)         &
                           + LossG3202(1:kproma,:)         &
                           + LossG6203(1:kproma,:)         &
                           + OHlossNMHC(1:kproma,:)        &
                           + OHlossHO2prodNMHC(1:kproma,:)

   redHO2prod(1:kproma,:) =  LossG2100(1:kproma,:)            &
                            + LossG2104(1:kproma,:)           &
                            + LossG2112(1:kproma,:)           &
                            + LossG3207(1:kproma,:)           &
                            + 0.94_dp * LossG6203(1:kproma,:) &
                            + HO2prodNMHCNOy(1:kproma,:)      &
                            + OHlossHO2prodNMHC(1:kproma,:)   &
                            + HO2prodNMHCphoto(1:kproma,:)    

   redHO2loss(1:kproma,:) =   LossG2106(1:kproma,:)         &
                            + LossG2107(1:kproma,:)         &
                            + LossG2109(1:kproma,:)         &
                            + 2._dp * LossG2110(1:kproma,:) &
                            + LossG3201(1:kproma,:)         &
                            + LossG3203(1:kproma,:)         & 
                            + LossG6204(1:kproma,:)         &
                            + LossG7201(1:kproma,:)         &
                            + HO2lossNMHC(1:kproma,:)

   ! correction term: residual of production and loss rates of reduced HOx reaction system
   resOH(1:kproma,:)  = redOHloss(1:kproma,:)  - redOHprod(1:kproma,:)
   resHO2(1:kproma,:) = redHO2loss(1:kproma,:) - redHO2prod(1:kproma,:)

   ! correction term for H
   resH(1:kproma,:) =   LossG2100(1:kproma,:)   &
                      - LossG2103(1:kproma,:)  &
                      - LossG2105(1:kproma,:)  &
                      - LossG4110(1:kproma,:)  &
                      - LossJ4101b(1:kproma,:)


   ! P_OH, L_OH, P_HO2, L_HO2 after Rieger et at. 2017
   P_OH(1:kproma,:) =   LossG2106(1:kproma,:)         &
                      + LossG2107(1:kproma,:)         &
                      + LossG3201(1:kproma,:)         &
                      - LossG2109(1:kproma,:)        

   L_OH(1:kproma,:) =   LossG2103(1:kproma,:)         &
                      + LossG2104(1:kproma,:)         &
                      + 2._dp * LossG2105(1:kproma,:) &
                      + LossG2109(1:kproma,:)         &
                      + 2._dp * LossG2112(1:kproma,:) &
                      + LossG3200(1:kproma,:)         &
                      + LossG3202(1:kproma,:)         &
                      + LossG4110(1:kproma,:)         &
                      + 2._dp * LossG4101(1:kproma,:) &
                      + 2._dp * LossG6203(1:kproma,:) &
                      + OHlossNMHC(1:kproma,:)        &
                      + OHlossHO2prodNMHC(1:kproma,:)

   P_HO2(1:kproma,:) =   LossG2103(1:kproma,:)           &
                       + LossG2104(1:kproma,:)           &
                       + 2._dp * LossG2105(1:kproma,:)   & 
                       + 2._dp * LossG2112(1:kproma,:)   &
                       + LossG4110(1:kproma,:)           & 
                       + 1.88_dp * LossG6203(1:kproma,:) &  ! 2*0.94
                       + OHlossHO2prodNMHC(1:kproma,:)   &
                       - LossG2109(1:kproma,:)

   L_HO2(1:kproma,:) =   LossG2106(1:kproma,:)         &
                       + LossG2107(1:kproma,:)         &
                       + LossG2109(1:kproma,:)         &
                       + 4._dp * LossG2110(1:kproma,:) & 
                       + LossG3201(1:kproma,:)         &
                       + LossG3203(1:kproma,:)         &
                       + 2._dp * LossG6204(1:kproma,:) &
                       + 2._dp * LossG7201(1:kproma,:) &
                       + HO2lossNMHC(1:kproma,:)     


   ! denominator (Grewe et at. 2017, equation (27) and (28))
   Par_fra_nsum(1:kproma,:) =   L_OH(1:kproma,:) * L_HO2(1:kproma,:) &
                              - P_OH(1:kproma,:) * P_HO2(1:kproma,:)




   ! loop over sectors
   do i=1,N_tag_species



    ! Calculate fractions: Tagged_Specie/Total
    do  k= 1, nlev
     do  j = 1, kproma

         !  CO
         if ( abs(co_tag_sum(j,k))  .gt.  limit_chem ) then    
             frac_co(j,k,i)= CO_tag2(j,k,i)/CO_tag_sum(j,k)
         else
             frac_co(j,k,i) = 0._dp
         endif


        ! PAN
         if (  abs(pan_tag_sum(j,k)) .gt.  limit_chem  ) then
              frac_pan(j,k,i)= pan_tag2(j,k,i)/pan_tag_sum(j,k)
         else 
              frac_pan(j,k,i) = 0._dp
         endif

         if ( abs(frac_pan(j,k,i)) .lt. limit_pan) then
              frac_pan(j,k,i)=sign(limit_pan, frac_pan(j,k,i))
         endif


        ! NMHC
         if ( abs(nmhc_tag_sum(j,k))  .gt.  limit_chem ) then
              frac_nmhc(j,k,i)= NMHC_tag2(j,k,i)/NMHC_tag_sum(j,k)
         else
              frac_nmhc(j,k,i)= 0._dp
         endif


         ! NOy
         if ( abs(NOy_tag_sum(j,k)) .gt. limit_chem ) then
            frac_noy(j,k,i)= NOy_tag2(j,k,i)/NOy_tag_sum(j,k)
         else 
            frac_noy(j,k,i) = 0._dp
         endif


         ! O3
         if ( abs(O3_tag_sum(j,k)) .gt. limit_chem ) then
            frac_o3(j,k,i)= O3_tag2(j,k,i)/O3_tag_sum(j,k)
         else 
            frac_o3(j,k,i) = 0._dp
         endif

     enddo
    enddo



   ! A_tag, B_tag after Rieger et al. 2017
    A_tag(1:kproma,:,i) =   LossG2106(1:kproma,:) * frac_o3(1:kproma,:,i)           &
                          + LossG2107(1:kproma,:) * frac_o3(1:kproma,:,i)           &
                          + 4._dp * LossG2111(1:kproma,:) * frac_o3(1:kproma,:,i)   &
                          + LossG3201(1:kproma,:) * frac_noy(1:kproma,:,i)          &
                          + 2._dp * LossJ3200(1:kproma,:) * frac_noy(1:kproma,:,i)  &
                          + 2._dp * LossJ3201(1:kproma,:) * frac_noy(1:kproma,:,i)  &
                          - LossG2103(1:kproma,:) * frac_o3(1:kproma,:,i)           &
                          - LossG2104(1:kproma,:) * frac_o3(1:kproma,:,i)           &
                          - LossG3200(1:kproma,:) * frac_noy(1:kproma,:,i)          &
                          - LossG3202(1:kproma,:) * frac_noy(1:kproma,:,i)          &
                          - LossG4110(1:kproma,:) * frac_co(1:kproma,:,i)           &
                          - OHlossNMHC(1:kproma,:)* frac_nmhc(1:kproma,:,i)         &
                          - OHlossHO2prodNMHC(1:kproma,:) * frac_nmhc(1:kproma,:,i) &
                          + 2._dp * resOH(1:kproma,:)/DBLE(N_tag_species) 

    B_tag(1:kproma,:,i) =   LossG2103(1:kproma,:) * frac_o3(1:kproma,:,i)                   &
                          + LossG2104(1:kproma,:) * frac_o3(1:kproma,:,i)                   &
                          + 2._dp * LossG3207(1:kproma,:) * frac_noy(1:kproma,:,i)          &
                          + LossG4110(1:kproma,:) * frac_co(1:kproma,:,i)                   & 
                          + 2._dp * LossJ4101b(1:kproma,:) * frac_nmhc(1:kproma,:,i)        &
                          + HO2prodNMHCNOy(1:kproma,:) * frac_nmhc(1:kproma,:,i)            &
                          + HO2prodNMHCNOy(1:kproma,:) * frac_noy(1:kproma,:,i)             &
                          + OHlossHO2prodNMHC(1:kproma,:) * frac_nmhc(1:kproma,:,i)         &
                          + 2._dp * HO2prodNMHCphoto(1:kproma,:) * frac_nmhc(1:kproma,:,i)  &
                          - LossG2106(1:kproma,:) * frac_o3(1:kproma,:,i)                   &
                          - LossG2107(1:kproma,:) * frac_o3(1:kproma,:,i)                   &
                          - LossG3201(1:kproma,:) * frac_noy(1:kproma,:,i)                  &
                          - LossG3203(1:kproma,:) * frac_noy(1:kproma,:,i)                  &
                          - HO2lossNMHC(1:kproma,:) * frac_nmhc(1:kproma,:,i)               &
                          + 2._dp * resH(1:kproma,:)/DBLE(N_tag_species)                    &
                          + 2._dp * resHO2(1:kproma,:)/DBLE(N_tag_species) 


   ! numerator (Grewe et at. 2017, equation (27) and (28))
    Numerator_oh(1:kproma,:,i) =   A_tag(1:kproma,:,i) * L_HO2(1:kproma,:) &
                                 + B_tag(1:kproma,:,i) * P_OH(1:kproma,:)

    Numerator_ho2(1:kproma,:,i) =   A_tag(1:kproma,:,i) * P_HO2(1:kproma,:) &
                                  + B_tag(1:kproma,:,i) * L_OH(1:kproma,:)


    ! Calculated OH and HO2 fractions, special care for singularities 
    do k=1,nlev
      do j=1,kproma

        if ( abs(Par_fra_nsum(j,k)) .gt. 1.e-45_dp ) then
          frac_oh(j,k,i)  = Numerator_oh(j,k,i)/Par_fra_nsum(j,k) 
          frac_ho2(j,k,i) = Numerator_ho2(j,k,i)/Par_fra_nsum(j,k) 
        else  
          frac_oh(j,k,i)  = 0._dp
          frac_ho2(j,k,i) = 0._dp
        end if

      enddo
    enddo


    ! calculate OH and HO2 contributions
    OH_tag(1:kproma,:,i)  = OH_tag_sum(1:kproma,:)  * frac_oh(1:kproma,:,i)
    HO2_tag(1:kproma,:,i) = HO2_tag_sum(1:kproma,:) * frac_ho2(1:kproma,:,i)

    
    ! tagged CH4 loss (needed for calculation of CH4 lifetime change)
    ch4loss_tag(1:kproma,:,i) = frac_oh(1:kproma,:,i) * ch4loss(1:kproma,:)

   enddo  !i_tag_species



  !write (*,*) 'src_select case, i_sp' , i_species
  
   SELECT CASE (i_species)
     CASE (0,1,2)
        oh_tag(1:kproma,:,:)=0._dp
       ho2_tag(1:kproma,:,:)=0._dp
     CASE (3)
      !  take calculated oh ho2 values
    CASE DEFAULT
       status=1
   END SELECT

   
   
!------------------------------------------------------
! 3. Determine prod and loss terms for advected species
!------

  SELECT CASE (i_species)
   CASE(0)     ! NOy O3 and Prod only
     do i=1,N_tag_species
        P_o3(1:kproma,:,i)=(o3prod_ro2(1:kproma,:)+o3prod_ho2(1:kproma,:))   & 
                  *(pan_tag2(1:kproma,:,i)+noy_tag2(1:kproma,:,i))/(pan_tag_sum(1:kproma,:)+noy_tag_sum(1:kproma,:))

        D_o3(1:kproma,:,i)=o3loss(1:kproma,:)*o3_tag2(1:kproma,:,i)/o3_tag_sum(1:kproma,:)
        P_pan(1:kproma,:,i)=0._dp                ! Not emitted no chemistry
        D_pan(1:kproma,:,i)=0._dp
        P_nmhc(1:kproma,:,i)=0._dp               ! NMHC might be emitted but does not affect chemistry
        D_nmhc(1:kproma,:,i)=coprod(1:kproma,:)*nmhc_tag2(1:kproma,:,i) /nmhc_tag_sum(1:kproma,:)
        P_noy(1:kproma,:,i)=0._dp                ! NOy only affected by wash-out
        D_noy(1:kproma,:,i)=0._dp
        P_co(1:kproma,:,i)=coprod(1:kproma,:)*nmhc_tag2(1:kproma,:,i)/nmhc_tag_sum(1:kproma,:)  ! 
        D_co(1:kproma,:,i)=coloss(1:kproma,:)*co_tag2(1:kproma,:,i)/co_tag_sum(1:kproma,:)
     enddo
     P_nmhc(1:kproma,:,ich4)=P_nmhc(1:kproma,:,ich4)+ch4loss(1:kproma,:)       !
     P_o3(1:kproma,:,istr)=P_o3(1:kproma,:,istr)+o3prod_o2(1:kproma,:)
   CASE(1)     ! NOy O3 
     do i=1,N_tag_species

        P_o3(1:kproma,:,i)=(o3prod_ro2(1:kproma,:)+o3prod_ho2(1:kproma,:))*(pan_tag2(1:kproma,:,i)+noy_tag2(1:kproma,:,i)) &
                                        /(pan_tag_sum(1:kproma,:)+noy_tag_sum(1:kproma,:)) 

        D_o3(1:kproma,:,i)=o3loss_no(1:kproma,:)*(pan_tag2(1:kproma,:,i)+noy_tag2(1:kproma,:,i))&
                    /(pan_tag_sum(1:kproma,:)+noy_tag_sum(1:kproma,:))   &
                   +o3loss_ua(1:kproma,:)*o3_tag2(1:kproma,:,i)/o3_tag_sum(1:kproma,:)


        P_pan(1:kproma,:,i)=0._dp                ! Not emitted no chemistry
        D_pan(1:kproma,:,i)=0._dp
        P_nmhc(1:kproma,:,i)=0._dp               ! NMHC might be emitted but does not affect chemistry



        D_nmhc(1:kproma,:,i)=coprod(1:kproma,:)*nmhc_tagging(1:kproma,:,i) / nmhc_tagging_sum(1:kproma,:)
        P_noy(1:kproma,:,i)=0._dp                ! NOy only affected by wash-out
        D_noy(1:kproma,:,i)=0._dp
        P_co(1:kproma,:,i)=coprod(1:kproma,:)*nmhc_tag2(1:kproma,:,i)/nmhc_tag_sum(1:kproma,:)  ! CH4->CO misssing for i=ich4
        D_co(1:kproma,:,i)=coloss(1:kproma,:)*co_tag2(1:kproma,:,i)/co_tag_sum(1:kproma,:)
     enddo
     P_nmhc(1:kproma,:,ich4)=P_nmhc(1:kproma,:,ich4)+ch4loss(1:kproma,:) 
     P_o3(1:kproma,:,istr)=P_o3(1:kproma,:,istr)+o3prod_o2(1:kproma,:)
   CASE(2) ! = CASE 2 or 3, since 2 = 3 with oh, ho2=0.    nevertheless do it separatey
   


   CASE(3)     ! i_species   = full method

        do i=1,N_tag_species     


! eq. (17) from Grewe et al. (2010) with Pxy =  o3prod_ho2(1:kproma,:) / ho2(1:kproma,:)/ noy(1:kproma,:)

         P_o3(1:kproma,:,i)  =  .5_dp*o3prod_ho2(1:kproma,:)* &
                                ( frac_ho2(1:kproma,:,i) + frac_noy(1:kproma,:,i)) &

                                + .5_dp*o3prod_ro2(1:kproma,:)* &   
                               (frac_nmhc(1:kproma,:,i)+ &
                                frac_noy(1:kproma,:,i)) 

      ! write (*,*) 'P_O3', MINVAL(P_o3(1:kproma,:,1:N_tag_species)),MAXVAL(P_o3(1:kproma,:,1:N_tag_species))
!------------------------------------------------------------------------------------------------------

        PO3_ho2(1:kproma,:,i) =  .5_dp*o3prod_ho2(1:kproma,:)*(frac_ho2(1:kproma,:,i) + frac_noy(1:kproma,:,i)) 

      ! write (*,*) 'P_O3_ho2', MINVAL(PO3_ho2(1:kproma,:,1:N_tag_species)),MAXVAL(PO3_ho2(1:kproma,:,1:N_tag_species))

        PO3_ro2(1:kproma,:,i) = .5_dp*o3prod_ro2(1:kproma,:)*(frac_nmhc(1:kproma,:,i)+ frac_noy(1:kproma,:,i)) 

     !   write (*,*) 'P_O3_ro2', MINVAL(PO3_ro2(1:kproma,:,1:N_tag_species)),MAXVAL(PO3_ro2(1:kproma,:,1:N_tag_species))


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         D_o3(1:kproma,:,i)  =  .5_dp*o3loss_oh(1:kproma,:)* &
                                ( frac_oh(1:kproma,:,i)+ frac_o3(1:kproma,:,i) ) + &
                                                                                                             
                                .5_dp*o3loss_ho2(1:kproma,:)*  &
                                ( frac_ho2(1:kproma,:,i) + frac_o3(1:kproma,:,i) ) + &
                                                                                                             
                                .5_dp*o3loss_no(1:kproma,:)* &
                               ( frac_noy(1:kproma,:,i) + &
                                frac_o3(1:kproma,:,i)) + &
                                                                                                              
                                .5_dp*o3loss_ro(1:kproma,:)* &
                               (frac_nmhc(1:kproma,:,i) + &
                                frac_o3(1:kproma,:,i)) + &
                                                                                                           
                                o3loss_xo(1:kproma,:) * frac_o3(1:kproma,:,i)


        !   write (*,*) 'D_O3', MINVAL(D_o3(1:kproma,:,1:N_tag_species)),MAXVAL(D_o3(1:kproma,:,1:N_tag_species))
!------------------------------------------------------------------------------------------------------

         Do3_oh(1:kproma,:,i) =  .5_dp*o3loss_oh(1:kproma,:)* &
                                ( frac_oh(1:kproma,:,i) + frac_o3(1:kproma,:,i) ) 

          !  write (*,*) 'D_O3_oh', MINVAL(Do3_oh(1:kproma,:,1:N_tag_species)),MAXVAL(Do3_oh(1:kproma,:,1:N_tag_species))


         Do3_ho2(1:kproma,:,i) = .5_dp*o3loss_ho2(1:kproma,:)*  &
                                ( frac_ho2(1:kproma,:,i) + frac_o3(1:kproma,:,i) ) 
      
         !  write (*,*) 'D_O3_ho2', MINVAL(Do3_ho2(1:kproma,:,1:N_tag_species)),MAXVAL(Do3_ho2(1:kproma,:,1:N_tag_species))


         Do3_no(1:kproma,:,i) = .5_dp*o3loss_no(1:kproma,:)* &
                                   ( frac_noy(1:kproma,:,i) + &
                                    frac_o3(1:kproma,:,i)) 


       !    write (*,*) 'D_O3_no', MINVAL(Do3_no(1:kproma,:,1:N_tag_species)),MAXVAL(Do3_no(1:kproma,:,1:N_tag_species))


         Do3_ro(1:kproma,:,i) = .5_dp*o3loss_ro(1:kproma,:)* &
                                  ( frac_nmhc(1:kproma,:,i) + &
                                    frac_o3(1:kproma,:,i)) 

       !   write (*,*) 'D_O3_ro', MINVAL(Do3_ro(1:kproma,:,1:N_tag_species)),MAXVAL(Do3_ro(1:kproma,:,1:N_tag_species))
    

        Do3_xo(1:kproma,:,i)  =  o3loss_xo(1:kproma,:) * frac_o3(1:kproma,:,i)

        !   write (*,*) 'D_o3_xo', MINVAL(Do3_xo(1:kproma,:,1:N_tag_species)),MAXVAL(Do3_xo(1:kproma,:,1:N_tag_species))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !write (*,*) 'panprod_src', MINVAL(panprod(1:kproma,:)), MAXVAL(panprod(1:kproma,:))
       !write (*,*) 'panloss_src', MINVAL(panloss(1:kproma,:)), MAXVAL(panloss(1:kproma,:))


          P_pan(1:kproma,:,i)  =  .5_dp*panprod(1:kproma,:)*(frac_nmhc(1:kproma,:,i) + frac_noy(1:kproma,:,i))              
          D_pan(1:kproma,:,i)  =  panloss(1:kproma,:)*frac_pan(1:kproma,:,i)

          P_co(1:kproma,:,i)   =   coprod(1:kproma,:)*frac_nmhc(1:kproma,:,i)
          D_co(1:kproma,:,i)   =   coloss(1:kproma,:)*frac_co(1:kproma,:,i)
          
          P_noy(1:kproma,:,i)  =   D_pan(1:kproma,:,i)
          D_noy(1:kproma,:,i)  =   P_pan(1:kproma,:,i)
          
          P_nmhc(1:kproma,:,i) =   2._dp*D_pan(1:kproma,:,i)
          D_nmhc(1:kproma,:,i) =   2._dp*P_pan(1:kproma,:,i) + coprod(1:kproma,:)*frac_nmhc(1:kproma,:,i)
       enddo


         P_nmhc(1:kproma,:,ich4)=P_nmhc(1:kproma,:,ich4)+ch4loss(1:kproma,:) 
         P_o3(1:kproma,:,istr)=P_o3(1:kproma,:,istr)+o3prod_o2(1:kproma,:)

         !write (*,*) 'P_pan_src', MINVAL(P_pan(1:kproma,:,i)), MAXVAL(P_pan(1:kproma,:,i))
         !write (*,*) 'D_pan_src', MINVAL(D_pan(1:kproma,:,i)), MAXVAL(D_pan(1:kproma,:,i))



! WRITE(*,*) 'vgvg in chemistry: P_O3, D_O3', P_o3(4,5,1),D_o3(4,5,1)
! WRITE(*,*) 'vgvg in chemistry: OH,HO2',OH_tag(4,5,1),HO2_tag(4,5,1)
! WRITE(*,*) 'vgvg in chemistry: Pco,Dco,co,',P_CO(4,5,1),D_CO(4,5,1),CO_tag(4,5,1)
    
  CASE DEFAULT
    status=1
  END SELECT    
   
!------
! Integration of  tagged species 
!     Problem for species with tag=0. and loss=0.
!     P and D given in mol/mol/s
!       dx=P-Dx  ->  xneu=(xalt+P*dt)/(1+D/xneu*dt)
!       (xneu-xalt)/dt=(xalt+P*dt-xalt-D/xneu*xalt*dt)/(1+D/xneu*dt)/dt
!                     = P-Dx   !!
!------
   
   do i=1,N_tag_species
      do3(1:kproma,:,i)  =P_o3(1:kproma,:,i)  -D_o3(1:kproma,:,i)
      dco(1:kproma,:,i)  =P_co(1:kproma,:,i)  -D_co(1:kproma,:,i)
      dnoy(1:kproma,:,i) =P_noy(1:kproma,:,i) -D_noy(1:kproma,:,i)
      dnmhc(1:kproma,:,i)=P_nmhc(1:kproma,:,i)-D_nmhc(1:kproma,:,i)
      ! dpan(1:kproma,:,i) = -dnoy(1:kproma,:,i)   !=P_pan(1:kproma,:,i)-D_pan(1:kproma,:,i)=D_noy(1:kproma,:,i)-P_noy(1:kproma,:,i)
      dpan(1:kproma,:,i) = P_pan(1:kproma,:,i)-D_pan(1:kproma,:,i)

   enddo 



  END SUBROUTINE tagging_chemistry    

! ----------------------------------------------------------------------

  SUBROUTINE tagging_read_nml(status, iou)

    ! READ TAGGING NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER(I4), INTENT(OUT) :: status ! error status
    INTEGER(I4), INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/  i_integrate,  &       ! 1: Euler backward   2: euler forward
                                           !  besser ueber CPL tp_4_o3_orig,   &     ! 1: WMO; 2: PV
                     i_tracer_init, &      ! How to initialize  0: uniformly  1: via File       
                     i_species, &          ! 0: NOy+O3 prod only
                                           ! 1: NOy+O3
                                           ! 2: NOy, CO, PAN, NMHC, O3
                                           ! 3: NOy, CO, PAN, NMHC, O3, HO2 and OH
                     i_advect_scaling, &   ! fix numerical diffusion caused by advection
                     l_adv_err_diag, &     ! diagnose impact of this diffusion
                     i_tagging_ctrl, &
                     category_in !op_mm_20171610+
                     


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tagging_read_nml'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER(I4)                 :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! 
      SELECT CASE(i_tracer_init)
       CASE (0)
         WRITE(*,*) '... Tagging tracers initialized from rerun files'
       CASE (1) 
         WRITE(*,*) '... Tagging tracers initialized uniformly for trop and strat' 
         WRITE(*,*) '...     only for new start'
       CASE (-1)  
         WRITE(*,*) '... Tagging tracers initialized uniformly for trop and strat' 
         WRITE(*,*) '...    Data on re-run files ignored' 
       CASE (2)
         WRITE(*,*) '... Tagging tracers initialized from file'
         WRITE(*,*) '...     only for new start'
       CASE (-2) 
         WRITE(*,*) '... Tagging tracers initialized from file'
         WRITE(*,*) '...    Data on re-run files ignored' 
       CASE DEFAULT
         WRITE(*,*) '... i_tracer_init out of range  ',i_tracer_init
         RETURN
      END SELECT 



      SELECT CASE(i_species)
       CASE(0)     ! NOy and O3 only and prod only -> Grewe 2004
         WRITE(*,*) '... use only NOy and O3 tagged tracers and Prod terms considered only'
       CASE(1)     ! NOy and O3 only
         WRITE(*,*) '... use only NOy and O3 tagged tracers '
       CASE(2)     ! medium
         WRITE(*,*) '... use NOy, PAN, NMHC, CO  and O3 as tagged tracers'
       CASE(3)     ! full
         WRITE(*,*) '... use NOy, PAN, NMHC, CO  and O3 as tagged tracers, and HO2 and OH'
         WRITE(*,*) '                                                    as steady state' 
       CASE DEFAULT
         WRITE(*,*) '... Incorrect value for i_species  ',i_species
         RETURN
       END SELECT


      SELECT CASE(i_advect_scaling)
       CASE (0)
         WRITE(*,*) '... No scaling of tagged tracers to sum to 100% of chemical tracer'
         WRITE(*,*) '... which could arise from numerical diffusion of the advection scheme'
       CASE (1)
         WRITE(*,*) '... Scaling of tagged tracers to sum to 100% of chemical tracer enabled'
         WRITE(*,*) '... Method: Linear Error analysis still to be done' 
       CASE DEFAULT 
         WRITE(*,*) '... i_advect_scaling is of range   ',i_advect_scaling
         RETURN
       END SELECT
     
     IF (l_adv_err_diag) then
         WRITE(*,*) '... Diagnosis of numerical advection error switched on'
     ELSE 
         WRITE(*,*) '... Diagnosis of numerical advection error switched off'
     ENDIF
    !

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

    END SUBROUTINE tagging_read_nml

    
    
! op_mm_20170818+

 SUBROUTINE parse_catstr(status, strlen, str, category, n)

    !AUTHOR:
    !       Mariano Mertens, DLR, 2017
    !
    !modified version from parse_datatstr (offemis)


    USE messy_main_tools,         ONLY: strcrack, str2num
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    INTEGER,  DIMENSION(:), POINTER :: category ! INTENT(OUT)
    ! LOCAL
    CHARACTER(LEN=strlen),               POINTER     :: sl1(:)
    INTEGER :: n
    INTEGER :: i


    status = 0 ! NO ERROR
    n = 0

    NULLIFY(sl1)


    !the different entries are seperated by a ;
    !Each entry then consists of a category_number
    CALL strcrack(str, ';', sl1, n) 
    ! -> n strings with tracer[_subname][,scaling];
    
    IF (n == 0) THEN
       status = -1
       RETURN
    END IF

    ALLOCATE(category(n))
    category(:) = 0

!loop over catagories 
    do i=1,n   
        CALL str2num(sl1(i),category(i), status)
        IF (status /= 0) RETURN
    end do 
    
    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)


    status = 0 ! NO ERROR
    
  END SUBROUTINE parse_catstr
!-------------------------------------------------------------------------
! op_mm_20170818-

!*************************************************************************
END MODULE messy_tagging
!*************************************************************************
