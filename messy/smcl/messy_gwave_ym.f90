MODULE messy_gwave_ym
!!===============================================================================
!!       THE SPECTRAL NONLINEAR GRAVITY WAVE PARAMETERIZATION
!!===============================================================================
! 
! This code calculates the momentum deposition by small-scale gravity waves,
! based on a specified GW spectrum and background atmosphere.   
! 
! (0) DEVELOPEMENT HISTORY:
! 
! The code has been developed by Erdal Yigit during his visit as a PhD student 
! at the Max Planck Institute for Solar System Research, Katlenburg-Lindau, in 
! collaboration with Dr. Alexander S. Medvedev. 
!
! (I) PUBLICATION HISTORY:
!
! --(A) Development of the parameterization: 
! 
! Yigit, E., A. D. Aylward, and A. S. Medvedev (2008), 
! Parameterization of the effects of vertically propagating gravity waves for 
! thermosphere general circulation models: Sensitivity study, 
! J. Geophys. Res., 113, D19106, doi:10.1029/2008JD010135. 
! 
! --(B) Implementation and validation of the scheme:
!
! Yigit, E., A. S. Medvedev, A. D. Aylward, P. Hartogh, and M. J. Harris (2009), 
! Modeling the effects of gravity wave momentum deposition on the circulation 
! above the turbopause, 
! J. Geophys. Res., 114, D07101, doi:10.1029/2008JD011132. 
!
! (II) THE CODE:
! 
! Basic inputs are, for example, background wind and temperature. The scheme
! calcualtes the gravity wave drag.  
! 
! (III) FURTHER DETAILS:
!
! -- Sept 16th, 2010.
! The version -v02 has been provided to Andreas Baumgaertner at Max Planck
! Institute for Chemistry to use in the coupled troposphere-stratosphere chemistry
! climate model-CMAT2.
!
! CONTACT:
! 
! Erdal Yigit, Ph.D.
! erdal@umich.edu
! http://aoss.engin.umich.edu/people/erdal
! University of Michigan, Ann Arbor
! Sept, 2010 
! 
!!==============================================================================


  USE messy_main_constants_mem, ONLY: prcn=>dp

  REAL(prcn), PUBLIC            :: p0_gw = 110.*1.e2  ! 15 km, mz_ab_20100916

  PUBLIC::EYGwave
  PUBLIC::SplineOn

CONTAINS

  SUBROUTINE EYGwave(v_eddy, rho, vy1d, vx1d, temp1d, h, &
       ut_gwd, vt_gwd, pres, eden, ht_dim)
    
    INTEGER,     INTENT(IN)   :: ht_dim
    REAL(prcn),  INTENT(IN)   :: vy1d(ht_dim)
    REAL(prcn),  INTENT(IN)   :: vx1d(ht_dim)
    REAL(prcn),  INTENT(IN)   :: temp1d(ht_dim)
    REAL(prcn),  INTENT(IN)   :: h(ht_dim)
    REAL(prcn),  INTENT(IN)   :: pres(ht_dim)
    REAL(prcn),  INTENT(IN)   :: rho(ht_dim)
    REAL(prcn),  INTENT(IN)   :: v_eddy(ht_dim)
    REAL(prcn),  INTENT(IN)   :: eden(ht_dim)
      
    REAL(prcn),  INTENT(OUT)  :: ut_gwd(ht_dim)
    REAL(prcn),  INTENT(OUT)  :: vt_gwd(ht_dim)

    REAL(prcn)                :: temp1d_m(ht_dim)
    REAL(prcn)                :: u_source_m(ht_dim)
    REAL(prcn)                :: pres_m(ht_dim)
    REAL(prcn)                :: rho_m(ht_dim)
    REAL(prcn)                :: h_m(ht_dim)
    REAL(prcn)                :: eden_m(ht_dim)
    
    INTEGER, PARAMETER        :: nh = 30

    REAL(prcn)                :: drag(ht_dim,nh)

    REAL(prcn), PARAMETER     :: g = 9.8067
    REAL(prcn), PARAMETER     :: pi = 3.1415926536

    REAL(prcn), SAVE          :: phasespeed(nh) 
    REAL(prcn)                :: uw_mom(nh) = 0.
    REAL(prcn)                :: brunt(ht_dim)
    REAL(prcn)                :: theta(ht_dim), thetap(ht_dim)
    INTEGER                   :: c_lev(nh)
    REAL(prcn)                :: dz(ht_dim)

    REAL(prcn)                :: tau(ht_dim, nh)     
    REAL(prcn)                :: beta(ht_dim, nh)           
    REAL(prcn)                :: beta_non(ht_dim,nh)
    REAL(prcn)                :: fac1, fac2

    REAL(prcn), ALLOCATABLE, SAVE :: alpha(:)   

    INTEGER                   :: n, i, j, init_zon = 1

    REAL(prcn)               :: c_int(ht_dim,nh) ! intrinsic phase speed, ci - u
    REAL(prcn)               :: sign_c(ht_dim,nh)
    REAL(prcn)               :: flux(ht_dim,nh)
    REAL(prcn)               :: up(ht_dim,nh)
    REAL(prcn)               :: upSq(ht_dim,nh)

    REAL(prcn), PARAMETER    :: flux0     = 0.00026
    REAL(prcn), PARAMETER    :: cw        = 35.
    REAL(prcn), PARAMETER    :: max_cp_y  = 80.
    REAL(prcn), PARAMETER    :: kx        = 2.*pi/300.e3 ! horizontal wave number (k_h)

    ! Total horizontal wind variance at the source level
    REAL(prcn)               :: rms=0.
    INTEGER                  :: sgn(nh)

    REAL(prcn)               :: m_vis(ht_dim), v_in(ht_dim)

    REAL(prcn)               :: sigma(ht_dim,nh)
    REAL(prcn)               :: sigmaSq(ht_dim,nh)
    REAL(prcn)               :: alpha_ins(ht_dim,nh)

    REAL(prcn)               :: max_ht = 300.e3
    REAL(prcn), PARAMETER    :: S2  =1.414213 ! 2.**0.5 !1.414213 
    REAL(prcn), PARAMETER    :: S2P = 2.506 ! (2.*pi)**0.5 !2.506 !
    REAL(prcn), PARAMETER    :: EF = 1._prcn 
    INTEGER, PARAMETER       :: SLEV = 1

    REAL(prcn)               :: u_source(ht_dim), gwd(ht_dim)
    REAL(prcn)               :: xv, yv
    
    DO n = 1, ht_dim
       
       IF(n==1) THEN
          u_source(n) = SQRT(vy1d(n)*vy1d(n) + vx1d(n)*vx1d(n))

          IF(u_source(n)==0.) u_source(n) = 1.

          yv = vy1d(n)/u_source(n)
          xv = vx1d(n)/u_source(n)
       ENDIF
       
       IF(n.GE.2) THEN
          u_source(n) = vy1d(n)*yv + vx1d(n)*xv
       ENDIF
    ENDDO 

    DO n = 1, ht_dim-1
       
       temp1d_m(n)   = 0.5*(temp1d(n)   + temp1d(n+1))
       u_source_m(n) = 0.5*(u_source(n) + u_source(n+1))
       pres_m(n)     = 0.5*(pres(n)     + pres(n+1))
       rho_m(n)      = 0.5*(rho(n)      + rho(n+1))
       h_m(n)        = 0.5*(h(n)        + h(n+1))
       eden_m(n)     = 0.5*(eden(n)     + eden(n+1))
       
    ENDDO 

    IF(init_zon == 1) THEN
       WRITE(6,*) "EY Gravity Wave Spectrum"
       WRITE(6,*) " i      Phasespeed     Variance      U_prime       Sigma_sq_tot    Hwhm"
       WRITE(6,*) "          [ms-1]       [m2s-2]        [m2s-2]         [m2s-2]           "
       WRITE(6,*)

       SPECTRUM : DO i = 1, nh

          IF(i==1)    phasespeed(i) = -80.
          IF(i.GE.2)  phasespeed(i) = phasespeed(i-1)*(40.**(-1./((nh/2.)-1.)))
          IF(i==(nh/2+1))   phasespeed(i) = 2.
          IF(i.GE.(nh/2+2)) phasespeed(i) = phasespeed(i-1)*(40.**(1./((nh/2.)-1.)))

          sgn(i)    = (phasespeed(i)-u_source_m(SLEV))/(ABS(phasespeed(i)-u_source_m(SLEV)))  
          uw_mom(i) = sgn(i)*(flux0*EXP(-((phasespeed(i)-u_source_m(SLEV))/cw)**2))
          
          upSq(1,i)   = ABS(uw_mom(i))*0.02/kx/ABS(phasespeed(i)-u_source_m(1))
          IF(ABS(upSq(1,i)).GT.0.) up(1,i) = SQRT(upSq(1,i))

          rms   = rms + 1./nh*upSq(1,i)
          IF (i == nh) rms = SQRT(rms) 
          
          WRITE(6,"(I3,2x,F20.16,3x,E14.3,1x,2F14.8,F8.2)") &
               i, phasespeed(i), uw_mom(i), up(1,i),  rms, cw
         
       ENDDO SPECTRUM
       WRITE(6,*) "-------------------------------------------------------"

       ALLOCATE(alpha(ht_dim))

       DO n=1, ht_dim-1
          alpha(n) = 3.e-6*(1.2+TANH(-7.*LOG(pres_m(n)/100000.)-50.)*0.833)
          write(6,"(I3,2x,2E8.2)") n, alpha(n), pres_m(n)
       ENDDO
    ENDIF

    init_zon = 0

    c_lev(:) = -1

    tau(:,:)   = 1.
    gwd(:)     = 0.

    drag(:,:)  = 0.
    sigmaSq(:,:) = 0.
    beta(:,:)    = 0.

    brunt(:)= 0.02
    beta_non = 0.   
    sigma = 0. 
    alpha_ins =  0. 
    upSq = 0.

    HT_LOOP : DO n = 1, ht_dim-2
       
       theta(n) = temp1d_m(n)*(100000./pres_m(n))**0.286
       
       IF (n==1) brunt(n) = 0.02
       
       IF(n.GE.2) THEN 
          dz(n)     = h_m(n) - h_m(n-1)
          thetap(n) = (theta(n) - theta(n-1))/dz(n)
          brunt(n)  = (ABS((g/theta(n))*thetap(n)))**0.5
       ENDIF
             
       m_vis(n) = 4.5E-05*(temp1d_m(n)/1000.0)**0.71/rho_m(n)  
       v_in(n)  = 7.22e-17*temp1d_m(n)**0.37*eden_m(n)       


       PS_LOOP : DO i = 1, nh
          
          sgn(i)    = (phasespeed(i)-u_source_m(SLEV))/(ABS(phasespeed(i)-u_source_m(SLEV)))  
          uw_mom(i) = sgn(i)*(flux0*EXP(-((phasespeed(i)-u_source_m(SLEV))/cw)**2))
                    
          IF(h(n) > max_ht) CYCLE
          
          c_int(n,i)  = phasespeed(i) - u_source_m(n)
          IF (c_int(n,i).NE.0) sign_c(n,i) = c_int(n,i)/ABS(c_int(n,i))
          IF (c_int(n,i).EQ.0) sign_c(n,i) = 1.


          CRITICAL : IF(c_lev(i).LT.0.AND.n.GE.2.AND.(sign_c(n,i).NE.sign_c(n-1,i).AND.c_int(n,i).NE.0.)) THEN
             c_lev(i)  = n
             flux(n,i) = 0.
             
           ELSE IF (c_lev(i).LT.0) THEN
          
             fac1      = 2.*brunt(n)*brunt(n)*brunt(n)/(kx*c_int(n,i)**4)
             fac2      = 2.*brunt(n)/(kx*c_int(n,i)*c_int(n,i))
             beta(n,i) = fac1*(m_vis(n)+v_eddy(n)) + fac2*(v_in(n)+alpha(n))


             IF(n==1) THEN
                tau(n,i)  = 1.
                flux(n,i) = uw_mom(i)
                upSq(n,i) = ABS(flux(n,i))*brunt(n)/kx/ABS(c_int(n,i))
                IF(ABS(upSq(n,i)).GT.0.) up(n,i) = SQRT(upSq(n,i))
             ENDIF 

             UPWARD : IF(n.GE.2) THEN 

                DO j = 1, nh
                   c_int(n,j) = phasespeed(j) - u_source_m(n)
                   IF (ABS(c_int(n,i)).GE.ABS(c_int(n,j))) THEN
                      sigmaSq(n,i) = sigmaSq(n,i) + upSq(n-1,j)
                   ENDIF
                ENDDO 
                
                IF(sigmaSq(n,i).GE.1E-30) THEN
                   sigma(n,i)     = SQRT(sigmaSq(n,i))
                   alpha_ins(n,i) = ABS(c_int(n,i))/(2.**0.5)/sigma(n,i)
                   beta_non(n,i)  = ((2.*pi)**0.5)*brunt(n)/sigma(n,i)*EXP(-1.*alpha_ins(n,i)*alpha_ins(n,i))
                ENDIF
       
                beta(n,i) = beta(n,i) + beta_non(n,i)
             
                tau(n,i)  = tau(n-1,i)*EXP(-1.*dz(n)*(beta(n,i)+beta(n-1,i))*0.5)
                flux(n,i) = uw_mom(i)*rho(1)/rho(n)*tau(n,i)
                upSq(n,i) = ABS(flux(n,i))*brunt(n)/kx/ABS(c_int(n,i))
                IF(ABS(upSq(n,i)).GT.0.) up(n,i) = SQRT(upSq(n,i))
                
                drag(n,i) = beta(n,i)*flux(n,i)*EF

                IF(alpha_ins(n,i).LT.0.75) drag(n,i) = 0.

                gwd(n)    = gwd(n) + drag(n,i)
             ENDIF UPWARD
          ENDIF CRITICAL
       ENDDO PS_LOOP
    ENDDO HT_LOOP

   DO n=2, ht_dim-2
      gwd(n) = 0.5*(gwd(n) + gwd(n-1))
   ENDDO

   SMOOTH : DO n = 2, ht_dim-1
      gwd(n)  = (gwd(n-1)  + 2.*gwd(n)  + gwd(n+1))*0.25
   ENDDO SMOOTH

   ut_gwd(ht_dim-1:ht_dim) = 0.
   vt_gwd(ht_dim-1:ht_dim) = 0.
   BACK_PROJECT : DO n=1, ht_dim-2
      ut_gwd(n) = yv * gwd(n)
      vt_gwd(n) = xv * gwd(n)
   ENDDO BACK_PROJECT
   
   RETURN

 END SUBROUTINE EYGwave
 
   SUBROUTINE SplineOn(cmat_field, interp_field, field_flag, zero_winds,       &
       placemat, placegwave, ht_dim, lev, status)

     USE messy_main_tools, ONLY: Spline1d, Splint1d

    INTEGER,    INTENT(IN)  :: ht_dim
    INTEGER,    INTENT(IN)  :: lev
    INTEGER,    INTENT(IN)  :: field_flag
    INTEGER,    INTENT(IN)  :: zero_winds
    REAL(prcn), INTENT(IN)  :: cmat_field(ht_dim) 
    REAL(prcn), INTENT(IN)  :: placemat(ht_dim)
    REAL(prcn), INTENT(IN)  :: placegwave(lev)

    REAL(prcn), INTENT(OUT) :: interp_field(lev)
    INTEGER,    INTENT(OUT) :: status


    ! arrays for the interpolation
    REAL(prcn) :: splineno(ht_dim),splineno11(ht_dim),splineno1(lev)

    INTEGER    :: k

    ! first zonal winds...(bottom up)
    DO k=1,ht_dim
       splineno(k)=cmat_field(ht_dim+1-k)
    ENDDO
    
    CALL Spline1d(placemat,splineno,ht_dim,0._prcn,0._prcn,splineno11,.TRUE.)
    
    DO k=1,lev
       CALL Splint1d(placemat,splineno,splineno11,ht_dim,placegwave(k),       &
            splineno1(k),status)
    ENDDO
    
    ! We make code handle any altitude range. mjh
    DO k=1,lev

       IF(placemat(ht_dim) >= placegwave(k)) THEN
          interp_field(lev+1-k) = splineno1(k)
       ELSE

          ! No gradient or zero below lower b
          interp_field(lev+1-k) = interp_field(lev+1-k +1)
          
          ! It's a velocity field and zero_winds is on
          IF(zero_winds == 1 .AND. field_flag == 1) interp_field(lev+1-k) = 0.

          ! It's a height field
          IF(field_flag == 2) interp_field(lev+1-k) = interp_field(lev+1-k +1) &
               - (interp_field(lev+1-k +2) - interp_field(lev+1-k +1))

       ENDIF
    ENDDO

    
  END SUBROUTINE SplineOn
 ! ****************************************************************

  SUBROUTINE ym_read_nml_ctrl(status, iou, modstr)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    CHARACTER(LEN=*), INTENT(IN) :: modstr

    NAMELIST /CTRL_YM/ p0_gw  

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='hines_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL_YM', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_YM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_YM', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST
    ! ...

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE ym_read_nml_ctrl

END MODULE messy_gwave_ym
