MODULE messy_gwave_hlm

  USE messy_main_constants_mem, ONLY: prcn=>dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: HLMGwave
  !PRIVATE :: PT_N

CONTAINS
  
  !r==================================================================
  !r=               Hybrid-Lindzen-Matsuno (HLM)                     =
  !r==================================================================
  !r 
  !r This is a gravity wave drag parameterization, based on
  !r Meyer [1999] JGR, 28,181-28,196.
  !r by Erdal Yigit
  !r It's a hybrid between the Matsuno transmission parameterization, and a
  !r Holton/Lindzen scheme. It has the advantage of iterating a few times such
  !r that GWs break according to the eddy diff of earlier iteration. 
  !r Wave-wave interactions are neglected, instead superposition of effects
  !r of all waves are evaluated.
  !r
  !r In CMAT2 this parameterization will be called inside Earth_Planet.f90
  !r
  !r Inputs
  !r --------------------------------------------------------------
  !r vy1d         : 1d array of east-west velocity (m*s-1)
  !r vx1d         : 1d array of north-south velocity (m*s-1)
  !r scht1d       : 1d array of scale height (m)
  !r temp1d       : 1d array of temperature (K)
  !r cp1d         : 1d array of specific heat capacity (K*kg-1*s-1)
  !r pres         : 1d array of pressure (Nm-2)
  !r direction    : switch for shifting the spectrum for EW
  !r                (1-EW 2-NS)
  !r ht_dim       : size of vertical 1d arrays
  !r
  !r Outputs
  !r -------------------------------------------------------------
  !r
  !r com_drag     : 1d array of drag (m*s-1*day-1))
  !r com_v_eddy   : 1d array of eddy diffusion (m* s-2)
  !r
  !r

  
  SUBROUTINE HLMGwave(vy1d, scht1d, temp1d, cp1d, h, com_drag, com_v_eddy, &
       pres, direction, com_gwh, ht_dim)

    IMPLICIT NONE

    INTEGER,      INTENT(IN)    :: ht_dim
    REAL(prcn),   INTENT(IN)    :: vy1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: scht1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: temp1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: cp1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: h(ht_dim)     
    REAL(prcn),   INTENT(IN)    :: pres(ht_dim)
    INTEGER,      INTENT(IN)    :: direction       ! Direction   

    REAL(prcn),   INTENT(OUT)   :: com_drag(ht_dim)
    REAL(prcn),   INTENT(OUT)   :: com_v_eddy(ht_dim)
    REAL(prcn),   INTENT(OUT)   :: com_gwh(ht_dim)
    


    ! Define spectrum
    ! The spectrum is currently gaussian
 
    INTEGER,    PARAMETER   :: no_waves = 15       ! Number of waves in spectrum
    INTEGER                 :: offset              ! How much the peak is moved from 0 ms-1
    REAL(prcn), PARAMETER   :: center   = 0.       ! Centre of spectrum ms-1
    REAL(prcn), PARAMETER   :: max_val  = 0.5e-4   ! Maximum momentum flux m2s-2
    REAL(prcn), PARAMETER   :: max_v    = 60.      ! How far out to sample Gaussian ms-1
    REAL(prcn), PARAMETER   :: fwhm     = 40.      ! Full width half maximum ms-1
    REAL(prcn), PARAMETER   :: kx = 1./100e3       ! Horiz wave # (denom. is w.length in m)
    REAL(prcn), PARAMETER   :: eff = 30.e-2        ! Intermittency factor
   

    ! Number of iterations

    INTEGER      :: num_iters = 1

!!$    ! Some constants
!!$
!!$    REAL(prcn), PARAMETER :: g = 9.8067       
!!$    REAL(prcn), PARAMETER :: pi = 3.1415926536 
   

    ! Declare some wave variables
    
    REAL(prcn)   :: phasespeed(no_waves) = 0.          
    REAL(prcn)   :: phasespeed_zon(no_waves) = 0.          
    REAL(prcn)   :: phasespeed_merid(no_waves) = 0.          


    REAL(prcn)   :: uw_mom(no_waves) = 0.
    REAL(prcn)   :: uw_mom_zon(no_waves) = 0.
    REAL(prcn)   :: uw_mom_merid(no_waves) = 0.

    REAL(prcn)   :: bruntv(ht_dim)              ! Brunt vaisala frequency, N
    ! REAL(prcn)   :: richardson(ht_dim)          ! Richardson number
    REAL(prcn)   :: u_tilde(no_waves) = 0.      ! Initial wave flux
    REAL(prcn)   :: u_tilde_zon(no_waves) = 0.  ! Initial wave flux
    REAL(prcn)   :: u_tilde_merid(no_waves) = 0.  ! Initial wave flux

    INTEGER      :: break_lev(no_waves)         ! little trick Matt has introduced
    REAL(prcn)   :: mfd_ls(ht_dim,no_waves)     ! mom. flux divergence due to wave saturation
    REAL(prcn)   :: mfd_trans(ht_dim,no_waves)  ! mom. flux divergence due to transmission
    REAL(prcn)   :: gwh(ht_dim,no_waves)        ! wave heating

    REAL(prcn)   :: com_mfd_ls(ht_dim)          ! composite value for lindzen mfd
    REAL(prcn)   :: com_mfd_trans(ht_dim)       ! composite value for matsuno mfd
    REAL(prcn)   :: c_i(ht_dim,no_waves)        ! intrinsic phase speed (u-c)
    REAL(prcn)   :: eddy_f(ht_dim)              ! 3*H*c_i^4*N^(-3) term 
    REAL(prcn)   :: sigma = 0.

    REAL(PRCN)   :: z_c(no_waves)
    INTEGER      :: c_level(no_waves)

    REAL(PRCN)   :: t_level(no_waves)
    
    ! Some local variables

    INTEGER      :: iter
    
    ! Some atmospheric variables

    REAL(prcn)   :: dz(ht_dim)                   !height increment
    REAL(prcn)   :: thetap(ht_dim)               !potential temperature
    REAL(prcn)   :: diff_coeff(ht_dim)           !molecular+eddy diff coeff ! mz_ab_20100901
    REAL(prcn)   :: tau(ht_dim)                  !transmissivity
    REAL(prcn)   :: z_b(no_waves)=0.             !evaluated breaking ht according the formula
    REAL(prcn)   :: z_h(no_waves)                !The first breaking ht
    REAL(prcn)   :: v_eddy(ht_dim,no_waves)      !eddy diffusion 
    REAL(prcn)   :: vy1d_del(ht_dim)             !gradient of the background zonal mean wind 
    REAL(prcn)   :: taup(ht_dim)                 !gradient of the transmissivity
    REAL(prcn)   :: fac
    
     
    INTEGER      :: i    
    INTEGER      :: nd    
    INTEGER      :: nu 
    INTEGER      :: nuu 
    INTEGER      :: n

    ! 0 - No off
    ! 1 - print out input data for every call then stop
    ! 2 - print out input and output then stop

    INTEGER      :: debug = 0

    ! Initilisation flag

    INTEGER      :: initialise_zon = 1   
    INTEGER      :: initialise_merid = 1

    ! Maximum height to calculate to

    REAL(prcn)   :: max_ht = 140.e3
    
    ! Height to damp gwaves above
    
    REAL(prcn)   :: turb_ht = 98.e3
    INTEGER      :: ht_100

    diff_coeff = 0.
    gwh = 0. 
    mfd_trans = 0.

    ! Set up the launch spectrum

   
    IF(initialise_zon == 1 .OR. initialise_merid == 1) THEN

       sigma    = fwhm/(2.*SQRT(2.*LOG(2.)))
       
       WRITE(6,*) "HLM Gravity Wave Spectrum"
       WRITE(6,*) "i        u(ms-1)        mom(m2s-2)        utilde(ms-1)        fwhm"
       WRITE(6,*)


       SPECTRUM : DO i = 1, no_waves

          IF (direction == 1 .AND. initialise_zon ==1) THEN 
             
             !peak position
             offset = 0

             phasespeed_zon(i) = 2.*max_v*(((i-1.)/(no_waves-1)) - 0.5) + offset

             uw_mom_zon(i) = max_val*exp((-((phasespeed_zon(i)  - (center &
                  + offset)))**2.)/(2.*(sigma**2.)))

              u_tilde_zon(i) = ((uw_mom_zon(i)*2.*0.02)/kx)**(1./3.)

            !  u_tilde_zon =  (/0.1,0.8,1.5,2.0,1.5,0.8,0.1/)

             WRITE(6,"(I2,1x,F15.7,3x,E15.7,3x,F15.7,3x,F14.7)") &
                  i, phasespeed_zon(i), uw_mom_zon(i), u_tilde_zon(i), fwhm 
             


          ENDIF


          IF (direction == 2 .AND. initialise_merid == 1) THEN 
             
             phasespeed_merid(i) = 2.*max_v*(((i-1.)/(no_waves-1)) - 0.5) 
             
             uw_mom_merid(i) = max_val*exp((-((phasespeed_merid(i)  - center))**2.) &
                  /(2.*(sigma**2.)))
             
              u_tilde_merid(i) = ((uw_mom_merid(i)*2.*0.02)/kx)**(1./3.)
             ! u_tilde_merid = (/0.1,0.8,1.5,2.0,1.5,0.8,0.1/)             

             WRITE(6,"(I2,1x,F15.7,3x,E15.7,3x,F15.7,3x,F14.7)") &
               i, phasespeed_merid(i), uw_mom_merid(i), u_tilde_merid(i), fwhm 

          ENDIF

          
       ENDDO SPECTRUM
       
       
       IF (direction == 1) THEN 
          initialise_zon = 0
       ENDIF

       IF (direction ==2) THEN 
          initialise_merid = 0
       ENDIF
    ENDIF


    IF (direction == 1) THEN
       phasespeed(:) = phasespeed_zon(:) 
       u_tilde(:) = u_tilde_zon(:) 
       uw_mom(:) = uw_mom_zon(:)
    ELSE
       phasespeed (:) = phasespeed_merid(:) 
       u_tilde(:) = u_tilde_merid(:) 
       uw_mom(:) = uw_mom_merid(:)
    ENDIF
    

    CALL PT_N(temp1d, pres, h, thetap, dz, bruntv, ht_dim)

     com_v_eddy(:) = 0.

    ITERATE_LOOP : DO iter = 1,num_iters
       
       PS_LOOP : DO i=1, no_waves

          break_lev(i) = -1
          c_level(i) = -1
          t_level(i) = -1
          
          HT_LOOP : DO n = 1, ht_dim-1

             ! Don't bother above 140km

             IF(h(n) > max_ht) CYCLE
             
             nd  = n - 1
             nu  = n + 1
             nuu = n + 2 
             
             ! Molecular diffusion + eddy diffusion
             
             diff_coeff(n)  = (200.*exp((h(n)-110000.)/scht1d(n))) + com_v_eddy(n)
             diff_coeff(nu) = (200.*exp((h(nu)-110000.)/scht1d(nu))) + com_v_eddy(nu)

             CALC_BREAK_HT : IF(break_lev(i) < 0) THEN
                
                ! Transmissivity
                
                IF(n == 1) THEN
                   
                   tau(n) = exp((-dz(n)/kx)*((( diff_coeff(n) *bruntv(n)**3)/ &
                        (phasespeed(i)-vy1d(n))**4)+ &
                        ((diff_coeff(nu)*bruntv(n+1)**3)/(phasespeed(i)-vy1d(n+1))**4)))  
                   
                   tau(nu) = exp((-dz(nu)/kx)*(((diff_coeff(nu)*bruntv(nu)**3)/ & 
                        (phasespeed(i)-vy1d(nu))**4)+&
                        ((diff_coeff(nuu)*bruntv(nuu)**3)/(phasespeed(i)-vy1d(nuu))**4)))  
                   
                ELSE
                   
                   tau(n) = tau(n-1)*exp((-dz(nd)/kx)*(((diff_coeff(nd)*bruntv(nd)**3)/ &
                        (phasespeed(i)-vy1d(nd))**4)+ &
                        ((diff_coeff(n)   *bruntv(n)**3)/(phasespeed(i)-vy1d(n))**4))) 
                   
                   tau(nu) = tau(n)*exp((-dz(n)/kx)*(((diff_coeff(n)*bruntv(n)**3)/ &
                        (phasespeed(i)-vy1d(n))**4)+ &
                        ((diff_coeff(nu)    *bruntv(nu)**3)/(phasespeed(i)-vy1d(nu))**4)))


                   ! IF tau is close to zero then put it to a small number
                   ! so that it doesn't cuase any numerical errors.
                   
                   IF (tau(n).LT.1.e-5) tau(n) = 1.e-6
                    
                ENDIF
                
                ! Gradient of the transmissivity
                
                taup(n) = ( (tau(nu)-tau(n)) /dz(n) )
                
                ! Determine the breaking height
                ! -----------------------------

                IF(tau(n) /= 0. .AND. u_tilde(i) /= 0.) THEN 
                   z_b(i) = h(1) + (3*scht1d(n)*log( (abs(vy1d(n)-phasespeed(i)))/ &
                        (tau(n)**(2./3.)*u_tilde(i)) ))
                ENDIF

                IF(z_b(i).LT.h(n)) THEN
                   
                   z_h(i) = z_b(i)
                   
                   ! Calculates the eddy diffusion at the breaking level 
                   
                   break_lev(i) = n
                   
                ENDIF
                
             ENDIF CALC_BREAK_HT
             
             ! Determine the critical heights 
             ! ------------------------------
             
             CR_HT : IF(c_level(i) < 0. .AND. n > 1.) THEN 
                
                c_i(n,i) = vy1d(n) - phasespeed(i)
                
                IF((phasespeed(i).GT.vy1d(n-1).AND.phasespeed(i).LT.vy1d(n)).OR.&
                     (phasespeed(i).LT.vy1d(n-1).AND.phasespeed(i).GT.vy1d(n))) THEN 
                   
                   z_c(i) = h(n)  
                   
                   c_level(i) = n           
                   
                ENDIF
                
             ENDIF CR_HT
             
             c_i(n,i) = vy1d(n) - phasespeed(i)
             
             vy1d_del(n) = (vy1d(n+1)-vy1d(n))/dz(n)
             
          ENDDO HT_LOOP
          
          vy1d_del(ht_dim) = 0.


          CALC_EDDY_DRAG : IF ( (break_lev(i) > 0. .AND. c_level(i) > 0).OR.&
               (break_lev(i) > 0. .AND. c_level(i) < 0.)) THEN      
             
             ABOVE_BREAK : DO n = break_lev(i), ht_dim-1                  

                IF(h(n) > max_ht) CYCLE

                ! Eddy Diffusion variables             
                
                eddy_f(n)   = eff*kx*((c_i(n,i))**4)*(1/(bruntv(n))**3)
                
                IF(c_level(i) /= -1 .AND. n >= c_level(i)) eddy_f(n)=0.
                
                
                ! Finally Eddy Diffusion
                v_eddy(n,i) = eddy_f(n) * (1/(2*scht1d(n)) - (1.5*vy1d_del(n)*(1/c_i(n,i))))
                
                ! and... Wave Saturation acceleration
                mfd_ls(n,i) = -1.*bruntv(n)**2*v_eddy(n,i)*(1/(c_i(n,i)))
                
                ! wave heating
                gwh(n,i)    = -1.*mfd_ls(n,i)*c_i(n,i)/cp1d(n) 

                ! Momentum flux divergence according to Matsuno
                IF(tau(n).GT.0.) mfd_trans(n,i) = -1.*uw_mom(i)*taup(n)


             ENDDO ABOVE_BREAK
             
             ! Exponential decay below breaking height                         
             EXP_DEC : DO n = break_lev(i)-1, 1, -1      
                ! the exponential decay factor
                fac = EXP(-5.*(ABS  (h(n) - z_b(i)) )/(scht1d(n)))
                ! Exp decay is applied to every waves' mfd seperately 
                mfd_ls(n,i)     = mfd_ls(break_lev(i),i)*fac
                v_eddy(n,i)     = v_eddy(break_lev(i),i)*fac
             ENDDO EXP_DEC
             
             ! Curently exp decay above the breaking height
                            
             ht_100 = 0
             EXP_DEC_ABOVE_BREAK : DO n = break_lev(i)+1, ht_dim, 1      
                
                ! Don't bother above 140km                              
                IF(h(n) > max_ht) CYCLE
                
                ! Don't bother if break_lev=-1
                IF(n<1) CYCLE

                IF(ht_100 == 0 .AND. h(n) > turb_ht) ht_100 = n

                v_eddy(n,i) = v_eddy(break_lev(i),i) * EXP(-(ABS( h(n) - z_b(i) ))/scht1d(n))
                
                mfd_ls(n,i) = mfd_ls(break_lev(i),i) * EXP(-(ABS( h(n) - z_b(i) ))/scht1d(n))
                gwh(n,i)    = gwh(break_lev(i),i) * EXP(-(ABS( h(n) - z_b(i) ))/scht1d(n))

             ENDDO EXP_DEC_ABOVE_BREAK
             
             IF (ht_100 < 1) ht_100 = 1

             ! Force Decay above 100km

             EXP_DEC_ABOVE_100 : DO n = ht_100, ht_dim, 1

                IF(h(n) > turb_ht) THEN
                   
                   v_eddy(n,i) = v_eddy(n,i) * EXP(-1.5*(ABS( h(n) - turb_ht))/scht1d(n))
                   mfd_ls(n,i) = mfd_ls(n,i) * EXP(-1.5*(ABS( h(n) - turb_ht))/scht1d(n))
                   gwh(n,i)    = gwh(n,i)* EXP(-1.5*(ABS(h(n) - turb_ht))/scht1d(n)) 
                ENDIF


             ENDDO EXP_DEC_ABOVE_100
             
          ELSE
             
             !Set things to zero in any other case
             
             OTHERS : DO n = 1, ht_dim
                
                v_eddy(n,i) = 0.
                mfd_trans(n,i) = 0.
                mfd_ls(n,i) = 0.
                gwh(n,i) = 0.
             ENDDO OTHERS
            
             
          ENDIF CALC_EDDY_DRAG
          
          DEB : IF (debug==2) THEN

             WRITE(6,*)''
             WRITE(6,'(A2,A12,A16,A12,1x,A14,1x,A10,1x,A10,2x,A10,3x,A12)') &
                  "N", "HEIGHT","PHASESPEED","CR_HEIGHT","BR_HEIGHT","TAU","EDDY",&
                  "LS_DRAG","TRANS_DRAG"
             DO n = 1, ht_dim
                WRITE(6,'(I2,1x,F12.2,F12.1,3x,F12.4,3x,F12.4,F12.7,F12.4,F12.4,F12.7)')    &
                     n, h(n), phasespeed(i), z_c(i), z_h(i), tau(n), v_eddy(n,i), &
                     mfd_ls(n,i)*86400., mfd_trans(n,i)*86400.
             ENDDO
          ENDIF DEB
          
       ENDDO PS_LOOP

       ! Initial values for the composite terms         

       com_v_eddy(:)    = 0.
       com_mfd_ls(:)    = 0.
       com_mfd_trans(:) = 0.
       com_drag(:)      = 0.     
       com_gwh(:)       = 0.
       
       ! Summate to get total drag and eddy
       TOTAL : DO i = 1, no_waves
          
          DO n = 2, ht_dim-1
             !contains the exp decay below the z_b
             
             com_v_eddy(n) = com_v_eddy(n) + ((v_eddy(n-1,i) + v_eddy(n,i) + v_eddy(n+1,i))/3.)
             com_mfd_ls(n) = com_mfd_ls(n) + ((mfd_ls(n-1,i) + mfd_ls(n,i) + mfd_ls(n+1,i))/3.) 
             com_gwh(n)    = com_gwh(n) + ((gwh(n-1,i) + gwh(n,i) + gwh(n+1,i))/3.)

             !no exp decay has been applied to mfd_trans 
             
             com_mfd_trans(n) = com_mfd_trans(n) + &
                  ((mfd_trans(n-1,i) + mfd_trans(n,i)+ mfd_trans(n+1,i))/3.)
             
          ENDDO
       ENDDO TOTAL
       
       
       DO n = 1, ht_dim-1
          
          ! Total drag, matsuno +  lindzen
          
          com_drag(n) = com_drag(n) + com_mfd_trans(n) + com_mfd_ls(n)
          
          DEB1 : IF (debug==2) THEN      
             
             WRITE(6,*) n, iter, h(n), com_v_eddy(n) , com_mfd_ls(n)*86400., & 
                  com_mfd_trans(n)*86400., com_drag(n)*86400.
             

          ENDIF DEB1
       ENDDO       

       ! 2nd iteration the transmissivity includes the com_v_eddy(n)
       ! Calculate the transmissivity with com_v_eddy(n)
       
       
    ENDDO ITERATE_LOOP
    
    IF (debug==2) THEN
       STOP    
    ENDIF
   

    ! This is for debugging, it prints out the input data and checks whether 
    ! they are alright
    
    DEB2 : IF (debug == 1) THEN
       
       INPUT : DO i = 1, no_waves
          
          WRITE(6,*) "NEXT WAVE"
          WRITE(6,*) "PHASESPEED = ",phasespeed(i)
          WRITE(6,*)
          WRITE(6,'(A2,1x,A10,1x,A12,1x,A12,1x,A12,1x,A9)') &
               "N","HEIGHT","Z_WIND","INT_P_SPEED","SCALE_HT","TEMP"

          DO n = 1,ht_dim
             
             WRITE(6,'(I2,F12.3,F12.3,F12.3,F12.3,F12.3,F12.3,F12.3,F12.4,F15.7,F15.7)') &
                  n, h(n), vy1d(n), c_i(n,i), scht1d(n), temp1d(n), & 
                  cp1d(n), pres(n), bruntv(n), vy1d_del(n), tau(n) 
          ENDDO

         

       ENDDO INPUT

       OUTPUT : DO i = 1, no_waves
          WRITE(6,*) "NEXT WAVE", phasespeed(i)
          DO n = 1,ht_dim
             WRITE(6,'(I2,1x,F13.2,F13.6,F13.6,F13.6,F13.6,F13.6)') &
                  n, h(n), v_eddy(n,i), mfd_ls(n,i)*86400, mfd_trans(n,i),86400.
          ENDDO
       ENDDO OUTPUT
       
       WRITE(6,'(5A15)') "N","HEIGHT","TOTAL DRAG","TOTAL EDDY","WAVE HEATING"

       FINAL : DO n = 1, ht_dim
          WRITE(6,'(I2,1x,4F15.3)') &
               n, h(n), com_drag(n)*86400, com_v_eddy(n), com_gwh(n)*cp1d(n)*86400.    
       ENDDO FINAL
    ENDIF DEB2
    

    RETURN

  END SUBROUTINE HLMGwave


  SUBROUTINE PT_N(temp1d, pres, z, thetap, dz, bf, ht_dim)

  !e ****************************************************************
  !e  Calculate Buoyancy frequency and Potential temperature variables
  !e  Erdal Yigit
  !e  Max Planck Instituet fuer Sonnensystemforschung
  !e  Nov 2007
  !e ****************************************************************

    IMPLICIT NONE

    INTEGER,    INTENT(IN)  :: ht_dim
    REAL(prcn), INTENT(IN)  :: temp1d(ht_dim), pres(ht_dim), z(ht_dim)
    REAL(prcn), INTENT(OUT) :: thetap(ht_dim), dz(ht_dim), bf(ht_dim)
    REAL(prcn)              :: theta(ht_dim)
    REAL(prcn), PARAMETER   :: g = 9.8067
    INTEGER                 :: s_debug = 0 ! if 1 print out variables. 
    INTEGER                 :: n

    thetap(:) = 0.
    theta(:)  = 0.  

    DO n = 1, ht_dim-2
       !  potential temperature
       theta(n)    = temp1d(n)*(100000./pres(n))**0.286
       theta(n+1)  = temp1d(n+1)*(100000./pres(n+1))**0.286
       theta(n+2)  = temp1d(n+2)*(100000./pres(n+2))**0.286
       
       dz(n) = z(n+1) - z(n)
    ENDDO
    
    dz(ht_dim-1) = z(ht_dim)-z(ht_dim-1)

    DO n = ht_dim-2, ht_dim
       theta(n)    = temp1d(n)*(100000./pres(n))**0.286
    ENDDO

    !  Calculate the potential temperature gradient
    thetap(1) = (theta(2)-theta(1)) / dz(1)
    
    PTEMP : DO n = 2, ht_dim-1
       thetap(n) = (theta(n+1)-theta(n-1)) / (dz(n)+dz(n-1))
    ENDDO PTEMP
    
    !  set the upper boundary gradient of the pot. temp to the one below
    thetap(ht_dim) = thetap(ht_dim-1)

    ! Buoyancy frequency
    bf(:) = 0.
    
    BRUNT : DO n = 1, ht_dim
       !   Brunt Vaisala Frequency
       bf(n) = ((g/theta(n))*(thetap(n)))**(0.5)
       
       IF (bf(n).LE.0.) bf(n) = 0.02

       IF (s_debug.EQ.1) THEN 
          WRITE(6,'(4F15.6)') z(n), bf(n), theta(n), thetap(n)
       ENDIF
    ENDDO BRUNT
    


  END SUBROUTINE PT_N
    

END MODULE messy_gwave_hlm

