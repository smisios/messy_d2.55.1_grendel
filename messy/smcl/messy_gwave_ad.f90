MODULE messy_gwave_ad

  USE messy_main_constants_mem, ONLY: prcn=>dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: AlexDGwave


CONTAINS


  !r--------------------------------------------------------------
  !r              ALEXANDER AND DUNKERTON 1999
  !r ------------------------------------------------------------
  !r This is a spectral parameterization based on the works of
  !r Lindzen and Holton.
  !r
  !r ERDAL YIGIT, Golders Green
  !r
  !r Inputs
  !r --------------------------------------------------------------
  !r vy           : 1d array of east-west velocity (m*s-1)
  !r vx1d         : 1d array of north-south velocity (m*s-1)
  !r scht1d       : 1d array of scale height (m)
  !r temp1d       : 1d array of temperature (K)
  !r cp1d         : 1d array of specific heat capacity (K*kg-1*s-1)
  !r pres         : 1d array of pressure (Nm-2)
  !r dir          : switch for shifting the spectrum for EW
  !r                (1-EW 2-NS)
  !r ht_dim       : size of vertical 1d arrays
  !r
  !r Outputs
  !r -------------------------------------------------------------
  !r
  !r tot_drag     : 1d array of drag (m*s-1*day-1))
  !r tot_dif      : 1d array of eddy diffusion (m* s-2)
  !r
  !r


 SUBROUTINE AlexDGwave(flag, u, scht1d, temp1d, cp1d, z, gwf, eddy, pres, &
       rho, bf, dir, gwh, ht_dim)

   USE messy_main_constants_mem, ONLY: g, pi

   IMPLICIT NONE ! op_pj_20110202

    INTEGER,      INTENT(IN)    :: ht_dim
    REAL(prcn),   INTENT(IN)    :: u(ht_dim)
    REAL(prcn),   INTENT(IN)    :: scht1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: temp1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: cp1d(ht_dim)
    REAL(prcn),   INTENT(IN)    :: z(ht_dim)
    REAL(prcn),   INTENT(IN)    :: pres(ht_dim)
    REAL(prcn),   INTENT(IN)    :: rho(ht_dim)    
    INTEGER,      INTENT(IN)    :: flag 
    INTEGER,      INTENT(IN)    :: dir

    REAL(prcn),   INTENT(OUT)   :: gwf(ht_dim) 
    REAL(prcn),   INTENT(OUT)   :: eddy(ht_dim)  
    REAL(prcn),   INTENT(INOUT) :: bf(ht_dim)
    REAL(prcn),   INTENT(OUT)   :: gwh(ht_dim) 

    ! constants

!!$    REAL(prcn), PARAMETER    :: g = 9.8067
!!$    REAL(prcn), PARAMETER    :: pi = 3.1415926536
!!$    REAL(prcn), PARAMETER    :: R_star = 287.7         ! gas constant in J kg-1 K-1
    REAL(prcn), PARAMETER    :: Bw = 0.018
    REAL(prcn), PARAMETER    :: cw = 30.
    REAL(prcn)               :: Bt, off 
    REAL(prcn)               :: c_0
    REAL(prcn), PARAMETER    :: k = 2*pi/(100.e3)

    INTEGER, PARAMETER       :: nc = 121
    INTEGER, PARAMETER       :: iz0 = 1
    REAL(prcn), PARAMETER    :: dc = 1.0

    REAL(prcn)         :: theta(ht_dim), dz(ht_dim),Hb(ht_dim)
    REAL(prcn)         :: thetap(ht_dim)       

    INTEGER            :: in_zon = 1
    INTEGER            :: in_mer = 1
    REAL(prcn)         :: c_0_zon =0.
    REAL(prcn)         :: c_0_mer=0.
    REAL(prcn)         :: Bt_zon=0.
    REAL(prcn)         :: Bt_mer=0.
    REAL(prcn)         :: off_z=0.
    REAL(prcn)         :: off_m=0.

    
    ! Vartiables to prevent calculation above altcutoff (km). mjh
    INTEGER            :: maxalt=0
    REAL(prcn)         :: altcutoff=130.e3

    INTEGER            :: n, i
    
    REAL(prcn)         :: fm,fe,Bsum
    REAL(prcn)         :: Nb,rb,ub,ubm,ubp,alp2,k2,fac,Foc,u0,N0,r0,H0,omc

    REAL(prcn)         :: B0(nc),c0(nc)
    REAL(prcn)         :: c,ci,cmax,c0j,B0j,cj0,L2,j0
    REAL(prcn)         :: eps,dfac,rbh, &
                          eddy2(ht_dim), gwf2(ht_dim), gwh2(ht_dim)
    
    INTEGER            :: rflg,flg0,rmsk,cmsk,icm(2),ict(2),sgn
    INTEGER            :: msk(nc)
    
    IF(dir==1.AND.in_zon==1) THEN 
       c_0_zon = 0.
       Bt_zon  = 9.0e-6
       off_z   = 0.
       WRITE(6,*)'****Zonal net mom flux & shift****************'
       WRITE(6,*)'STARTING VALUES : ', Bt_zon, c_0_zon, off_z
    ENDIF

    IF(dir==2.AND.in_mer==1) THEN
       c_0_mer = 0.
       Bt_mer  = 0.2*Bt_zon
       off_m   = 0. 
       WRITE(6,*)'****Meridional net mom flux********************'
       WRITE(6,*)'STARTING VALUES : ', Bt_mer, c_0_mer, off_m
    ENDIF

    IF(dir==1) in_zon=0
    IF(dir==2) in_mer=0
    IF(dir==1) THEN
       c_0 =c_0_zon
       Bt  = Bt_zon
       off = off_z
    ENDIF

    IF(dir==2) THEN
       c_0 = c_0_mer 
       Bt  = Bt_mer
       off = off_m
    ENDIF

    thetap(:) = 0.
    theta(:)  = 0.  
        
    DO n = 1, ht_dim-2
       
       ! potential temperature
       
       theta(n)    = temp1d(n)*(100000./pres(n))**0.286
       theta(n+1)   = temp1d(n+1)*(100000./pres(n+1))**0.286
       theta(n+2)  = temp1d(n+2)*(100000./pres(n+2))**0.286
       
       dz(n) = z(n+1) - z(n)
       
     !  WRITE(6,*) n, dz(n)
    
    ENDDO
    
    dz(ht_dim-1) = z(ht_dim)-z(ht_dim-1)

    DO n = ht_dim-2, ht_dim
       theta(n)    = temp1d(n)*(100000./pres(n))**0.286
    ENDDO

    ! Calculate the potential temperature gradient

    thetap(1) = (theta(2)-theta(1)) / dz(1)
    
    PTEMP : DO n = 2, ht_dim-1
       
       thetap(n) = (theta(n+1)-theta(n-1)) / (dz(n)+dz(n-1))
       
    ENDDO  PTEMP
    
    ! set the upper boundary gradient of the pot. temp to the one below
    thetap(ht_dim) = thetap(ht_dim-1)
    
    ! Buoyancy frequency
    
    bf(:) = 0.
    
    BRUNT : DO n = 1, ht_dim
       
       ! Brunt Vaisala Frequency
       bf(n) = ((g/theta(n))*(thetap(n)))**(0.5)
       
       IF (bf(n).LT.0.) bf(n) = 0.02
       
       !  WRITE(6,*) z(n), bf(n), theta(n), thetap(n)
       
    ENDDO BRUNT


    ! WRITE(6,*) 'Passed 2' 

    !   Lower boundary set up
    u0   = u(iz0)
    N0   = bf(iz0)
    r0   = rho(iz0)
    H0   = -(z(iz0+1)-z(iz0))/LOG(rho(iz0+1)/r0)
    flg0 = 1
    
    do n = iz0, ht_dim
       gwf(n) = 0.
    enddo

    !cc   = u0*(1-flag)
    cmax = (nc-1)*dc/2.

    fac  = k/2./N0
    k2   = k*k
    alp2 = 1./(4.e6*(H0*H0))
    omc  = sqrt((N0*N0*k2)/(k2+alp2))
    rflg = 1
    Bsum = 0.
    rmsk = 0
    L2   = alog(2.)

    !   Loop through phase speeds
    PHASE: DO i = 1, nc
       
       msk(i) = 1
       cmsk   = 1
       c0j    = ((i-1)*dc-cmax)+off
       c0(i)  = c0j

      !  WRITE(6,*) c0(i)

       ci     = c0j-u0
       c      = c0j*flag+ci*(1-flag)
       if (ci.eq.0) cmsk = 0
       sgn  = cmsk*ci/(abs(ci)+(1-cmsk))
       B0j  = sgn*(Bw*exp(-L2*((c-c_0)/cw)**2))
       ! B0j  = sgn*(Bw*exp(-1.*((c-c_0)/cw)**2))
       B0(i)=B0j


       if (sgn.le.0) j0 = i
       if ((c0j.gt.u0).and.(rflg.gt.0)) rflg=0
       if (c0j.eq.u0) msk(i)=0
       if ((abs(c0j-u0)*k).lt.omc) then
          rmsk = 1
          icm(rflg+1) = i
          rflg = 0
          
       else
          rmsk = 0
       endif
       msk(i) = msk(i)*rmsk
       if (msk(i).eq.1) then
          Foc = B0j/(c0j-u0)**3
          if ((Foc-fac).ge.0.0) msk(i) = 0
       endif
       Bsum = Bsum+abs(B0j)
    enddo PHASE

      ! STOP


    IF (Bsum.eq.0.0) THEN
       print *,'Zero flux input'
       return
    ENDIF
    
    cj0 = (j0-1)*dc-cmax
    if (cj0.eq.u0) flg0 = 0
    
    !     Intermittency factor
    eps  = Bt/(r0*Bsum)
    !     Conversion factor for /km/s to /m/day
    dfac = 24.*3600./1.e3
    !      print *,'Intermittency=',eps,' dfac=',dfac  
    
    maxalt=0

    gwh(:) = 0._prcn
    eddy(:) = 0._prcn

    DO n  = iz0+1, ht_dim-1

       ! Don't bother above altcutoff. mjh
       IF(z(n).GE.altcutoff .AND. maxalt==0) THEN
          gwf(n) = 0.
          gwh(n) = 0.
          eddy(n) = 0.
          maxalt = n
          CYCLE
       ENDIF

       fm = 0.
       fe = 0.
      
       Hb(n) = -dz(n)/LOG(rho(n)/rho(n-1))
       rb   = rho(n)

       rbh = sqrt(rb*rho(n-1))
       Nb  = bf(n)
       ub  = u(n)
       ubm = u(n-1)
       ubp = u(n+1)
       fac = rb/r0*k/2./Nb
       alp2 = 1./(4e6*Hb(n)*Hb(n))
       omc = sqrt((Nb*Nb*k2)/(k2+alp2))
       rflg=1
       
       !       Phase speed loop                                                                           
       do i = icm(2),icm(1)
          if (msk(i).eq.1) then
             B0j = B0(i)
             c0j = c0(i)
             if ((c0j.gt.ub).and.(rflg.gt.0)) then
                ict(rflg+1) = j0+flg0
                rflg = 0
             endif
             if ((abs(c0j-ub)*k).lt.omc) then
                rmsk        = 1
                ict(rflg+1) = i
                rflg        = 0
             else
                rmsk = 0
             endif
          endif
          msk(i) = msk(i)*rmsk

          if (((c0j.LE.ubm).AND.(c0j.GE.ubp)).OR. &
               ((c0j.GE.ubm).AND.(c0j.LE.ubp)).OR.(c0j.EQ.ub)) msk(i) = 0

          if (msk(i).eq.1) then
             Foc = B0j/(c0j-ub)**3
             if((Foc-fac.ge.0).or.((c0j-u0)*(c0j-ub).le.0))msk(i) = 0
             fm = fm+(1-msk(i))*B0j
             fe = fe+(1-msk(i))*(c0j-ub)*B0j
          endif
       enddo
       icm(1) = min(ict(1),icm(1))
       icm(2) = max(ict(2),icm(2))
       
            
       !  Compute the force and eddy diffusion coeff.  
       !  (Eddy diffusion here is simply energy dissipation rate over  
       !  buoyancy frequency squared.)                                    

       gwf(n)   =(dfac*r0/rbh)*fm*eps/dz(n)
       gwf(n-1) =(gwf(n-1)+gwf(n))/2.

       gwh(n)   = ((dfac*r0/rbh)*fe*eps/dz(n))/cp1d(n)
       gwh(n-1) = (gwh(n-1)+gwh(n))/2.

       eddy(n) = 0.5*(dfac*r0/rbh)*fe*eps/(dz(n)*Nb*Nb)
       eddy(n-1) = (eddy(n) + eddy(n-1))/2.

       !       WRITE(6,'(I3,5F12.3)') n, ub, ubm, ubp, dz(n), dzb

    enddo

    ! Smooth and add background eddy diffusion
    DO n = maxalt-1, 2 ,-1
       ! Smooth eddy diffusion and gwave drag
       eddy2(n) = 0.25*(eddy(n-1) + 2.*(eddy(n)) + eddy(n+1))
       gwf2(n)  = 0.25*(gwf(n-1) + 2.*(gwf(n)) + gwf(n+1))
       gwh2(n)  = 0.25*(gwh(n-1) + 2.*(gwh(n)) + gwh(n+1))
    ENDDO

    DO n = maxalt-1, 2 ,-1
       eddy(n) = eddy2(n)
       ! Add background eddy diffusion
       IF(z(n) < 100.e3) eddy(n) = eddy(n) + 1.*EXP(-1.*(100.e3-z(n))/30.e3)
       gwf(n) = gwf2(n)
       gwh(n) = gwh2(n)
    ENDDO

    RETURN
  END SUBROUTINE AlexDGwave

END MODULE MESSY_GWAVE_AD
