module messy_clamssedi_hetero_shi
  
contains

  !****************************************************************************
  ! part from from hetero_shi / subroutine liquid  (jug, 05/2015)
  ! Numbers for one air parcel:
  ! input:
  !   hno3 : HNO3 in [m^3/m^3] 
  !   h2o  : H2O  in [m^3/m^3] 
  ! output: parthno3 is used for the following calculation:
  !   HNO3(g) = HNO3(total) * parthno3
  !   HNO3(c) = HNO3(total) * (1 - parthno3)
  !****************************************************************************  
  subroutine calc_parthno3(press, temp, hno3, h2o, aer_h2so4, parthno3)
    
    use messy_clams_global,       only: prec, k_boltzmann, pi
    use messy_main_constants_mem, only: R_gas
    
    implicit none
    
    real(kind=prec) :: press                   ! in hPa
    real(kind=prec), intent(in) :: temp        ! in K
    real(kind=prec), intent(in) :: hno3, h2o   ! volume mixing ratio
    real(kind=prec), intent(in) :: aer_h2so4   ! 
    real(kind=prec), intent(out) :: parthno3 
    real(kind=prec) :: chno3, ch2o             ! concentration in cm^-3
    real(kind=prec) :: tnd                     ! air molecule density in cm^-3

    real(kind=prec)  :: ppbh2so4      != 0.2         ! eq. volume mixing ratio [ppbv]
    real(kind=prec)  :: ms, msb, mn, mnb
    
    
    real(kind=prec), parameter :: ctoa = 7.336E+21 
    
    
    real(kind=prec) :: lnp, nsul,  ph2o, tice, t, &
         pn0, a, b, c, xs, xn, x, tt, &
         hns, hnn, c0, d, alpha, beta, xx, phi, vpn
    
    
    real(kind=prec), dimension(10) :: hnnarr = (/ 0.1457341419E+02, 0.6159942877E-01, &
         -0.1148954003E+01, 0.6916933619E+00, -0.9886352167E-01, 0.5157939752E-02, &
         0.1234728041E+00, -0.1155743833E+00, 0.1101132109E-01, &
         0.9791396551E-02/)  
    
    real(kind=prec), dimension(10) :: hnsarr = (/ 0.1446998319E+02, 0.6387958319E-01, &
         -0.3295968441E+01, 0.1778224331E+01, -0.2232437780E+00, 0.8648607708E-02, &
         0.5366954950E+00, -0.3351643646E+00, 0.2651532224E-01, &
         0.1575503129E-01/)  
    
    real(kind=prec), dimension(7) ::  mnbarr = (/ 0.1985306995E+03, -0.1194830927E+05, &
         -0.2846953556E+02, -0.3913646516E+02, 0.8328785711E+02, 0.6358493057E+04, &
         -0.1764983192E+05/)  
    
    real(kind=prec), dimension(7) ::  msbarr = (/ 0.4700412633E+02, -0.6969007247E+04, &
         -0.4618369994E+01, -0.2166107102E+02, 0.5181581907E+02, 0.2724212325E+04, &
         -0.1573208296E+05/)  
    
    ppbh2so4 = aer_h2so4

    ! convert k_boltzmann (m^3 -> cm^3)
    ! total number desity: tnd in cm^-3
    tnd = press * 100 / ( k_boltzmann*1.0e6 * temp )
    ! determine concentration of h2o and hno3 [molec/ccm]
    ch2o = h2o * tnd
    chno3 = hno3 * tnd
    
    ! ======================================================================
    ! PREVENT PH2O OUTSIDE LIMITS OF MODEL
    ph2o = ch2o / ctoa * temp 
    ph2o = max(2.0D-8, ph2o) 
    ph2o = min(2.0D-5, ph2o) 
    
    ! ======================================================================
    ! PREVENT T<185 K AND T<TICE-3 K, AS EXPECTED BY MODEL
    tice = 2668.7 / (10.431 - (log(ph2o) + log(760.0))/log(10.0)) 
    if (temp <= tice - 3) then 
       t = tice - 3.0 
    else if (temp > 240.0) then 
       t = 240.0 
    else 
       t = temp 
    endif
    t = max(185.D0, t) 
    
    ! PARTIAL PRESSURE IN ATM FOR H2O & HNO3
    nsul = ppbh2so4 * 1.E-9 * press * 100./ R_gas / t  !moles/m3 sulphate when pure liq 
    pn0 = chno3 / ctoa * t 

    lnp = log(ph2o)  ! CALCULATE MSB, MOLALITY OF H2SO4 IN BINARY
    a = msbarr(5) + msbarr(7)/t 
    b = msbarr(4) + msbarr(6)/t 
    c = msbarr(1) + msbarr(2)/t + msbarr(3)*log(t) - lnp 
    xs = ((-b) - sqrt(abs(b**2 - 4.0*a*c)))/(2.0*a) 
    msb = 55.51*xs/(1.0 - xs) 
    if (t <= 210.0) then 
       ! DON'T CALCULATE HNO3 UPTAKE FOR T>210K (TOO SMALL)
       ! CALCULATE MNB, MOLALITY OF HNO3 IN BINARY
       a = mnbarr(5) + mnbarr(7)/t 
       b = mnbarr(4) + mnbarr(6)/t 
       c = mnbarr(1) + mnbarr(2)/t + mnbarr(3)*log(t) - lnp 
       xn = ((-b) - sqrt(abs(b**2 - 4.0*a*c)))/(2.0*a) 
       mnb = 55.51*xn/(1.0 - xn) 
       ! CALCULATE HNS, HENRY'S LAW CONSTANT OF HNO3 IN H2SO4-H2O
       x = lnp + 18.4 
       tt = 1.E+4*(1.0/t - 1.0/230.0) 
       hns = exp(hnsarr(1)+hnsarr(2)*tt**2+(hnsarr(3)+hnsarr(4)*tt+hnsarr(5)*&
            tt**2+hnsarr(6)*tt**3)*x+(hnsarr(7)+hnsarr(8)*tt+hnsarr(9)*tt**2)*&
            x**2+(hnsarr(10)*tt)*x**3) 
       ! CALCULATE HNN, HENRY'S LAW CONSTANT OF HNO3 IN HNO3-H2O
       hnn = exp(hnnarr(1)+hnnarr(2)*tt**2+(hnnarr(3)+hnnarr(4)*tt+hnnarr(5)*&
            tt**2+hnnarr(6)*tt**3)*x+(hnnarr(7)+hnnarr(8)*tt+hnnarr(9)*tt**2)*&
            x**2+(hnnarr(10)*tt)*x**3) 
       ! NOW ANALYTICAL EXPRESSION FOR COMPOSITION OF HNO3-H2SO4-H2O AEROSOL
       ! Factor 1E-5 to convert Pa to bar
       c0 = R_gas * 1E-5 * t * nsul 
       d = (c0*hnn*mnb*msb**2)/(mnb - msb) 
       c = msb*((-2.0*c0*hnn*mnb) + c0*hns*msb + mnb*msb - hnn*msb*pn0)/(mnb-msb)
       b = (c0*hnn*mnb**2 - c0*hns*mnb*msb - 2.0*mnb**2*msb + mnb*msb**2 + &
            hnn*mnb*msb*pn0 - hns*msb**2*pn0)/(mnb**2 - mnb*msb) 
       alpha = (-2.0*b**3) + 9.0*b*c - 27.0*d 
       beta = sqrt(4.0*(b**2 - 3.0*c)**3 - alpha**2) 
       xx = beta/alpha 
       
       phi = atan(xx) 
       
       if (phi < 0.) phi = phi + pi 
       
       ms = -1.0/3.0*(b + 2.0*sqrt(b**2 - 3.0*c)*cos((1.0*pi + phi)/3.0)) 
       mn = mnb*(1.0 - ms/msb) 
       ! let mn not be negative! (jug, 26.06.2002)
       mn = max(mn,0.d0)
       ! ======================
       ! NOW THE OUTPUT VARIABLES
       vpn = mn/(hnn*mn/(mn + ms) + hns*ms/(mn + ms)) 
       
       parthno3 = 1.0 - (pn0 - vpn)/pn0 
       
       ! jug  line added to stay in valid range of parthno3 (J.U. Grooss, 26.10.98)
       parthno3 = min(parthno3, 1.D0) 
    else 
       parthno3 = 1.0
    endif
    
  end subroutine calc_parthno3
  
end module messy_clamssedi_hetero_shi
