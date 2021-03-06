MODULE messy_clamscirrus 
! **************************************************************************
! MODULE FOR CLaMS CIRRUS MODULE
! **************************************************************************

  USE messy_clams_global, ONLY: DP, dnparts_max

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamscirrus'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  PUBLIC :: cirrus
  PUBLIC :: clamscirrus_read_nml
  PRIVATE :: allocate_arrays, deallocate_arrays
  PRIVATE :: freeze_out_at_100percent, freeze_out_with_supersaturation
  PRIVATE :: densice_fit_Ni
  PRIVATE :: compute_mean_free_path, compute_slip_correction
  PRIVATE :: compute_terminal_velocity

!----------- 
CONTAINS
!----------- 

!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!
!! cirrus.f90 --- Liest SH ein und berechnet Saettigungsdampfdruck ueber Eis.
!! Hieraus ergibt sich das zugehoerige Wasserdampfsaettigungsmischungs-        
!! verhaeltnis, das mit dem momentanen verglichen wird. Hieraus resultiert die
!! Menge Wasser, die ausfriert und damit der Gasphase verlorengeht.
!! Author          : Gebhard Guenther
!! Created On      : Tue Nov 11 12:21:05 2003
!! Last Modified By: c.hoppe
!! Last Modified On: Fri Jan 20 17:08:04 2012
!! Update Count    : 945
!! Status          : Unknown, Use with caution!

subroutine cirrus (lat_array,lon_array,theta_array,temp_array,pressure_array,pv_array, &
                   h2o_array, h2o_100_array, iwc, iwc_100, clwc, delta_t, dnparts_local)

  use messy_clamscirrus_global

  implicit none

  real(DP) :: lat_array(dnparts_max),lon_array(dnparts_max),theta_array(dnparts_max)
  real(DP) :: temp_array(dnparts_max),pressure_array(dnparts_max),pv_array(dnparts_max)
  real(DP) :: h2o_array(dnparts_max)
  real(DP) :: h2o_100_array(dnparts_max), iwc_100(dnparts_max)
  real(DP) :: clwc(dnparts_max), iwc(dnparts_max)
  real(DP) :: delta_t ![s]
  integer  :: dnparts_local

  integer :: itime

  pi=4*atan(1.)

  !   Starten der Stopuhr

  call system_clock (count=start)

  !!   use_sh: Use Specific Humidity [kg/kg] instead of Volume Mixing Ratio (For
  !!           compatibility reasons only, supply SH instead of H2O

  if (rank==0) then
     write (*,*)
     write (*,*) "START of CIRRUS"
     write (*,*)
  end if

  ! Temperatur-Offset:
  delta_temp = 0.

  call allocate_arrays
  ! Initialize arrays:
  densice_fit = 0.
  t_rad = 0.
  t_massendichte = 0.
  t_masse = 0.
  t_vol = 0.
  mean_free_path = 0.
  slip_correction = 0.
  terminal_velocity = 0.
  sed_laenge = 0.

  if (rank==0) then
     print*,'Achtung: Temperatur-Offset:', delta_temp
     print*,'Achtung: Zeitschritt :', delta_t
     if (troposphere) print*,'Cirrus is restricted to the Troposphere.'
     if (freeze_out_type==1) then
        write (*,*) 'Ausfrieren auf Schwellwert'
     else
        write (*,*) 'Ausfrieren auf ISMR'
     endif
  endif
     
  ! use trajectory for dehydration 

  !! Achtung: Temperatur-Offset
  temp_array = temp_array+delta_temp  
  
  ! Berechnungen fuer Standardausfrieren bei 100%
     
  call freeze_out_at_100percent (temp_array, pressure_array, &
       theta_array, pv_array, &
       iwc_100, h2o_100_array, iwc, delta_t, dnparts_local)
     
  ! Berechnungen fuer Martinas Parametrisierung
  
  call freeze_out_with_supersaturation (temp_array, pressure_array, &
       theta_array, pv_array, &
       h2o_array, iwc, clwc, delta_t, dnparts_local)

  ! Aufraeumen
  call deallocate_arrays

  !   Stoppen der Stopuhr
  call system_clock (COUNT_RATE = clock(2),COUNT_MAX  = clock(3))
  ! if (rank==0) then
  !    print *, "System clock runs at ", clock(2)," ticks per second"
  !    print *, " "
  ! endif
  call system_clock ( COUNT=finish )
  if ( finish .ge. start ) then
     seconds = float(finish-start)/float(clock(2))
  else
     seconds = 0.0
  end if
  if (rank==0) write(6,333) seconds
  
333 format("This job has taken ",f10.3, " seconds to execute!")

end subroutine cirrus

  !*****************************************************************************
  ! Allokieren der globalen Felder
  !*****************************************************************************
  subroutine allocate_arrays

    USE messy_clamscirrus_global

    implicit none

    allocate(densice_fit(dnparts_max))
    allocate(t_rad(dnparts_max))

    allocate(t_massendichte(dnparts_max),t_masse(dnparts_max),t_vol(dnparts_max))
    allocate(mean_free_path(dnparts_max))
    allocate(slip_correction(dnparts_max))
    allocate(terminal_velocity(dnparts_max))
    allocate(sed_laenge(dnparts_max))
    
  end subroutine allocate_arrays

  !*****************************************************************************
  ! Deallokieren der globalen Felder
  !*****************************************************************************
  subroutine deallocate_arrays

    USE messy_clamscirrus_global
   
    implicit none

    deallocate(densice_fit)
    deallocate(t_rad)

    deallocate(t_massendichte, t_masse, t_vol)
    deallocate(mean_free_path, slip_correction, terminal_velocity, sed_laenge)
    
  end subroutine deallocate_arrays

  !*****************************************************************************
  !  Berechnungen fuer Standardausfrieren bei 100%
  !*****************************************************************************
  subroutine freeze_out_at_100percent (temp_array, pressure_array, &
       theta_array, pv_array, &
       iwc_100, sh_100_array, iwc, delta_t, dnparts_local)

    USE messy_clamscirrus_global

    implicit none

    real(DP) :: temp_array(dnparts_max),pressure_array(dnparts_max)
    real(DP) :: theta_array(dnparts_max), pv_array(dnparts_max)
    real(DP) :: iwc_100(dnparts_max), sh_100_array(dnparts_max)
    real(DP) :: iwc(dnparts_max)
    real(DP) :: delta_t
    integer  :: dnparts_local

    real(DP), dimension(:), allocatable :: wasser_mv_100
    real(DP), dimension(:), allocatable :: e_sat_ice
    real(DP), dimension(:), allocatable :: ismr
    integer :: i
    real(DP) :: iwc_temp,clwc_temp,iwc_0,delta_iwc !Hilfsvariablen fuer die Verdunstung
    real(DP) :: uebersaettigung,sh_to_mv
!    real(DP) :: uebersaettigung_ueber_wasser

    allocate(wasser_mv_100(dnparts_local))
    allocate(e_sat_ice(dnparts_local))
    allocate(ismr(dnparts_local))

    !! Umrechnungsfaktor SH nach Volumenmischungsverhaeltnis: Verhaeltnis der
    !!      molaren Massen (28.9644/18.015)

    if (use_sh) then
       sh_to_mv=masse_luft/masse_h2o
    else
       sh_to_mv=1.
    endif

    wasser_mv_100=0.  !-1.e30
    e_sat_ice=0.      !-1.e30
    ismr=0.           !-1.e30


    ! Uebersaettigung fest eingestellt: 130% laut Alex Mangold, 110% sind aber im globalen
    ! Mittel wahrscheinlicher. 
     
    uebersaettigung=1.0 
!    uebersaettigung_ueber_wasser=1.0 
     

    do i=1, dnparts_local

!!!!! ACHTUNG: Erweiterung von Ines:
    if ( ((troposphere) .and. (theta_array(i) < 380.) .and. (pv_array(i) < 2.)) .or. &
         (.not. troposphere)) then

       iwc_0 = iwc_100(i) 
       delta_iwc = 0.     

       if (abs((sh_100_array(i)-mdi)/mdi) > eps  .and. &
           abs((temp_array(i)-mdi)/mdi) > eps  .and. &
           abs((pressure_array(i)-mdi)/mdi) > eps ) then
          
          !! Umrechnung von spez. Feuchte in Mischungsverhaeltnis bei use_sh=1
          
!!$!         if (use_sh) then
!!$!            wasser_mv_100(i)=(28.9644/18.015)*sh_100_array(i)
!!$!         else
!!$!            wasser_mv_100(i)=sh_100_array(i)  
!!$!         endif
          
          wasser_mv_100(i)=sh_100_array(i)*sh_to_mv

          !! Berechnung des Saettigungsdampfdrucks ueber Eis in hPa -> "/100"
          
          e_sat_ice(i)=(10**((-2663.5/temp_array(i))+12.537))/100.
          
          !! Berechnung des zugehoerigen Wasserdampfmischungsverhaeltnisses
          !! nach Wiederhold
          
          if((pressure_array(i)-e_sat_ice(i)) .eq. 0.) then   
             ! .or. pressure_array(i) .le. e_sat_ice(i)
             !!! For the upper stratosphere
             ismr(i)=-1.
          else
             ismr(i)=(e_sat_ice(i)/(pressure_array(i)-e_sat_ice(i)))
          endif
       endif
       
       !! Ueberpruefen, ob aktuelles Wasser ueber dem zulaessigen Wert liegt
       !! Fuer Wasser, Umrechnung in spez. Feuchte
       
!!$!       if (wasser_mv(i) >= uebersaettigung_ueber_wasser*wsmr(i)) then
!!$!          clwc(i)=clwc(i)+(wasser_mv(i)-uebersaettigung_ueber_wasser*wsmr(i))&
!!$!               &/sh_to_mv
!!$!          wasser_mv(i)=uebersaettigung_ueber_wasser*wsmr(i)
!!$!       endif
!!$!    

       !! Fuer Eis, Umrechnung in spez. Feuchte
       if ((wasser_mv_100(i) >= uebersaettigung*ismr(i)) .and. (ismr(i) .gt. 0)) then
          
          select case (freeze_out_type)
             
          case (1)
             
             !! Ausfrieren auf den Schwellwert
             delta_iwc       = (wasser_mv_100(i) - uebersaettigung*ismr(i))/sh_to_mv
             iwc_100(i)=iwc_100(i)+(wasser_mv_100(i)-uebersaettigung*ismr(i))/sh_to_mv
             wasser_mv_100(i)=uebersaettigung*ismr(i)

          case default
             
             !! Ausfrieren auf ISMR
             delta_iwc = (wasser_mv_100(i)-ismr(i))/sh_to_mv
             iwc_100(i)=iwc_100(i)+(wasser_mv_100(i)-ismr(i))/sh_to_mv
             wasser_mv_100(i)=ismr(i)

          end select

      endif

       !! Sedimentieren: Umrechnung von iwc in Vol.-MV, mit Masse_Wasser pro
       !!  mol, Ergebnis Masse in kg/m^3
     
       if (abs((temp_array(i)-mdi)/mdi) > eps) then
          t_massendichte(i)=(iwc_100(i)*sh_to_mv*(pressure_array(i)*100.*masse_h2o)&
               &/(k_boltzmann*temp_array(i)*avogadro))!*1.e-3 ! [kg/m^3]
        
          !! Hier werden die Teilchenzahldichten explizit ueber verschiedene Fits
          !! berechnet. 
        
! ############################################ FP - 130111

          select case (type_ice_fit)
          case (0)
             !! Alte Parametrisierung der Eisteilchenzahldichte (im klimatolog. Lauf)
             densice_fit(i) = (10**((-11.7563+3.*0.159595) &
                               +(0.0547315+3.*0.00077228)*temp_array(i)))*1.e6 
          case (1)
          !! Berechnung nach den neuesten Parametrisierungen 20110708
!!$!         Ni_median cm^-3
!!$!         a0= -28.7694 ; - 15 flights
!!$!         a1= 0.998214
!!$!         a2= 18.2202
!!$!         Ni_median_T = 10^(a0*a1^Tempgrid + a2) 

!!$!         densice_fit(i)=(10**(-28.7694*0.998214**(temp_array(i))+18.2202))*1.e6
             densice_fit(i)=densice_fit_Ni(-28.7694_DP,0.998214_DP,18.2202_DP,temp_array(i),1.e6_DP)
          end select

!!$!         Ni_max
!!$!         a0= -1.73736e+14 ; - 15 flights
!!$!         a1= 0.837851
!!$!         a2= 0.586513
!!$!         Ni_max_T = 10^(a0*a1^Tempgrid + a2)
!!$!         densice_fit(i)=(10**(-1.73736e+14*0.837851**(temp_array(i))+0.586513))*1.e6
!!$!         densice_fit(i)=densice_fit_Ni(-1.73736e+14,0.837851,0.586513,temp_array(i),1.e6)

! ############################################ FP - 130111

          t_masse(i)=t_massendichte(i)/densice_fit(i)
        
          !! Hier werden die Teilchenzahldichten explizit gesetzt.
        
!!$!     t_masse(i)=t_massendichte(i)/densice ! [kg]
        
          t_vol(i)=t_masse(i)/rhoice ! [m^3]
          t_rad(i)=((3.*t_vol(i))/(4.*pi))**(1./3.) ![m]
          
          !print*,minval(pressure_array),minval(temp_array)

          !! Berechnung der Radien direkt aus Martinas Parametrisierung

!!$!         t_rad(i)=10**(-60777.5*0.938722**temp_array(i)+1.36289)*1.e-6
!!$!         t_vol(i)=(4.*pi/3.)*t_rad(i)**3
!!$!         t_masse(i)=t_vol(i)*rhoice
!!$!         densice_fit(i)=t_massendichte(i)/t_masse(i)

       endif
     
       if (abs((pressure_array(i)-mdi)/mdi) > eps  .and. &
           abs((temp_array(i)-mdi)/mdi) > eps) then
!!$!         mean_free_path(i)=(2.*viscosity)/(pressure_array(i)*100.*sqrt(8.*masse_luft/&
!!$!              &(pi*universal_gas_constant*temp_array(i))))
          mean_free_path(i)=compute_mean_free_path(pressure_array(i),temp_array(i),100._DP)
       else
          mean_free_path(i)=mdi
       endif
       
       if (t_rad(i) > 0. .and. mean_free_path(i) > 0. ) then
          
!!$!         slip_correction(i)=1.+(mean_free_path(i)/t_rad(i))*(1.257+.4&
!!$!              &*exp(-1.1*t_rad(i)/mean_free_path(i)))
          
          slip_correction(i)=compute_slip_correction(mean_free_path(i),t_rad(i))
          
!!$!         terminal_velocity(i) = -(4.*t_rad(i)*t_rad(i)&
!!$!              & *rhoice*9.81*slip_correction(i))/(18.*viscosity) ! [m/s]

          terminal_velocity(i) = compute_terminal_velocity(t_rad(i),slip_correction(i))
          
          sed_laenge(i)=terminal_velocity(i)*delta_t ![m]
          
       end if
     
       !! Verminderung des Eiswassergehalts um die sedimentierte Komponente

! ############################################ FP - 130110
       if (characteristic_height_100 > 0) then
          if (.not. use_traj) then  
             iwc_100(i) = iwc_0 + delta_iwc*(1.-min(1.,-sed_laenge(i)/characteristic_height_100))
          else
             iwc_100(i)=iwc_100(i)*(1.-min(1.,-sed_laenge(i)/characteristic_height_100))
          endif

       else
          !! Komplettes Aussedimentieren
          iwc_100(i)=0.  

       endif
! ############################################ FP - 130110

       !! Evaporation von iwc, wenn wasser_mv kleiner 100% Saettigung
     
       !! Fuer Eis
     
       if ((wasser_mv_100(i) < ismr(i)) .and. (ismr(i) .gt. 0)) then
          iwc_temp=iwc(i)
          iwc_100(i)=iwc_100(i)-min(iwc_100(i),(ismr(i)-wasser_mv_100(i))/sh_to_mv)
          wasser_mv_100(i)=wasser_mv_100(i)+min(iwc_temp*sh_to_mv,(ismr(i)-wasser_mv_100(i)))
       endif
     
       !! Zurueckspeichern auf spez. Feuchte
       
       if (abs((wasser_mv_100(i)-mdi)/mdi) > eps) then
             sh_100_array(i)=wasser_mv_100(i)/sh_to_mv
       endif

!!!!! ACHTUNG: Erweiterung von Ines:
    endif

    enddo

    deallocate (wasser_mv_100)
    deallocate (e_sat_ice)
    deallocate (ismr)

  end subroutine freeze_out_at_100percent

  !*****************************************************************************
  !  Berechnungen fuer Martinas Parametrisierung
  !*****************************************************************************
  subroutine freeze_out_with_supersaturation (temp_array, pressure_array, &
       theta_array, pv_array, &
       sh_array, iwc, clwc, delta_t, dnparts_local)

    USE messy_clamscirrus_global

    implicit none

    real(DP) :: temp_array(dnparts_max),pressure_array(dnparts_max)
    real(DP) :: theta_array(dnparts_max), pv_array(dnparts_max)
    real(DP) :: sh_array(dnparts_max)
    real(DP) :: iwc(dnparts_max), clwc(dnparts_max)
    real(DP) :: delta_t
    integer  :: dnparts_local

    real(DP), dimension(:), allocatable :: wasser_mv
    real(DP), dimension(:), allocatable :: e_sat_water
    real(DP), dimension(:), allocatable :: e_sat_ice
    real(DP), dimension(:), allocatable :: ismr
    real(DP), dimension(:), allocatable :: wsmr
    integer :: i

    real(DP) :: iwc_temp,clwc_temp,iwc_0,delta_iwc !Hilfsvariablen fuer die Verdunstung
    real(DP) :: uebersaettigung,uebersaettigung_ueber_wasser,sh_to_mv

    allocate(wasser_mv(dnparts_local))
    allocate(e_sat_water(dnparts_local))     ! Neu fuer Wasserwolken
    allocate(e_sat_ice(dnparts_local))
    allocate(ismr(dnparts_local))
    allocate(wsmr(dnparts_local))            ! Neu fuer Wasserwolken

    !! Umrechnungsfaktor SH nach Volumenmischungsverhaeltnis: Verhaeltnis der
    !!      molaren Massen (28.9644/18.015)

    if (use_sh) then
       sh_to_mv=masse_luft/masse_h2o
    else
       sh_to_mv=1.
    endif

    wasser_mv=0.      !-1.e30
    e_sat_water=0.
    e_sat_ice=0.      !-1.e30
    ismr=0.
    wsmr=0.

    uebersaettigung_ueber_wasser=1.0 
    uebersaettigung=1.0

    do i=1, dnparts_local


!!!!! ACHTUNG: Erweiterung von Ines:
    if ( ((troposphere) .and. (theta_array(i) < 380.) .and. (pv_array(i) < 2)) .or. &
         (.not. troposphere)) then

       iwc_0 = iwc(i)
       delta_iwc = 0. 

       ! Uebersaettigung nach RH_ice_freeze_het=212.-0.4*temp (laut Martina Kraemer
       ! abgeleitet aus AIDA Messungen)

! ############################################ FP - 130110
       if (sat_crit > 0) then
          uebersaettigung = sat_crit
       else
          select case (int(sat_crit))
          case (-1)
             uebersaettigung = max(1.0,(212.-0.4*temp_array(i))/100.)
          case (-2)
             uebersaettigung = max(1.1,(200.-0.4*temp_array(i))/100.) ! Versuch
          case (-3)
             ! Uebersaettigung nach RH_ice_freeze_hom=238.-0.4*temp (laut Martina Kraemer
             ! abgeleitet nach Koop)
             uebersaettigung = max(1.1,(238.-0.4*temp_array(i))/100.)
          case default
             uebersaettigung = max(1.0,(212.-0.4*temp_array(i))/100.)
          end select
       endif
! ############################################ FP - 130110
      
       if (abs((sh_array(i)-mdi)/mdi) > eps  .and. &
           abs((temp_array(i)-mdi)/mdi) > eps  .and. &
           abs((pressure_array(i)-mdi)/mdi) > eps ) then

          !! Umrechnung von spez. Feuchte in Mischungsverhaeltnis wenn use_sh=1
          
             wasser_mv(i)=sh_array(i)*sh_to_mv

          ! Uebersaettigung ueber Wasser nach Wiederholt in hPa 
          ! (siehe auch Fishwiki-Formelsammlung): 
          ! esat_w = (1.0007+3.46*10-6*press)*611.21*exp(17.502*temp)/(240.97
          !          +temp))/100
          ! Achtung: Temperatur in C
          
          e_sat_water(i)=(1.0007+3.46e-6*pressure_array(i))*611.21&
               &*exp((17.502*(temp_array(i)-273.16))/(240.97+temp_array(i)-273.16))/100
          
          ! Nach Sonntag 1994 in hPa
          
!!$!          e_sat_water(i)=(exp(-6096.9585/temp_array(i) + 16.635794 - &
!!$!               & 2.711193*1.e-2 * temp_array(i) + &
!!$!               & 1.673952*1.e-5*temp_array(i)*temp_array(i) + &
!!$!               & 2.433502*log(temp_array(i))))
          
          !! Berechnung des zugehoerigen Wasserdampfmischungsverhaeltnisses 
          
          if((pressure_array(i)-e_sat_water(i)) .eq. 0.) then 
             wsmr(i)=-1.
          else
             wsmr(i)=e_sat_water(i)/(pressure_array(i)-e_sat_water(i))
          endif
          
          !! Berechnung des Saettigungsdampfdrucks ueber Eis in hPa -> "/100"
          
          e_sat_ice(i)=(10**((-2663.5/temp_array(i))+12.537))/100.
          
          !! Berechnung des zugehoerigen Wasserdampfmischungsverhaeltnisses
          !! nach Wiederhold
          
          if((pressure_array(i)-e_sat_ice(i)) .eq. 0.) then               
             ! .or. pressure_array(i) .le. e_sat_ice(i)
             !!! For the upper stratosphere
             ismr(i)=-1.
          else
             ismr(i)=(e_sat_ice(i)/(pressure_array(i)-e_sat_ice(i)))
          endif
       endif


       !! Ueberpruefen, ob aktuelles Wasser ueber dem zulaessigen Wert liegt
       !! Fuer Wasser, Umrechnung in spez. Feuchte
       
       if (wasser_mv(i) >= uebersaettigung_ueber_wasser*wsmr(i) .and.&
            & (temp_array(i) > 273.16) .and. (wsmr(i) .gt. 0)) then
          clwc(i)=clwc(i)+(wasser_mv(i)-uebersaettigung_ueber_wasser*wsmr(i))&
               &/sh_to_mv
          
!!$!           clwc(i)=0.
        
          wasser_mv(i)=uebersaettigung_ueber_wasser*wsmr(i)

       endif

       !! Fuer Eis, Umrechnung in spez. Feuchte
       if (wasser_mv(i) >= uebersaettigung*ismr(i) .and. (ismr(i) .gt. 0)) then !$
          !             & .and. (temp_array(i) <= 273.16)) then

          select case (freeze_out_type)
             
          case (1)
             
             !! Ausfrieren auf den Schwellwert
             delta_iwc       = (wasser_mv(i) - uebersaettigung*ismr(i))/sh_to_mv
             iwc(i)=iwc(i)+(wasser_mv(i)-uebersaettigung*ismr(i))/sh_to_mv
             wasser_mv(i)=uebersaettigung*ismr(i)

          case default
             
             !! Ausfrieren auf ISMR
             delta_iwc = (wasser_mv(i)-ismr(i))/sh_to_mv
             iwc(i)=iwc(i)+(wasser_mv(i)-ismr(i))/sh_to_mv
             wasser_mv(i)=ismr(i)

          end select
          
       endif

       !! Sedimentieren: Umrechnung von iwc in Vol.-MV, mit Masse_Wasser pro mol,
       !!  Ergebnis Masse in kg/m^3
       
       if (abs((temp_array(i)-mdi)/mdi) > eps ) then
          t_massendichte(i)=(iwc(i)*sh_to_mv*(pressure_array(i)*100.&
               &*masse_h2o)&
               &/(k_boltzmann*temp_array(i)*avogadro))!*1.e-3 ! [kg/m^3]
          
          !! Hier werden die Teilchenzahldichten explizit ueber verschiedene Fits
          !! berechnet. 
          
! ############################################ FP - 130111

          select case (type_ice_fit)

          case (0)
             !! Alte Parametrisierung der Eisteilchenzahldichte (im klimatolog. Lauf)
             densice_fit(i) = (10**((-11.7563+3.*0.159595) &
                              +(0.0547315+3.*0.00077228)*temp_array(i)))*1.e6 
          case (1)
             !! Berechnung nach den neuesten Parametrisierungen 20110708
!!$!          Ni_median cm^-3
!!$!          a0= -28.7694 ; - 15 flights
!!$!          a1= 0.998214
!!$!          a2= 18.2202
!!$!          Ni_median_T = 10^(a0*a1^Tempgrid + a2) 

!!$!          densice_fit(i)=(10**(-28.7694*0.998214**(temp_array(i))+18.2202))*1.e6
             densice_fit(i)=densice_fit_Ni(-28.7694_DP,0.998214_DP,18.2202_DP,temp_array(i),1.e6_DP)
          end select

!!$!         Ni_max
!!$!         a0= -1.73736e+14 ; - 15 flights
!!$!         a1= 0.837851
!!$!         a2= 0.586513
!!$!         Ni_max_T = 10^(a0*a1^Tempgrid + a2)
!!$!         densice_fit(i)=(10**(-1.73736e+14*0.837851**(temp_array(i))+0.586513))*1.e6
!!$!         densice_fit(i)=densice_fit_Ni(-1.73736e+14,0.837851,0.586513,temp_array(i),1.e6)

! ############################################ FP - 130111
          
          t_masse(i)=t_massendichte(i)/densice_fit(i)
          
          !! Hier werden die Teilchenzahldichten explizit gesetzt.
          
!!$!           t_masse(i)=t_massendichte(i)/densice ! [kg]
          
          t_vol(i)=t_masse(i)/rhoice ! [m^3]
          t_rad(i)=((3.*t_vol(i))/(4.*pi))**(1./3.) ![m]
          

          !--------------------------------------------------------------------
          ! -- check for const. particle size
          !t_rad(i)=10.e-6               !!! -- FP - 101020
          !--------------------------------------------------------------------


          !print*,minval(pressure_array),minval(temp_array)
          
          !! Berechnung der Radien direkt aus Martinas Parametrisierung

!!$!         t_rad(i)=10**(-60777.5*0.938722**temp_array(i)+1.36289)*1.e-6
!!$!         t_vol(i)=(4.*pi/3.)*t_rad(i)**3
!!$!         t_masse(i)=t_vol(i)*rhoice
!!$!         densice_fit(i)=t_massendichte(i)/t_masse(i)

       endif
       
       if (abs((pressure_array(i)-mdi)/mdi) > eps .and. &
           abs((temp_array(i)-mdi)/mdi) > eps ) then
!!$!          mean_free_path(i)=(2.*viscosity)/(pressure_array(i)*100.*sqrt(8.*masse_luft/&
!!$!               &(pi*universal_gas_constant*temp_array(i))))
          mean_free_path(i)=compute_mean_free_path(pressure_array(i),temp_array(i),100._DP)
       else
          mean_free_path(i)=mdi
       endif
       
       if (t_rad(i) > 0. .and. mean_free_path(i) > 0. ) then
          
!!$!          slip_correction(i)=1.+(mean_free_path(i)/t_rad(i))*(1.257+.4&
!!$!               &*exp(-1.1*t_rad(i)/mean_free_path(i)))

          slip_correction(i)=compute_slip_correction(mean_free_path(i),t_rad(i))
          
!!$!          terminal_velocity(i) = -(4.*t_rad(i)*t_rad(i)&
!!$!               & *rhoice*9.81*slip_correction(i))/(18.*viscosity) ! [m/s]

          terminal_velocity(i) = compute_terminal_velocity(t_rad(i),slip_correction(i))

          sed_laenge(i)=terminal_velocity(i)*delta_t ![m]
          
       end if
  
       !! Verminderung des Eiswassergehalts um die sedimentierte Komponente
     
! ############################################ FP - 130110
       if (characteristic_height > 0) then
          if (.not. use_traj) then  
             iwc(i) = iwc_0 + delta_iwc*(1.-min(1.,-sed_laenge(i)/characteristic_height))
          else
             iwc(i)=iwc(i)*(1.-min(1.,-sed_laenge(i)/characteristic_height))
          endif

       else
          !! Komplettes Aussedimentieren
          iwc(i)=0.  
       endif
! ############################################ FP - 130110

       !! Evaporation von iwc, wenn wasser_mv kleiner 100% Saettigung
       
       !! Fuer Eis
       
       if ((wasser_mv(i) < ismr(i)) .and. (ismr(i) .gt. 0.)) then
          iwc_temp=iwc(i)
          iwc(i)=iwc(i)-min(iwc(i),(ismr(i) - wasser_mv(i))/sh_to_mv)
          wasser_mv(i)=wasser_mv(i)+min(iwc_temp*sh_to_mv,(ismr(i) - wasser_mv(i)))
       endif
       
       !! Fuer Wasser
       
! ############################################ FP - 130110
!!!!! sollte auskommentiert sein ?!
!!$!! ------------------- comment out for calc. without CLWC ---------------------
!!$       if ((wasser_mv(i) < wsmr(i)) .and. (wsmr(i) .gt. 0.)) then
!!$          clwc_temp=clwc(i)
!!$          clwc(i)=clwc(i)-min(clwc(i),(wsmr(i) - wasser_mv(i))/sh_to_mv)
!!$          wasser_mv(i)=wasser_mv(i)+min(clwc_temp*sh_to_mv,(wsmr(i) - wasser_mv(i)))
!!$       endif
!!$!! ----------------------------------------------------------------
     
!!$!       do iii=1,nparts
!!$!          if (wasser_mv(iii) >= uebersaettigung*ismr(iii)) then
!!$!             iwc(iii)=iwc(iii)+(wasser_mv(iii)-uebersaettigung*ismr(iii))
!!$!             wasser_mv(iii)=uebersaettigung*ismr(iii)
!!$!          end if
!!$!
!!$!          if (wasser_mv(iii) < ismr(iii)) then
!!$!             iwc(iii)=iwc(iii)-min(iwc(iii),(ismr(iii) - wasser_mv(iii)))
!!$!             wasser_mv(iii)=wasser_mv(iii)+min(iwc(iii),(ismr(iii) - wasser_mv(iii)))
!!$!          end if
!!$!       end do
! ############################################ FP - 130110
     
       !! Zurueckspeichern auf spez. Feuchte
       
       if (abs((wasser_mv(i)-mdi)/mdi) > eps) then

          sh_array(i)=wasser_mv(i)/sh_to_mv

       endif

!!$! print*,minval(ismr),pressure_array(minloc(ismr)),maxval(ismr)&
!!$!      &,pressure_array(maxloc(ismr))
!!$! print*,minval(ismr),temp_array(minloc(ismr)),maxval(ismr)&
!!$!      &,temp_array(maxloc(ismr))
   

!!!!! ACHTUNG: Erweiterung von Ines:
    endif       

    enddo

    deallocate(wasser_mv)
    deallocate(e_sat_water)
    deallocate(e_sat_ice)
    deallocate (ismr)
    deallocate (wsmr)

  end subroutine freeze_out_with_supersaturation

! --+----+----+----+----+---- Computing Fit for Ice Particle Density +--

  real(DP) function densice_fit_Ni(a0,a1,a2,temp,scale_factor)

    implicit none
    real(DP), intent(in) :: a0,a1,a2,scale_factor
    real(DP), intent(in) :: temp

    densice_fit_Ni=(10**(a0*a1**(temp)+a2))*scale_factor

    return
  end function densice_fit_Ni

! --+----+----+----+----+----+----+----+--- Computing Mean Free Path +--

  real(DP) function compute_mean_free_path(pressure,temp,scale_factor)

    USE messy_clamscirrus_global

    implicit none
    real(DP), intent(in) :: scale_factor
    real(DP), intent(in) :: pressure,temp

    compute_mean_free_path=(2.*viscosity)/(pressure*scale_factor*sqrt(8.*masse_luft/&
         &(pi*universal_gas_constant*temp)))

    return
  end function compute_mean_free_path

! --+----+----+----+----+----+----+----+-- Computing Slip Correction +--

  real(DP) function compute_slip_correction(mean_free_path_length,radius)

    implicit none
    real(DP), intent(in) :: mean_free_path_length,radius

    compute_slip_correction = 1.+(mean_free_path_length/radius)*(1.257+.4&
         &*exp(-1.1*radius/mean_free_path_length))

    return
  end function compute_slip_correction

! --+----+----+----+----+----+- Computing Terminal Settling velocity +--

  real(DP) function compute_terminal_velocity(radius,slip_corr)

    USE messy_clamscirrus_global

    implicit none
    real(DP), intent(in) :: slip_corr,radius

    compute_terminal_velocity = -(4.*radius*radius&
         & *rhoice*9.81*slip_corr)/(18.*viscosity) ! [m/s]

    return
  end function compute_terminal_velocity
!---------------------------------------------------------------------------

!---------------------------------------------------------------------------
  SUBROUTINE clamscirrus_read_nml(status, iou)

    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_clamscirrus_global, ONLY: timestep_cirrus, type_ice_fit, sat_crit, &
                                        characteristic_height, characteristic_height_100, &
                                        troposphere, rank, freeze_out_type

    IMPLICIT NONE
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamscirrus_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output

    NAMELIST /CTRL/ timestep_cirrus, type_ice_fit, sat_crit, &
                    characteristic_height, characteristic_height_100, &
                    troposphere, freeze_out_type

    status = 1 !ERROR

    if (rank==0) then
       l_print = .true.
    else
       l_print = .false.
    endif

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr, l_print)
    
    status = 0 !NO ERROR

  END SUBROUTINE clamscirrus_read_nml

! **************************************************************************
END MODULE messy_clamscirrus
! **************************************************************************
