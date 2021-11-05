!******************************************************************************
!
! Modul clamssedi_compute
!
! Written by Gebhard Guenther, Nicole Thomas, Thomas Breuer
! IEK-7, Forschungszentrum Juelich
! Last Modified By: N.Thomas
! Last Modified On: Fri Apr 24 10:19:30 2015
!
! subroutine interpolate_background_chem (ipart)
! function phmhno3 (t, ph2o)
! function phmh2o (t)
! subroutine calc_nucrate(ipart)
! function pmkh2oliq (t)
! function pmkh2oice (t)
! subroutine compute_particle_growth(ipart, status)
! function d_eff_hno3(r,t,p)
! function d_eff_h2o(r,t,p)
! function d_hno3(t,p)
! function d_h2o(t,p)
! function rho_nat(t)
! function c_hno3(t)
! function c_h2o(t)
! subroutine compute_settling_velocity(ipart)
! 
!******************************************************************************

Module messy_clamssedi_compute

  use messy_clams_global, only: prec
  
  real(prec), dimension(3)   :: weight_down, weight_up
  real(prec)                 :: vert_down, vert_up

CONTAINS
  
  !****************************************************************************
  ! determine the current H2O- and HNO3-background
  !****************************************************************************
  subroutine interpolate_background_chem (ipart, ipart_local)
    
    use messy_clamssedi_global, only: particles, nairparcels50_lev, &
                                      nairparcels50_lev_init, airparcels, &
                                      triangles, r_coarse, r_high, & 
                                      lat_down, lat_up, lev_window
    
    use messy_clamssedi_triang,     only: set_particle_kart
    use messy_clams_global,         only: prec, radius_earth
    use messy_clams_tools_interpol, only: determ_weight, vert_interpol
    
    implicit none

    integer, intent(in)        :: ipart, ipart_local
    integer                    :: i, ilev, lev_ind_down, lev_ind_up
    real(prec), dimension(3)   :: particle_kart
    real(prec), dimension(3,3) :: vertex_coords
    integer, dimension(3)      :: ap_ids_down, ap_ids_up
    real(prec), dimension(3)   :: lev_arr
    real(prec)                 :: h2o_down, h2o_up, hno3_down, hno3_up
    real(prec)                 :: icebin_down, icebin_up, natbin_down, natbin_up
    real(prec)                 :: lev_nearest_down, lev_nearest_up
    real(prec)                 :: delta_lev, ratio

    lev_ind_down = particles%lev_down(ipart_local) 
    lev_ind_up = particles%lev_up(ipart_local)
    
    particles%airparcel_density_change(ipart) = &
         real(nairparcels50_lev(lev_ind_down))/real(nairparcels50_lev_init(lev_ind_down))
        
    call set_particle_kart(particles%lon(ipart), particles%lat(ipart), particle_kart)

    !!! only calc_type==3

    if (particles%lat(ipart) < lat_down .or. &
         particles%lat(ipart) > lat_up) then
       ratio = r_coarse/radius_earth
    else
       ratio = r_high/radius_earth
    end if
    
    ! airparcel indices of the triangle vertices on the level below and above the particle
    ap_ids_down = triangles%airparcel_indices(particles%tr_ind_down(ipart_local),particles%lev_down(ipart_local),:)
    ap_ids_up = triangles%airparcel_indices(particles%tr_ind_up(ipart_local),particles%lev_up(ipart_local),:)

    ! determine weights for the level below the particle ipart
    do i=1,3
       vertex_coords(:,i) = airparcels%coor(:,ap_ids_down(i))
    enddo
    lev_arr = airparcels%lev(ap_ids_down)

    ilev = airparcels%ilev(ap_ids_down(1))
    delta_lev = lev_window(2,ilev) - lev_window(1,ilev)

    call determ_weight (vertex_coords, lev_arr, particle_kart, particles%lev(ipart), &
                        weight_down, lev_nearest_down, delta_lev, ratio, 2)
    
    ! determine H20 and HNO3 for the level below the current particle
    h2o_down  = sum(weight_down(:) * airparcels%h2o(ap_ids_down))
    hno3_down = sum(weight_down(:) * airparcels%hno3(ap_ids_down))

    ! determine weights for the level above the current particle ipart
    do i=1,3
       vertex_coords(:,i) = airparcels%coor(:,ap_ids_up(i))
    enddo
    lev_arr = airparcels%lev(ap_ids_up)

    ilev = airparcels%ilev(ap_ids_up(1))
    delta_lev = lev_window(2,ilev) - lev_window(1,ilev)

    call determ_weight (vertex_coords, lev_arr, particle_kart, particles%lev(ipart), &
                        weight_up, lev_nearest_up, delta_lev, ratio, 2)

    ! determine H20 and HNO3 for the level above the current particle ipart
    h2o_up  = sum(weight_up(:) * airparcels%h2o(ap_ids_up))
    hno3_up = sum(weight_up(:) * airparcels%hno3(ap_ids_up))
    
    ! linear interpolation in lev-coordinates for H2O
    call vert_interpol(lev_nearest_down, lev_nearest_up, particles%lev(ipart), &
         h2o_down, h2o_up, particles%h2o(ipart), vert_up, vert_down, .true.)
    ! linear interpolation in lev-coordinates for HNO3
    call vert_interpol(lev_nearest_down, lev_nearest_up, particles%lev(ipart), &
         hno3_down, hno3_up, particles%hno3(ipart), vert_up, vert_down, .true.)
    
    ! linear interpolation in lev-coordinates for ice- and natbin
    ! (should be executed only onetime, once the particle is created)
    if (particles%particle_id(ipart) < 0 .AND. &
         particles%icebin(ipart) == 0 .AND. particles%natbin(ipart) == 0 ) then

       icebin_down = sum(weight_down(:) * airparcels%icebin(ap_ids_down))
       natbin_down = sum(weight_down(:) * airparcels%natbin(ap_ids_down))
       icebin_up = sum(weight_up(:) * airparcels%icebin(ap_ids_up))
       natbin_up = sum(weight_up(:) * airparcels%natbin(ap_ids_up))

       call vert_interpol (lev_nearest_down, lev_nearest_up, particles%lev(ipart), &
            icebin_down, icebin_up, particles%icebin(ipart), vert_up, vert_down, .true.)
       call vert_interpol (lev_nearest_down, lev_nearest_up, particles%lev(ipart), &
            natbin_down, natbin_up, particles%natbin(ipart), vert_up, vert_down, .true.)

       if (particles%icebin(ipart) < 0) then
          write(*,*) 'interpolate_background_chem: ipart, particles%icebin(ipart)= ', ipart, particles%icebin(ipart)
          write(*,*) 'interpolate_background_chem: weights=', weight_down(:),weight_up(:)
          write(*,*) 'interpolate_background_chem: icebins=', airparcels%icebin(ap_ids_down),airparcels%icebin(ap_ids_up)
       endif
    endif
    
  end subroutine interpolate_background_chem


  !****************************************************************************
  ! Berechnung des Dampfdruckes von HNO3 ueber NAT nach Hanson and Mauersberger
  ! Siehe hierzu auch hetero_ken_6.f90!
  !****************************************************************************
  function phmhno3 (t, ph2o)

    use messy_clams_global, only: prec

    implicit none

    real(prec) :: phmhno3
    real(prec) , intent(in) :: t    ! Temperatur [K]
    real(prec) , intent(in) :: ph2o ! Partialdampfdruck von H2O [atm]
    real(prec) :: dum1, dum2 

    dum1 = (-2.7836) - 0.00088*t
    dum2 = 38.9855 - 11397.0/t + 0.009179*t
    phmhno3 = 10.0**(dum1*log10(ph2o*760.0) + dum2)/760.0 ! [atm]

    return

  end function phmhno3

  !****************************************************************************
  ! Berechnung des Dampfdruckes von H2O ueber Eis nach Hanson and Mauersberger
  ! Siehe hierzu auch hetero_ken_6.f90!
  !****************************************************************************
  function phmh2o (t)

    use messy_clams_global, only: prec

    implicit none

    real(prec) :: phmh2o
    real(prec) , intent(in) :: t    ! Temperatur [K]

    phmh2o =  10**(10.431 - 2668.7/t)/760.0 ! [atm]

    return

   end function phmh2o
   
   !****************************************************************************
   ! Calculation of the nucleation rate
   ! nucrate1 = NAT; nucrate2 = ice
   !****************************************************************************
   subroutine calc_nucrate(ipart, ipart_local)
     
     use messy_clamssedi_global, only: particles, densnat_val, sice_table, &
                                       timestep_sedi, ntbins, xnice_table, &
                                       nbins, snat_table, xnnat_table, factor, &
                                       nat_tstep, ice_tstep
     use messy_clams_global,     only: prec, pi
     USE messy_clams_global,     ONLY: HOUR, MINUTE

     implicit none

     integer, intent(in) :: ipart, ipart_local
     integer             :: ibin, tbin
     real(prec)          :: nucrate1, nucrate2, sice_crit, snat_crit
     real(prec)          :: vwater, vice, ph2ov, dg, awice, awliq, x 
     real(prec)          :: dum_homo, xjnuc, xn
     real(prec)          :: icebin0, natbin0
     real(prec)          :: snatval, siceval, tempval

     snatval = 0.
     siceval = 0.
     tempval = 0.


     if (ipart <= 0) write(*,*) 'calc_nucrate warning 1: ipart=',ipart 
     if (ipart > 400000) write(*,*) 'calc_nucrate warning 2: ipart=',ipart 
     if (particles%icebin(ipart)< 0) write(*,*) 'calc_nucrate warning 3: ipart, icebin(ipart)=',ipart, particles%icebin(ipart)
     if (particles%icebin(ipart).ne. particles%icebin(ipart)) write(*,*) 'calc_nucrate warning 4: ipart, icebin(ipart)=',ipart,particles%icebin(ipart)

     icebin0 = floor(particles%icebin(ipart))
     natbin0 = floor(particles%natbin(ipart))

     nucrate1 = 0.
     nucrate2 = 0.
     
     ! Ice nucleation
     ! Every 24 hours
     if ( (ice_tstep .EQ. 24) .AND. HOUR == 12 .and. MINUTE == 0 ) then
        siceval = particles%sicemax(ipart)
        tempval = particles%tmin(ipart)

        if (particles%sicemax(ipart) < particles%sice(ipart)) &
             & siceval = particles%sice(ipart)
        if (particles%tmin(ipart) > particles%temperature(ipart)) &
             & tempval = particles%temperature(ipart)
        particles%sicemax(ipart) = 0.
     ! Every hour
    elseif (ice_tstep .NE. 24) then 
        siceval = particles%sice(ipart)
        tempval = particles%temperature(ipart)
    endif 

    if (siceval > 1.05 .AND. densnat_val <= -2) then
        
        ! Homogeneous ice nucleation resulting in an ice nuber density of 10 #/cm3
        ! Barrier for sice calculated following the parametrisation of Kärcher und Lohmann (2002)
        if (siceval > (2.583 - tempval/207.83)) then
           
           nucrate2 = 10
           particles%icebin(ipart) = 0.
           
        ! Heterogeneous ice nucleation
        elseif (siceval < maxval(sice_table)) then   
           tbin = floor(tempval) - 169.
           tbin = max(tbin, 1)
           tbin = min(tbin, ntbins-1)
           
           ibin = floor(particles%icebin(ipart))
           ibin = max(ibin, 1)
           
           ! active contact angle bins for nucleation
           do
              sice_crit = sice_table(ibin, tbin)
              
              if (siceval < sice_crit .or. ibin >= nbins ) exit
              
              nucrate2 = nucrate2 + xnice_table(ibin)
              ibin = ibin + 1
              
           enddo
           
           if (ibin > 1) particles%icebin(ipart) = real(ibin)
           
        endif
     endif

     ! NAT nucleation
     ! Every 24 hours
     if ( (nat_tstep .EQ. 24) .AND. HOUR == 12 .and. MINUTE == 0 ) then
        snatval = particles%snatmax(ipart)
        tempval = particles%tmin(ipart)
        if (particles%snatmax(ipart) < particles%snat(ipart)) &
             & snatval = particles%snat(ipart)
        if (particles%tmin(ipart) > particles%temperature(ipart)) &
             & tempval = particles%temperature(ipart)
        particles%snatmax(ipart) = 0.
        particles%tmin(ipart) = 0.
     ! Every hour
     elseif (nat_tstep .NE. 24) then  
        snatval = particles%snat(ipart)
        tempval = particles%temperature(ipart)
     endif

     if (snatval > 1.05 .AND. (densnat_val == -1 .OR. densnat_val == -3)) then
        
        tbin = floor(tempval) - 179.
        tbin = max(tbin, 1)
        tbin = min(tbin, ntbins-1)
        
        ibin = floor(particles%natbin(ipart))
        ibin = max(ibin,1)
        
        ! active contact angle bins for nucleation
        do
           snat_crit = snat_table(ibin, tbin)
           
           if (snatval < snat_crit .or. ibin >= nbins ) exit
           
           nucrate1 = nucrate1 + xnnat_table(ibin)
           ibin = ibin + 1
           
        enddo
        
        if (ibin > 1) particles%natbin(ipart) = real(ibin)
        
     endif
     
     ! nucleation rate, units: #particles cm^-3 day^-1
     if (nucrate2 > 1.e-3) then
        
        particles%density(ipart) = nucrate2 * (1./factor) &
             * particles%airparcel_density_change(ipart) * 1E6
        particles%class(ipart) = 21.
        if (nucrate2 == 10.) particles%class(ipart) = 20.
        particles%natbin(ipart) = natbin0
        
     elseif (nucrate1 > 1.e-5) then
        
        particles%density(ipart) = nucrate1 * (1./factor) &
             * particles%airparcel_density_change(ipart) * 1E6
        particles%class(ipart) = 11.
        particles%icebin(ipart) = icebin0
        
     else
        particles%density(ipart) = 0.
        particles%radius(ipart) = 0.
        particles%class(ipart) = 0.
        particles%natbin(ipart) = natbin0
        particles%icebin(ipart) = icebin0
     endif

     if (particles%icebin(ipart) < 0) then
          write(*,*) 'calc_nucrate: ipart, particles%icebin(ipart)= ', ipart, particles%icebin(ipart)
     endif
     
     particles%icebin_diff(ipart_local) = max(particles%icebin(ipart) - icebin0, 0.) * & 
          (1./factor) * particles%airparcel_density_change(ipart)
     particles%natbin_diff(ipart_local) = max(particles%natbin(ipart) - natbin0, 0.) * &
          (1./factor) * particles%airparcel_density_change(ipart)

  end subroutine calc_nucrate
  
  !****************************************************************************
  ! Berechnung des Dampfdruckes von H2O ueber Wasser nach Murphy and Koop
  ! (ph2oliq)
  !****************************************************************************
  function pmkh2oliq (t)
    
    use messy_clams_global, only: prec
    
    implicit none
     
    real(prec)             :: pmkh2oliq
    real(prec), intent(in) :: t    ! Temperatur [K]
    real                   :: pmkh2oliq_Pa
    
    pmkh2oliq_Pa = exp(54.842763 - 6763.22/t - 4.210*log(t) + 0.000367*t + &
         tanh(0.0415*(t - 218.8))*(53.878 - 1331.22/t - 9.44523*log(t) + &
         0.014025*t)) ! [Pa]
    pmkh2oliq    = pmkh2oliq_Pa * 9.8692e-6 ! [atm]
    
    return
    
  end function pmkh2oliq
  
  !****************************************************************************
  ! Berechnung des Dampfdruckes von H2O ueber Eis nach Murphy and Koop
  ! (ph2oice)
  !****************************************************************************
  function pmkh2oice (t)
    
    use messy_clams_global, only: prec
    
    implicit none
    
    real(prec)             :: pmkh2oice
    real(prec), intent(in) :: t    ! Temperatur [K]
    real                   :: pmkh2oice_Pa
    
    pmkh2oice_Pa = exp(9.550426 - 5723.265/t + 3.53068*log(t) - 0.00728332*t) ![Pa]
    pmkh2oice    = pmkh2oice_Pa * 9.8692e-6 ! [atm]
    
    return
    
  end function pmkh2oice
  
  !****************************************************************************
  ! calculation of the particle growth using the method of Carslaw et al.
  !****************************************************************************
  subroutine compute_particle_growth(ipart, ipart_local, status)
    
    use messy_clamssedi_global, only: particles, masse_nat, timestep_sedi, &
                                      triangles, densnat_val, airparcels, &
                                      hno3_minval, h2o_minval, lev_start, &
                                      factor
    use messy_clams_global,     only: prec, k_boltzmann, masse_luft, avogadro, &
                                      masse_h2o, rho_ice, pi, rank, ntasks
    use messy_main_constants_mem, only: R_gas
    
    implicit none
    
    integer, intent(in) :: ipart, ipart_local
    integer             :: status, k
    integer             :: ap_id
    real(prec)          :: density, density_m
    real(prec)          :: ch2o, chno3, ph2o, phno3
    real(prec)          :: h2o_ice, hno3_nat, max_h2o, max_hno3
    real(prec)          :: dum1, dum2, dum3
    real(prec)          :: growth, densnat, densice
    real(prec)          :: mass, volumen, volumen_old, radius_old
    logical             :: completely_molten
    real                :: amount_hno3, amount_h2o 
    real(prec)          :: eq_mix_ratio, number_molecules
    integer             :: max_nnat, part_nat_down, part_nat_up
    integer             :: max_nice, part_ice_up, part_ice_down
    
    ! number density of air parcels [#/m^3]. pressure [hPa].
    density = particles%pressure(ipart) * 100. / (k_boltzmann * particles%temperature(ipart))
    ! mass density of air parcels [kg/m^3]
    density_m = (particles%pressure(ipart) * 100. / &
         (k_boltzmann * particles%temperature(ipart))) * masse_luft / avogadro
    ! calculation of the concentrations [molec/m^3]
    ch2o = particles%h2o(ipart) * density_m * avogadro / masse_luft
    chno3 = particles%hno3(ipart) * density_m * avogadro / masse_luft
    ! calculate the partial pressure of H2O and HNO3 in [atm]
    ph2o = (particles%h2o(ipart) * particles%temperature(ipart) * density_m &
         & * R_gas / masse_luft) * 9.87e-6
    phno3 = (particles%hno3(ipart) * particles%temperature(ipart) * density_m &
         & * R_gas / masse_luft) * 9.87e-6
    
    ! calculate the differenc of pressure back into [Pa]
    ! check for nat (class == 1x) or ice (class == 2x)
    if ( (particles%class(ipart) >= 10.) .AND. (particles%class(ipart) < 20.) ) then
       dum1 = d_eff_hno3(particles%radius(ipart), particles%temperature(ipart), &
            particles%pressure(ipart)) * masse_nat
       dum2 = rho_nat(particles%temperature(ipart)) * R_gas * &
            particles%temperature(ipart)
       dum3 = (phno3 - phmhno3(particles%temperature(ipart), ph2o)) * 1.01325e5
       growth = (dum1 / dum2) * dum3
    elseif ( particles%class(ipart) >= 20. ) then
       dum1 = d_eff_h2o(particles%radius(ipart), particles%temperature(ipart), &
            particles%pressure(ipart)) * masse_h2o
       dum2 = rho_ice * R_gas * particles%temperature(ipart)
       dum3 = (ph2o - pmkh2oice(particles%temperature(ipart))) * 1.01325e5
       growth = (dum1 / dum2) * dum3
    else
       growth = 0.
    endif
    
    volumen_old = (4. * pi * particles%radius(ipart)**3) / 3.
    radius_old = particles%radius(ipart)
    
    completely_molten = .false.
    
    ! if conditions for growth are fulfilled: calculate new radius
    if (particles%radius(ipart)**2 + 2 * growth * timestep_sedi > 0) then
       particles%radius(ipart) = sqrt(particles%radius(ipart)**2 + 2 * growth * timestep_sedi)
       
    ! vaporize (Verdampfen)
    else
       particles%radius(ipart) = 0.
       completely_molten = .true.
    endif
    
    ! NAT
    if ( (particles%class(ipart) >= 10.) .AND. (particles%class(ipart) < 20.)  ) then
       
       ! NAT = 1*HNO3 and 3*H2O
       amount_hno3 = 1.
       amount_h2o  = 3.
       
       densnat = particles%density(ipart) ! 100.
       volumen = (4. * pi * particles%radius(ipart)**3) / 3.
       mass = (volumen - volumen_old) * rho_nat(particles%temperature(ipart)) * densnat
       number_molecules = mass / (masse_nat / avogadro)
       eq_mix_ratio = number_molecules * masse_luft / (density_m * avogadro)
       
       ! Das aequivalente Mischungsverhaeltnis darf nicht groesser als HNO3_NAT
       ! sein, wobei das zur Verfuegung stehende HNO3 gleichmaessig auf die
       ! Partikel verteilt, die sich innerhalb einer Triangel befinden.
       ! Das aequivalente Mischungsverhaeltnis darf 25% des gesamten HNO3-Wertes pro
       ! Zeitschritt nicht ueberschreiten.
       ! (Pressure hPa -> atm)

       hno3_nat = (phno3 - phmhno3(particles%temperature(ipart), ph2o)) &
            & / (particles%pressure(ipart) * 9.8692e-4)

       ! maximum number of nat particles in the triangle above or below the current particle
       part_nat_up = triangles%nparticles_nat(particles%tr_ind_up(ipart_local), &
                                              particles%lev_up(ipart_local))
       part_nat_down = triangles%nparticles_nat(particles%tr_ind_down(ipart_local), &
                                                particles%lev_down(ipart_local))
       max_nnat = max(part_nat_up, part_nat_down)
       
       ! Verfuegbares HNO3 aufgeteilt auf die Anzahl NAT Teilchen und korrigiert
       ! mit der c_box_density
       max_hno3 = hno3_nat / max_nnat / particles%airparcel_density_change(ipart)
       
       ! growth
       if ( (eq_mix_ratio > 0.) .AND. (eq_mix_ratio > max_hno3) ) then
          eq_mix_ratio = min(eq_mix_ratio, max_hno3, 0.25 * particles%hno3(ipart))
          number_molecules = eq_mix_ratio * density_m * avogadro / masse_luft
          mass = number_molecules * (masse_nat / avogadro)
          volumen = volumen_old + mass / (rho_nat(particles%temperature(ipart)) * densnat)
          particles%radius(ipart) = (3. * volumen/(4. * pi))**(1./3.)
       endif
       
       ! Wenn der NAT-Rock vollkommen aufschmilzt, wird das gesamte alte Volumen
       ! zurueck an die C-Boxen gegeben
       if (completely_molten .eqv. .true.) then
          mass = -volumen_old * rho_nat(particles%temperature(ipart)) * densnat
          number_molecules = mass / (masse_nat / avogadro)
          eq_mix_ratio = number_molecules * masse_luft / (density_m * avogadro)
       endif
       
       ! Beruecksichtigung der C-Box-Dichte beim Fluss
       eq_mix_ratio = eq_mix_ratio * particles%airparcel_density_change(ipart)

       particles%hno3(ipart) = particles%hno3(ipart) - amount_hno3 * eq_mix_ratio
       particles%h2o(ipart) = particles%h2o(ipart) - amount_h2o * eq_mix_ratio

       if (particles%hno3(ipart) < -1E-11 ) then
          write(*,*)'Error in compute: ipart: ',ipart, ' hno3 ',particles%hno3(ipart)
          write(*,*)' eq_mix_ratio', eq_mix_ratio,' radius ',particles%radius(ipart)
          write(*,*)' volumen, volumen_old:', volumen,volumen_old
          write(*,*) dum1,dum2,dum3,growth
          status = 1
          return
       endif
       
    ! ice   
    elseif ( particles%class(ipart) >= 20. ) then
       
       ! Eis = 0*HNO3 und 1*H2O
       amount_hno3 = 0.
       amount_h2o  = 1.
       
       densice = particles%density(ipart) !100.

       volumen = (4. * pi * particles%radius(ipart)**3) / 3.
       mass = (volumen - volumen_old) * rho_ice * densice
       number_molecules = mass / (masse_h2o / avogadro)
       eq_mix_ratio = number_molecules * masse_luft / (density_m * avogadro)
       
       ! Das aequivalente Mischungsverhaeltnis darf nicht groesser als H2O_ice
       ! sein, wobei das zur Verfuegung stehende H2O gleichmaessig auf die
       ! Partikel verteilt, die sich innerhalb einer Triangel befinden.
       ! Das aequivalente Mischungsverhaeltnis darf 25% des gesamten H2O-Wertes pro
       ! Zeitschritt nicht ueberschreiten.
       ! (Pressure hPa -> atm)

       h2o_ice = (ph2o - pmkh2oice(particles%temperature(ipart))) &
            & / (particles%pressure(ipart) * 9.8692e-4)
       
       ! maximum number of ice particles in the triangle above or below the current particle
       part_ice_up = triangles%nparticles_ice(particles%tr_ind_up(ipart_local), &
                                              particles%lev_up(ipart_local))
       part_ice_down = triangles%nparticles_ice(particles%tr_ind_down(ipart_local), &
                                                particles%lev_down(ipart_local))
       max_nice = max(part_ice_down, part_ice_up)
       
       ! Verfuegbares H2O aufgeteilt auf die Anzahl Eis Teilchen und korrigiert
       ! mit der c_box_density
       max_h2o = h2o_ice / max_nice / particles%airparcel_density_change(ipart)

       ! growth
       if ( (eq_mix_ratio > 0.) .AND. (eq_mix_ratio > max_h2o) ) then
          eq_mix_ratio = min(eq_mix_ratio, max_h2o, 0.5*particles%h2o(ipart))
          number_molecules = eq_mix_ratio * density_m * avogadro / masse_luft
          mass = number_molecules * (masse_h2o / avogadro)
          volumen = volumen_old + mass / (rho_ice * densice)
          particles%radius(ipart) = (3. * volumen/(4. * pi))**(1./3.)
       endif

       ! Wenn das Eis vollkommen aufschmilzt, wird das gesamte alte Volumen
       ! zurueck an die C-Boxen gegeben
       if (completely_molten .eqv. .true.) then
          mass = -volumen_old * rho_ice * densice
          number_molecules = mass / (masse_h2o / avogadro)
          eq_mix_ratio = number_molecules * masse_luft / (density_m * avogadro)
       endif
       
       ! Beruecksichtigung der C-Box-Dichte beim Fluss
       eq_mix_ratio = eq_mix_ratio * particles%airparcel_density_change(ipart)
       particles%hno3(ipart) = particles%hno3(ipart) - amount_hno3 * eq_mix_ratio
       particles%h2o(ipart) = particles%h2o(ipart) - amount_h2o * eq_mix_ratio
       
       ! NAT nucleation on preexisting ice particles following Luo et al. (2003)

       ! Ueberpruefe, ob die Bedingungen fuer NAT auf Eis gegeben sind:

       ph2o = (particles%h2o(ipart) * particles%temperature(ipart) * density_m &
            & * R_gas / masse_luft) * 9.87e-6
       phno3 = (particles%hno3(ipart) * particles%temperature(ipart) * density_m &
            & * R_gas / masse_luft) * 9.87e-6

       hno3_nat = (phno3 - phmhno3(particles%temperature(ipart), ph2o)) &
            & / (particles%pressure(ipart) * 9.8692e-4)

       if ((completely_molten .eqv. .true.) .and. (densnat_val == -3) .and. &
            (particles%lev(ipart) .gt. lev_start) .and. (hno3_nat > 0) .and. &
            (particles%particle_id(ipart) .gt. 0) ) then

          do k = 1, 3
             ! update airparcel k of the triangle below the current particle
             ap_id = triangles%airparcel_indices(particles%tr_ind_down(ipart_local), &
                  particles%lev_down(ipart_local), k)
             
             !!airparcels%h2o(ap_id) = max(h2o_minval, &
             !!     airparcels%h2o(ap_id) - amount_h2o * eq_mix_ratio * weight_down(k) * vert_down)
             airparcels%diff_h2o(ap_id) = airparcels%diff_h2o(ap_id) + &
                  amount_h2o * eq_mix_ratio * weight_down(k) * vert_down
             
             ! update airparcel k of the triangle above the current particle
             ap_id = triangles%airparcel_indices(particles%tr_ind_up(ipart_local), &
                  particles%lev_up(ipart_local), k)
             
             !!airparcels%h2o(ap_id) = max(h2o_minval, &
             !!     airparcels%h2o(ap_id) - amount_h2o * eq_mix_ratio * weight_up(k) * vert_up)
             airparcels%diff_h2o(ap_id) = airparcels%diff_h2o(ap_id) + amount_h2o * &
                  eq_mix_ratio * weight_up(k) * vert_up
          enddo

          ! NAT = 1*HNO3 und 3*H2O
          amount_hno3 = 1.
          amount_h2o  = 3.

          particles%density(ipart) = min(particles%density(ipart) * 0.5, (1.*1d6*(1./factor)))
          densnat = particles%density(ipart) !100.
          particles%class(ipart) = 12.

          ! Verfügbares HNO3 aufgeteilt auf die Anzahl NAT Teilchen und korrigiert
          ! mit der c_box_density
          max_hno3 = hno3_nat / max_nice / particles%airparcel_density_change(ipart)
          if  (particles%hno3(ipart) <= 0) write(*,*)'warning clamssedi_compute ',ipart, particles%hno3(ipart)
          eq_mix_ratio = min(max_hno3, 0.25 * particles%hno3(ipart))
          number_molecules = eq_mix_ratio * density_m * avogadro / masse_luft
          mass = number_molecules * (masse_nat / avogadro)
          volumen = mass / (rho_nat(particles%temperature(ipart)) * densnat)
          particles%radius(ipart) = (3. * volumen/(4. * pi))**(1./3.)

          completely_molten = .false.
                
          eq_mix_ratio = eq_mix_ratio * particles%airparcel_density_change(ipart)

          particles%hno3(ipart) = particles%hno3(ipart) - amount_hno3 * eq_mix_ratio
          particles%h2o(ipart) = particles%h2o(ipart) - amount_h2o * eq_mix_ratio

       endif
       
       if (particles%h2o(ipart) < -1.e-7 ) then
          write(*,*)'Error in comupte: ipart: ',ipart,' h2o ',particles%h2o(ipart)
          write(*,*)' eq_mix_ratio', eq_mix_ratio,' radius ',particles%radius(ipart)
          write(*,*)' volumen, volumen_old:', volumen,volumen_old
          write(*,*) particles%lev(ipart)
          status = 1 
          return
       endif

    else
       amount_hno3 = 0.
       amount_h2o = 0.
       volumen = 0.
       mass = 0.
       number_molecules = 0.
       eq_mix_ratio = 0.
    endif
    
    ! particles%radius(ipart) is already set to 0.
    if (completely_molten) then
       particles%density(ipart) = 0.
    endif
    
    ! only calc_type == 3
    do k = 1, 3
       ! update airparcel k of the triangle below the current particle
       ap_id = triangles%airparcel_indices(particles%tr_ind_down(ipart_local), &
                                           particles%lev_down(ipart_local), k)
    
       airparcels%diff_h2o(ap_id) = airparcels%diff_h2o(ap_id) + &
            amount_h2o * eq_mix_ratio * weight_down(k) * vert_down
       airparcels%diff_hno3(ap_id) = airparcels%diff_hno3(ap_id) + &
            amount_hno3 * eq_mix_ratio * weight_down(k) * vert_down
       
       if (particles%particle_id(ipart) < 0) then
          airparcels%icebin_diff(ap_id) = airparcels%icebin_diff(ap_id) + &
                                          particles%icebin_diff(ipart_local) * &
                                          weight_down(k) * vert_down
          
          airparcels%natbin_diff(ap_id) = airparcels%natbin_diff(ap_id) + &
                                          particles%natbin_diff(ipart_local) * &
                                          weight_down(k) * vert_down
       endif
       
       ! if the particle has vaporized, surrounding bins are set to 0.
       if (completely_molten .eqv. .true.) then
          airparcels%natbin(ap_id) = 0.
          airparcels%icebin(ap_id) = 0.
          airparcels%natbin_diff(ap_id) = -1000. * (ntasks+1)
          airparcels%icebin_diff(ap_id) = -1000. * (ntasks+1)
       endif
       
       ! update airparcel k of the triangle above the current particle
       ap_id = triangles%airparcel_indices(particles%tr_ind_up(ipart_local), &
                                           particles%lev_up(ipart_local), k)
       
       airparcels%diff_h2o(ap_id) = airparcels%diff_h2o(ap_id) + amount_h2o * &
                                    eq_mix_ratio * weight_up(k) * vert_up
       airparcels%diff_hno3(ap_id) = airparcels%diff_hno3(ap_id) + amount_hno3 * &
                                     eq_mix_ratio * weight_up(k) * vert_up
       
       if (particles%particle_id(ipart) < 0) then
          airparcels%icebin_diff(ap_id) = airparcels%icebin_diff(ap_id) + &
                                          particles%icebin_diff(ipart_local) * &
                                          weight_up(k) * vert_up
          
          airparcels%natbin_diff(ap_id) = airparcels%natbin_diff(ap_id) + &
                                          particles%natbin_diff(ipart_local) * &
                                          weight_up(k) * vert_up
       endif
       
       ! if the particle has vaporized, surrounding bins are set to 0.  
       if (completely_molten .eqv. .true.) then
          airparcels%natbin(ap_id) = 0.
          airparcels%icebin(ap_id) = 0.
          airparcels%natbin_diff(ap_id) = -1000. * (ntasks+1)
          airparcels%icebin_diff(ap_id) = -1000. * (ntasks+1)
       endif
              
    enddo
    
  end subroutine compute_particle_growth
   
   !****************************************************************************
   ! Berechnung der effektiven Diffusion (siehe Pruppacher&Klett)
   !   ! Die Einheiten von c_hno3, r und d_** muessen passend sein!!!
   !****************************************************************************
   function d_eff_hno3(r,t,p)
     
     use messy_clams_global, only: prec
     
     implicit none
     
     real(prec)             :: d_eff_hno3
     real(prec), intent(in) :: r, t, p
     
     d_eff_hno3 = d_hno3(t,p) / (1 + 4 * d_hno3(t,p) / (c_hno3(t) * r)) ! [m**2/s]
     
     return
     
   end function d_eff_hno3
   
   !****************************************************************************
   ! Berechnung der effektiven Diffusion (siehe Pruppacher&Klett)
   ! Die Einheiten von c_hno3, r und d_** muessen passend sein!!!
   !****************************************************************************
   function d_eff_h2o(r,t,p)
     
     use messy_clams_global, only: prec
          
     implicit none
     
     real(prec)             :: d_eff_h2o
     real(prec), intent(in) :: r,t,p
     
     d_eff_h2o = d_h2o(t,p) / (1 + 4 * d_h2o(t,p) / (c_h2o(t) * r)) ! [m**2/s]
     
     return
     
   end function d_eff_h2o


   !****************************************************************************
   ! Berechnung der Diffusion von HNO3 (siehe Doktorarbeit Meilinger bzw.
   ! Khosrawi). Achtung: Nur Naeherung! Die Einheit kommt aus d_h2o.
   !****************************************************************************
   function d_hno3(t,p)
     
     use messy_clamssedi_global, only: masse_nat
     use messy_clams_global,     only: masse_h2o, prec
     
     implicit none
     
     real(prec)             :: d_hno3
     real(prec), intent(in) :: t ! temperature [K]
     real(prec), intent(in) :: p ! pressure [hPa]
     
     d_hno3 = d_h2o(t,p) * sqrt(masse_h2o/masse_nat) ! [m**2/s]
     
     return
     
   end function d_hno3

   !****************************************************************************
   ! Berechnung der Diffusion von H2O (siehe Pruppacher &Klett)
   ! Gilt fuer den Temperaturbereich [-40,...,+40] C!
   ! Der Faktor 1.e-4 sorgt fuer die Umrechnung von cm**2/s in m**2/s
   !****************************************************************************
   function d_h2o(t,p)
     
     use messy_clamssedi_global, only: T_0, p_0
     use messy_clams_global, only: prec
     
     implicit none
     
     real(prec)             :: d_h2o
     real(prec), intent(in) :: t ! temperature [K]
     real(prec), intent(in) :: p ! pressure [hPa]
     
     d_h2o = 0.211 * ((t/T_0)**1.94) * (p_0/p) * 1.e-4 ! [m**2/s]
     
     return
     
   end function d_h2o
   
   !****************************************************************************
   ! Berechnung der NAT-Kristalldichte (nach K. Carslaw)
   ! Umrechnung der Einheiten von [g/cm**3] nach [kg/m**3]:
   ! 1.e3=1.e-3*1.e6
   !****************************************************************************
   function rho_nat(t)
     
     use messy_clams_global, only: prec
     
     implicit none
     
     real(prec)             :: rho_nat
     real(prec), intent(in) :: t
     
     rho_nat = (1.6252 - 2.3585e-5 * t) * 1.e3   ! [kg/m**3]
     
     return
     
   end function rho_nat
   
   !****************************************************************************
   ! Berechnung der mittleren Molekulargeschwindigkeit von HNO3 (Kuchling, S. 210)
   !****************************************************************************
   function c_hno3(t)
     
     use messy_clams_global,     only: pi, prec
     use messy_clamssedi_global, only: masse_hno3
     use messy_main_constants_mem, only: R_gas
    
     implicit none
     
     real(prec)             :: c_hno3
     real(prec), intent(in) :: t
     
     c_hno3 = sqrt((8 * R_gas * t) / (pi * masse_hno3)) ![m/s]
     
     return
     
   end function c_hno3
   
   
   !****************************************************************************
   ! Berechnung der mittleren Molekulargeschwindigkeit von H2O
   !****************************************************************************
   function c_h2o(t)
     
     use messy_clams_global, only: pi, masse_h2o, prec
     use messy_main_constants_mem, only: R_gas
     
     implicit none
     
     real(prec)             :: c_h2o
     real(prec), intent(in) :: t
     
     c_h2o = sqrt((8 * R_gas * t) / (pi * masse_h2o)) ![m/s]
     
     return
     
   end function c_h2o

   !****************************************************************************
   !
   !****************************************************************************
   subroutine compute_settling_velocity(ipart)

     use messy_clamssedi_global, only: particles, timestep_sedi
     use messy_clams_global,     only: prec, pi, viscosity, masse_luft, &
                                       rho_ice
     use messy_main_constants_mem, only: R_gas
     
     implicit none

     integer, intent(in) :: ipart
     real(prec)          :: mean_free_path ! [m]
     real(prec)          :: slip_correction, terminal_velocity

     ! Berechnung der mittleren freien Weglaenge (Seinfeld, S. 455)
     ! Angabe des Drucks in Pa (Faktor 100), Temperatur in K, alles andere
     ! siehe global.f90. Das Ergebnis hat die Dimension m.
     mean_free_path = (2. * viscosity) / (particles%pressure(ipart) * 100. * &
          sqrt(8. * masse_luft / (pi * R_gas * particles%temperature(ipart))))

     ! Berechnung der Fallgeschwindigkeit des Teilchens (Seinfeld, S. 466)
     ! Dichte : 1.6 g/cm^3 ???
     ! Genauer 1.6252 -2.3585e-5*T [g/cm^3]
     ! Dynamische Viskositaet : 1.72 g/(cm s)
     slip_correction = 1. + (mean_free_path / particles%radius(ipart)) * (1.257 + .4 &
          & * exp(-1.1 * particles%radius(ipart) / mean_free_path))
     
     ! Check for ice or NAT
     if ( (particles%class(ipart) >= 10.) .AND. (particles%class(ipart) < 20.) ) then
        terminal_velocity = -(4. * particles%radius(ipart)**2 &
             & * rho_nat(particles%temperature(ipart)) &
             & * 9.81 * slip_correction) / (18. * viscosity)
     elseif ( particles%class(ipart) >= 20. ) then
        terminal_velocity = -(4. * particles%radius(ipart)**2 &
             & * rho_ice &
             & * 9.81 * slip_correction) / (18. * viscosity)
     else
        terminal_velocity = 0
     endif
     
     particles%sedimentation(ipart) = particles%sedimentation(ipart) + &
          terminal_velocity * timestep_sedi
     particles%tsv(ipart) = terminal_velocity
     
   end subroutine compute_settling_velocity

End Module messy_clamssedi_compute
