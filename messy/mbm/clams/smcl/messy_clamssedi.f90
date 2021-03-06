MODULE messy_clamssedi

  ! **************************************************************************
  ! 
  ! Written by Thomas Breuer, Nicole Thomas
  ! IEK-7, Forschungszentrum Juelich
  ! 
  ! subroutine clamssedi(status)
  ! subroutine sedi_prepare(status)
  ! SUBROUTINE clamssedi_read_nml(status, iou)
  ! 
  ! **************************************************************************

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamssedi'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  PUBLIC :: clamssedi
  PUBLIC :: clamssedi_read_nml

!----------- 
CONTAINS
!-----------

!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!

subroutine clamssedi(status, lcoupled)
  
  !----------------------------------------------
  !
  !----------------------------------------------
  
  USE messy_clams_global,            ONLY: DP, rank, ntasks, pi, &
                                           pre_year, pre_month, pre_day, pre_sec
  USE messy_clamssedi_global,        ONLY: nparticles, particles, timestep_sedi, &
                                           part_ind_start, part_ind_end 
  use messy_clamssedi_compute,       only: interpolate_background_chem, &
                                           compute_particle_growth, &
                                           compute_settling_velocity
  USE messy_clamssedi_traj,          ONLY: compute_trajectory
  USE messy_clams_tools_dateconv,    ONLY: ymds2js
  USE messy_clams_tools_interpolreg, ONLY: interpolate_param_one_point
  USE messy_clams_global,            ONLY: YEAR, MONTH, DAY, HOUR, &
                                           MINUTE, SECOND
  
  implicit none

  integer  :: status
  logical  :: lcoupled
  
  integer  :: ipart, ipart_local
  real(DP) :: itime, ipasttime

  ! if (rank==0) THEN
  !   WRITE (*,*) 
  !   WRITE (*,*) 'START OF SEDI'
  !   WRITE (*,*) 
  ! endif
  
  do ipart = part_ind_start, part_ind_end

     ipart_local = ipart - part_ind_start + 1

     if (particles%radius(ipart) > 0. .and. particles%density(ipart) > 0.) then
        call interpolate_background_chem(ipart, ipart_local)
        call compute_particle_growth(ipart, ipart_local, status)
        if (status /= 0) return
        if (particles%radius(ipart) > 0. .and. particles%density(ipart) > 0.) then
           call compute_settling_velocity(ipart)
        endif
     endif
  enddo

  write (*,*) 'call compute_trajectory', rank
  call compute_trajectory (lcoupled)
  
  ! interpolate temp/press for 'local' particles
  itime = ymds2js (YEAR,MONTH,DAY,HOUR*3600+MINUTE*60+SECOND) + timestep_sedi
  ipasttime = ymds2js (pre_year,pre_month,pre_day,pre_sec)

  do ipart = part_ind_start, part_ind_end
     if (particles%radius(ipart) > 0. .and. particles%density(ipart) > 0.) then
        particles%lat(ipart) = particles%lat(ipart) * 180./pi
! op_pj_20160606+
!!$        particles%lon(ipart) = mod(particles%lon(ipart) * 180./pi, 360.)
        particles%lon(ipart) = mod(particles%lon(ipart) * 180._dp/pi, 360._dp)
! op_pj_20160606-

        call interpolate_param_one_point (status, particles%lon(ipart), &
                                          particles%lat(ipart), particles%lev(ipart), &
                                          itime, ipasttime, "TEMP", &
                                          particles%temperature(ipart))
        if (status/=0) return
        call interpolate_param_one_point (status, particles%lon(ipart), &
                                          particles%lat(ipart), particles%lev(ipart), &
                                          itime, ipasttime, "PRESS", &
                                          particles%pressure(ipart))
        if (status/=0) return
     endif
  enddo
  
  return
  
end subroutine clamssedi
!---------------------------------------------------------------------------
subroutine sedi_prepare(status)

  use messy_clamssedi_global,   only: nparticles, particles, densnat_val, &
                                      triangles, part_ind_start, part_ind_end
  use messy_clamssedi_triang,   only: get_triangles_for_particle
  use messy_clamssedi_compute,  only: interpolate_background_chem, & 
                                      pmkh2oice, phmhno3, calc_nucrate
  use messy_clams_global,       only: rank, ntasks, prec, &
                                      k_boltzmann, avogadro, masse_luft
  use messy_main_constants_mem, only: R_gas
  
  implicit none

  integer    :: status, ipart
  real(prec) :: density_m, ph2o, phno3
  integer    :: ipart_local

!!! Nicole: phmh2o, phmhno3 bereits vorhanden. Evtl eigenes Modul  
  
  do ipart = part_ind_start, part_ind_end
     
     ipart_local = ipart - part_ind_start + 1
     
     call get_triangles_for_particle(ipart, ipart_local, status)
     if (status /= 0) return 
     
     call interpolate_background_chem(ipart, ipart_local)
     
     ! density of air [kg/m^3]
     density_m = ( particles%pressure(ipart) * 100./ &
          (k_boltzmann * particles%temperature(ipart)) ) &
          & * masse_luft / avogadro
     
     ph2o = (particles%h2o(ipart) * particles%temperature(ipart) * density_m &
          & * R_gas / masse_luft) * 9.87e-6
     
     phno3 = (particles%hno3(ipart) * particles%temperature(ipart) * density_m &
          & * R_gas / masse_luft) * 9.87e-6
     
     ! determine the current supersaturation
     particles%sice(ipart) = ph2o / pmkh2oice(particles%temperature(ipart))
     particles%snat(ipart) = phno3 / phmhno3(particles%temperature(ipart), ph2o)

     particles%sicemax(ipart) = max(particles%sicemax(ipart), particles%sice(ipart))
     particles%snatmax(ipart) = max(particles%snatmax(ipart), particles%snat(ipart))
     particles%tmin(ipart) = min(particles%tmin(ipart), particles%temperature(ipart))
    
     ! calculate nucleation rate (only once for each particle)
     if (densnat_val < 0 .and. particles%density(ipart) <= 0.) then  
        call calc_nucrate(ipart, ipart_local)
     else if (densnat_val >= 0 .and. particles%particle_id(ipart) > 0) then
        particles%class(ipart) = 11.
     endif

     ! increase the counter of ice or nat particles in the triangle above 
     ! and below the current particle based on the class of the particle
     if (particles%class(ipart) >= 20.) then
        triangles%nparticles_ice(particles%tr_ind_down(ipart_local), particles%lev_down(ipart_local)) = &
             triangles%nparticles_ice(particles%tr_ind_down(ipart_local), particles%lev_down(ipart_local)) + 1
     
        triangles%nparticles_ice(particles%tr_ind_up(ipart_local), particles%lev_up(ipart_local)) = &
             triangles%nparticles_ice(particles%tr_ind_up(ipart_local), particles%lev_up(ipart_local)) + 1
     elseif ( (particles%class(ipart) >= 10.) .and. (particles%class(ipart) < 20.) ) then
        triangles%nparticles_nat(particles%tr_ind_down(ipart_local), particles%lev_down(ipart_local)) = &
             triangles%nparticles_nat(particles%tr_ind_down(ipart_local), particles%lev_down(ipart_local)) + 1

        triangles%nparticles_nat(particles%tr_ind_up(ipart_local), particles%lev_up(ipart_local)) = &
             triangles%nparticles_nat(particles%tr_ind_up(ipart_local), particles%lev_up(ipart_local)) + 1
     endif

  enddo
  
end subroutine sedi_prepare

!---------------------------------------------------------------------------
SUBROUTINE clamssedi_read_nml(status, iou)
  
  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  USE messy_clams_global, ONLY: met_freq
  USE messy_clamssedi_global
  
  IMPLICIT NONE
  !I/O
  INTEGER, INTENT(OUT) ::   status
  INTEGER, INTENT(IN)  ::   iou
  !LOCAL
  CHARACTER(LEN=*),PARAMETER :: substr='clamssedi_read_nml'
  LOGICAL              :: lex      ! file exists ?
  INTEGER              :: fstat    ! file status 
  
  NAMELIST /CTRL/ lat_min, lat_max, lev_start, lev_end, factor, nfine, &
       & densnat_val, flaggary, timestep_sedi,  &
       & nparticle_max, ice_nuc_table, nat_nuc_table, part_init_file, &
       & nat_tstep, ice_tstep

  status = 1 !ERROR
  CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, .false.)
  IF (.not.lex) RETURN    ! <modstr>.nml does not exist 
  
  READ(iou, NML=CTRL, IOSTAT=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, .false.)
  IF (fstat /= 0) RETURN  ! error while reading namelist
  
  CALL read_nml_close(substr, iou, modstr, .false.)
  
  status = 0 !NO ERROR
  

END SUBROUTINE clamssedi_read_nml
  
! **************************************************************************
END MODULE messy_clamssedi
! ************************************************************************** 
