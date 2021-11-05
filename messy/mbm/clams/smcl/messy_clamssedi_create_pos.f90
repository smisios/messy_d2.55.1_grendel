!******************************************************************************
!
! Modul clamssedi_create_pos
!
! Written by Gebhard Guenther, Thomas Breuer                                                  
! IEK-7, Forschungszentrum Juelich
! Last Modified By: N.Thomas
! Last Modified On: Mon Sep 21 10:09:13 2015
!
!  subroutine create_position
!  subroutine latlon_step(lon_in, lat_in, area, lon_out, lat_out)
!  subroutine lev_step(lev_min, lev_grid, levs_delta, nfine, lswitch, i, &
!                      count, countold, idummy, tdummy, is, levold, lev)
!
!******************************************************************************

Module messy_clamssedi_create_pos

CONTAINS


  !*********************************************************************                  
  !
  !*********************************************************************                  
  subroutine create_position (status)

    USE messy_clams_global,        ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

    USE messy_clams_global,        ONLY: DP, prec, rank, ntasks, mdi, &
                                         pi, radius_earth, &
                                         pre_year, pre_month, pre_day, pre_sec
    USE messy_clamssedi_global,    ONLY: nparticle_max, pos_lev_grid, pos_lev_delta, &
                                         pos_r_grid, pos_nlevs, lev_min, r_high, &
                                         nfine, factor,  &
                                         lat_min, lat_max, &
                                         nparticles, particles, nparticles_added, &
                                         hno3_background, h2o_background
    USE messy_clams_tools_dateconv,    ONLY: ymds2js
    USE messy_clams_tools_interpolreg, ONLY: interpolate_param_one_point
    
    implicit none

    integer :: status

    real(prec), parameter :: lon_start = 0.

    integer    :: nlevs_help, counter_lev       ! , totalcount
    integer    :: idummy, tdummy, countold, lswitch, n_thresh, is, ilev
    real(prec) :: ldummy
    real(prec) :: l_diff, l_diffmin, step_lev
    real(prec) :: levold, lev
    real(prec) :: r_high_act, area_high, p_hm, p_mk
    real(prec) :: temperature, pressure, partialdruck_hno3, partialdruck_h2o
    real(prec) :: lon, lat
    real(prec) :: lon_out, lat_out
    real(DP)   :: itime, ipasttime
    
    logical :: lat_range

    status = 0
    
    ! Initializations for the subroutine lev_step
    lswitch = 0
    ldummy = 999.99
    countold = 1
    levold = 1
    idummy = 0
    tdummy = 0

    nlevs_help = ubound(pos_lev_grid,1)

    ! totalcount = 0
    counter_lev = 1           ! added particles for one level 
    nparticles_added = 0      ! added particles in total for each rank
    
    ! Start of main loop
    do ilev=rank+1, nlevs_help, ntasks

       ! if (rank /= 0) totalcount = 0
       
       if (associated(pos_r_grid)) then
          r_high_act = pos_r_grid(ilev) / sqrt(factor)
       else
          r_high_act = r_high / sqrt(factor)
       endif
       area_high = (r_high_act / (radius_earth/1000.))**2
       
       l_diffmin = sqrt(area_high)*(180./pi)

       idummy = 0
       tdummy = ilev - 1
       
       ! make fine grid
       lon = lon_start
       lat = lat_max
       lat_range = .true.
       do while (lat_range)                     ! Set remaining fine-grid points
          call latlon_step(lon, lat, area_high, lon_out, lat_out)
          if (lat_out >= lat_min) then
             lon = lon_out
             lat = lat_out
             lev = pos_lev_grid(ilev)
             
             if (nfine > 0 ) then
                if (ldummy /= lat) lswitch=1
                
                call lev_step(pos_lev_grid, pos_lev_delta, nfine, &
                     lswitch, ilev, counter_lev, countold, idummy, tdummy, is, levold, lev)
                
                ldummy = lat
                lswitch = 0
                ! The following condition relates to the threshold problem                            
                ! at lat=lat_min, the meetingpoint of the coarse and fine grid (START
                ! and END of a run on a given theta_level.
                ! This condition prevents that the gridpoints of last latitude-zone
                ! of the fine grid has the same lev-shift as the gridpoints of the
                ! first latitude-zone of the coarse grid:
                ! Latitude:
                ! North pole    ....       lat_min         ......         South pole
                !--------------------------------------------------------------------
                !                            |START at lat_min
                !                            |with the coarse grid     -----> | Jump
                !                            |                                |  to                   
                !                            |                                |  the                  
                !                            |                                | North                 
                !  <---------------------------------------------------------<| Pole
                ! Continue with the fine grid|                                                        
                !   ------>   END at lat_min |
                
                n_thresh = mod(idummy, nfine)
                l_diff = abs(lat_min - lat)
                if (l_diff .lt. l_diffmin .and. n_thresh .eq. 1) then
                   step_lev  = pos_lev_delta(ilev)/real(nfine)
                   if (is .lt. nfine-1) lev = lev + step_lev
                   if (is .eq. nfine-1) lev = lev - (nfine-1) * step_lev
                endif
             endif

             ! interpolate temperature and pressure
             itime = ymds2js (YEAR,MONTH,DAY,HOUR*3600+MINUTE*60+SECOND)
             ipasttime = ymds2js (pre_year,pre_month,pre_day,pre_sec)
             
             call interpolate_param_one_point (status, lon, lat, lev, itime, ipasttime, &
                                         "TEMP", temperature)
             if (status/=0) return
             call interpolate_param_one_point (status, lon, lat, lev, itime, ipasttime, &
                                         "PRESS", pressure)
             if (status/=0) return

             partialdruck_hno3 = hno3_background * pressure
             partialdruck_h2o = h2o_background * pressure
             p_hm = exp((-2.7836 - 0.00088 * temperature)*log(h2o_background*pressure) + &
                  &90.8556 - 26242.6/temperature + 0.0213885 * temperature)
             p_mk = exp(9.550426 - 5723.265/temperature + 3.53068*log(temperature) - &
                  &0.00728332*temperature)
             
!!!!!in if hinzufuegen:::  .OR. partialdruck_h2o  >= p_mk             
             if (partialdruck_hno3 >= p_hm) then
                nparticles_added = nparticles_added + 1
                
                particles%lat(nparticles+nparticles_added) = lat
                particles%lon(nparticles+nparticles_added) = lon
                particles%lev(nparticles+nparticles_added) = lev
                particles%pressure(nparticles+nparticles_added) = pressure
                particles%temperature(nparticles+nparticles_added) = temperature
                
                counter_lev = counter_lev + 1
             endif
             
          else
             lat_range=.false.
          endif
       
       enddo


    enddo

  end subroutine create_position


!******************************************************************************
!
!******************************************************************************
  subroutine latlon_step(lon_in, lat_in, area, lon_out, lat_out)

    USE messy_clams_global,        ONLY: prec

    ! Input: 0 <= lon_in < 360, -90 <= lat_in <= 90
    !        area - area on a unit sphere
    ! Output: lon_out, lat_out
    ! Purpose: Starting from (lon_in, lat_in), the next point
    ! of a lon-lat grid with equal areas for the grid mashes will
    ! be determined
    !
    ! (lon_in, lat_in)
    !       *------------------------>* (lon_out, lat_out) along lat=const
    !       |                         . 
    !       |         area            .
    !       V                         .
    !       *.........................*
    ! (lon_out, lat_out) after a circle with lat=const was filled

    USE messy_clams_global,        ONLY: pi

    implicit none

    real(prec), intent(in)                            :: lon_in, lat_in, area
    real(prec), intent(out)                           :: lon_out, lat_out

    real(prec)                                        :: lon, lat, &
                                                         delta_lat, delta_lon, &
                                                         new_lon

    !transfer to rad
    lon=lon_in*(pi/180.)
    lat=lat_in*(pi/180.)

    !determine delta_lon
    if (abs(abs(lat) - pi/2.) >= 1.e-4) then ! avoid division by 0
       delta_lon=sqrt(area)/cos(lat)
    else
       delta_lon=3.*pi                       ! near the poles, make
    endif                                    ! a maximal lon-step

    ! determine delta_lat
    delta_lat=sqrt(area)

    ! carry out the lon-step 
    new_lon=lon+delta_lon

    ! check if the lon-edge (2*pi) was reached
    if (new_lon <= (2.*pi)) then
       ! no -> do not change lat
       lon_out=lon_in+delta_lon*(180./pi)
       lat_out=lat_in
    else
       ! yes -> carry out the lat-step, set lon=0 
       lon_out=0.
       lat_out=lat_in-delta_lat*(180./pi)
    endif
  end subroutine latlon_step

!******************************************************************************
!
!******************************************************************************
  subroutine lev_step (lev_grid, levs_delta, nfine, &
                         lswitch, i, count, &
                         countold, idummy, tdummy, is, &
                         levold, lev)

    USE messy_clams_global,        ONLY: prec

    implicit none

    !Input:  lev-levels, as indicated in bounds_vert.inp 
    !Output: shifted lev-levels, shared on nfine sublevels around the original
    !         lev value:
    !
    !            INPUT              -->                          OUTPUT
    ! lev_max ----------------------------------------------------------------
    !                            lev_min(i)+(2*nfine-1)*delta/2 --------------
    !                                                                     .
    !                                                                     .
    !                                                                     .
    ! lev     --------------                                            .
    !                                                                     .
    !                                   lev_min(i)+3*delta/2    --------------
    !                                   lev_min(i)+  delta/2    --------------
    ! lev_min ----------------------------------------------------------------
    !                         
    ! As explained in the subroutine lonlat-step on a given lev-level a lon-lat
    ! grid is rastered, following latitude-zones.
    ! For the i-th lat-zone (lat=const) is:
    !
    ! lon=0°                                                 lon=<360°
    !   *------------------------------------------------------>*
    ! count=count[lat(i-1)]+1                                count=count[lat(i)]
    ! countold=count
    !
    ! Count counts the gridpoints, countold saves the gridpoint, after which
    ! the latitude changes and their difference (cdummy) gives the number of
    ! gridpoints in one latitude-zone. idummy counts the lat-zones on a given
    ! lev-level and lswitch marks the transition of the latitude-zones.
    !  The vertical rastering is performed similar to the horizontal rastering
    ! using lev, levold and tdummy. So results a shifting factor 'is' which
    ! permutes from north to south and with the altitude (lev).

    real(prec), dimension(:), pointer :: lev_grid, levs_delta

    integer, intent(in)   :: nfine, lswitch, i, count
    integer, intent(inout):: tdummy
    integer, intent(out)  :: countold, idummy, is
    real(prec), intent(out)     :: levold, lev

    integer                     :: cdummy
    real(prec)                  :: delta, lev_shift, lev_min_act

    if (levold.ne.lev_grid(i)) then
       idummy = 0
       tdummy = tdummy + 1
    endif

    lev_min_act = lev_grid(1) + sum(levs_delta(1:i-1))

    if (lswitch .ne. 0) then
       idummy = idummy+1
       countold = count
       levold = lev_grid(i)
    endif

    cdummy = count - countold
    is = mod(cdummy+(idummy-1) + (tdummy-1), nfine)
    delta = levs_delta(i) / real(nfine)
    lev_shift = delta*is + delta/2.
    lev = lev_min_act + lev_shift

  end subroutine lev_step


End Module messy_clamssedi_create_pos
