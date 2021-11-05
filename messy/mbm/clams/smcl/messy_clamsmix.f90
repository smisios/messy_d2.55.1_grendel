MODULE messy_clamsmix
! **************************************************************************
! MODULE FOR CLaMS MIXING MODULE 'MIX'
!
! This module contains the following subroutines:
!   SUBROUTINE mix (status, my_nparts, levelrange, &
!                   lat, lat_old, lat_new, lon, lon_old, lon_new, &
!                   zeta, zeta_new, specarr, specarr_new)
!  SUBROUTINE clamsmix_read_nml(status, iou)
!  SUBROUTINE clamsmix_set_config (status)
!  SUBROUTINE clamsmix_nparts_per_level (irange, dnparts_level, nparts_level)
!  SUBROUTINE clamsmix_set_level (levelrange, nparts_level)
!  SUBROUTINE set_adapt_par (a_p, ilev, grid_switch)
!  SUBROUTINE clamsmix_set_bounds 
!
! **************************************************************************

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamsmix'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'


!----------- 
CONTAINS
!----------- 

  !****************************************************************************
  !  PROGRAM MIX -- IMPLICIT 2d/3d-MIXING BETWEEN AIR PARCELS 
  !****************************************************************************
  SUBROUTINE mix (status, my_nparts, &
                  lat, lat_old, lat_new, lon, lon_old, lon_new, &
                  zeta, zeta_new, state, state_vert, specarr, specarr_new, &
                  theta, theta_old, bvf, bvf_old)

    use messy_clamsmix_global
    use messy_clams_global,         ONLY: rank, mdi, prec, dp,&
                                          species_type
    use messy_clamsmix_lib_io,      ONLY: find_index, get_ap_s, save_ap_s
    use messy_clamsmix_lib_mix,     ONLY: adapt_grid
    use messy_clamsmix_lib_triang,  only: qhull
    use messy_clamsmix_ap_m_access, only: ap_struct, object, free, n, &
                                          reset_state, reset_state_vert

    implicit none
    
    integer    :: status, my_nparts

!!!!! alle Felder sind mit dnparts_max_shuffle dimensioniert !!!
    real(prec),         dimension(:), pointer :: lat, lat_old, lat_new    
    real(prec),         dimension(:), pointer :: lon, lon_old, lon_new
    real(prec),         dimension(:), pointer :: zeta, zeta_new
    real(prec),         dimension(:), pointer :: state, state_vert
    TYPE(species_type), DIMENSION(:), POINTER :: specarr, specarr_new
    real(prec),         dimension(:), pointer :: theta, theta_old
    real(prec),         dimension(:), pointer :: bvf, bvf_old

    type(ap_struct)                          :: ap_s
    type(ap),dimension(:),allocatable, target:: ap_s_0
    integer, dimension(:), pointer           :: layer_ind

    REAL(DP)     :: seconds  
    
    integer      :: finish, start, counts_per_sec
    integer      :: step_adapt, &
                    nparts,  nparts_max, &
                    ilev, j, i, writepos, counter, n_old

    logical   :: found

    status = 0 ! no error

    !ctrl_out = .true.

    if(rank==0)  then
       CALL system_clock (count=start)
       write (*,*) 
       write (*,*) 'START OF MIX'
       write (*,*) 
    endif

    if (dates30) then
       write (*,*)
       write (*,*) 'ACHTUNG: 30-Tage-Monate'
       write (*,*)
    endif

    nparts_max = size(lat)

    lat_new  = mdi
    lon_new  = mdi
    zeta_new = mdi
    do i = 1, size(specarr_new)
       specarr_new(i)%values  = mdi
    enddo

    !write (*,*) 'in mix: rank, my_nparts, nparts_max:',rank, my_nparts, nparts_max

    ! write names of species
    if (rank==0) then
       print *, 'specarr:'
       do i = 1, size(specarr)
          write (*,*) trim(specarr(i)%name)
       end do
    endif
  
    if (rank==0) write (*,*) 'switch_mixing=',switch_mixing    
    !--------------------------------------------------------------------
    ! Begin of mixing loop
    !--------------------------------------------------------------------
    if (switch_mixing==1 .or. switch_mixing==2) then           ! implicit mixing part
 
       writepos = 1
     
       do ilev = levelrange(1,rank), levelrange(2,rank)   

          !write (*,*) 'rank, ilev=', rank, ilev
          
          ! alle Punkte zwischen lev_min_act und lev_max_act ermitteln
          ! ! read layer
          found = find_index(status, lev_min_act(ilev), lev_max_act(ilev), &
                             layer_ind, lat, zeta)
          if (status/=0) then
             write (*,*) 'Error in find_index !!!'
             return
          endif
          
          nparts = 0
          if (found) nparts=size(layer_ind)
          
          if (nparts==0) then
             print *, 'Not enough APs in the layer ',ilev
             cycle
          else
             if (ctrl_out) print *, 'rank, # of APs in the layer: ', rank,nparts
             !print *, 'rank, ilev, # of APs: ', rank,ilev,nparts
             
             allocate(ap_s_0(1:nparts))
             do j=1, nparts
                nullify(ap_s_0(j)%ind)
             enddo
             
            call get_ap_s (status, switch_mixing, ap_s_0, layer_ind, lat, lon,  &
                  zeta, lat_old, lon_old, specarr, theta, bvf, theta_old, bvf_old)
             if (status/=0) then
                write (*,*) 'Error in get_ap_s!!!'
                return
             endif
             
             deallocate(layer_ind)
             
             ! mark points between l_min_act and l_max_act
             ! -> alle Punkte markieren, fuer die ilev=levelind(i)
             do j=1, nparts
                if (l_min_act(ilev)<=ap_s_0(j)%lev .and. &
                     ap_s_0(j)%lev<l_max_act(ilev)) then
                   ap_s_0(j)%subset = .true.
                else
                   ap_s_0(j)%subset = .false.
                endif
             enddo
!           write (*,*) 'rank, ilev, nparts, count: ', &
!                rank, ilev, nparts, count(ap_s_0%subset)
          
             if (count(ap_s_0%subset)==0) then
                print *, 'Not enough APs in the layer'
                cycle
             endif
          
             do j=1, nparts
                nullify(ap_s_0(j)%ind)
             enddo
             
             ap_s = object(ap_s_0) 
             call reset_state (ap_s)
             call reset_state_vert(ap_s)
             
             ! begin of implicit mixing        
             if ((grid_switch >0 .and. ilev > 1  .and. ilev < l_nlevs) .or. &
                  (grid_switch == 0)) then

                if (ctrl_out) print *, 'Implicit mixing'
                call qhull (status,ap_s,n,3,'o')
                if (status/=0) then
                   write (*,*) 'Error in qhull !!!'
                   return
                endif

                if (switch_mixing == 2) then
                   call adapt_grid (status, ap_s,adapt_par,4,l_min_act(ilev),ilev)              
                else ! switch_mixing == 1
                   call adapt_grid (status, ap_s,adapt_par,3,l_min_act(ilev),ilev)              
                endif

                if (status/=0) then
                   write (*,*) 'Error in adapt_grid !!!'
                   return
                endif

                do step_adapt=1, adapt_par%no_steps
                   if (ctrl_out) print *, '# Elimination loop: ', step_adapt
                   n_old = n
                   call qhull (status,ap_s,n,3,'a')
                   if (status/=0) then
                      write (*,*) 'Error in qhull !!!'
                      return
                   endif
                   call adapt_grid (status, ap_s,adapt_par,1,l_min_act(ilev),ilev)
                   if (status/=0) then
                      write (*,*) 'Error in adapt_grid !!!'
                      return
                   endif
                   if (ctrl_out) write (*,*) 'n_old, n, rel. Abweichung:',n_old,n,abs(n_old-n)/REAL(n_old)
                   if (abs(n_old-n)/REAL(n_old) <= adapt_par%r_dev) exit
                enddo

             else 
                if (ctrl_out) print *, 'No mixing !'
             endif
             ! end of implicit mixing
             
             if (writepos-1+n > nparts_max) then
                print *,'Too many points'
                status = 99
                return
             endif
          
             call save_ap_s (status, ap_s, lat_new, lon_new, zeta_new, &
                             state, state_vert, specarr_new, &
                             writepos, counter, &
                             l_min_act(ilev),lev_max, &
                             adapt_par%lat_down, adapt_par%lat_up, &
                             adapt_par%lev_grid, adapt_par%lev_delta, &
                             adapt_par%r_grid,  &
                             check_subset=.true.)   
             if (status/=0) then
                write (*,*) 'Error in save_ap_s !!!'
                return
             endif
             
             writepos = writepos + counter

             call free(ap_s)    
             deallocate(ap_s_0)
             
          
          endif
     
       enddo

    else        ! no mixing part

       if (rank==0)  print *, 'No mixing !!!'

       found = find_index(status, lev_min, lev_max, layer_ind, lat, zeta)
       if (status/=0) then
          write (*,*) 'Error in find_index !!!'
          return
       endif
       
       nparts=size(layer_ind)
       if (ctrl_out) print *, 'rank, # of APs: ', rank, nparts
       
       allocate(ap_s_0(1:nparts))
       do j=1, nparts
          nullify(ap_s_0(j)%ind)
       enddo
          
       call get_ap_s (status, switch_mixing,ap_s_0, layer_ind, lat, lon,  &
            zeta, lat_old, lon_old, specarr, theta, bvf, theta_old, bvf_old)
       if (status/=0) then
          write (*,*) 'Error in get_ap_s !!!'
          return
       endif
       
       deallocate(layer_ind)
       
       ap_s = object(ap_s_0) 
       call reset_state (ap_s)
       call reset_state_vert(ap_s)
       
       writepos = 1
       
       call save_ap_s (status, ap_s, lat_new, lon_new, zeta_new, &
                       state, state_vert, specarr_new,  &
                       writepos, counter, &
                       lev_min,lev_max, &
                       adapt_par%lat_down, adapt_par%lat_up, &
                       adapt_par%lev_grid, adapt_par%lev_delta, &
                       adapt_par%r_grid,  &
                       check_subset=.false.)   
       if (status/=0) then
          write (*,*) 'Error in save_ap_s !!!'
          return
       endif
          
       writepos = writepos + counter
          
       call free(ap_s)    
       deallocate(ap_s_0)
             
    endif

    my_nparts = writepos - 1
    !write (*,*) 'nach mixing: rank, my_nparts:',rank, my_nparts


    !--------------------------------------------------------------------
    ! End of mixing loop
    !--------------------------------------------------------------------
    
    if(rank==0)then
       print *, 'Normal termination of mix'
        
       CALL system_clock (COUNT_RATE=counts_per_sec)
       CALL system_clock (count=finish)
       IF (finish > start) THEN
          seconds = float(finish-start) / float(counts_per_sec)
       ELSE
          seconds = 0.0
       ENDIF
       
       !WRITE (*,*)
       !WRITE (*,*) 'System clock runs at ', counts_per_sec, 'ticks per second'
       WRITE (*,*)
       WRITE (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.' 
    endif
    
  END SUBROUTINE mix

  !****************************************************************************
  !
  !****************************************************************************
  SUBROUTINE clamsmix_read_nml(status, iou)

    USE messy_main_tools,       ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_clams_global,     ONLY: rank, ldiagout
    USE messy_clamstraj_global, ONLY: idtsec2
    USE messy_clamsmix_global
    
    IMPLICIT NONE
    
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamsmix_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output
    
    NAMELIST /CTRL/ switch_mixing, lexp,&
                    fac_limit_outside, fac_limit_inside, &
                    fac_limit_lev_down, fac_limit_lev_up, &
                    fac_eliminate, &
                    fac_bvf_min, vert_mix_param, &
                    delta_lev, no_steps, r_dev, &
                    grid_switch, nintervals, timestep_mix, &
                    nmixspec, mixspec

    status = 1 !ERROR
    
    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif

    ! Read namelist variables:
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr, l_print)
        
    status = 0 !NO ERROR
    
  END SUBROUTINE clamsmix_read_nml
  
   
  !****************************************************************************
  ! Set nlevs, l_nlevs, lev_min, lev_max
  !****************************************************************************
  SUBROUTINE clamsmix_set_config (status)

    USE messy_clams_global,    ONLY: rank
    USE messy_clamsmix_global
    
    IMPLICIT NONE

    integer :: status

    status = 0 ! no error
    
    ! get nlevs, l_nlevs, lev_min, lev_max
    if (adapt_par%nlevs == 0 .OR. adapt_par%nlevs == 1) then     ! 2d
       !  isentropic configuration
       nlevs=1
       l_nlevs=1
       lev_min=adapt_par%lev_min-5.
       lev_max=adapt_par%lev_max+5.
    else
       nlevs = adapt_par%nlevs 
       lev_min=adapt_par%lev_min
       lev_max=adapt_par%lev_max
       select case (grid_switch)
       case (0)
          if(rank==0) print*,'Original levels'
          l_nlevs = nlevs * nintervals   
       case (1)
          if(rank==0) print*, 'Shifted levels'
          l_nlevs = nlevs * nintervals + 1 
       case (2)
          if(rank==0) print*, 'Shifted levels (old)'
          if (nintervals==1) then
             l_nlevs = nlevs + 1
          else
             if (rank==0) print*, 'Number of intervals must be 1 !!!' 
             status = -1
             return
          endif
       case default
          if (rank==0) print*, 'Wrong grid_shift' 
          status = -1
          return
       end select
    endif
    
    !write (*,*) 'in clamsmix_set_config: nlevs, l_nlevs=',nlevs, l_nlevs

  END SUBROUTINE clamsmix_set_config

  !****************************************************************************
  ! Get number of particles per level
  ! nparts_level (1:l_nlevs)  : total number of particles on each level (for all ranks)
  ! dnparts_level(1:l_nlevs,i): number of particles on each level on rank i
  !****************************************************************************
  SUBROUTINE clamsmix_nparts_per_level (irange, dnparts_level, nparts_level)

    USE messy_clams_global,    ONLY: ntasks
    USE messy_clamsmix_global, ONLY: l_nlevs

    IMPLICIT NONE

    integer,    dimension(:,:,:),pointer :: irange   ! (2,nlev,0:ntasks-1)
    integer,    dimension(:,:),  pointer :: dnparts_level ! (l_nlevs,0:ntasks-1)
    integer,    dimension(:),    pointer :: nparts_level  ! (l_nlevs)

    integer :: ilev, irank

    allocate (dnparts_level(l_nlevs,0:ntasks-1))
    allocate (nparts_level(l_nlevs))

    nparts_level  = 0
    dnparts_level = 0
    do ilev = 1, l_nlevs
       do irank = 0, ntasks-1
          if (irange(1,ilev,irank) > 0) then
             dnparts_level(ilev,irank) = irange(2,ilev,irank)-irange(1,ilev,irank)+1
          endif
       enddo
       nparts_level(ilev) = sum (dnparts_level(ilev,:))
    enddo

  END SUBROUTINE clamsmix_nparts_per_level

  !****************************************************************************
  ! Set start- and endlevel for each task
  !****************************************************************************
  SUBROUTINE clamsmix_set_level (levelrange, nparts_level)

    USE messy_clams_global,    only: rank, ntasks, nparts
    USE messy_clamsmix_global, only: l_nlevs

    IMPLICIT NONE

    integer, dimension(:,:), pointer :: levelrange   ! (2,0:ntasks-1)
    integer, dimension(:),   pointer :: nparts_level ! (l_nlevs)

    integer :: irank, i
    integer :: startlevel, endlevel, my_nparts, all_nparts

    allocate (levelrange(2,0:ntasks-1))
    levelrange = -1
    

    ! Jede CPU bearbeitet gleich viele Level:
    ! if (ntasks>=l_nlevs) then
    !    do irank = 0, l_nlevs-1
    !       levelrange(1,irank) = irank + 1
    !       levelrange(2,irank) = irank + 1
    !    enddo
    ! else
    !    do irank = 0, ntasks-1
    !       levelrange(1,irank) =  irank   *(l_nlevs/ntasks) + 1
    !       levelrange(2,irank) = (irank+1)*(l_nlevs/ntasks) +  &
    !                          ((irank+1)/ntasks)*modulo(l_nlevs,ntasks)
    !    enddo
    ! endif


    ! Wenn mehr CPUs als Level genutzt werden, bearbeiten die ersten nlevs CPUs
    ! jeweils ein Level, die weiteren bleiben ungenutzt
    if (ntasks>=l_nlevs) then
       do irank = 0, l_nlevs-1
          levelrange(1,irank) = irank + 1
          levelrange(2,irank) = irank + 1
       enddo
     
    ! Anz. Punkte moeglichst gleichmaessig auf CPUs verteilen
    else

       startlevel = 1
       all_nparts = 0

       ! Ermittle die Level, die CPU irank bearbeiten soll:
       do irank = 0, ntasks-1

          endlevel = startlevel
          my_nparts = nparts_level(startlevel)
          
          ! Es wird solange jeweils das naechste Level hinzugefuegt, bis die Anzahl der Punkte 
          ! moeglichst nah an nparts/ntasks liegt:
          do i = startlevel+1, l_nlevs
             if (abs(nparts/ntasks-my_nparts)<abs(nparts/ntasks-(my_nparts+nparts_level(i)))) then
                exit
             else
                endlevel = i
                my_nparts = my_nparts+nparts_level(i)
             endif
          enddo

          ! Es wird zusaetzlich ueberprueft, ob die Gesamtzahl der bisher zugeteilten Punkte
          ! moeglichst nah an (nparts/ntasks)*(irank+1) liegt.
          ! Sonst kann es passieren, dass bei allen CPUs die Anzahl unter (nparts/ntasks) bleibt,
          ! und fuer die letzte CPU zuviele Punkte uebrig bleiben.
          all_nparts = all_nparts + my_nparts
          if (endlevel<l_nlevs) then
             if (abs((irank+1)*(nparts/ntasks)-all_nparts)>  &
                  abs((irank+1)*(nparts/ntasks)-(all_nparts+nparts_level(endlevel+1)))) then
                endlevel = endlevel+1
                my_nparts = my_nparts+nparts_level(endlevel)
                all_nparts = all_nparts+nparts_level(endlevel)
             endif
          endif
      
          ! Falls noch nicht alle Level zugeteilt sind, macht die letzte CPU den Rest
          if (irank == ntasks-1  .AND. endlevel<l_nlevs) endlevel=l_nlevs

          ! erstes und letztes zu bearbeitendes Level merken
          levelrange(1,irank) = startlevel
          levelrange(2,irank) = endlevel

          if (endlevel == l_nlevs) exit

          startlevel = endlevel + 1

       enddo
    endif
    
  END SUBROUTINE clamsmix_set_level

  !*************************************************************************
  !
  !*************************************************************************
  SUBROUTINE set_adapt_par (a_p, ilev, grid_switch)
 
    USE messy_clams_global,    ONLY: prec
    USE messy_clamsmix_global, ONLY: adapt_set, ctrl_out, nlevs

    IMPLICIT NONE 

    type(adapt_set)                      :: a_p
    integer, intent(in)                  :: ilev, grid_switch
    
    real(prec), parameter                :: r_earth=6378169.0, &  ! in m
                                            limit_max=0.4
    real(prec)                           :: r_mean_c, r_mean_h, phi_rot, hh

    ! define param for calc of Delta_lev
    hh = (2./7.)*(1./7000.)  
    ! Pot. Temp:      lev/lev_0=(p_0/p)^(R/c_p), R/c_p=2/7
    ! Standard atmos: dz/H = ln(p_0/p)             , H=7 km
    ! With lev=lev_0+dlev => (R/c_p)
    !      ln(1+dlev/lev_0)=(R/c_p)*(1/H)
    
    ! define mean distances on a unit sphere
    if (associated(a_p%r_grid)) then
       if (grid_switch==0) then
          r_mean_h = 1000. * a_p%r_grid(ilev) / r_earth
          r_mean_c = 1000. * a_p%r_grid(ilev)*(a_p%r_mean_c/a_p%r_mean_h) / r_earth
       else
          if (ilev==1) then
             r_mean_h = 1000. * a_p%r_grid(1) / r_earth
             r_mean_c = 1000. * a_p%r_grid(1)*(a_p%r_mean_c/a_p%r_mean_h) / r_earth
          elseif (ilev>nlevs) then
             r_mean_h = 1000. * a_p%r_grid(nlevs) / r_earth
             r_mean_c = 1000. * a_p%r_grid(nlevs)*(a_p%r_mean_c/a_p%r_mean_h) / r_earth
          else
             r_mean_h = 1000. * 0.5*(a_p%r_grid(ilev)+a_p%r_grid(ilev-1)) / r_earth
             r_mean_c = 1000. * 0.5*(a_p%r_grid(ilev)+a_p%r_grid(ilev-1)) & 
                  *(a_p%r_mean_c/a_p%r_mean_h) / r_earth
          endif
       endif
    else
       r_mean_h = 1000. * a_p%r_mean_h / r_earth
       r_mean_c = 1000. * a_p%r_mean_c / r_earth
    endif
    
    !if (ctrl_out) write (*,*) 'r_mean_h, r_mean_c:',r_mean_h,r_mean_c
    
!    use Ljapunov exponent for definition of fac_min and fac_max
!     phi_rot = 0.5*atan(2./(a_p%lexp*time_step))
!     a_p%fac_min=sqrt(1.-a_p%lexp*time_step*tan(phi_rot))
!     a_p%fac_max=sqrt(1.+a_p%lexp*time_step*(1./tan(phi_rot)))
    a_p%fac_min(ilev) = exp(-a_p%timestep * a_p%lexp / (24.*3600.))
    a_p%fac_max(ilev) = exp(a_p%timestep * a_p%lexp / (24.*3600.))
    
    ! define r_min/r_max/r_limit
    !a_p%r_min_c(ilev) = a_p%fac_min(ilev) * r_mean_c
    a_p%r_min_c(ilev) = a_p%fac_min(ilev) * r_mean_c * a_p%fac_eliminate
    a_p%r_max_c(ilev) = a_p%fac_max(ilev) * r_mean_c
    a_p%r_lim_c_outside(ilev) = a_p%fac_limit_outside * a_p%r_max_c(ilev)
    a_p%r_lim_c_inside(ilev)  = a_p%fac_limit_inside * a_p%r_max_c(ilev)
    
    !a_p%r_min_h(ilev) = a_p%fac_min(ilev) * r_mean_h
    a_p%r_min_h(ilev) = a_p%fac_min(ilev) * r_mean_h * a_p%fac_eliminate
    a_p%r_max_h(ilev) = a_p%fac_max(ilev) * r_mean_h
    a_p%r_lim_h_outside(ilev) = a_p%fac_limit_outside * a_p%r_max_h(ilev)
    a_p%r_lim_h_inside(ilev)  = a_p%fac_limit_inside * a_p%r_max_h(ilev)
  
    ! redefine r_limit for runs within a reduced lat-range
    if (a_p%lat_down > -89. .or. a_p%lat_up < 89.) then 
       if (a_p%r_lim_c_outside(ilev) > limit_max) then
          a_p%r_lim_c_outside(ilev) = limit_max
       endif
       if (a_p%r_lim_c_inside(ilev) > limit_max) then
          a_p%r_lim_c_inside(ilev) = limit_max
       endif
       if (a_p%r_lim_h_outside(ilev) > limit_max) then 
          a_p%r_lim_h_outside(ilev) = limit_max
       endif
       if (a_p%r_lim_h_inside(ilev) > limit_max) then 
          a_p%r_lim_h_inside(ilev) = limit_max
       endif
    endif
    
  END SUBROUTINE set_adapt_par

  !*************************************************************************
  ! Set boundaries for each level 
  !*************************************************************************
  SUBROUTINE clamsmix_set_bounds 

    USE messy_clams_global,     only: rank
    use messy_clamsmix_global
    use messy_clamsmix_lib_mix, only: lev2levdelta
    
    IMPLICIT NONE

    real(prec)   :: l_delta_act_before
    
    integer      :: lev_int, l_lev_int, lev_int_before 
    integer      :: ilev, i, irank

    allocate (lev_min_act(l_nlevs))
    allocate (lev_max_act(l_nlevs))
    allocate (lev_delta_act(l_nlevs))
    allocate (l_min_act(l_nlevs))
    allocate (l_max_act(l_nlevs))
    allocate (l_delta_act(l_nlevs))

    ALLOCATE(adapt_par%fac_min(l_nlevs))
    ALLOCATE(adapt_par%fac_max(l_nlevs))
    ALLOCATE(adapt_par%r_min_c(l_nlevs))
    ALLOCATE(adapt_par%r_max_c(l_nlevs))
    ALLOCATE(adapt_par%r_min_h(l_nlevs))
    ALLOCATE(adapt_par%r_max_h(l_nlevs))
    ALLOCATE(adapt_par%r_lim_c_outside(l_nlevs))
    ALLOCATE(adapt_par%r_lim_h_outside(l_nlevs))
    ALLOCATE(adapt_par%r_lim_c_inside(l_nlevs))
    ALLOCATE(adapt_par%r_lim_h_inside(l_nlevs))

    do ilev = 1, l_nlevs
        
       if (nlevs == 1) then    ! 2d isentropic configuration
           
          !  define parameters for grid adaptation
          call set_adapt_par(adapt_par, ilev, grid_switch)
          
          print *, '2d isentropic'
          
          l_delta_act(ilev) = 10.
          l_min_act(ilev)   = adapt_par%lev_min-l_delta_act(ilev)/2.
          l_max_act(ilev)   = adapt_par%lev_max+l_delta_act(ilev)/2.
          lev_delta_act(ilev) = l_delta_act(ilev)
          lev_min_act(ilev)   = l_min_act(ilev)
          lev_max_act(ilev)   = l_max_act(ilev)
          
       else
          
          if (grid_switch==0) then  ! original levels
             
             ! Nummer des grossen Intervalls: ((ilev-1)/nintervals)+1
             ! Nummer des kleineren Intervalls (innerhalb des grossen):
             !  MOD((ilev-1),nintervals)+1
             
             lev_int = ((ilev-1)/nintervals)+1
             l_lev_int = MOD((ilev-1),nintervals)+1
             
             !--------------------------------------------------------
             !  define parameters for grid adaptation
             !--------------------------------------------------------
             
             call set_adapt_par(adapt_par, lev_int, grid_switch)
             
             !--------------------------------------------------------
             ! determine l_min_act and l_max_act
             !--------------------------------------------------------
             
             l_delta_act(ilev) = adapt_par%lev_delta(lev_int) / nintervals
             
             l_min_act(ilev) = adapt_par%lev_min 
             do i = 1, lev_int-1
                l_min_act(ilev) = l_min_act(ilev) + adapt_par%lev_delta(i)
             enddo
             l_min_act(ilev) = l_min_act(ilev) + (l_lev_int-1)*l_delta_act(ilev)
             
             l_max_act(ilev) = l_min_act(ilev) + l_delta_act(ilev)
             
             !--------------------------------------------------------
             ! determine lev_min_act and lev_max_act
             !--------------------------------------------------------
             
             if (nintervals==1) then
                lev_min_act(ilev) = l_min_act(ilev)
                lev_max_act(ilev) = l_max_act(ilev)
                lev_delta_act(ilev) = lev_max_act(ilev) - lev_min_act(ilev)
                
             else
                lev_delta_act(ilev) = lev2levdelta(adapt_par%lev_grid, &
                     adapt_par%lev_delta, &
                     (l_min_act(ilev)+l_max_act(ilev))/2.)
                
                lev_min_act(ilev) = (l_min_act(ilev)+l_max_act(ilev))/2. - &
                                         lev_delta_act(ilev)/2.
                lev_max_act(ilev) = (l_min_act(ilev)+l_max_act(ilev))/2. + &
                                         lev_delta_act(ilev)/2.
                if (lev_min_act(ilev) < adapt_par%lev_min) &
                    lev_min_act(ilev)=adapt_par%lev_min
                if (lev_max_act(ilev) > adapt_par%lev_max) &
                    lev_max_act(ilev)=adapt_par%lev_max
                 
             endif

          elseif (grid_switch==1) then    ! shifted levels
              
             lev_int_before = ((ilev-2)/nintervals)+1
             if (lev_int_before==0) lev_int_before=1
             lev_int   = ((ilev-1)/nintervals)+1
             l_lev_int = MOD((ilev-1),nintervals)+1
              
             !--------------------------------------------------------
             !  define parameters for grid adaptation
             !--------------------------------------------------------
             call set_adapt_par(adapt_par, lev_int, grid_switch)

             !--------------------------------------------------------
             ! determine l_min_act and l_max_act
             !--------------------------------------------------------

             l_delta_act_before = adapt_par%lev_delta(lev_int_before) / nintervals
             
             if (lev_int > nlevs) then
                l_delta_act(ilev) = adapt_par%lev_delta(nlevs) / nintervals
             else
                l_delta_act(ilev) = adapt_par%lev_delta(lev_int) / nintervals
             endif
             
             l_min_act(ilev) = adapt_par%lev_min 
             do i = 1, lev_int-1
                l_min_act(ilev) = l_min_act(ilev) + adapt_par%lev_delta(i)
             enddo
             l_min_act(ilev) = l_min_act(ilev) + (l_lev_int-1)*l_delta_act(ilev)
             l_max_act(ilev) = l_min_act(ilev)
             
             l_min_act(ilev) = l_min_act(ilev) - l_delta_act_before/2.
             if (l_min_act(ilev) < adapt_par%lev_min) &
                 l_min_act(ilev) = adapt_par%lev_min
             l_max_act(ilev) = l_max_act(ilev) + l_delta_act(ilev)/2.
             if (l_max_act(ilev) > adapt_par%lev_max) &
                 l_max_act(ilev) = adapt_par%lev_max
             l_delta_act(ilev) = l_max_act(ilev) - l_min_act(ilev)
               
             !--------------------------------------------------------
             ! determine lev_min_act and lev_max_act
             !--------------------------------------------------------
             
             if (nintervals==1) then
                lev_min_act(ilev) = l_min_act(ilev)
                lev_max_act(ilev) = l_max_act(ilev)
                lev_delta_act(ilev) = lev_max_act(ilev) - lev_min_act(ilev)
             else
                lev_delta_act(ilev) = lev2levdelta(adapt_par%lev_grid, &
                     adapt_par%lev_delta, &
                     (l_min_act(ilev)+l_max_act(ilev))/2.)
                
                lev_min_act(ilev) = (l_min_act(ilev)+l_max_act(ilev))/2. - &
                     lev_delta_act(ilev)/2.
                lev_max_act(ilev) = (l_min_act(ilev)+l_max_act(ilev))/2. + &
                     lev_delta_act(ilev)/2.
                if (lev_min_act(ilev) < adapt_par%lev_min) &
                    lev_min_act(ilev)=adapt_par%lev_min
                if (lev_max_act(ilev) > adapt_par%lev_max) &
                    lev_max_act(ilev)=adapt_par%lev_max
                 
             endif
              
          else  ! grid_switch==2   ! shifted levels (old version)
              
             !--------------------------------------------------------
             !  define parameters for grid adaptation
             !--------------------------------------------------------
             call set_adapt_par(adapt_par, ilev, grid_switch)
             
             !--------------------------------------------------------
             ! determine l_min_act and l_max_act
             !--------------------------------------------------------
             if (ilev==1) then
                l_min_act(ilev) = adapt_par%lev_min
                l_max_act(ilev) = adapt_par%lev_min + adapt_par%lev_delta(1)/2.
             elseif (ilev==nlevs+1) then
                l_min_act(ilev) = adapt_par%lev_grid(nlevs)
                l_max_act(ilev) = adapt_par%lev_max
             else
                l_min_act(ilev) = adapt_par%lev_grid(ilev-1)
                l_max_act(ilev) = adapt_par%lev_grid(ilev)
             endif
             l_delta_act(ilev) = l_max_act(ilev) - l_min_act(ilev)
             
             !--------------------------------------------------------
             ! determine lev_min_act and lev_max_act
             !--------------------------------------------------------
             lev_min_act(ilev) = l_min_act(ilev)
             lev_max_act(ilev) = l_max_act(ilev)
             lev_delta_act(ilev) = lev_max_act(ilev) - lev_min_act(ilev)
             
          endif
          
       endif
       
       ! write (*,*) &
       !      'ilev, l_min_act, l_max_act, l_delta_act: ', &
       !       ilev, l_min_act(ilev), l_max_act(ilev), l_delta_act(ilev)
       
    end do
    
    if (rank==0) then
       do ilev = 1, l_nlevs
          write (*,'(A,I5,2F12.4)') 'ilev, Boundaries (K): ', &
               ilev, lev_min_act(ilev), lev_max_act(ilev) !, lev_delta_act(ilev)
          !write (*,'(A,2F12.4)') &
          !     'l_min/max: ', &
          !      l_min_act(ilev), l_max_act(ilev)!, l_delta_act(ilev)
       enddo
    endif

!!$          write(*,*) 'adapt_par%fac_min', adapt_par%fac_min
!!$          write(*,*) 'adapt_par%fac_max', adapt_par%fac_max
!!$          write(*,*) 'adapt_par%r_min_c', adapt_par%r_min_c
!!$          write(*,*) 'adapt_par%r_max_c', adapt_par%r_max_c
!!$          write(*,*) 'adapt_par%r_lim_c_inside', adapt_par%r_lim_c_inside
!!$          write(*,*) 'adapt_par%r_lim_c_outside', adapt_par%r_lim_c_outside
!!$          write(*,*) 'adapt_par%r_min_h', adapt_par%r_min_h
!!$          write(*,*) 'adapt_par%r_max_h', adapt_par%r_max_h
!!$          write(*,*) 'adapt_par%r_lim_h_inside', adapt_par%r_lim_h_inside
!!$          write(*,*) 'adapt_par%r_lim_h_outside', adapt_par%r_lim_h_outsid e
    
  END SUBROUTINE clamsmix_set_bounds


END MODULE messy_clamsmix
