!****************************************************************************** 
!
! Module compute_traj
! 
! Written by Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jun 29 09:19:16 2016
!
! subroutine compute_trajectory
!
!****************************************************************************** 

Module messy_clamssedi_traj

CONTAINS


  !************************************************************************
  !
  !************************************************************************
  subroutine compute_trajectory (lcoupled)

    USE messy_clamssedi_global,     ONLY: particles, nparticle_max, nparticles, &
                                          udt_sedi, vdt_sedi, wdt_sedi, & 
                                          leveldt_sedi, dlevdzdt_sedi, &
                                          timestep_sedi, part_ind_start, part_ind_end
    USE messy_clams_global,         ONLY: prec, dp, rank, ntasks, &
                                          mdi, eps, pi, nx, ny, nz, &
                                          longrid_rad, latgrid_rad, levelgrid, &
                                          loglev, &
                                          ufut, vfut, wfut, levelfut, dlevdzfut, &
                                          fut_year, fut_month, fut_day, fut_sec
    USE messy_clams_tools_dateconv, ONLY: datsec
    USE messy_clamstraj_calc3d,     ONLY: rungek3d
    USE messy_clams_global,         ONLY: delta_time, &
                                          YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

    IMPLICIT NONE

    logical :: lcoupled    

    ! winds at time t
    real(DP), dimension(:,:,:),pointer :: u0, v0, w0, dlevdz0, level0   
    ! winds at time t+delta_t/2    
    real(DP), dimension(:,:,:),pointer :: umid, vmid, wmid,dlevdzmid, levelmid 
    
    integer    :: iyear, imonth, iday, isec
    integer    :: k, j, i
    logical    :: calc(nparticle_max)
    real(prec) :: dt
    real(prec) :: pslat
    real(dp)   :: future, past
    
    pslat = pi * .4  

    calc = .TRUE.

    dt = real(timestep_sedi)

    ! current time
    iyear    = YEAR
    imonth   = MONTH
    iday     = DAY
    isec     = HOUR*3600+MINUTE*60+SECOND

    ! the function datsec returns the date in seconds from 0,1,1,1960               
    past   = datsec(isec,iday,imonth,iyear)
    if (.not. lcoupled) then
       future = datsec(fut_sec,fut_day,fut_month,fut_year)
    else
       future = past + delta_time
    endif

!!!!!
    if (rank==0) then
       write (*,*) 'in compute_trajectory:'
       write (*,*) 'current time:',iyear,imonth,iday,isec
       write (*,*) 'future time:',fut_year,fut_month,fut_day,fut_sec
    endif

    allocate (u0(nx,ny,nz))
    allocate (v0(nx,ny,nz))
    allocate (w0(nx,ny,nz))
    allocate (level0(nx,ny,nz))
    allocate (dlevdz0(nx,ny,nz))
    allocate (umid(nx,ny,nz))
    allocate (vmid(nx,ny,nz))
    allocate (wmid(nx,ny,nz))
    allocate (levelmid(nx,ny,nz))
    allocate (dlevdzmid(nx,ny,nz))


    DO k=1,nz
       DO j=1,ny      
          DO i=1,nx 
             IF (ABS((udt_sedi(i,j,k)-MDI)/MDI)<=eps .OR. ABS((vdt_sedi(i,j,k)-MDI)/MDI)<=eps &
                  .OR. ABS((wdt_sedi(i,j,k)-MDI)/MDI)<=eps) THEN
                u0(i,j,k) = MDI
                v0(i,j,k) = MDI
                w0(i,j,k) = MDI
                level0(i,j,k) = MDI
                dlevdz0(i,j,k) = MDI
                umid(i,j,k) = MDI 
                vmid(i,j,k) = MDI
                wmid(i,j,k) = MDI
                levelmid(i,j,k) = MDI
                dlevdzmid(i,j,k) = MDI
                udt_sedi(i,j,k)  = MDI 
                vdt_sedi(i,j,k)  = MDI 
                wdt_sedi(i,j,k)  = MDI 
                leveldt_sedi(i,j,k) = MDI
                dlevdzdt_sedi(i,j,k) = MDI
             ELSEIF (ABS((ufut(i,j,k)-MDI)/MDI)<=eps .OR. ABS((vfut(i,j,k)-MDI)/MDI)<=eps &
                  .OR. ABS((wfut(i,j,k)-MDI)/MDI)<=eps) THEN
                u0(i,j,k) = MDI
                v0(i,j,k) = MDI
                w0(i,j,k) = MDI
                level0(i,j,k) = MDI
                dlevdz0(i,j,k) = MDI
                umid(i,j,k) = MDI 
                vmid(i,j,k) = MDI
                wmid(i,j,k) = MDI
                levelmid(i,j,k) = MDI
                dlevdzmid(i,j,k) = MDI
                udt_sedi(i,j,k)  = MDI 
                vdt_sedi(i,j,k)  = MDI 
                wdt_sedi(i,j,k)  = MDI 
                leveldt_sedi(i,j,k) = MDI
                dlevdzdt_sedi(i,j,k) = MDI
             ELSE
                u0(i,j,k) = udt_sedi(i,j,k)
                v0(i,j,k) = vdt_sedi(i,j,k)
                w0(i,j,k) = wdt_sedi(i,j,k)
                level0(i,j,k) = leveldt_sedi(i,j,k)
                dlevdz0(i,j,k) = dlevdzdt_sedi(i,j,k)
                if (.not. lcoupled) then
                   umid(i,j,k) = &
                     (0.5*dt*ufut(i,j,k) + &
                     (future-past-0.5*dt)*u0(i,j,k))/(future-past)
                   vmid(i,j,k) = & 
                     (0.5*dt*vfut(i,j,k) +  &
                     (future-past-0.5*dt)*v0(i,j,k))/(future-past)
                   wmid(i,j,k) = &
                     (0.5*dt*wfut(i,j,k) + &
                     (future-past-0.5*dt)*w0(i,j,k))/(future-past)
                   levelmid(i,j,k) = &
                     (0.5*dt*levelfut(i,j,k) + &
                     (future-past-0.5*dt)*level0(i,j,k))/(future-past)
                   dlevdzmid(i,j,k) = &
                     (0.5*dt*dlevdzfut(i,j,k) + &
                     (future-past-0.5*dt)*dlevdz0(i,j,k))/(future-past)
                   udt_sedi(i,j,k) = &
                     (    dt*ufut(i,j,k) + &
                     (future-past-    dt)*u0(i,j,k))/(future-past)
                   vdt_sedi(i,j,k) = &
                     (    dt*vfut(i,j,k) + &
                     (future-past-    dt)*v0(i,j,k))/(future-past)
                   wdt_sedi(i,j,k) = &
                     (    dt*wfut(i,j,k) + &
                     (future-past-    dt)*w0(i,j,k))/(future-past)
                   leveldt_sedi(i,j,k) = &
                     (    dt*levelfut(i,j,k) + &
                     (future-past-dt)*level0(i,j,k))/(future-past)
                   dlevdzdt_sedi(i,j,k) = &
                     (    dt*dlevdzfut(i,j,k) + &
                     (future-past-dt)*dlevdz0(i,j,k))/(future-past)
                endif
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    if (lcoupled) then
       umid = u0
       vmid = v0
       wmid = w0
       levelmid = level0
       dlevdzmid = dlevdz0
       udt_sedi = u0
       vdt_sedi = v0
       wdt_sedi = w0
       leveldt_sedi = level0
       dlevdzdt_sedi = dlevdz0
    endif

    !-------------------------------------------------------------
    ! Runge kutta scheme to advect particles in three dimensions
    !-------------------------------------------------------------
    

!!!! Testen:: anstelle der kompletten Felder (nparticle_max) nur welche der Laenge 
!!!!          nparticles uebergeben.
!!!!          Noch besser waeren evtl sogar Teilfelder von part_ind_start bis part_ind_end.
    CALL rungek3d (particles%lon, particles%lat, particles%lev, calc, particles%tsv, &
                   u0, v0, w0, level0, dlevdz0,  &
                   umid, vmid, wmid, levelmid, dlevdzmid,   &
                   udt_sedi, vdt_sedi, wdt_sedi, leveldt_sedi, dlevdzdt_sedi,   &    
                   nx, ny, nz,    &
                   longrid_rad, latgrid_rad, levelgrid,   &
                   nparticle_max, part_ind_start, part_ind_end, dt, pslat, loglev, .true.)  
    
    deallocate (u0, v0, w0, level0, dlevdz0)
    deallocate (umid, vmid, wmid, levelmid, dlevdzmid)

  end subroutine compute_trajectory

End Module messy_clamssedi_traj
