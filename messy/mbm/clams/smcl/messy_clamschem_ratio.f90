Module messy_clamschem_ratio

Contains


!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine ratio(ip, mratio) 

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_asad_mod,                 only: speci, tnd, y, lvmr

  use messy_clamschem_global,   ONLY: ntraj, iodumpo

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer, intent(in) :: ip 
  real(kind=DP), dimension(ntraj),intent(out) :: mratio 
!-----------------------------------------------
  if (iodumpo) then 
     write (*, *) 'ratio: ip,=', ip
     write (6, *) 'species ', speci(ip) 
     write (6, *) 'TND(1)= ', tnd(1)
     write (6, *) 'Y(1,ip)= ', y(1,ip)
  endif
  if (lvmr) then 
     mratio = y(:ntraj,ip)/tnd(:ntraj)
  else 
     mratio = y(:ntraj,ip)
  endif
  if (iodumpo) write (6, *) 'MRATIO= ', ip, mratio(1)
  return  
end subroutine ratio 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine aero(values) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_asad_mod, only: tnd, lvmr

  use messy_clamschem_global,   ONLY: ntraj, iodumpo, &
                                      con, wt, ar, densaero, aerh2so4, &
                                      nhetspec, nhetpar
  use messy_clamschem_globalhet,only: denssat_old, densnat_old, densice_old,  &
                                      vliq_save


  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  real(kind=DP), intent(out) :: values(ntraj,nhetspec+nhetpar) 
!-----------------------------------------------

  ! values(:,1:nhetspec) => SPECARR 
  ! values(:,1:9) => SPECARR 

  if (lvmr) then 
     values(:ntraj,1) = con(:ntraj,1)/tnd(:ntraj) 
     values(:ntraj,2) = con(:ntraj,2)/tnd(:ntraj) 
     values(:ntraj,3) = con(:ntraj,3)/tnd(:ntraj) 
     values(:ntraj,4) = con(:ntraj,4)/tnd(:ntraj) 
  else 
     values(:ntraj,1) = con(:ntraj,1) 
     values(:ntraj,2) = con(:ntraj,2) 
     values(:ntraj,3) = con(:ntraj,3) 
     values(:ntraj,4) = con(:ntraj,4) 
  endif

  values(:ntraj,5)  = densaero(:ntraj)
  values(:ntraj,6)  = denssat_old(:ntraj) 
  values(:ntraj,7)  = densnat_old(:ntraj) 
  values(:ntraj,8)  = densice_old(:ntraj) 
  values(:ntraj,9)  = aerh2so4(:ntraj) 
  
  
  ! values(:,nhetspec+1:nhetspec+nhetpar) => HETPAR 
  ! values(:,10:20) => HETPAR 

  values(:ntraj,nhetspec+1) = wt(:ntraj,1)*100. 
  values(:ntraj,nhetspec+2) = wt(:ntraj,2)*100. 
  values(:ntraj,nhetspec+3) = wt(:ntraj,3)*100. 
  values(:ntraj,nhetspec+4) = wt(:ntraj,4)*100. 
  values(:ntraj,nhetspec+5) = wt(:ntraj,5)*100. 
  values(:ntraj,nhetspec+6) = wt(:ntraj,6)*100. 

  values(:ntraj,nhetspec+7) = ar(:ntraj,1) 
  values(:ntraj,nhetspec+8) = ar(:ntraj,2) 
  values(:ntraj,nhetspec+9) = ar(:ntraj,3) 
  values(:ntraj,nhetspec+10) = ar(:ntraj,4) 

  values(:ntraj,nhetspec+11) = vliq_save(:ntraj) 
 
  if (iodumpo) write (6, *) 'aero value', values(1,:) 
  return  

end subroutine aero 

End Module messy_clamschem_ratio
