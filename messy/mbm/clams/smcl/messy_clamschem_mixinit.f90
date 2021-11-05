Module messy_clamschem_mixinit

Contains

subroutine mixinit(status, chemspecarr) 

  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------

  ! MESSy Main
  use messy_main_constants_mem, only: DP

  ! ASAD
  use messy_clamschem_asad_mod, only: speci, ctype, moffam, nodd, &
                                      madvtr, advt, family, &
                                      jpfm, jpif, jpsp, jpna, jpco, &
                                      fn2, fo2, fco2, fh2, tnd
  use messy_clamschem_asad_mod_clams, only: jpctr, jpspec
  USE messy_clamschem_asad_totnud,    ONLY: ASAD_TOTNUD

  ! CLaMS
  use messy_clams_global,       only: rank, species_type, nchemspec
  USE messy_clamschem_global,   ONLY: ntraj, jpnl_offset, jpnl_count, ipart, &
                                      iodump, nfamily, fnames, ftrindex , &
                                      ftr, ftr_ini, missing_index, &
                                      slt


  implicit none

  !-----------------------------------------------
  !   A r g u m e n t s
  !-----------------------------------------------
  INTEGER :: status
  TYPE(species_type), DIMENSION(:) :: chemspecarr

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: start, end
  integer :: is, if
  real(DP), allocatable :: values(:) 
!  real(DP), parameter ::zboltz = 1.3806e-23 * 1.0e6
  INTEGER :: i


  status = 0 ! no error
 
  allocate(values(ntraj))
 
  start = jpnl_offset(ipart)
  end = start + jpnl_count(ipart) - 1  
  !rite (*,*) 'mixinit: start,end=',start,end

  ! set jpco species to zero and constants to constant values
  do is = 1, jpspec 
     values = 0.0
     if (ctype(is)==jpco) then ! initialize constant mixing ratios
        select case (trim(speci(is)))
        case('N2')   
           values = fn2
        case('O2')   
           values = fo2
        case('CO2')  
           values = fco2
        case('H2')   
           values = fh2
        end select

        DO i = 1, nchemspec ! write in chemspecarr
           IF (TRIM(speci(is)) == chemspecarr(i)%name) THEN
              chemspecarr(i)%values(start:end) = values
           END IF
        END DO
        
     endif
  end do

!!!!! Familien fuer Ausgabe sichern ?!?
  do if = 1, nfamily
     if (rank==0 .and. iodump) write (*,*) 'if,fnames(if):',if,fnames(if)
  enddo

  if (iodump) write (6, *) 'ftr array set to zero' 
  ftr = 0.0 
     
!!!!! Evtl. vom init-file einlesen ?!?
  slt(:ntraj) = 0.
 
  ! calculate total number density TND of first timestep
  call asad_totnud (ntraj)


!!!!! Werden nirgendwo gespeichert ?!?
  values = tnd(:ntraj)

  !-----------------------------------------------------------------------------
  !   Now check to see if a member of a family and initialise  FTR array
  !
  !   Now identify the families within ftr and generate ftrindex
  !
  do is = 1, jpspec 

     if (ctype(is)==jpfm .or. ctype(is)==jpif) then 

        DO i = 1, nchemspec  
           IF (TRIM(speci(is)) == chemspecarr(i)%name) THEN

              !  changed 27.3.98 mk/rm ( nodd multiplied for correct total chlorine )
              ftr(:ntraj,moffam(is)) = ftr(:ntraj,moffam(is)) + &
                                       chemspecarr(i)%values(start:end)*nodd(is)

              if (rank==0 .and. iodump) then 
                 write (*,*) speci(is), 'identified with family ', family(is) 
                 write (*,*) speci(is),'  values(1)= ', chemspecarr(i)%values(start) 
              endif

           END IF
        END DO


!!$        write(*,*) 'Family concept not yet implemented for MESSy!'
!!$        status = 131
!!$        return
!!$
     else if (ctype(is) == jpsp) then 

        DO i = 1, nchemspec  

           IF (TRIM(speci(is)) == chemspecarr(i)%name) THEN

              ftr(:ntraj,madvtr(is)) = chemspecarr(i)%values(start:end)

              if (rank==0 .and. iodump) then
                 write (*,*) speci(is), 'identified with tracer ', advt(madvtr(is)) 
                 write (*,*) speci(is),'  values(1)= ', chemspecarr(i)%values(start) 
              endif
              
           END IF
        END DO

     else if (ctype(is)==jpna .or. ctype(is)==jpco) then 
        if (rank==0 .and. iodump) write (*,*) speci(is), &
             'identified with steady State or Constant ' 
     endif

  end do
 

!!!!! Familien fuer Ausgabe sichern ?!?


! Moeglichgeit der Initialisierung falls Anfangswert der Trajektorie noch
! ungueltig -- 
  if (any(missing_index) ) then 
     allocate(ftr_ini(ntraj,jpctr))
     ftr_ini=ftr
  end if

  deallocate(values)
  return  
end subroutine mixinit

end Module messy_clamschem_mixinit
