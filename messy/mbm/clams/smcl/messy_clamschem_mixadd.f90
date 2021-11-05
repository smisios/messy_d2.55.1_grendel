Module messy_clamschem_mixadd

Contains

subroutine mixadd(specarr, &
                  BRATES, TRATES, JRATES, HRATES, &
                  BCONST, TCONST, JCONST, HCONST, HETPAR) 
  !
  ! Add data points into the new netCDF dataset
  !
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  ! MESSy Main
  use messy_main_constants_mem, only: DP

  ! ASAD
  use messy_clamschem_asad_mod_clams, only: jpctr, jpspec, &
                                            jpbk, jphk, jppj, jptk, &
                                            lhet, lphotol
  use messy_clamschem_asad_mod,       only: speci, tnd, emr, majors, nemit, nlemit

  USE messy_clams_global,       ONLY: rank, species_type, nchemspec

  USE messy_clamschem_global,   ONLY: rc_type, missing_value, ntraj,  &
                                      jpnl_count, jpnl_offset, ipart, &
                                      iodump, iodumpo, &
                                      fnames, ftrindex, nfamily, &
                                      ftr, missing_index, &
                                      const, rates, hetparam, emrates, &
                                      shindex, nhetspec, nhetpar, hetvar

  USE messy_clamschem_ratio,    ONLY: ratio, aero
  use messy_clamschem_drates,   only: get_bconst, get_tconst, &
                                      get_jconst, get_hconst, &
                                      get_brates, get_trates, &
                                      get_jrates, get_hrates

 
  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  TYPE(species_type), DIMENSION(:) :: specarr

  TYPE(species_type), DIMENSION(:) :: HETPAR

  TYPE(rc_type),      DIMENSION(:), POINTER :: BRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: TRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: JRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: HRATES

  TYPE(rc_type),      DIMENSION(:), POINTER :: BCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: TCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: JCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: HCONST


  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  INTEGER :: ispec
  integer :: start, end
  integer :: spid, tnid, rcode, vaid 
  integer ::  ip, if, ibr, itr, ijr, ihr, is, j , ihet
  real(DP), dimension(ntraj) :: values
  real(DP), dimension(ntraj,nhetspec+nhetpar) :: aerosol 
  character :: varname*20
  !-----------------------------------------------

  if (all(missing_index)) then  !Barbara
     write(*,*)"rank= ",rank,"ntraj=",ntraj, 'mixadd warning:  no valid output data'
     ! return 
  endif

  !if (iodump) write (6, *) 'Extract total number density' 
  values=tnd(:ntraj)
  where(missing_index) values=missing_value  
  ! TND to MESSy output?

  start = jpnl_offset(ipart)
  end = start + jpnl_count(ipart) - 1 
  if (iodump) write (*,*) 'mixadd: rank,start,end=',rank,start,end

  do ip = 1, jpspec 
     call ratio (ip, values)
     where(missing_index) values=missing_value   
!if (ip == 1) write(*,*) 'PE', rank, 'values1', values, 'ntraj', ntraj

     ! Write to  MESSy array => SPECARR(1:nchemspec-nhetspec)
     DO ispec = 1, nchemspec
        IF (TRIM(speci(ip)) == SPECARR(ispec)%name) THEN
           SPECARR(ispec)%values(start:end) = values
           if (rank==0 .and. iodump) &
                write(*,*) ip, 'speci mixadd ', speci(ip), SPECARR(ispec)%values(1)
        END IF
     END DO
  end do

!!!!!
  ! Do the same for families (not yet implemented in MESSy)    

  if (rates) then 
     if (rank==0 .and. iodump) write (6, *) 'Extract bimolecular rates' 
     do ibr = 1, jpbk 
        call get_brates (ibr, values) 
        if (rank==0 .and. iodumpo) write (6, *) ibr, values(1)
        where(missing_index) values=missing_value
        ! Write to MESSy array
        BRATES(ibr)%values(start:end) = values
     end do
     if (rank==0 .and. iodump) write (6, *) 'Extract trimolecular rates' 
     do itr = 1, jptk 
        call get_trates (itr, values) 
        where(missing_index) values=missing_value 
        ! Write to MESSy array
        TRATES(itr)%values(start:end) = values
     end do
 
     if (rank==0 .and. iodump) write (6, *) 'Extract photolysis rates' 
     do ijr = 1, jppj 
        if (.not. lphotol) cycle
        call get_jrates (ijr, values) 
        if (rank==0 .and. iodumpo) write (6, *) 'R(', ijr, ')= ', values(1)
        where(missing_index) values=missing_value
        ! Write to MESSy array
        JRATES(ijr)%values(start:end) = values
     end do
     if (lhet) then 
        if (rank==0 .and. iodump) write (6, *) 'Extract heterogeneous rates' 
        do ihr = 1, jphk 
           call get_hrates (ihr, values)
           where(missing_index) values=missing_value
           ! Write to MESSy array
           HRATES(ihr)%values(start:end) = values
       end do
     endif
  endif

  if (const) then 
     if (rank==0 .and. iodump) write (6, *) 'Extract bimolecular rate constants' 
     do ibr = 1, jpbk 
        call get_bconst (ibr, values) 
        if (iodumpo) write (6, *) ibr, values(1)
        where(missing_index) values=missing_value
        ! Write to MESSy array
        BCONST(ibr)%values(start:end) = values
     end do
       
     if (rank==0 .and. iodump) write (6, *) 'Extract trimolecular rate constants' 
     do itr = 1, jptk 
        call get_tconst (itr, values) 
        where(missing_index) values=missing_value
        ! Write to MESSy array
        TCONST(itr)%values(start:end) = values
     end do
     if (rank==0 .and. iodump) write (6, *) 'Extract photolysis rate constants' 
     do ijr = 1, jppj 
        if (.not. lphotol) cycle
        call get_jconst (ijr, values) 
        if (iodumpo) write (6, *) 'J(', ijr, ')= ', values(1)
        where(missing_index) values=missing_value
        ! Write to MESSy array
        JCONST(ijr)%values(start:end) = values
    end do
     if (lhet) then 
        if (rank==0 .and. iodump) write (6, *) &
             'Extract heterogeneous rates constants' 
        do ihr = 1, jphk 
           call get_hconst (ihr, values) 
           where(missing_index) values=missing_value 
           ! Write to MESSy array
           HCONST(ihr)%values(start:end) = values
        end do
     endif
  endif
     
  if (lhet) then 

     ! heterogeneous chemistry output:
     if (rank==0 .and. iodump) write (6, *) 'Extract condensed phase info' 
  
     call aero (aerosol)
        
     ! Write to MESSy array => SPECARR(nchemspec-nhetspec+1:nchemspec)
     if (rank==0) write (*,*) 'mixadd: nhetspec=',nhetspec
     do ihet = 1, nhetspec
        values = aerosol(:ntraj,ihet)
        where(missing_index) values = missing_value
        do ispec = nchemspec-nhetspec+1, nchemspec
           IF (TRIM(hetvar(ihet)%name) == SPECARR(ispec)%name) THEN
              SPECARR(ispec)%values(start:end) = values
              if (rank==0 .and. iodump) &
                   write(*,*) 'mixadd: ', SPECARR(ispec)%name, SPECARR(ispec)%values(1)
           END IF
        enddo
     end do
     
     ! Write to MESSy array => HETPAR
     if (hetparam) then
        if (rank==0) write (*,*) 'mixadd: nhetpar=',nhetpar
        do ihet = nhetspec+1, nhetspec+nhetpar
           values = aerosol(:ntraj,ihet)
           where(missing_index) values = missing_value
           do ispec = 1, nhetpar
              IF (TRIM(hetvar(ihet)%name) == HETPAR(ispec)%name) THEN
                 HETPAR(ispec)%values(start:end) = values
                 if (rank==0 .and. iodump) &
                      write(*,*) 'mixadd: ', HETPAR(ispec)%name, HETPAR(ispec)%values(1)
              END IF
           enddo
        enddo
     endif

     ! HCl, HNO3, H2O, HBr => sum of gasphase and condensed phase
 
     ! add condensed phase to normal mixing ratios
     ! temporary fix, should be replaced by initialization of condensed 
     ! phase and gas phase. J.U. Grooss, 17.05.99
     !
     ! condensed phase species (HCl, HONO2, H2O, HBr) cannot be member 
     ! of any family. majors() contains also the representation number 
     ! of the tracers in the species array...

     do is = 1, 4

        ip = majors(shindex(is))
        call ratio (ip, values)
        values = values + aerosol(:ntraj,is)
        where(missing_index) values=missing_value

        ! Write to  MESSy array => SPECARR
        DO ispec = 1, nchemspec
           IF (TRIM(speci(ip)) == SPECARR(ispec)%name) THEN
              SPECARR(ispec)%values(start:end) = values
              if (rank==0 .and. iodump) &
                   write(*,*) ip, 'speci mixadd ', speci(ip), SPECARR(ispec)%values(1)
           END IF
        END DO

     enddo

  endif  ! if (lhet)
 
  if (emrates) then 
     do j = 1, nemit  
        varname='em_'//trim(speci(nlemit(j)))
        values=emr(:ntraj,nlemit(j))
        ! Write to MESSy array...
     end do
  endif

  if (iodumpo) write (6, *) 'MIXADD complete ' 
  return  
end subroutine mixadd

End Module messy_clamschem_mixadd
