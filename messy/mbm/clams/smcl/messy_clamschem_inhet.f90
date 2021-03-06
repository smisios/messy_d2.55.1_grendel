Module messy_clamschem_inhet

Contains

subroutine inhet (chemspecarr)
  !
  !     subroutine for initialization of heterogeneous chemistry.
  !     version for het. chmistry only on liquid ternary solution particles
  !     numhet needs to be set to 33 for this version
  !

  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  use messy_clams_global,         only: rank, species_type
  use messy_clamschem_global,     only: ntraj, iodump, ftr, &
                                        con, chindex, shindex, &
                                        densaero, aerh2so4, parth, wt, ar
  use messy_clamschem_globalhet,  only: numhet
  use messy_clamschem_inhet_ken,  only: inhet_ken
  use messy_clamschem_asad_mod_clams, only: jpctr, jphk, jphkp1
  use messy_clamschem_asad_mod,       only: advt, sph

  implicit none

  TYPE(species_type), DIMENSION(:) :: chemspecarr

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: nmatch, i, ij, ik, ip1, ip2, iq1, iq2, ir1, ir2, is1, is2, isp, ish
  logical :: found, lrec, lpro 
  character, dimension(numhet,4) :: chet*10 
  character, dimension(8)        :: shet*10 
  !-----------------------------------------------
  !
 
  ! reactions on NAT
  data (chet(1,i),i=1,4)/ 'HCl', 'ClONO2', 'Cl2', 'HONO2'/  
  data (chet(2,i),i=1,4)/ 'H2O', 'ClONO2', 'HOCl', 'HONO2'/  
  data (chet(3,i),i=1,4)/ 'HOCl', 'HCl', 'Cl2', 'H2O'/  
  data (chet(4,i),i=1,4)/ 'N2O5', 'HCl', 'ClNO2', 'HONO2'/  
  data (chet(5,i),i=1,4)/ 'N2O5', 'H2O', 'HONO2', 'HONO2'/  
  data (chet(6,i),i=1,4)/ 'ClONO2', 'HBr', 'BrCl', 'HONO2'/  
  data (chet(7,i),i=1,4)/ 'BrONO2', 'HCl', 'BrCl', 'HONO2'/  
  data (chet(8,i),i=1,4)/ 'HBr', 'HOCl', 'BrCl', 'H2O'/  
  data (chet(9,i),i=1,4)/ 'HOBr', 'HCl', 'BrCl', 'H2O'/  
  data (chet(10,i),i=1,4)/ 'HOBr', 'HBr', 'Br2', 'H2O'/  
  data (chet(11,i),i=1,4)/ 'BrONO2', 'H2O', 'HOBr', 'HONO2'/  
  
  ! reactions on ICE
  data (chet(12,i),i=1,4)/ 'HCl', 'ClONO2', 'Cl2', 'HONO2'/  
  data (chet(13,i),i=1,4)/ 'H2O', 'ClONO2', 'HOCl', 'HONO2'/  
  data (chet(14,i),i=1,4)/ 'HOCl', 'HCl', 'Cl2', 'H2O'/  
  data (chet(15,i),i=1,4)/ 'N2O5', 'HCl', 'ClNO2', 'HONO2'/  
  data (chet(16,i),i=1,4)/ 'N2O5', 'H2O', 'HONO2', 'HONO2'/  
  data (chet(17,i),i=1,4)/ 'ClONO2', 'HBr', 'BrCl', 'HONO2'/  
  data (chet(18,i),i=1,4)/ 'BrONO2', 'HCl', 'BrCl', 'HONO2'/  
  data (chet(19,i),i=1,4)/ 'HBr', 'HOCl', 'BrCl', 'H2O'/  
  data (chet(20,i),i=1,4)/ 'HOBr', 'HCl', 'BrCl', 'H2O'/  
  data (chet(21,i),i=1,4)/ 'HOBr', 'HBr', 'Br2', 'H2O'/  
  data (chet(22,i),i=1,4)/ 'BrONO2', 'H2O', 'HOBr', 'HONO2'/  
  
  ! reactions on liquid aerosol
  data (chet(23,i),i=1,4)/ 'HOCl', 'HCl', 'Cl2', 'H2O'/  
  data (chet(24,i),i=1,4)/ 'HCl', 'ClONO2', 'Cl2', 'HONO2'/  
  data (chet(25,i),i=1,4)/ 'H2O', 'ClONO2', 'HOCl', 'HONO2'/  
  data (chet(26,i),i=1,4)/ 'N2O5', 'H2O', 'HONO2', 'HONO2'/  
  data (chet(27,i),i=1,4)/ 'HOBr', 'HCl', 'BrCl', 'H2O'/  
  data (chet(28,i),i=1,4)/ 'HOBr', 'HBr', 'Br2', 'H2O'/  
  data (chet(29,i),i=1,4)/ 'HBr', 'HOCl', 'BrCl', 'H2O'/  
  data (chet(30,i),i=1,4)/ 'BrONO2', 'H2O', 'HOBr', 'HONO2'/  
  
  ! reactions on SAT
  data (chet(31,i),i=1,4)/ 'HCl', 'ClONO2', 'Cl2', 'HONO2'/  
  data (chet(32,i),i=1,4)/ 'H2O', 'ClONO2', 'HOCl', 'HONO2'/  
  data (chet(33,i),i=1,4)/ 'N2O5', 'H2O', 'HONO2', 'HONO2'/  
  
  data (shet(i),i=1,8)/ 'HCl', 'HONO2', 'H2O', 'HBr', 'HOCl', 'ClONO2', &
       'HOBr', 'BrONO2'/  
 
  allocate (shindex(8))
  allocate (chindex(jphkp1,4))

  allocate (con(ntraj,4))
  allocate (densaero(ntraj))
  allocate (aerh2so4(ntraj))

  allocate (parth(ntraj,2,4))
  allocate (wt(ntraj,6))
  allocate (ar(ntraj,5))

!!$  allocate (sedinucl(ntraj))

 
  do ij = 1, jphk 
     found = .FALSE. 
  !jug: numhet changed to 33; now multiple matches are possible for different 
  !     substrates here, the last match will be used for the comparison
  !     default value for chindex is numhet+1 (currently 34) if the reaction
  !     takes not place on all 4 substrates
  !
     chindex(ij,:) = numhet + 1 
     
     !     nmatch: number of matches
     nmatch = 0 
 
     do ik = 1, numhet 
        ir1 = index(sph(ij,1),chet(ik,1)) 
        ir2 = index(sph(ij,2),chet(ik,2)) 
        is1 = index(sph(ij,2),chet(ik,1)) 
        is2 = index(sph(ij,1),chet(ik,2)) 
        ip1 = index(sph(ij,3),chet(ik,3)) 
        ip2 = index(sph(ij,4),chet(ik,4)) 
        iq1 = index(sph(ij,4),chet(ik,3)) 
        iq2 = index(sph(ij,3),chet(ik,4)) 
        lrec = ir1*ir2>0 .or. is1*is2>0 
        lpro = ip1*ip2>0 .or. iq1*iq2>0 
        if (.not.(lrec .and. lpro)) cycle  
        nmatch = nmatch + 1 
        found = .TRUE. 
        chindex(ij,nmatch) = ik 
     end do
     if (nmatch == 0) then 
        if (iodump .and. rank==0) write (6, *) 'no match for ', (sph(ij,i),i=1,5), ij 
     else 
        if (iodump .and. rank==0) write (6, *) nmatch, ' matches for ', (sph(ij,i),i=1,5),&
              ij, (chindex(ij,i),i=1,nmatch) 
     endif
  end do 

  !     ---- now match up scale species using the just computed chindex...
  !     scale_index determination left out:  not needed for hetero6
  if (iodump .and. rank==0) write (*,*) 'match up soluble species' 
  do ish = 1, 8 
     found = .FALSE. 
     do isp = 1, jpctr 
        if (trim(advt(isp)) .ne.trim(shet(ish))) cycle  
        shindex(ish) = isp 
        found = .TRUE. 
     end do
     if (.not.found) then 
        if (iodump .and. rank==0) write (*,*) 'no match for ', shet(ish) 
        shindex(ish) = -1 
     else 
        if (iodump .and. rank==0) write (*,*) '   match for ', shet(ish) 
     endif
  end do

  con=0. ! initialize concentrations of condensed phase concentrations

  !
  !     Now initialize the hetero procedure.
  !
  call inhet_ken (chemspecarr)

  return  
end subroutine inhet 

end Module messy_clamschem_inhet
