Module messy_clamschem_inhet_ken

Contains

subroutine inhet_ken (chemspecarr)

  ! Routine for initializing heterogeneous chemistry module written by 
  ! Ken Carslaw. Adopted into CLaMS Chemistry program by J.-U.Grooss 1998-1999
  !
  ! 19.03.2002 Initial state of aerosol taken from initial.nc file (jug)

  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  use messy_clams_global,         only: PREC, rank, species_type
  use messy_clamschem_global,     only: ntraj, iodump, ftr, missing_index, &
                                        densaero, aerh2so4, wt, ar, parth, shindex
  use messy_clamschem_cirrus_clim,only: densice_cirrus,rhice_freeze_clim
  use messy_clamschem_hetero_shi, only: liquid, phase
  use messy_clamschem_globalhet,  only: liq_sdist_sigma, sat_meltallowed, param_nat_HR, &
                                        saturation_criteria, transform, gamma, &
                                        allice, allnat, cice_from_clim, &
                                        teold, prold, told, pressold, &
                                        densnat_old, densice_old, denssat_old, &
                                        cnat_old, cice_old, vliq_save, &
                                        astate, lstate, log_state, &
                                        kstate, laststate, liqtest, &
                                        mixedpop, natcore, &
                                        liquids, satmelting, &
                                        densnat, densice, denssat, &
                                        aer_h2so4_default, densaero_default, &
                                        cnatinit, ciceinit
  use messy_clamschem_asad_mod,    only: p, t, tnd
  USE messy_clamschem_asad_totnud, ONLY: ASAD_TOTNUD

  implicit none      

  TYPE(species_type), DIMENSION(:) :: chemspecarr
 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------

  integer    :: itr 
  real(PREC) :: temp, press, chcl, chocl, ch2o, cclno3, cbrno3, chbr, chobr, &
         chno3, parthcl, parthcl0, parthno3, parthno30, parth2o, parth2o0, &
         parthbr, parthbr0, aliq, asat, anat, aice, asatliq, wts, wtn, wtcl, &
         wtbr, wthocl, wthobr, sulreduction,dens,snat

  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: pi = 3.1415926 
  !real(PREC), parameter :: dum1 = 0.0 

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer    :: k, j
  real(PREC) :: dum, vliq,radius
  real(PREC) :: help_prec
  !-----------------------------------------------
  !
  !
  !******************************************

!!!!! Jetzt aus clamschem.nml einlesen:

!!$  ! read in heterogeneous chemistry and micrphysics parameters
!!$  open(1, file='data/hetinit6.dat', status='OLD', position='asis') 
!!$  read (1, *) 
!!$  read (1, *) liq_sdist_sigma, densaero_default, aer_h2so4_default, ciceinit, cnatinit 
!!$  allice = (ciceinit == 99.0)
!!$  allnat = (cnatinit == 99.0)
!!$  read (1, *) 
!!$  read (1, *) k 
!!$  sat_meltallowed = (k == 1)
!!$  read (1, *) 
!!$  read (1, *) k 
!!$  param_nat_HR = (k == 1)
!!$  read (1, *) 
!!$  do j = 1, 4 
!!$     read (1, *) (transform(j,k),k=1,4) 
!!$  end do
!!$  read (1, *) 
!!$  do j = 1, 4 
!!$     read (1, *) (saturation_criteria(j,k),k=1,4) 
!!$  end do
!!$  read (1, *) 
!!$  read (1, *) (dum,gamma(j),j=1,33) 
!!$  close(1) 

  if (rank==0 .and. iodump) then
     write (*,*) 'Matrizen in inhet_ken:'
     write (*,*) 'transform:',transform
     write (*,*) 'saturation:',saturation_criteria
  endif


  ! allocate all hetero arrays

  allocate(teold(ntraj), prold(ntraj), densnat_old(ntraj), densice_old(ntraj),  &
           denssat_old(ntraj), cnat_old(ntraj),cice_old(ntraj),vliq_save(ntraj), &
           astate(ntraj),lstate(ntraj), log_state(4,ntraj) )

  cice_from_clim = (ciceinit==0.)
  if (cice_from_clim) then 
     write(*,*)'****************************************************'
     write(*,*)'Saturation threshold from Kraemer et al. climatology'
!!!!! ACHTUNG: DOUBLE <=> REAL
!!$     dum=rhice_freeze_clim(0.d0)
!!$     radius=densice_cirrus(0.d0,-100.d0)
     dum=rhice_freeze_clim(0.d0)
     radius=densice_cirrus(0.d0,-100.d0)
     write(*,'(a30,f5.1,a3)')'Constant Ice particle radius:',radius*1E4,' um'
     write(*,*)'****************************************************'
  endif
  densaero(:ntraj)= densaero_default  ! initial guess, may be overwritten by get_astate
  aerh2so4(:ntraj)=aer_h2so4_default  ! initial guess, may be overwritten by get_astate
  call get_astate (chemspecarr)

  lstate(:ntraj) = 0 
  
  ! set initial surface areas to zero (they will be calculated later)
  aice=0.0
  anat=0.0
  aliq=0.0
  asat=0.0
  asatliq=0.0

  ! calculate total number density --> common /cmphys/
  call asad_totnud (ntraj)

  ! loop over air parcels
  do itr=1,ntraj

     if (t(itr) < 150. .or. p(itr) <= 0.) missing_index(itr)=.true.
     if (missing_index(itr)) cycle

     ! TRANSFER INPUT VARIABLES TO INTERNAL PARAMETERS
     temp   = t(itr)
     press  = p(itr)/100.
     chcl   = ftr(itr,shindex(1))*tnd(itr)    ! HCl in molec/ccm
     chno3  = ftr(itr,shindex(2))*tnd(itr)    ! HNO3 in molec/ccm
     ch2o   = ftr(itr,shindex(3))*tnd(itr)    ! H2O  in molec/ccm
     chbr   = ftr(itr,shindex(4))*tnd(itr)    ! HBr  in molec/ccm
     chocl  = ftr(itr,shindex(5))*tnd(itr)    ! HOCl in molec/ccm
     cclno3 = ftr(itr,shindex(6))*tnd(itr)    ! ClONO2 in molec/ccm
     chobr  = ftr(itr,shindex(7))*tnd(itr)    ! HOBr in molec/ccm
     cbrno3 = ftr(itr,shindex(8))*tnd(itr)    ! BrONO2 in molec/ccm
     parthcl  = 1.0                   ! (old) gas fraction of HCl (0<parthcl<1)
     parthno3 = 1.0                   ! (old) gas fraction of HNO3 (0<part<1)
     parth2o  = 1.0                   ! (old) gas fraction of H2O (0<part<1)
     parthbr  = 1.0                   ! (old) gas fraction of HBr (0<part<1) 
     
     ! Store old values of het and output variables
     parthcl0 = parthcl 
     parthno30 = parthno3 
     parth2o0 = parth2o 
     parthbr0 = parthbr 


     ! read in initial particle size parameters and matrix for
     ! transformation of phases
     pressold = press 
     told = temp 
     kstate = astate(itr) 
     laststate = kstate 
     satmelting = .false. 
     liqtest = .false. 
 
 
     parth2o = 1.0 
     parthno3 = 1.0 
     ! prevent unrealistic values of quantities that are divided
     chcl = abs(chcl) + 1.0 
     chobr = abs(chobr) + 1.0 
     chocl = abs(chocl) + 1.0 
     chbr = abs(chbr) + 1.0 
     !
     liquids =(kstate == 1)
     denssat=denssat_old(itr)
     densnat=densnat_old(itr)
     densice=densice_old(itr)
     dens=densaero(itr)
     aer_h2so4_default=aerh2so4(itr)

     mixedpop=(kstate >= 3 .or. (kstate == 2 .and. denssat < dens*0.999))


     if (kstate==1 .or. kstate==2 .and. mixedpop) then 
        sulreduction = 1.0 - (densnat + denssat + densice)/dens
        call liquid (press, temp, chcl, chbr, chocl, chobr, chno3, ch2o, parthcl, &
             parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
             aer_h2so4_default, vliq, satmelting, sulreduction) 
 
     endif
     if (transform(1,2)==1.0 .and. saturation_criteria(1,2)==1.0) then 
        if (sat_meltallowed) then 
           write (6, *) &
      'SAT MELTING UPON COOLING NOT ALLOWED WHEN TRANSFORM(1,2)=1. RESET TO 0' 
           sat_meltallowed = .false. 
        endif
     endif
     !
     ! DETERMINE NEW PARTICLE PHASE
!!!!! ACHTUNG: DOUBLE <=> REAL
!!$     call phase (press, temp, 0.0d0, chno3, ch2o, chcl, aer_h2so4_default, parthno3, &
!!$       parth2o, transform, saturation_criteria, densnat, densice, denssat, ciceinit, cnatinit, &
!!$       dens, satmelting, liquids, sat_meltallowed,snat)
     help_prec = 0.0
     if (temp < 150. .or. press <= 0.) write(*,*) 'warning inhet_ken: temp, press =',temp, press

     call phase (press, temp, help_prec, chno3, ch2o, chcl, aer_h2so4_default, parthno3, &
       parth2o, transform, saturation_criteria, densnat, densice, denssat, ciceinit, cnatinit, &
       dens, satmelting, liquids, sat_meltallowed,snat)

 
     ! ================================================================================
     ! transfer results of output array for transfer to ASAD:
     parth(itr,1,1) = parthcl 
     parth(itr,1,2) = parthno3 
     parth(itr,1,3) = parth2o 
     parth(itr,1,4) = parthbr 
     parth(itr,2,1) = parthcl0 
     parth(itr,2,2) = parthno30 
     parth(itr,2,3) = parth2o0 
     parth(itr,2,4) = parthbr0 
     wt(itr,1)      = wts      ! Weight fraction of Sulphur
     wt(itr,2)      = wtn      ! Weight fraction of Nitrogen
     wt(itr,3)      = wtcl     ! Weight fraction of Chlorine
     wt(itr,4)      = wtbr     ! Weight fraction of Bromine
     wt(itr,5)      = wthocl   ! Weight fraction of HOCl
     wt(itr,6)      = wthobr   ! Weight fraction of HOBr
     ar(itr,1)      = aice     ! surface area ice
     ar(itr,2)      = anat     ! surface area NAT 
     ar(itr,3)      = aliq     ! surface area liq
     ar(itr,4)      = asat     ! surface area SAT
     ar(itr,5)      = asatliq  ! surface area SATLIQ

     astate(itr)    = kstate 
     lstate(itr)    = laststate 
     teold(itr)     = Told
     Prold(itr)     = pressold
     densnat_old(itr) = densnat
     densice_old(itr) = densice
     denssat_old(itr) = denssat      
     ! variables cnatinit and ciceinit are set for the first 
     ! air parcel and not changed within inhet:
     cnat_old(itr)    = cnatinit
     cice_old(itr)    = ciceinit
     log_state(1,itr) = mixedpop
     log_state(2,itr) = natcore
     log_state(3,itr) = satmelting
     log_state(4,itr) = liquids
  enddo
  return  
end subroutine inhet_ken


subroutine get_astate  (chemspecarr)

  use messy_clams_global,         only: species_type, nchemspec
  use messy_clamschem_globalhet
  use messy_clamschem_global,     only: rank, ntraj, jpnl_offset, ipart, iodump, &
                                        densaero, aerh2so4

  implicit none

  TYPE(species_type), DIMENSION(:) :: chemspecarr

  integer           :: i 
  integer           :: start, end



  denssat_old(:)=0.0
  densnat_old(:)=0.0
  densice_old(:)=0.0
  

  start = jpnl_offset(ipart)
  !end = start + jpnl_count(ipart) - 1
  end = start + ntraj - 1

  !ju_jg_20201118 ---  densaero and aerh2so4 should not be set to zero
  !                    therefore zero values in the initialisation should
  !                    not be transferred
  DO i = 1, nchemspec  
     SELECT CASE (TRIM(chemspecarr(i)%name))
     CASE ("N(aero)")
        where(chemspecarr(i)%values(start:end) > 0.) &
             densaero(1:ntraj) = chemspecarr(i)%values(start:end)
     CASE ("N(SAT)")
        where(chemspecarr(i)%values(start:end) > 0.) &
             denssat_old(1:ntraj) = chemspecarr(i)%values(start:end)
     CASE ("N(NAT)")
        where(chemspecarr(i)%values(start:end) > 0.) &
             densnat_old(:ntraj) = chemspecarr(i)%values(start:end)
     CASE ("N(ice)")
        where(chemspecarr(i)%values(start:end) > 0.) &
             densice_old(:ntraj) = chemspecarr(i)%values(start:end)
     CASE ("aer_H2SO4")
        where(chemspecarr(i)%values(start:end) > 0.) &
             aerh2so4(1:ntraj) = chemspecarr(i)%values(start:end)
     END SELECT
  ENDDO

  ! Vorbesetzung:
  astate(:ntraj) = 1          !set initial particle phase to liquid 

  where(denssat_old /= 0.0) astate = 2
  where(densnat_old /= 0.0) astate = 3
  where(densice_old /= 0.0) astate = 4
  
  if (iodump .and. rank==0) then
     write(*,*)'get_astate: ',count(astate == 1),'liq air parcels' ,rank 
     write(*,*)'get_astate: ',count(astate == 2),'SAT air parcels' ,rank
     write(*,*)'get_astate: ',count(astate == 3),'NAT air parcels' ,rank
     write(*,*)'get_astate: ',count(astate == 4),'ice air parcels' ,rank
  endif 

  return 

end subroutine get_astate


End Module messy_clamschem_inhet_ken
