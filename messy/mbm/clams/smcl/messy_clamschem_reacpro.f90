Module messy_clamschem_reacpro

Contains

subroutine breac(ibr, lab, name, reactants, products, lpro, lrec) 

  use messy_clamschem_asad_mod, only: spb, jpspb
  use messy_clamschem_global,   only: iodump, iodumpo, slenb

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: ibr 
  character , intent(in) :: lab 
  character , intent(out) :: name*30 

  integer , intent(out) :: lpro 
  integer , intent(out) :: lrec 
  character , intent(out) :: reactants*30 
  character , intent(out) :: products*30 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: il, ib 
!-----------------------------------------------
  if (iodumpo) write (6, *) 'spb=', (spb(ibr,il),il=1,jpspb) 
  do ib = 1, jpspb - 1 
     slenb(ib) = index(spb(ibr,ib),' ') - 1 
     if(iodumpo) write (6, *) ib, slenb(ib), spb(ibr,ib) 
  end do
  name = lab//'('//spb(ibr,1)(1:slenb(1))//'+'//spb(ibr,2)(1:slenb(2))//')' 
  reactants = spb(ibr,1)(1:slenb(1))//'+'//spb(ibr,2)(1:slenb(2)) 
  write (products, *) (spb(ibr,il)(1:slenb(il)),'+',il=3,jpspb - 2), spb(&
       ibr,jpspb-1)(1:slenb(jpspb-1)) 
  lpro = index(products(2:30),' ') + 1 
  lrec = index(reactants(2:30),' ') + 1 
  if (iodump) then 
     write (6, *) 'name= ', name 
     write (6, *) 'reactants= ', reactants 
     write (6, *) 'products= ', products, lpro 
  endif
  return  
end subroutine breac 


 
subroutine treac(itr, lab, name, reactants, products, lpro, lrec, therm_flag) 

  use messy_clamschem_asad_mod, only: spt, jpspt
  use messy_clamschem_global,   only: iodump, iodumpo, slent

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: itr 
  integer , intent(out) :: lpro 
  integer , intent(out) :: lrec 
  integer , intent(out) :: therm_flag 
  character , intent(in) :: lab 
  character , intent(out) :: name*30 
  character , intent(out) :: reactants*30 
  character , intent(out) :: products*30 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: indm, il, it 
!-----------------------------------------------
 
  if (iodumpo) write (6, *) (spt(itr,il),il=1,jpspt) 
  
  do it = 1, jpspt - 1 
     slent(it) = index(spt(itr,it),' ') - 1 
     if(iodumpo) write (6, *) it, slent(it), spt(itr,it) 
  end do
  
  indm = index(spt(itr,2)(1:slent(2)),'m') 
  if (indm == 1) then 
     if (iodumpo) write (6, *) 'Thermal reaction detected' 
     name = lab//'('//spt(itr,1)(1:slent(1))//'+M)' 
     reactants = spt(itr,1)(1:slent(1))//'+M' 
     therm_flag = 1 
     write (products, *) spt(itr,3)(1:slent(3)), '+', spt(itr,4) 
  else 
     name = lab//'('//spt(itr,1)(1:slent(1))//'+'//spt(itr,2)(1:slent(2))//'+M)' 
     reactants = spt(itr,1)(1:slent(1))//'+'//spt(itr,2)(1:slent(2))//'+M' 
     therm_flag = 0 
     write (products, *) spt(itr,3)(1:slent(3)) 
  endif
  
  lpro = index(products(2:30),' ') + 1 
  lrec = index(reactants(2:30),' ') + 1 
  if (iodump) then 
     write (6, *) 'name= ', name 
     write (6, *) 'reactants= ', reactants 
     write (6, *) 'products= ', products, lpro 
  endif
  return  
end subroutine treac 


 
subroutine jreac(ijr, lab, name, reactants, products, lpro, lrec) 

  use messy_clamschem_asad_mod, only: spj, jpspj
  use messy_clamschem_global,   only: iodump, iodumpo, slenj

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: ijr 
  integer , intent(out) :: lpro 
  integer , intent(out) :: lrec 
  character , intent(in) :: lab 
  character , intent(out) :: name*30 
  character , intent(out) :: reactants*30 
  character , intent(out) :: products*30 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: il, ij 
!-----------------------------------------------
  if (iodumpo) write (6, *) (spj(ijr,il),il=1,jpspj) 
  do ij = 1, jpspj - 1 
     slenj(ij) = index(spj(ijr,ij),' ') - 1 
     if(iodumpo) write (6, *) ij, slenj(ij), spj(ijr,ij) 
  end do
  name = lab//'('//spj(ijr,1)(1:slenj(1))//'+hv)' 
  reactants = spj(ijr,1)(1:slenj(1))//'+hv' 
  write (products, *) (spj(ijr,il)(1:slenj(il)),'+',il=3,jpspj - 2), &
       spj(ijr,jpspj-1)(1:slenj(jpspj-1)) 
  lrec = index(reactants(2:30),' ') + 1 
  lpro = index(products(2:30),' ') + 1 
  if (iodump) then 
     write (6, *) 'name= ', name 
     write (6, *) 'reactants= ', reactants 
     write (6, *) 'products= ', products, lpro 
  endif
  return  
end subroutine jreac 


 
subroutine hreac(ihr, lab, name, reactants, products, lpro, lrec) 

  use messy_clamschem_asad_mod, only: sph, jpsph
  use messy_clamschem_global,   only: iodump, iodumpo, slenh

  implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
  integer , intent(in) :: ihr 
  integer , intent(out) :: lpro 
  integer , intent(out) :: lrec 
  character , intent(in) :: lab 
  character , intent(out) :: name*30 
  character , intent(out) :: reactants*30 
  character , intent(out) :: products*30 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer :: il, ih 
!-----------------------------------------------
  if (iodumpo) write (6, *) (sph(ihr,il),il=1,jpsph) 
  do ih = 1, jpsph - 1 
     slenh(ih) = index(sph(ihr,ih),' ') - 1 
     if(iodumpo) write (6, *) ih, slenh(ih), sph(ihr,ih) 
  end do

  name = lab//'('//sph(ihr,1)(1:slenh(1))//'+'//sph(ihr,2)(1:slenh(2))//')' 
  reactants = sph(ihr,1)(1:slenh(1))//'+'//sph(ihr,2)(1:slenh(2)) 
  write (products, *) (sph(ihr,il)(1:slenh(il)),'+',il=3,jpsph - 2), sph(&
       ihr,jpsph-1)(1:slenh(jpsph-1)) 
  lpro = index(products(2:30),' ') + 1 
  lrec = index(reactants(2:30),' ') + 1 
  if (iodump) then 
     write (6, *) 'name= ', name 
     write (6, *) 'reactants= ', reactants 
     write (6, *) 'products= ', products, lpro 
  endif
  return  
end subroutine hreac 


subroutine check_dup_reac (status, allnames, name, n_names)

  ! check duplicate reactants

  implicit none

  CHARACTER(30), dimension(:), pointer :: allnames
  CHARACTER(30)                        :: name
  INTEGER                              :: status, n_reac, n_names

  integer :: lnam, ia, i, k, i_char
  character(30) :: namea, nameb
  logical       :: found

  status = 0 ! no error

  ia=ichar('a')

  lnam = index(name,' ') - 1 
  namea = name(1:lnam)//'a' 
!!$  if (n_names==1) allnames(n_names) = name
  allnames(n_names) = name
  do i = 1, n_names-1
     if (name==allnames(i) .or. namea==allnames(i)) then
        allnames(i) = namea
        do i_char = ia, ia+26
           nameb = name(1:lnam)//char(i_char)
           found = .false.
           do k = 1, n_names-1
              if (nameb==allnames(k)) found = .true.
              if (found .and. i_char==ia+26) status = -1  ! error
           enddo
           if (.not. found) exit
        enddo
        allnames(n_names) = nameb
!!$     else
!!$        allnames(n_names) = name
     endif
  enddo

end subroutine check_dup_reac

End Module messy_clamschem_reacpro
