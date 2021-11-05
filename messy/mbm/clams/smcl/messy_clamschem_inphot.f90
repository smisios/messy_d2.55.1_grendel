Module messy_clamschem_inphot

contains

subroutine inphot 
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:09:34  11/20/98  
!...Switches: -ymhf               
!     ------
!     link between photolysis code by David Lary originally from Mainz Box model
!     (FACSIMILE) and ASAD code.
!
!     11.08.1997: updated with an alternative  name array (cphota) that
!     alternativly used species names in the photolysis file ratj.d,
!     e.g. Cl+O2 instead of ClOO. In the alternative list, only the first two
!     reactands are checked for a match.
!     (Jens-Uwe Grooss)
!     ------
!     27.11.97, updated to use rk rather than pj
!     to hold photolsis rates: necessary for ASAD 2
!     (Rolf Mueller)
!
!     28.09.98 added switch to call the Roeth photolysis table (iphot=2) instead
!     of the Lary code (iphot=1)                    (J. Ankenbrand /J.U. Grooss)
!
!     29.09.98 data statements for cphot/cphota replaced with regular assignment
!     (J.U. Grooss)
!     28.04.99 added paths Cl2O2->ClO+ClO and ClONO2->ClO+NO2 
!     for simple chemistry case  (jug)
!     
!     09.04.08 changed string comparisons from "index" to compare strings to
!      avoid ambiguous interpretation (like F11 and F113)  (jug)
!
!     =================================

  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_asad_mod_clams, only: jppj
  use messy_clamschem_asad_mod,       only: spj

  USE messy_clams_global,       ONLY: rank
  USE messy_clamschem_global,   ONLY: ipart, jpdim, ip_messy, iodump, &
                                      iphot, cindex
  USE messy_cmn_photol_mem
  USE messy_dissoc,             ONLY: davg 
  USE messy_clams_tools_utils,  ONLY: quick_sort_index

  implicit none 
 
!!$ jpdim and ip_messy moved to messy_clamschem_global.f90 
!!$  !-----------------------------------------------
!!$  !   L o c a l   P a r a m e t e r s
!!$  !-----------------------------------------------
!!$  integer, parameter :: jpdim = 48  ! number of implemented photolysis rates in dissoc
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
!!$  integer , dimension(jpdim) :: ip_messy
  integer , dimension(jpdim) :: indarr, indarra, ix, ip_dissoc
  integer :: j, i, ij, ik
  logical :: lr,lp1,lp2,lq1,lq2,lra,lp1a,lp2a,lq1a,lq2a,lunk1,lunk2,lunk3,lunk4
  logical :: found 
  character, dimension(jpdim,3) :: cphot*10, cphota*10 
  !-----------------------------------------------
  !

!!$#ifndef USE_MPI
!!$  if (ipart > 1) return
!!$#endif
  
  cphota = ' '

  cphot(1,:) = (/'O2     ','O(3P)  ','O(3P)  '/)   
  cphota(1,:)= (/'O2     ','O3     ','O3     '/) 
  cphot(2,:) = (/'O3     ','O2     ','O(3P)  '/) 
  cphot(3,:) = (/'O3     ','O2     ','O(1D)  '/) 
  cphot(4,:) = (/'H2O2   ','OH     ','OH     '/) 
  cphot(5,:) = (/'Cl2    ','Cl     ','Cl     '/) 
  cphot(6,:) = (/'Cl2O2  ','Cl     ','Cl     '/) 
  cphota(6,:) = (/'Cl2O2  ','ClO    ','ClO    '/) 
  cphot(7,:) = (/'HOCl   ','Cl     ','OH     '/) 
  cphot(8,:) = (/'ClNO2  ','Cl     ','NO2    '/) 
  cphot(9,:) = (/'ClONO2 ','Cl     ','NO3    '/) 
  cphota(9,:) = (/'ClONO2 ','ClO    ','NO2    '/) 
  cphot(10,:) = (/'HONO2  ','OH     ','NO2    '/) 
  cphot(11,:) = (/'NO2    ','O(3P)  ','NO     '/) 
  ! the following line is needed for simple chmistry without O(3P):
  cphota(11,:) = (/'NO2    ','O3     ','NO     '/)   
  cphot(12,:) = (/'N2O5   ','NO2    ','NO3    '/) 
  cphot(13,:) = (/'HO2NO2 ','HO2    ','NO2    '/) 
  cphot(14,:) = (/'NO3    ','O(3P)  ','NO2    '/) 
  cphot(15,:) = (/'NO3    ','NO     ','O2     '/) 
  cphot(16,:) = (/'HCHO   ','CO     ','H2     '/) 
  cphot(17,:) = (/'HCHO   ','HCO    ','H      '/) 
  cphota(17,:) = (/'HCHO   ','HO2    ','CO     '/) 
  cphot(18,:) = (/'MeOOH  ','OH     ','MeO    '/) 
  cphota(18,:) = (/'MeOOH  ','OH     ','HCHO   '/) 
  cphot(19,:) = (/'BrONO2 ','BrO    ','NO2    '/) 
  cphot(20,:) = (/'BrCl   ','Cl     ','Br     '/) 
  cphot(21,:) = (/'OClO   ','O(3P)  ','ClO    '/) 
  cphot(22,:) = (/'H2O    ','OH     ','H      '/) 
  cphot(23,:) = (/'HCl    ','Cl     ','H      '/) 
  cphot(24,:) = (/'NO     ','N      ','O(3P)  '/) 
  cphot(25,:) = (/'N2O    ','N2     ','O(1D)  '/) 
  cphot(26,:) = (/'MeOCl  ','Cl     ','MeO    '/) 
  cphot(27,:) = (/'HOBr   ','Br     ','OH     '/) 
  cphot(28,:) = (/'Br2    ','Br     ','Br     '/) 
  cphot(29,:) = (/'MeO2NO2','MeOO   ','NO2    '/) 
  cphot(30,:) = (/'BrO    ','Br     ','O(3P)  '/) 
  cphot(31,:) = (/'ClONO2 ','ClO    ','NO2    '/) 
  cphot(32,:) = (/'BrONO2 ','Br     ','NO3    '/) 
  cphot(33,:) = (/'HO2NO2 ','OH     ','NO3    '/) 
  cphot(34,:) = (/'F11    ','Cl     ','Cl2    '/) 
  cphot(35,:) = (/'F12    ','Cl     ','Cl     '/) 
  cphot(36,:) = (/'F22    ','Cl     ','?      '/) 
  cphot(37,:) = (/'F113   ','Cl     ','Cl2    '/) 
  cphot(38,:) = (/'MeCl   ','MeOO   ','Cl     '/) 
  cphota(38,:) = (/'MeCl   ','CH3    ','HCl    '/) 
  cphot(39,:) = (/'CCl4   ','Cl2    ','Cl2    '/) 
  cphot(40,:) = (/'CH3Br  ','Br     ','MeOO   '/) 
  cphot(41,:) = (/'H1211  ','Br     ','Cl     '/) 
  cphot(42,:) = (/'H1301  ','Br     ','?      '/) 
  cphot(43,:) = (/'CHBr2Cl','Br2    ','Cl     '/) 
  cphot(44,:) = (/'CH2BrCl','Br     ','Cl     '/) 
  cphot(45,:) = (/'CHBrCl2','Br     ','Cl2    '/) 
  cphot(46,:) = (/'CH2Br2 ','Br2    ','?      '/) 
  cphot(47,:) = (/'H2402  ','Br2    ','?      '/) 
  cphot(48,:) = (/'CHBr3  ','Br2    ','Br     '/) 


!!!!! => moved to messy_clamschem_initialize 
!!$  ! corresponding MESSy photolysis indices from messy_cmn_photol_mem
!!$  !
!!$  ! *** must be in the same order than the above list ***
!!$  ip_messy = (/ip_O2, ip_O3P, ip_O1D, ip_H2O2, ip_Cl2, ip_Cl2O2, &
!!$       ip_HOCl, ip_ClNO2, ip_ClNO3, ip_HNO3, ip_NO2, ip_N2O5, &
!!$       ip_HO2NO2, ip_NO2O, ip_NOO2, ip_CHOH, ip_COH2, ip_CH3OOH, &
!!$       ip_BrONO2, ip_BrCl, ip_OClO, ip_H2O, ip_HCl, ip_NO, &
!!$       ip_N2O, ip_CH3OCl, ip_HOBr, ip_Br2, ip_MEO2NO2, ip_BrO, &
!!$       ip_ClONO2, ip_BrNO3, ip_OHNO3, ip_CFCl3, ip_CF2CL2, ip_CHF2Cl, &
!!$       ip_F113, ip_CH3Cl, ip_CCl4/)
 
  !----------------------------------------------Photolysis-switch ----
  ! Photolysis - Code from Roeth (iphot=2 -- currently disabled)
  ! Photolysis - Code from Lary (iphot=1)
  iphot=1

  ! setup array ip_dissoc that should contain the indices 1..jpdim in the order
  ! of the messy indices
  call  quick_sort_index(ip_messy,ix)
  do i=1,jpdim
     ip_dissoc(ix(i)) = i
  enddo
  !write(*,*)'inphot: ip_messy ',ip_messy
  !write(*,*)'inphot: ip_dissoc ',ip_dissoc

  if(rank==0 .and. ipart==1 .and. .not. davg) write (6, *) 'Lary - Photolysis called'
  if(rank==0 .and. ipart==1 .and. davg) write (6, *) 'Lary - Photolysis called for diurnal average'
  
  if (iphot == 2) then 
     if(rank==0 .and. ipart==1) write (6, *) 'Roeth - Photolysis called' 
     ! cphot-index for roeth-list
     indarr(:) = -1
     indarr(1:30) = (/ 1, 3, 2, 6,18,24,30,33,35,16, &
                       8,13,17,10, 9,75,74,90,71,72, &
                      21, 5,29,7,11,-1,110,69,109,70/)
     ! cphota-index for roeth list
     indarra(:) = -1
     indarra(1:18) = (/-1,-1,-1,-1,-1,24,-1,-1,35,-1, &
                        8,-1,-1,-1,-1,-1,74,90/)
  endif

  do ij = 1, jppj 
     found = .FALSE. 
     ! check if unknown products are present in data file
     lunk1= (index(spj(ij,3),'?') > 0) 
     lunk2= (index(spj(ij,4),'?') > 0)

     do ik = 1, jpdim 
        ! check if unknown products are present in program list
        lunk3= (index(cphot(ik,2),'?') > 0) 
        lunk4= (index(cphot(ik,3),'?') > 0) 

        lr = (trim(spj(ij,1)) == trim(cphot(ik,1)))
        lp1 = (trim(spj(ij,3)) == trim(cphot(ik,2)))
        lp2 = (trim(spj(ij,4)) == trim(cphot(ik,3)))
        lq1 = (trim(spj(ij,4)) == trim(cphot(ik,2)))
        lq2 = (trim(spj(ij,3)) == trim(cphot(ik,3)))
           
        if (lunk1) lp1=.true. 
        if (lunk1) lq2=.true. 
        if (lunk2) lp2=.true. 
        if (lunk2) lq1=.true. 
        if (lunk3) lp1=.true. 
        if (lunk3) lq1=.true. 
        if (lunk4) lp2=.true. 
        if (lunk4) lq2=.true. 

 
        if ( (lr.and.lp1.and.lp2 ) .or. (lr.and.lq1.and.lq2)) then 
           if (iphot/=2) then
!              cindex(ij) = ik
              cindex(ij) = ip_dissoc(ik)
           else 
              cindex(ij)= indarr(ik)
           endif
           if (iodump .and. rank==0) write (6, *) 'match found between ASAD ', ij, &
                'and Photolysis', cindex(ij),lp1,lp2,lq1,lq2 
           found = .TRUE. 
        endif
     enddo
 
     if (.not.found) then 
        do ik = 1, jpdim 
           ! look for alternative matching
           lra = (trim(spj(ij,1)) == trim(cphota(ik,1)))
           lp1a = (trim(spj(ij,3)) == trim(cphota(ik,2)))
           lp2a = (trim(spj(ij,4)) == trim(cphota(ik,3)))
           lq1a = (trim(spj(ij,4)) == trim(cphota(ik,2)))
           lq2a = (trim(spj(ij,3)) == trim(cphota(ik,3)))

           if ( (lra.and.lp1a.and.lp2a ) .or. (lra.and.lq1a.and.lq2a)) then 
              if (iphot/=2) then
                 cindex(ij) = ip_dissoc(ik)
              else   
                 cindex(ij) = indarra(ik)
              endif
              if (iodump .and. rank==0) write (6, *) &
                   'alternative match found between ASAD ', ij, &
                   'and Photolysis',cindex(ij),lp1a,lp2a,lq1a,lq2a 
              found = .TRUE. 
           endif
        enddo
        if (.not. found) then
           write (*,*)'**********   WARNING   **********'
           write (*,*)'ASAD/inphot: No match found for ',spj(ij,1)
           write (*,*)' ASAD ', ij,(spj(ij,i),i=1,5)
           write (*,*)' setting '//trim(spj(ij,1))//' photolysis rate to zero'
           write (*,*)'*********************************'
           cindex(ij) = -1 
        endif
     endif
     if (iodump .and.found .and. rank==0) write (*,*) ij,cindex(ij),(spj(ij,i),i=1,5)

  end do 

  return
end subroutine inphot 




end Module messy_clamschem_inphot
 
 
 
