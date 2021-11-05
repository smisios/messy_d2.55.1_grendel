Module messy_clamschem_hetero_shi

Contains



!********************************************************************************!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2010
! Kenneth S. Carslaw, Jens-Uwe Grooss, Rolf Mueller, Tobias Wegner
! Forschungszentrum Juelich GmbH
! Last Modified By: Nicole  Thomas
! Last Modified On: Fri Jun 17 12:01:42 2016
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version. This program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more
! details. You should have received a copy of the GNU General Public
! License along with this program; if not, see <https://www.gnu.org/licenses/>.
!
! ------------------------------------------------------------------------------
!
! Subroutine of CLaMS module chem
!
! calculates the microphysics of particle formation and heterogeneous 
! reaction rates for single air parcels.
! Originally written in fortran 77 format by Ken Carslaw
! translated into fortran 90 and adapted to CLaMS by Jens-Uwe Grooss 
! Incorporation of the parameterisation of chemistry on liquid aerosols
! of Shi et al. (Kinetic model for reaction of ClONO2 with H2O and HCl 
! and HOCl with HCl} in sulfuric acid solutions, J. Grophys. Res., 106,
! 24259-24274, doi:10.1029/2000JD000181, 2001) by Tobias Wegner.
!
! -----------------------------------------------------------------------------


subroutine hetero_shi(press, t, chno3, chcl, ch2o, chocl, chbr, chobr, &
     cclno3, cbrno3, parthno3, parthcl, parthbr, parth2o, wts, wtn, wtcl, &
     wtbr, wthocl, wthobr, aice, anat, aliq, asat, asatliq, n, khet1&
     &,terminal_velocity,mean_free_path,vliq,snat) 
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  ! correction of typo in hetsolid jug/tw 03/2009
  !
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,        ONLY: prec
  use messy_clamschem_globalhet

  implicit none
       
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) :: press, t, chno3, chcl, ch2o, chocl, chbr, chobr, cclno3, &
       cbrno3, parthno3, parthcl, parthbr, parth2o, wts, wtn, wtcl, wtbr, &
       wthocl, wthobr, aice, anat, aliq, asat, asatliq, n,snat
  real(PREC) :: khet1(numhet+1) 
  real(PREC) :: terminal_velocity(2),mean_free_path
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: pi = 3.14159265358979324 
  real(PREC), parameter :: dum1 = 0.0 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: ndrop, denschange, tflag, sulreduction, vliq, vtot
  
  !-----------------------------------------------
  !
  ! EXPLANATION OF LOGICAL VARIABLES:
  ! SATMELTING MEANS THAT SAT IS IN THE PROCESS OF MELTING
  ! LIQUIDS MEAN THAT LIQUID PARTICLES EXIST IN ISOLATION,
  !     AND NOT AS AN EXTERNAL LIQUID/SOLID MIXTURE.
  ! SAT_MELTALLOWED MEANS THAT SAT MELTING IS CONSIDERED AS A PROCESS
  ! PARAM_NAT_HR IS A SWITCH BETWEEN HET RATES OF HANSON AND RAVI OR ABBATT AND MOLINA
  !     (SEE HETSOLID)
  ! ALLICE/ALLNAT ARE TRUE WHEN NO EXTERNAL SOLID/LIQUID MIXTURES ARE ALLOWED
  !     (THESE ARE SET AT START OF PROGRAM IN HETINIT FILE)
  ! MIXEDPOP IS A RUN-TIME LOGICAL VARIABLE THAT CHECKS WHETHER A SOLID/LIQUID
  ! EXTERNAL MIXTURE EXISTS.
  !
 
  ! =====================================================================
  ! ====== SUBSEQUENT CALLS TO HETERO ===================================
  ! =====================================================================

  ! set reaction rates to zero
  khet1 = 0.0
  ! PREVENT UNREALISTIC VALUES OF QUANTITIES THAT ARE DIVIDED
  chcl = abs(chcl) + 1.0 
  chobr = abs(chobr) + 1.0 
  chocl = abs(chocl) + 1.0 
  chbr = abs(chbr) + 1.0 
  ch2o = max(ch2o,1.0d0)
  chno3 = max(chno3,1.0d0)
  ! CALCULATE SOME QUANTITIES USED IN FOLLOWING SUBROUTINES
  denschange = (told/t)*(press/pressold) 
  ! CHANGE PARTICLE NUMBER DENSITY WITH PRESSURE AND TEMPERATURE
  n = n*denschange 
  denssat = denssat*denschange 
  densnat = densnat*denschange 
  densice = densice*denschange 
  ! CHANGE ICE AND NAT NUMBER DENSITIES THAT CAN FORM
  ciceinit = ciceinit*denschange 
  cnatinit = cnatinit*denschange 
  ! LIMIT NAT AND ICE NUMBER DENSITIES TO MAXIMUM OF N
  cnat = min(cnatinit,n)   !don't let cnat>n. must be done every step 
  cice = min(ciceinit,n)   !don't let cice>n 
  if (allice) cice = n 
  if (allnat) cnat = n 
  !
  tflag = t - told 
  pressold = press 
  told = t 
  !
  satmelting = .false. 
  asatliq=0.0   !this line added by JUG, 07/99
  !
  !
  ! IF PURE LIQUID OR LIQUID/SAT MIXED POPULATION.
  ! UPDATE LIQUID COMPOSITION ETC BEFORE DETERMINING PHASE CHANGES.
  ! THE IMPORTANT PARAMETER CALCULATED IS PARTHNO3, WHICH DETERMINES
  ! WHETHER NAT WILL FORM FROM A LIQUID OR LIQUID/SAT MIXED POPULATION.
  !
  if (kstate==1 .or. kstate==2 .and. mixedpop) then 
     sulreduction = 1.0 - (densnat + denssat + densice)/n 
     ! SULREDUCTION IS FRACTIONAL REDUCTION IN LIQUID NUMBER DENSITY DUE TO
     ! EXISTENCE OF OTHER PARTICLES
     call liquid (press, t, chcl, chbr, chocl, chobr, chno3, ch2o, parthcl, &
          parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          aer_h2so4_default, vliq, satmelting, sulreduction) 

  endif
  !
  ! DETERMINE NEW PARTICLE PHASE
  call phase (press, t, tflag, chno3, ch2o, chcl, aer_h2so4_default, parthno3, &
       parth2o, transform, saturation_criteria, densnat, densice, denssat, cice, cnat, n, &
       satmelting, liquids, sat_meltallowed,snat) 
  !---
  !
  !
  !---
  if (satmelting .and. .not.mixedpop) then 
     ! SATMELTING MEANS THAT SAT IS IN EQM WITH LIQUID LAYER.
     ! CALCULATE SAT/LIQUID LAYER COMOSITION AND VOLUME.
     ! DO THIS ONLY IF PURE SAT EXISTS (NOT WITH MIXED POPULATION).
     ! SAT CAN STILL MELT WITH MIXED POP., BUT ONLY COMPLETE MELTING CONSIDERED.
     !
     call satliquid (press, t, chcl, chno3, ch2o, parthcl, parthno3, &
          parth2o, wts, wtcl, wtn, aer_h2so4_default, vliq, vtot, satmelting) 
     !
     ! CALCULATE HETEROGENEOUS REACTIONS ON LIQUIDS (THE MELTING SAT)
     !
     ndrop = n 
     call hetliquid (khet1, t, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          ch2o, chcl, chbr, chocl, chobr, cclno3, vtot, asatliq, ndrop, &
          parthbr) 
 
     !
     ! NOTE THAT IT IS VTOT IN THIS CALL, SINCE THIS IS THE VOLUME OF THE
     ! SAT/LIQUID MIXTURE FROM WHICH THE SURFACE AREA WILL BE CALCULATED.
     ! SET OTHER PARTICLE AREAS AND HETEROGENEOUS REACTIONS TO ZERO
     asat = 0.0 
     anat = 0.0 
     aice = 0.0 
     khet1(:22) = 0. 
     khet1(31:33) = 0. 

  else if (liquids) then 
     ! CALCULATE PURE LIQUID COMPOSITION (UNLESS WAS ALREADY PURE LIQUID, IN
     ! CASE, LIQUID() WAS ALREADY CALLED EARLIER.
     if (laststate /= 1) then 
        sulreduction = 1.0 

        ! jug: calculate the volume vliq with temperature t and the 
        ! remainig parameters with temperature t as solubilities etc
        ! are interchanged though the common blocks molality and hstar
        call liquid (press, t, chcl, chbr, chocl, chobr, chno3, ch2o, &
             parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, &
             wthocl, wthobr, aer_h2so4_default, vliq, satmelting, sulreduction) 


     endif
     ! CALCULATE HETEROGENEOUS REACTIONS ON LIQUIDS
     !
     ndrop = n 
     call hetliquid (khet1, t, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          ch2o, chcl, chbr, chocl, chobr, cclno3, vliq, aliq, ndrop, &
          parthbr) 
     !
     ! SET SOLID PARTICLE RATE COEFFS ETC. TO ZERO
     anat = 0.0 
     aice = 0.0 
     asat = 0.0 
     asatliq = 0.0 
     khet1(:22) = 0. 
     khet1(31:33) = 0. 
     !---
     !
     !
     !---
  else if (mixedpop) then 
     ! MIXED SOLID/LIQUID POPULATION
     ! FIRST CALCULATE HETEROGENEOUS REACTIONS ON LIQUIDS
     sulreduction = 1.0 - (densnat + denssat + densice)/n 
     !
     call liquid (press, t, chcl, chbr, chocl, chobr, chno3, ch2o, parthcl, &
          parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          aer_h2so4_default, vliq, satmelting, sulreduction) 
     !
     !
 
     ! jug test printout warning if parthno3 is invalid -> see subroutine liquid
     !         IF (parthno3 .gt. 1.0) then
     !            write(*,*)'!!!!! parthno3=',parthno3
     !            parthno3=1.0
     !         endif
     ndrop = n*sulreduction 
     call hetliquid (khet1, t, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          ch2o, chcl, chbr, chocl, chobr, cclno3, vliq, aliq, ndrop, &
          parthbr) 
     !
     ! NOW CALCULATE HETEROGENEOUS REACTIONS ON SOLIDS
     !
     ! jug added the following lines (29.03.2004)
     if (parthno3 >= 1.0) densnat=0.
     if (parth2o >= 1.0) densice=0.
     call hetsolid (khet1, press, chcl, chbr, chocl, chobr, cclno3, &
          cbrno3, chno3, ch2o, aer_h2so4_default, parthno3, parth2o, n, densnat, &
          densice, denssat, anat, aice, asat, t, terminal_velocity, &
          mean_free_path)  

  else 
     ! ONLY SOLIDS EXIST.
     ! CALCULATE HETEROGENEOUS REACTIONS ON SOLIDS
     ! SET PARTHCL TO 1. PARTHNO3 AND PARTH2O ARE SET IN PHASE().
     parthcl = 1.0 
     ! SET WT-FRAC IN LIQUID TO 0. THESE ARE NOT USED, BUT MIGHT BE PLOTTED.
     wts = 0.0 
     wtn = 0.0 
     wtcl = 0.0 
     wtbr = 0.0 
     wthocl = 0.0 
     wthobr = 0.0 
     !
     call hetsolid (khet1, press, chcl, chbr, chocl, chobr, cclno3, &
          cbrno3, chno3, ch2o, aer_h2so4_default, parthno3, parth2o, n, densnat, &
          densice, denssat, anat, aice, asat, t, terminal_velocity, &
          mean_free_path) 
     !
     ! SET LIQUID PARTICLE RATE COEFFS TO ZERO
     khet1(23:30) = 0.0 
     aliq = 0.0 
     asatliq = 0.0 

     !
  endif 
  !
  return  
end subroutine hetero_shi 



! ========================================================================
subroutine phase(press, t, tflag, chno3, ch2o, chcl, aer_h2so4_default, parthno3, &
     parth2o, transform, saturation_criteria, densnat, densice, denssat, cice, cnat, n, &
     satmelting, liquids, sat_meltallowed, snat) 
       
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,         ONLY: prec
  USE messy_clamschem_cirrus_clim,ONLY: densice_cirrus,rhice_freeze_clim
  use messy_clamschem_globalhet,  ONLY: natcore, mixedpop, liqtest, t_lt_tnat, &
                                        kstate, laststate, cice_from_clim

  implicit none

  !
  !
  ! **************************************************************
  ! CONTROLS THE PHASE CHANGES BETWEEN LIQUIDS, SAT, MELTING-SAT,
  ! NAT AND ICE
  ! 1=LIQUIDS WITHOUT ANY SOLIDS PRESENT AT ALL
  ! 2=SAT EITHER IN ISOLATION OR AS SAT/LIQUID EXTERNAL MIXTURE
  !     (I.E., SAT AND LIQUID IN SEPARATE PARTICLES)
  ! 3=NAT EITHER IN ISOLATION OR AS NAT/LIQUID OR NAT/LIQUID/SAT
  !     EXTERNAL MIXTURE
  ! 4=ICE EITHER IN ISOLATION OR AS ICE/LIQUID, ICE/NAT/LIQUID OR
  !     ICE/NAT/SAT/LIQUID EXTERNAL MIXTURE
  !
  ! THUS, KSTATE REFERS TO THE 'HIGHEST' PHASE THAT EXISTS
  !
  ! **************************************************************
  !
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
       
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC)  :: press 
  real(PREC)  :: t 
  real(PREC) , intent(in) :: tflag 
  real(PREC)  :: chno3 
  real(PREC)  :: ch2o 
  real(PREC)  :: chcl 
  real(PREC)  :: aer_h2so4_default 
  real(PREC)  :: parthno3 
  real(PREC)  :: parth2o 
  real(PREC) , intent(inout) :: densnat 
  real(PREC) , intent(inout) :: densice 
  real(PREC) , intent(inout) :: denssat 
  real(PREC) , intent(in) :: cice 
  real(PREC) , intent(in) :: cnat 
  real(PREC) , intent(in) :: n 
  logical  :: satmelting 
  logical  :: liquids 
  logical , intent(in) :: sat_meltallowed 
  real(PREC) , intent(in) :: transform(4,4) 
  real(PREC) , intent(in) :: saturation_criteria(4,4) 
  real(PREC) :: mod_criteria(4,4) 

!  real(PREC),external ::phmhno3, phmh2o, densice_cirrus,rhice_freeze_clim

  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: ctomb = 7.2427E+18 
  real(PREC), parameter :: ctoa = 7.336E+21 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: k 
  real(PREC) , dimension(4) :: satpresent 
  real(PREC) :: ph2o, ph2oice, svapice, phno3, phno3nat, svapnat, tsat, parthcl, &
       wts, wtcl, wtn, vliq, vtot, sulreduction, chbr, chocl, chobr, parthbr, &
       wtbr, wthocl, wthobr, totsolid, dum, rhice_freeze,snat 
  logical :: sat, nat, ice, oldliquids, oldsat, oldnat, oldice, tincrease, &
       tdecrease 
  !-----------------------------------------------
  satmelting = .false. 
  mixedpop = .false. 
  snat=0.
  !
  !
  !
  ! CHECK TEMPERATURE TENDENCY
  tincrease = (tflag >= 0)
  tdecrease = (tflag <= 0)
  ! STORE EXISTING PARTICLE STATE
  laststate = kstate 
  oldliquids = (laststate == 1)
  oldsat = (laststate == 2)
  oldnat = (laststate == 3)
  oldice = (laststate == 4)
  ! HANSON & MAUERSBERGER SATURATION WITH RESPECT TO ICE
  ph2o = ch2o*t/ctoa          !water partial pressure (atm) 

  ! filter out unrealistically low or zero values
  if (ph2o < 1E-9) then
     write(*,*) ' hetero_shi warning: ph2o ',ph2o,' corrected'
     ph2o = 1E-9
  endif

  ph2oice = phmh2o(t)         !vapour pressure over ice (atm) 
  svapice = ph2o/ph2oice      !saturation ratio wrt ice 
  satpresent(4) = svapice 
  ! HANSON & MAUERSBERGER SATURATION WITH RESPECT TO NAT
  phno3 = chno3*t/ctoa        !HNO3 partial pressure (assuming all HNO3 in gas) 

  ! filter out unrealistically low or zero values
  if (phno3 < 1E-12)then
     write(*,*) ' hetero_shi warning: phno3 ',phno3,' corrected'
     phno3 = 1E-12
  endif

  phno3nat = phmhno3(t,ph2o) 
  t_lt_tnat = ( phno3nat < phno3)
  svapnat = phno3/phno3nat 
  satpresent(3) = svapnat 
  !
  ! DECIDE WHETHER T>TSAT
  ! TSAT FROM TABAZADEH ET AL., GRL, 21, 1619-, 1994.
  tsat = 3236.0/(11.502 - log10(ph2o*760.0)) 
  if (t >= tsat) then          !only liquids can exist 
     denssat = 0.0 
     densnat = 0.0 
     natcore = .false. 
     densice = 0.0 
     kstate = 1 
     liquids = .true. 
     return  
  endif
  ! SET SAT SATURATION. THIS CANNOT BE CALCULATED, SO IS SET
  ! TO A DUMMY SUCH THAT SAT WILL FORM BELOW TSAT IF
  ! TRANSFORM(1,2) = 1 AND SATURATION_CRITERIA (1,2) = 1.
  if (t <= tsat) satpresent(2) = 99.0 
  if (t<=tsat .and. transform(1,2)==1.0 .and. saturation_criteria(1,2)==1.0 .and. &
       kstate==1) kstate = 2 
  !
  if (tdecrease) then 
     mod_criteria=saturation_criteria
     if (cice_from_clim) then
          rhice_freeze = rhice_freeze_clim(T)
          mod_criteria(1:3,4) = rhice_freeze 
     endif

     satpresent(3) = svapnat*parthno3 
     satpresent(4) = svapice*parth2o 
     ! IE. USE ACTUAL AMOUNT OF HNO3 IN THE GAS PHASE, SINCE SOME IS IN LIQUID
     ! WHEN TEMPS ARE INCREASING NEED TO USE TOTAL POSSIBLE HNO3 PARTIAL PRESSURE
     !
     ! TEST FOR PHASE CHANGES UPON COOLING
     do k = kstate + 1, 4 
     !write(*,'(a6,2i2,4f10.3)') 'test+:',kstate,k,transform(kstate,k),&
     !                 satpresent(kstate),transform(kstate,k)*satpresent(kstate),&
     !                 saturation_criteria(kstate,k)
        if (transform(kstate,k)*satpresent(k) < mod_criteria(kstate,k)) cycle  
        kstate = k 
     end do
  else if (tincrease) then 
     ! TEST FOR PHASE CHANGES UPON WARMING
     do k = kstate - 1, 1, -1 
     !write(*,'(a6,2i2,4f10.3)') 'test-:',kstate,k,transform(kstate,k),&
     !                 satpresent(kstate),transform(kstate,k)*satpresent(kstate),&
     !                 saturation_criteria(kstate,k)
        if (transform(kstate,k)*satpresent(kstate) > saturation_criteria(kstate,k)) &
             cycle  
        kstate = k 
     end do
!  jug 02/2012: 
!  previeous case neither tdecrease or tincrease 
!  should never happen due to "<=" and ">=" but left in...
  else 
     kstate = laststate 
  endif

  !-------------------------------------------------------------------
  ! fix to strange "Ice melting upon cooling" (sb, 2.8.02)
  ! found at ph2o=4.7810088049058184E-007, T=192.3991394042969, T decreasing very slowly
  ! tflag=-1.3732910156250000E-004
  if (tdecrease .and. svapice < 1. .and. kstate==4) then
     !if (tflag < -0.1) write(*,*) 'info: Ice melting with decreasing temperature ',t,tflag,svapice
     ! TEST FOR PHASE CHANGES UPON WARMING
     do k = kstate - 1, 1, -1
        if (transform(kstate,k)*satpresent(kstate) > saturation_criteria(kstate,k)) &
             cycle
        kstate = k
     end do
  endif

  !-------------------------------------------------------------------
  ! fix to strange "NAT melting upon cooling" (jug)
  ! found at phno3=2.14755E-10, T=192.853, T decreasing very slowly
  if (tdecrease .and. svapnat < 1. .and. kstate==3) then 
     !if (tflag < -0.001) write(*,*) 'info: NAT melting with decreasing temperature'
     ! TEST FOR PHASE CHANGES UPON WARMING
     do k = kstate - 1, 1, -1 
        if (transform(kstate,k)*satpresent(kstate) > saturation_criteria(kstate,k)) &
             cycle  
        kstate = k 
     end do
  endif 
  !-------------------------------------------------------------------

  ! SET LOGICAL VARIABLES FOR CURRENT PHASE
  liquids = (kstate == 1)
  sat = (kstate == 2)
  nat = (kstate == 3)
  ice = (kstate == 4)

  ! jug 02/2012 save NAT supersaturation
  snat = svapnat * parthno3
  !
  ! DETERMINE THE PARTICLE CHARACTERISTICS
  !
  if (liquids) then 
     densnat = 0.0 
     natcore = .false. 
     densice = 0.0 
     denssat = 0.0 
     return  
  endif
  !
  if (sat) then 
     if (oldice) denssat = denssat + densice 
     if (oldnat) denssat = denssat + densnat 
     if (oldliquids) denssat = n 
     densnat = 0.0 
     natcore = .false. 
     densice = 0.0 
     parthno3 = 1.0 
     parth2o = 1.0 
     ! SEE IF SAT IS MELTING
     if (sat_meltallowed) then 
        ! CHECKMELT RETURNS LOGICAL SATMELTING = YES OR NO,
        ! AND LOGICAL LIQUIDS = YES OR NO
        call checkmelt (press, t, ch2o, chno3, chcl, aer_h2so4_default, liquids, &
             satmelting) 
        ! SAT MELTING -> NAT IF NAT ALLOWED TO FORM FROM LIQUID
        ! FIRST, GET PARTHNO3 FOR LIQUID OR LIQUID/SAT WHICH
        ! IS USED TO CALCULATE THE NEW SVAPNAT
        if (satmelting .and. .not.mixedpop) then 
           call satliquid (press, t, chcl, chno3, ch2o, parthcl, parthno3, &
                parth2o, wts, wtcl, wtn, aer_h2so4_default, vliq, vtot, satmelting) 
        else if (liquids) then 
           sulreduction = 1.0 
           ! SULREDUCTION IS FRACTIONAL REDUCTION IN LIQUID NUMBER DENSITY DUE
           ! TO EXISTENCE OF OTHER PARTICLES

           call liquid (press, t, chcl, chbr, chocl, chobr, chno3, ch2o, &
                parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, &
                wthocl, wthobr, aer_h2so4_default, vliq, satmelting, sulreduction) 
        endif
        satpresent(3) = svapnat*parthno3 
        if ((liquids .or. satmelting) .and. transform(1,3)*satpresent(3)>=&
             saturation_criteria(1,3)) then 
           kstate = 3 
           nat = .true. 
           satmelting = .false. 
           liquids = .false. 
        else if (liquids) then 
           kstate = 1 
           satmelting = .false. 
           return  
        else 
           ! MUST STILL BE PURE SAT
           go to 99 
        endif
     endif
  endif
  !
  if (nat) then 
     parthno3 = 1.0/svapnat 
     parth2o = 1.0 
     natcore = .true. 
     if (oldsat) then 
        densnat = min(denssat,cnat) 
        densice = 0.0 
        denssat = max(denssat - cnat,0.D0) 
        go to 99 
     else if (oldliquids) then 
        densnat = cnat 
        densice = 0.0 
        denssat = 0.0 
        go to 99 
     else if (oldice) then 
        densnat = densnat + densice 
        densice = 0.0 
        ! denssat unchanged
        go to 99 
     endif
  endif
  !
  if (ice) then 
     ! NEXT 2 LINES ASSUME THAT NAT ALWAYS FORMS WHEN THERE IS ICE
     natcore = .true. 
     parthno3 = 1.0/svapnat 
     parth2o = 1.0/svapice 
     if (oldliquids) then 
        densice = cice 
        densnat = 0.0 
        denssat = 0.0 
        natcore = .true. 
        go to 99 
     else if (oldnat) then 
        natcore = .true. 
        densice = min(cice,densnat) 
        densnat = max(densnat - densice,0.D0) 
        parthno3 = 1.0/svapnat 
        go to 99 
        ! ASSUME DENSSAT UNCHANGED (=> ICE NUCLEATES ON LARGE NAT PARTICLES)
     else if (oldsat) then 
        densice = min(cice,denssat) 
        denssat = max(denssat - cice,0.D0) 
        go to 99 
     endif
  endif
  ! ELIMINATE VERY SMALL NUMBER DENSITIES CAUSED BY NUMERICAL ROUNDING
99 continue 
  if (cice_from_clim) then
     ! use formula from Martina Kraemers climatology
     densice =  densice_cirrus(ch2o,parth2o)
  endif
  !jug, 03/2002 minimum density set to 1E-6 instead of 1E-4
  if (denssat <= 1E-6) denssat = 0.0 
  if (densnat <= 1E-6) densnat = 0.0 

  ! CHECK WHETHER SOME LIQUID PARTICLES REMAIN. METHOD ACCOUNTS FOR
  ! SMALL PROPAGATION OF ROUNDING ERRORS IN DENSNAT...ETC
  totsolid = denssat + densnat + densice 
  dum = abs(totsolid/n - 1.0) 
  ! IF NUMBER OF LIQUID PARTICLES < 1E-4 * TOTAL NUMBER OF SOLID PARTICLES
  ! THEN ASSUME ALL PARTICLES ARE SOLID
  if (dum > 1.E-6) mixedpop = .true. 
  !
  return  
end subroutine phase 


!
!
!
!
!
!
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
!
subroutine checkmelt(press, t, ch2o, chno3, chcl, aer_h2so4_default, liquids, satmelting)
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: press 
  real(PREC) , intent(in) :: t 
  real(PREC)  :: ch2o 
  real(PREC)  :: chno3 
  real(PREC)  :: chcl 
  real(PREC)  :: aer_h2so4_default 
  logical , intent(out) :: liquids 
  logical  :: satmelting 
  !-----------------------------------------------
  !   C o m m o n   B l o c k s
  !-----------------------------------------------
! common /molality/ ms, mn 
! real(PREC)   ms, mn 
  common /molality/ ms, mn, mst, wtst, msb, wtsb
  real(PREC)   ms, mn, mst, wtst, msb, wtsb
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: ctoa = 7.336E+21 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i 
  real(PREC) , dimension(10) :: q, p 
  real(PREC) :: nsulsatmelt, phno3, ph2o, pwfit, pnfit, tmeltapprox, temp, chbr, &
       chocl, chobr, parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, &
       wtbr, wthocl, wthobr, vliq, tmelt, tupper, tlower 
  !-----------------------------------------------
  !
  !
  !
  ! ********************************************************************************
  ! WHEN SAT EXISTS, THIS ROUTINE CHECKS TO SEE WHETHER THE SAT MELTING TEMPERATURE
  ! HAS BEEN REACHED (TMELTAPPROX).  IF NOT, RETURNS TO PHASE.  IF YES, THEN CALCULATES
  ! THE UPPER SAT MELTING TEMP (TUPPER).  IF T STILL > TUPPER, THEN RETURNS TO PHASE.
  ! IF T<TUPPER, THEN CALCULATES LOWER MELTING TEMP (TLOWER). IF TLOWER<T<TUPPER, THEN
  ! SETS SATMELTING LOGICAL VARIABLE TO .TRUE. AND RETURNS TO PHASE.  IF T<TLOWER
  ! THEN SAT HAS COMPLETELY MELTED, SO SET LIQUIDS LOGICAL VARIABLE TO .TRUE.
  ! AND RETURN TO PHASE.
  ! ********************************************************************************
  !
  !
  !
  DATA P/ 40.195593, 10.086927, -9.1598623E-003, 0.37268462, -0.61415295, &
       -3.11643687E-002, 47.436130, -3.571783E-004, 5.07800820E-006, &
       -8.2090956E-005/  
  DATA Q/ 182.213, 3.2854, -0.2332034, 0.376993, -3.52459E-002, &
       1.766402E-003, 1.257800, -0.43670, -3.388164E-005, 4.18979691E-002/  
  ! ========================================================
  phno3 = chno3*t/ctoa 
  ph2o = ch2o*t/ctoa 
  pwfit = ph2o*1013.0*1e4 
  pnfit = phno3*1013.0*1e7 
  tmeltapprox = q(1) + q(2)*pwfit + q(3)*pwfit**2 + (q(4)*pnfit+q(5)*pnfit**2+&
       q(6)*pnfit**3+q(9)*pnfit**4)*(1.0 + q(7)*pwfit+q(8)*pwfit**2+q(10)&
       *pwfit**3) 
  if (t - tmeltapprox >= 1.0) then 
     satmelting = .false.           !WELL ABOVE MELTING TEMPERATURE 
     return  
  else 
     satmelting = .true.            !CLOSE TO MELTING, SO FIND REAL UPPER MELTING TEMPERATURE 
     temp = tmeltapprox + 1.0 
     do i = 1, 200 
        temp = temp - 0.1 
        nsulsatmelt = 1.E-15 
        call liquid (press, temp, chcl, chbr, chocl, chobr, chno3, ch2o, &
             parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, &
             wthocl, wthobr, aer_h2so4_default, vliq, satmelting, nsulsatmelt) 
 
        tmelt = p(1) + p(2)*mn + p(3)*mn**3 + p(5)*mn**0.5 + (p(4)*ms**3 + &
             p(6)*ms**4+p(7)*ms**0.5)*(1.0 + p(8)*mn**2+p(9)*mn**4+p(10)*mn**3)
        if (tmelt > temp) cycle  
        exit                        !FOUND MELTING TEMPERATURE 
     end do
     tupper = temp 
  endif
  if (tupper <= t) then 
     satmelting = .false.           ! REAL MELTING TEMP STILL LOWER, SO RETURN 
     return  
  else                              !SAT REALLY IS MELTING! 
     !...SO FIND LOWER MELTING TEMPERATURE ALSO
     temp = tupper 
     do i = 1, 100 
        temp = temp - 0.1 
        nsulsatmelt = press*100.0*aer_h2so4_default*1.E-9/8.314/temp 
        call liquid (press, temp, chcl, chbr, chocl, chobr, chno3, ch2o, &
             parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, &
             wthocl, wthobr, aer_h2so4_default, vliq, satmelting, nsulsatmelt) 
        
        tmelt = p(1) + p(2)*mn + p(3)*mn**3 + p(5)*mn**0.5 + (p(4)*ms**3+ &
             p(6)*ms**4+p(7)*ms**0.5)*(1.0 + p(8)*mn**2+p(9)*mn**4+p(10)*mn**3)
        if (tmelt > temp) cycle  
        exit  
     end do
     tlower = temp 
     if (t <= tlower) then 
        satmelting = .false. ! WE ARE BELOW MELTING TEMP, FULLY MELTED => LIQUID ONLY 
        liquids = .true. 
        return  
     else 
        satmelting = .true. 
     endif
  endif
  return  
end subroutine checkmelt 


!
!
!
!
!
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
!
subroutine satliquid(press, t, chcl, chno3, ch2o, parthcl, parthno3, &
     parth2o, wts, wtcl, wtn, aer_h2so4_default, vliq, vtot, satmelting) 
  !
  !
  ! **********************************************************************
  ! CALCULATES THE COMPOSITION OF LIQUID HNO3/H2SO4/H2O LAYER ON A MELTING
  ! SAT PARTICLE
  ! **********************************************************************
  !
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC), intent(in)  :: press 
  real(PREC), intent(in)  :: t 
  real(PREC), intent(in)   :: chcl 
  real(PREC), intent(in)   :: chno3 
  real(PREC), intent(in)   :: ch2o 
  real(PREC), intent(out)  :: parthcl 
  real(PREC), intent(inout)  :: parthno3 
  real(PREC), intent(inout)  :: parth2o 
  real(PREC)  :: wtn 
  real(PREC)  :: wts 
  real(PREC), intent(out)  :: wtcl 
  real(PREC), intent(in)  :: aer_h2so4_default 
  real(PREC), intent(out)  :: vliq 
  real(PREC), intent(out) :: vtot 
  logical, intent(in)  :: satmelting 
  !-----------------------------------------------
  !   C o m m o n   B l o c k s
  !-----------------------------------------------
! common /molality/ ms, mn 
! real(PREC)   ms, mn 
  common /molality/ ms, mn, mst, wtst, msb, wtsb
  real(PREC)   ms, mn, mst, wtst, msb, wtsb
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: pi = 3.14159265358979324 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i 
  real(PREC) :: ns, nsliq, nssat 
  real(PREC) , dimension(10) :: p 
  real(PREC) :: fraction, chbr, chocl, chobr, parthbr, wtbr, wthocl, wthobr,&
       tsat, vsat 
  !-----------------------------------------------
  ! =========================================================
  DATA P/ 40.195593, 10.086927, -9.1598623E-003, 0.37268462, -0.61415295, &
       -3.11643687E-002, 47.436130, -3.571783E-004, 5.07800820E-006, &
       -8.2090956E-005/  
  ! =========================================================
  !
  fraction = 1.0 
  do i = 1, 19 
     fraction = fraction - 0.05 
     ns = press*100.0*aer_h2so4_default*1.E-9/8.314/t 
     nsliq = ns*fraction 
     nssat = (1.0 - fraction)*ns 
     call liquid (press, t, chcl, chbr, chocl, chobr, chno3, ch2o, parthcl, &
          parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, wthocl, wthobr, &
          aer_h2so4_default, vliq, satmelting, nsliq) 
     !
     tsat = p(1) + p(2)*mn + p(3)*mn**3 + p(5)*mn**0.5 + (p(4)*ms**3+p(6)*&
          ms**4+p(7)*ms**0.5)*(1.0 + p(8)*mn**2+p(9)*mn**4+p(10)*mn**3) 
     if (tsat > t) cycle  
     exit  
  end do
  vsat = nssat*170.0/1.6*1.E+6 
  vtot = vliq + vsat 
  !
  return  
end subroutine satliquid 


!
!
!
!
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
!
subroutine liquid(press, temp, chcl, chbr, chocl, chobr, chno3, ch2o, &
     parthcl, parthbr, parthno3, parth2o, wts, wtn, wtcl, wtbr, wthocl, &
     wthobr, aer_h2so4_default, vliq, satmelting, nsulreduced) 

  !
  ! *******************************************************************
  ! CALCULATES COMPOSITION OF LIQUID HNO3/H2SO4/H2O AEROSOL
  ! USING THE ANALYTIC PROGRAM OF CARSLAW ET AL., GRL, 22, 1877-, 1995,
  ! AND SOLUBILITIES OF HCL, HBR, HOCL, HOBR.
  !
  ! ALSO CALCULATES COMPOSITION OF LIQUID AEROSOLS IN EQUILIBRIUM WITH
  ! NAT OR ICE OR SAT
  ! *******************************************************************
  !
  !
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  use messy_clamschem_globalhet, ONLY: natcore, mixedpop, liqtest, t_lt_tnat, &
                                       kstate, laststate
  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: press 
  real(PREC) , intent(in) :: temp 
  real(PREC) , intent(in) :: chcl 
  real(PREC) , intent(in) :: chbr 
  real(PREC) , intent(in) :: chocl 
  real(PREC) , intent(in) :: chobr 
  real(PREC) , intent(in) :: chno3 
  real(PREC) , intent(in) :: ch2o
  real(PREC) , intent(out) :: parthcl 
  real(PREC) , intent(out) :: parthbr 
  real(PREC) , intent(inout) :: parthno3 
  real(PREC) , intent(inout) :: parth2o 
  real(PREC)  :: wtn
  real(PREC)  :: wts
  real(PREC) , intent(out) :: wtcl 
  real(PREC) , intent(out) :: wtbr 
  real(PREC) , intent(out) :: wthocl 
  real(PREC) , intent(out) :: wthobr 
  real(PREC) , intent(in) :: aer_h2so4_default 
  real(PREC) , intent(out) :: vliq 
  real(PREC) , intent(in) :: nsulreduced 
  logical , intent(in) :: satmelting 
  !-----------------------------------------------
  !   C o m m o n   B l o c k s
  !-----------------------------------------------
  common /molality/ ms, mn, mst, wtst, msb, wtsb
  real(PREC)   ms, mn, mst, wtst, msb, wtsb
  common /hstar/ hhcl, hhbr, hhocl, hhobr, hclono2, mhcl, phcl 
  real(PREC)   hhcl, hhbr, hhocl, hhobr, hclono2, mhcl, phcl
!  real(PREC),external :: density, density_ter
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: ctomb = 7.2427E+18 
  real(PREC), parameter :: ctoa = 7.336E+21 
  real(PREC), parameter :: r = 8.314E-5 ! gas constant in bar m^3  K^-1 mol^-1
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) , dimension(10) :: hnnarr, hnsarr 
  real(PREC) , dimension(7) :: mnbarr, msbarr 
  real(PREC) :: mnb, lnp, nsul,  mbr, mcl,  mhocl, mhobr, ph2o, tice, t, &
       phcl0, pn0, phbr0, phobr0, phocl0, a, b, c, xs, xn, wtnb, x, tt, &
       hns, hnn, c0, d, alpha, beta, xx, phi, pi, dd, vpn, wn, ws, phbr, &
       phocl, phobr, xmf, shocl, mh2so4, sclono2

  save hnnarr, hnsarr, mnbarr, msbarr, /hstar/, /molality/ 
  ! =====================================================================
  ! DATA ARRAYS FOR HNO3-H2SO4-H2O COMPOSITION CALCULATION
  DATA hnnarr/ 0.1457341419E+02, 0.6159942877E-01, -0.1148954003E+01, &
       0.6916933619E+00, -0.9886352167E-01, 0.5157939752E-02, &
       0.1234728041E+00, -0.1155743833E+00, 0.1101132109E-01, &
       0.9791396551E-02/  
 
  DATA hnsarr/ 0.1446998319E+02, 0.6387958319E-01, -0.3295968441E+01, &
       0.1778224331E+01, -0.2232437780E+00, 0.8648607708E-02, &
       0.5366954950E+00, -0.3351643646E+00, 0.2651532224E-01, &
       0.1575503129E-01/  
 
  DATA mnbarr/ 0.1985306995E+03, -0.1194830927E+05, -0.2846953556E+02, &
       -0.3913646516E+02, 0.8328785711E+02, 0.6358493057E+04, &
       -0.1764983192E+05/  
 
  DATA msbarr/ 0.4700412633E+02, -0.6969007247E+04, -0.4618369994E+01, &
       -0.2166107102E+02, 0.5181581907E+02, 0.2724212325E+04, &
       -0.1573208296E+05/  
  ! ======================================================================
  !
  ! PREVENT PH2O OUTSIDE LIMITS OF MODEL
  !
  ph2o = ch2o/ctoa*temp 
  ph2o = max(2.0D-8,ph2o) 
  ph2o = min(2.0D-5,ph2o) 
  ! ======================================================================
  !
  ! PREVENT T<185 K AND T<TICE-3 K, AS EXPECTED BY MODEL
  !
  tice = 2668.7/(10.431 - (log(ph2o) + log(760.0))/log(10.0)) 
  if (temp <= tice - 3) then 
     t = tice - 3.0 
  else if (temp > 240.0) then 
     t = 240.0 
  else 
     t = temp 
  endif
  t = max(185.D0,t) 
  !
  ! PARTIAL PRESSURE IN ATM FOR H2O & HNO3
  nsul = aer_h2so4_default*press*1.E-9*1.E+2/8.314/t  !moles/m3 sulphate when pure liq 
  ! NOTE: NSULREDUCED SERVES 2 PURPOSES. WHEN LIQUID() IS CALLED FROM SATMELTING
  ! ROUTINES, NSULSATMELT IS THE ACTUAL REDUCED NSUL. WITH A MIXED POPULATION
  ! (MIXEDPOP=.TRUE.) NSULREDUCED IS A REDUCTION 'FACTOR' CALCULATED IN THE MAIN PROG.
  phcl0 = chcl/ctoa*t 
  pn0 = chno3/ctoa*t 
  phbr0 = chbr/ctoa*t 
  phobr0 = chobr/ctoa*t 
  phocl0 = chocl/ctoa*t 
  ! CHECK TO SEE IF LIQUID IS PURE, OR IN EQM WITH ICE, NAT OR SAT
  if (.not.liqtest) then 
     if (satmelting) nsul = nsulreduced 
     if (mixedpop) nsul = aer_h2so4_default*press*1.E-9*1.e2/8.314/t*nsulreduced 
     if (natcore) pn0 = pn0*parthno3         !use nat partial pressure when nat 
     if (kstate == 4) ph2o = ph2o*parth2o    !use ice partial pressure when ice exists 
  endif
  lnp = log(ph2o) 
  ! CALCULATE MSB, MOLALITY OF H2SO4 IN BINARY
  a = msbarr(5) + msbarr(7)/t 
  b = msbarr(4) + msbarr(6)/t 
  c = msbarr(1) + msbarr(2)/t + msbarr(3)*log(t) - lnp 
  xs = ((-b) - sqrt(abs(b**2 - 4.0*a*c)))/(2.0*a) 
  msb = 55.51*xs/(1.0 - xs) 
  wtsb = msb*98.0/(msb*98.0 + 1000.0) 
  if (t <= 210.0) then 
     ! DON'T CALCULATE HNO3 UPTAKE FOR T>210K (TOO SMALL)
     ! CALCULATE MNB, MOLALITY OF HNO3 IN BINARY
     a = mnbarr(5) + mnbarr(7)/t 
     b = mnbarr(4) + mnbarr(6)/t 
     c = mnbarr(1) + mnbarr(2)/t + mnbarr(3)*log(t) - lnp 
     xn = ((-b) - sqrt(abs(b**2 - 4.0*a*c)))/(2.0*a) 
     mnb = 55.51*xn/(1.0 - xn) 
     wtnb = mnb*63.0/(mnb*63.0 + 1000.0) 
     ! CALCULATE HNS, HENRY'S LAW CONSTANT OF HNO3 IN H2SO4-H2O
     x = lnp + 18.4 
     tt = 1.E+4*(1.0/t - 1.0/230.0) 
     hns = exp(hnsarr(1)+hnsarr(2)*tt**2+(hnsarr(3)+hnsarr(4)*tt+hnsarr(5)*&
          tt**2+hnsarr(6)*tt**3)*x+(hnsarr(7)+hnsarr(8)*tt+hnsarr(9)*tt**2)*&
          x**2+(hnsarr(10)*tt)*x**3) 
     ! CALCULATE HNN, HENRY'S LAW CONSTANT OF HNO3 IN HNO3-H2O
     hnn = exp(hnnarr(1)+hnnarr(2)*tt**2+(hnnarr(3)+hnnarr(4)*tt+hnnarr(5)*&
          tt**2+hnnarr(6)*tt**3)*x+(hnnarr(7)+hnnarr(8)*tt+hnnarr(9)*tt**2)*&
          x**2+(hnnarr(10)*tt)*x**3) 
     ! NOW ANALYTICAL EXPRESSION FOR COMPOSITION OF HNO3-H2SO4-H2O AEROSOL
     c0 = r*t*nsul 
     d = (c0*hnn*mnb*msb**2)/(mnb - msb) 
     c = msb*((-2.0*c0*hnn*mnb) + c0*hns*msb + mnb*msb - hnn*msb*pn0)/(mnb-msb)
     b = (c0*hnn*mnb**2 - c0*hns*mnb*msb - 2.0*mnb**2*msb + mnb*msb**2 + &
          hnn*mnb*msb*pn0 - hns*msb**2*pn0)/(mnb**2 - mnb*msb) 
     alpha = (-2.0*b**3) + 9.0*b*c - 27.0*d 
     beta = sqrt(4.0*(b**2 - 3.0*c)**3 - alpha**2) 
     xx = beta/alpha 
!!!!! datan: verursacht bei Compilation mit ifc (ohne -r8) einen Compilerfehler
!     phi = datan(xx) 
     phi = atan(xx) 
     pi = 2.0*asin(1.0) 
     if (phi < 0.) phi = phi + pi 
     dd = sqrt(1.0 - 3.0*c/b**2) 
     ms = -1.0/3.0*(b + 2.0*sqrt(b**2 - 3.0*c)*cos((1.0*pi + phi)/3.0)) 
     mn = mnb*(1.0 - ms/msb) 
     !let mn not be negative! (jug, 26.06.2002)
     mn=max(mn,0.d0)
     ! ======================
     ! NOW THE OUTPUT VARIABLES
     wts = ms*98.12/(1000.0 + 98.12*ms + mn*63.0) 
     wtn = mn*63.0/(1000.0 + 98.12*ms + mn*63.0) 
     vpn = mn/(hnn*mn/(mn + ms) + hns*ms/(mn + ms)) 
     if (.not.natcore) parthno3 = 1.0 - (pn0 - vpn)/pn0 
     ! DON'T CALCULATE PARTHNO3 IF NAT EXISTS, SINCE IT IS THEN
     ! FIXED BY THE NAT VAPOUR PRESSURE. PARTHNO3 IS THEN CALCULATED IN PHASE().
     
     !jug  line added to stay in valid range of parthno3 (J.U. Grooss, 26.10.98)
     parthno3 = min(parthno3,1.D0) 
     
     vliq = nsul*98.12/(wts*density_ter(wts,wtn,t))*1.E-6 
     ! VLIQ IN UM3/CM3    !jug...should be cm^3/cm^3
     ! ======================

mst = ms
wtst = wts

msb = msb
wtsb = wtsb

  else 
     ms = msb 
     mst = msb
     mn = 0.0 
     wts = wtsb
     wtsb = wtsb
     wtst = wtsb
     wtn = 0.0 

     parthno3 = 1.0 
     vliq = nsul*98.12/(wts*density_ter(wts,wtn,t))*1.E-6 
  endif
  ! EQUATIONS ASSUME NO H2O REMOVAL, SO SET PARTH2O=1
  if (.not.mixedpop) parth2o = 1.0 
  !
  !
  !
  ! NOW THE SOLUBILITIES OF HCL,HBR,HOCL AND HOBR
  ! =====================================================================
  ! THE SOLUBILITY OF HCL AND HBR
  ! HHCl (MOL/KG/ATM) taken form Shi et al., JGR 2001
  ! HHBr original code taken from Luo et al. JGR 1995
  ! CALCULATED CONCENTRATIONS ASSUME THAT HCL AND HBR ARE TRACE
  ! COMPONENTS OF THE AEROSOL.
  wn = wtn 
  ws = wtsb
  xmf = (ws*100.)/((ws*100.)+(100.-(ws*100.))*98./18.)
  hhcl = (0.094-0.61*xmf+1.2*xmf**2.)*exp(-8.68+(8515.-10718.*xmf**(0.7))/t) 
  mcl = (1.0/r/t*phcl0)/(nsul/mst + 1.0/r/t/hhcl) 
  wtcl = mcl*36.5/(1000.0 + 98.12*mst + 63.03*mn) 
  phcl = mcl/hhcl 
  mhcl = hhcl*phcl
  parthcl = 1.0 - (phcl0 - phcl)/phcl0 
  !
  hhbr = exp((-(17.83 + 1.02*wn - 1.08*wtst + 3.9*sqrt(wn) + 4.38*sqrt(wtst) - &
       8.87*wn**2 - 17.0*wn*wtst + 3.73*wtst**2)) - 1.0/t*(-8220.50 - 362.76*wn+&
       658.93*wtst - 914.0*sqrt(wn) - 955.3*sqrt(wtst) + 9976.6*wn**2 + &
       19778.5*wn*wtst + 7680.0*wtst**2) - log(wn + 0.410*wtst) - &
       log(80.918/(1000.0 + 98.076*mst + 63.012*mn)))*1.013e3 
  mbr = (1.0/r/t*phbr0)/(nsul/mst + 1.0/r/t/hhbr) 
  wtbr = mbr*80.918/(1000.0 + 98.12*mst + 63.03*mn) 
  phbr = mbr/hhbr 
  parthbr = 1.0 - (phbr0 - phbr)/phbr0 
  ! =====================================================================
  ! THE SOLUBILITY OF HOCL AND HOBR
  ! H* HOCL (MOL/KG/ATM) FROM SHI ET AL JGR 2001. Assuming H depends on molarity
  ! of the HNO3/H2SO4 solution. Former expression depended on molality
  shocl = 0.0776+59.18/t
  mh2so4 = density(msb,t)*ws*100./9.8
  hhocl = 1.91E-6*exp(5862.4/t)*exp(-shocl*mh2so4)
  mhocl = (1.0/r/t*phocl0)/(nsul/mst + 1.0/r/t/hhocl) 
  wthocl = mhocl*52.46/(1000.0 + 98.076*mst + 63.012*mn) 
  !
  ! SOLUBILITY OF HOCL IS LOW ENOUGH TO IGNORE GAS PHASE REMOVAL
  !
  ! SOLUBILITY OF HOBR IS NOT KNOWN FOR ALL STRATOSPHERIC CONDITIONS.
  ! LIMITED DATA (HANSON AND RAVISHANKARA 1995 AT 210K, 60 WT% H2SO4
  ! INDICATE THAT HHOBR=APPROX 18*HHOCL.  FOR HOBR AN EFFECTIVE HENRY'S
  ! LAW CONSTANT =18*HHOCL IS USED.
  hhobr = 18.0*hhocl 
  mhobr = (1.0/r/t*phobr0)/(nsul/mst + 1.0/r/t/hhobr) 
  wthobr = mhobr*96.91/(1000.0 + 98.076*mst + 63.012*mn) 
  phobr = phobr0 
  !
  !HClONO2 FROM SHI ET AL JGR 2001
  !
  sclono2 = 0.306+24/t
  hclono2 = 1.6E-6*exp(4710/t)*exp(-sclono2*mh2so4)
  !
  return  
end subroutine liquid 


!
!
!
!
!
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
!
subroutine hetliquid(khet1, t, wts, wtn, wtcl, wtbr, wthocl, &
     wthobr, ch2o, chcl, chbr, chocl, chobr, cclno3, vliq, aliq, ndrop, &
     parthbr) 
  !
  !
  !
  ! ***********************************************************************
  ! ROUTINE TO CALCULATE UPTAKE COEFFICIENTS (GAMMA VALUES).
  ! GAMMA VALUES ARE INDICATED BY VARIABLES WITH PREFIX 'G', FOR
  ! EXAMPLE GHOCLHCL IS THE GAMMA VALUE OF HOCL DUE TO REACTION WITH
  ! HCL IN THE DROPLETS.
  ! FROM THE GAMMA VALUES, SECOND ORDER RATE CONSTANTS ARE CALCULATED.
  ! THESE ARE CALLED KHET1(...) AND HAVE UNITS CM3 MOLECULE-1 S-1. FOR
  ! EXAMPLE, THE LOSS OF CLNO3 AND HCL DUE TO THE HETEROGENEOUS REACTION
  ! CLNO3+HCL -> CL2+HNO3 IS D(CLNO3)/DT (UNITS MOLECULE CM-3 S-1) =
  ! -RCLNO3HCL.[CLNO3].[HCL], WHERE [CLNO3] AND [HCL] ARE THE
  ! ****TOTAL**** AMOUNTS OF THESE SPECIES IN UNITS MOLECULE CM-3.
  ! ***********************************************************************
  !
  ! 06/07/1999 : gamma and liq_sdist_sigma are now defined in messy_clamchem_hetpar
  !              instead of being
  !              an argument of this subroutine (J.-U. Grooss)

  !
  !
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  use messy_clamschem_globalhet, ONLY: numhet, gamma, liq_sdist_sigma
  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC), intent(in)  :: t 
  real(PREC), intent(in)  :: wtn 
  real(PREC), intent(in)  :: wts
  real(PREC) , intent(in) :: wtcl 
  real(PREC) , intent(in) :: wtbr 
  real(PREC) , intent(in) :: wthocl 
  real(PREC) , intent(in) :: wthobr 
  real(PREC) , intent(in) :: ch2o 
  real(PREC) , intent(in) :: chcl 
  real(PREC) , intent(in) :: chbr 
  real(PREC) , intent(in) :: chocl 
  real(PREC) , intent(in) :: chobr 
  real(PREC) , intent(in) :: cclno3 
  real(PREC) , intent(in) :: vliq 
  real(PREC) , intent(out) :: aliq 
  real(PREC) , intent(in) :: ndrop 
!  real(PREC) , intent(in) :: liq_sdist_sigma 
  real(PREC) , intent(in) :: parthbr 
  real(PREC) , intent(out) :: khet1(numhet+1) 
!  real(PREC) , intent(in) :: gamma(numhet) 
  !-----------------------------------------------
  !   C o m m o n   B l o c k s
  !----------------------------------------------
  common /molality/ ms, mn, mst, wtst, msb, wtsb
  real(PREC)   ms, mn, mst, wtst, msb, wtsb
  common /hstar/ hhcl, hhbr, hhocl, hhobr, hclono2, mhcl, phcl 
  real(PREC)  hhcl, hhbr, hhocl, hhobr, hclono2, mhcl, phcl 

!  real(PREC),external :: density, dhocl, dclono2, density_ter
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: ksur = 576.0 
  real(PREC), parameter :: k2hoclhcl = 1.e5 
  real(PREC), parameter :: k2hobrhcl = 1.e5 
  real(PREC), parameter :: k2hbrhobr = 1.e7 
  real(PREC), parameter :: k2hbrhocl = 1.e6 
  real(PREC), parameter :: ctoa = 7.336E+21 
  real(PREC), parameter :: alpha = 1.0 
  real(PREC), parameter :: rho = 2.E+3 
  real(PREC), parameter :: pi = 3.14159265358979324 
  real(PREC), parameter :: e = 1.0 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: k1hoclhcl, k1hobrhcl, k1hbrhobr, k1hbrhocl, ph2o, rmode, rmean, &
       dl, dens, chclliq, choclliq, chbrliq, chobrliq, factor, q, fq, cbar, &
       ghoclhcl, ah2o, g0, gs, p, gcalc, adivl, f, ge, gclno3hcl, gclno3h2o, &
       gn2o5h2o, ghobrhcl, ghbrhobr, ghbrhocl, grxn, gbrno3h2o, ah, fhocl,&
       & lhocl, grxnhocl, fhcl, xmf, cbar2, kh, kh2o, khydr, gbh2o, khcl,&
       & lclono2, fclono2, grxnclono2, gbhcl, gs_strich, gbhcl_strich, gb, pclono2 
  !-----------------------------------------------
  !
  khet1(23:30) = 0. 
  !
  ph2o = ch2o/ctoa*t                                  ! ATM 
  pclono2 = cclno3/ctoa*t
  ! NOW THE LIQUID PARAMETERS NEEDED FOR CALCULATING REACTION RATES
  rmode = (3.0*vliq/(4.0*pi*ndrop)*exp((-9.0/2.0*log(liq_sdist_sigma)**2)))**(1.0/3.0)
  !CM 
  rmean = rmode*exp(0.5*log(liq_sdist_sigma)**2)                ! CM  
  aliq = 4.0*pi*rmode**2*ndrop*exp(2.0*log(liq_sdist_sigma)**2) ! CM^2
  !
  ! THE LIQUID PHASE DIFFUSION CONSTANT (CM2/S)
  ! DL = DIFFUSION CONSTANT OF HOCL IN H2SO4/HNO3 FROM Shi ET AL., JGR 2001
  ! Diffusion COefficient HBR AND HOBR taken from Huthwelker et al.,
  ! J. ATMOS. CHEM., 21, 81-95, 1995
  dl = dhocl(wtsb,t)
  ! THE SOLUTION DENSITY (G/CM3)
  dens = density_ter(wtst,wtn,t)  
 ! LIQUID PHASE CONCENTRATIONS (MOL DM-3 ATM-1)
  chclliq = wtcl*dens/0.036461 
  choclliq = wthocl*dens/0.05246 
  chbrliq = wtbr*dens/0.08091 
  chobrliq = wthobr*dens/0.09691 
  ! CONVERT EFFECTIVE HENRY'S LAW CONSTANTS TO MOL/DM3/ATM
  factor = dens/(1.0 + mst*0.098076 + mn*0.063012)
  hhobr = hhobr*factor 
  hhbr = hhbr*factor 
  ! FIRST ORDER RATE CONSTANTS FOR REACTION IN LIQUID - DM3 MOL-1 S-1
  k1hbrhocl = k2hbrhocl*choclliq 
  k1hbrhobr = k2hbrhobr*chobrliq 
  k1hobrhcl = k2hobrhcl*chclliq 
  ! =====================================================================
  ! ========================== HOCL + HCL ===============================
  ! ======================= ON LIQUID AEROSOL ===========================
  ! Taken from Shi et al., JGR 2001
  ! First order rate constants depend on acid activity
  ! Gamma value is calculated later to account for HCl depletion in the bulk
  ! due to fast surface reaction with ClONO2
  
  ah = exp(60.51-0.095*(wtsb*100.)+0.0077*(wtsb*100.)**2-1.61E-5*(wtsb*100.)**3-(1.76&
       &+2.52E-4*(wtsb*100.)**2)*sqrt(t)+(-805.89+253.05*(wtsb*100.)**(0.076))/sqrt(t))
  k1hoclhcl = 1.25E9*ah*dhocl(wtsb,t)*mhcl
  lhocl = sqrt(dhocl(wtsb,t)/k1hoclhcl)
  fhocl = 1./tanh(rmean/lhocl)-lhocl/rmean
  cbar = 2009.0*sqrt(t)                      ![2009=sqrt(8.r/pi.mhocl)]
  grxnhocl = fhocl*4.*hhocl*0.082*t*sqrt(dhocl(wtsb,t)*k1hoclhcl)/cbar 
  
  ! =====================================================================
  ! ====================== CLONO2+HCL -> CL2+HNO3 =======================
  ! ====================== CLONO2+H2O -> HOCL+HNO3 ======================
  ! ========================= ON LIQUID AEROSOL =========================
  ! =====================================================================
  ! Parameters taken directly from Shi et al., JGR 2001, 106, 24259
  
  xmf = (wtsb*100.)/((wtsb*100.)+(100.-(wtsb*100.))*98./18.)
  ah2o = exp((-69.775*xmf-18253.7*xmf**2.+31072.2*xmf**3.-25668.8*xmf**4.)*(1./t&
       &-26.9033/t**2.)) 
  cbar2 = 1474.*sqrt(t)
  kh2o = 1.95E10*exp(-2800./t)
  kh = 1.22E12*exp(-6200./t)
  khydr = kh2o*ah2o+kh*ah*ah2o
  gbh2o = 4.*hclono2*0.082*t*sqrt(dclono2(wtsb,t)*khydr)/cbar2
  khcl = 7.9E11*ah*dclono2(wtsb,t)*mhcl
  lclono2 = sqrt((dclono2(wtsb,t)/(khydr+khcl)))
  fclono2 = 1./tanh(rmean/lclono2)-lclono2/rmean
  grxnclono2 = fclono2*gbh2o*sqrt(1.+khcl/khydr)
  gbhcl = grxnclono2*khcl/(khcl+khydr)
  gs = 66.12*exp(-1374/t)*hclono2*mhcl
  fhcl = 1./(1.+0.612*(gs+gbhcl)*pclono2/phcl)
  gs_strich = fhcl*gs
  gbhcl_strich = fhcl*gbhcl
  gb = gbhcl_strich+grxnclono2*khydr/(khcl+khydr)
  ge = 1./(1.+1./(gs_strich+gb))
  gclno3hcl = ge*(gs_strich+gbhcl_strich)/(gs_strich+gb)
  gclno3h2o = ge-gclno3hcl 
  khet1(24) = 0.25*gclno3hcl*1474.0*sqrt(t)*aliq/(max(chcl,cclno3) + e)*&
       gamma(24) 
  khet1(25) = 0.25*gclno3h2o*1474.0*sqrt(t)*aliq/ch2o*gamma(25) 
  !=====================================================================
  ! Gamma HOCl + HCl is calculated at this point to account for HCl depletion
  ! in the bulk
  ghoclhcl = 1./(1.+1./(grxnhocl*fhcl)) 
  khet1(23) = 0.25*ghoclhcl*cbar*aliq/(max(chcl,chocl) + e)*gamma(23)

  ! TAKEN DIRECTLY FROM HANSON AND RAVISHANKARA, J. PHYS. CHEM.,
  ! 98, 5728, 1994, EXCEPT FOR THE FUNCTION F - SEE FOLLOWING COMMENTS -
  !
  ! THE FUNCTION F: THE FORM OF F USED BY HANSON AND RAVISHANKARA CAN
  ! 'EXPLODE' UNDER CERTAIN CONDITIONS. IT HAS BEEN REPLACED HERE BY
  ! A STABLE FUNCTION THAT IS ACCURATE WITHIN ABOUT 4%. THIS IS ALSO
  ! THE CASE FOR OTHER REACTIONS THAT FOLLOW
  !
  ! =====================================================================
  ! ======================== N2O5 + H2O =================================
  ! ====================== ON LIQUID AEROSOL ============================
  ! =====================================================================
  !
  gn2o5h2o = 0.1 
  cbar = 1400.1*sqrt(t) 
  khet1(26) = 0.25*gn2o5h2o*cbar*aliq/ch2o*gamma(26) 
  ! ====================================================================
  ! =========================== HOBR + HCL =============================
  ! ======================= ON LIQUID AEROSOL ==========================
  ! ====================================================================
  !
  q = rmean*sqrt(k1hobrhcl/dl) 
  fq = (q + 0.312*q**2)/(3.0 + q + 0.312*q**2) 
  cbar = 1477.0*sqrt(t)                      ![1477=sqrt(8.r/pi.mhobr)] 
  ghobrhcl = fq/(fq + cbar/(4.0*hhobr*0.082*t*sqrt(k1hobrhcl*dl))) 
  khet1(27) = 0.25*ghobrhcl*cbar*aliq/(max(chcl,chobr) + e)*gamma(27) 
  ! =====================================================================
  ! ========================= HBR + HOCL/HOBR ===========================
  ! ======================= ON LIQUID AEROSOL ===========================
  ! =====================================================================
  !
  q = rmean*sqrt(k1hbrhobr/dl) 
  fq = (q + 0.312*q**2)/(3.0 + q + 0.312*q**2) 
  cbar = 1616.0*sqrt(t)                      ![1616=sqrt(8.r/pi.mhbr)] 
  ghbrhobr = fq/(fq + cbar/(4.0*hhbr*0.082*t*sqrt(k1hbrhobr*dl))) 
  ! RHBRHOBR IS MULTIPLIED BY PARTHBR TO ACCOUNT FOR GAS PHASE REMOVAL OF
  ! HBR TO THE DROPLETS.  THIS ENABLES TOTAL HBR TO BE USED IN RATE
  ! CALCULATION. NOTE THIS IS NOT NECESSARY FOR REACTIONS INVOLVING HCL
  ! SINCE CHCL APPEARS IN RATE EXPRESSIONS RATHER THAN HHCL.
  khet1(28) = 0.25*ghbrhobr*cbar*aliq*parthbr/(max(chobr,chbr) + e)*gamma(28) 
  !
  q = rmean*sqrt(k1hbrhocl/dl) 
  fq = (q + 0.312*q**2)/(3.0 + q + 0.312*q**2) 
  ghbrhocl = fq/(fq + cbar/(4.0*hhbr*0.082*t*sqrt(k1hbrhocl*dl))) 
  khet1(29) = 0.25*ghbrhocl*cbar*aliq*parthbr/(max(chocl,chbr) + e)*gamma(29) 
  ! =====================================================================
  ! ======================== BRONO2 + H2O ===============================
  ! ====================== ON LIQUID AEROSOL ============================
  ! =====================================================================
  ! FROM HANSON ET AL., JGR, 101, 9063-9069, 1996.
  ! =====================================================================
  !
  ah2o = 1013.25*ph2o/10**(9.217 - 2190.0/(t - 12.70)) 
  grxn = 211.0*ah2o**1.37 
  gbrno3h2o = (0.84*grxn)/(grxn + 0.84) 
  cbar = 1221.4*sqrt(t) 
  khet1(30) = 0.25*gbrno3h2o*cbar*aliq/ch2o*gamma(30) 
!
 
  return  
end subroutine hetliquid 


!
!
!
!
!
!
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
! ========================================================================
!
subroutine hetsolid(khet1, press, chcl, chbr, chocl, chobr, cclno3, &
     cbrno3, chno3, ch2o, aer_h2so4_default, parthno3, parth2o, n, densnat, &
     densice, denssat, anat, aice, asat, t,terminal_velocity,mean_free_path)  
  ! *******************************************************************
  ! HETEROGENEOUS RATES ON NAT, SAT AND ICE
  ! *******************************************************************
  !
  ! 06/07/1999 : gamma and param_nat_HR are now defined in module hetpar instead of being
  !              an argument of this subroutine (J.-U. Grooss)


  !
  !
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  use messy_clamschem_globalhet, ONLY: natcore, mixedpop, liqtest, t_lt_tnat, &
                                       kstate, laststate, &
                                       numhet, gamma, param_nat_HR


  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: press 
  real(PREC) , intent(in) :: chcl 
  real(PREC) , intent(in) :: chbr 
  real(PREC) , intent(in) :: chocl 
  real(PREC) , intent(in) :: chobr 
  real(PREC) , intent(in) :: cclno3 
  real(PREC) , intent(in) :: cbrno3 
  real(PREC) , intent(in) :: chno3 
  real(PREC) , intent(in) :: ch2o 
  real(PREC) , intent(in) :: aer_h2so4_default 
  real(PREC) , intent(in) :: parthno3 
  real(PREC) , intent(in) :: parth2o 
  real(PREC) , intent(in) :: n 
  real(PREC) , intent(in) :: densnat 
  real(PREC) , intent(in) :: densice 
  real(PREC) , intent(in) :: denssat 
  real(PREC) , intent(out) :: anat 
  real(PREC) , intent(out) :: aice 
  real(PREC) , intent(out) :: asat 
  real(PREC) , intent(in) :: t
  real(PREC) , intent(out) :: khet1(numhet+1) 
  real(PREC) , intent(out) :: terminal_velocity(2) 
  real(PREC) , intent(in) :: mean_free_path

  ! real(PREC) , intent(in) :: gamma(numhet) 
  ! logical, intent(in) :: param_nat_HR

!  real(PREC),external :: phmh2o
  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  real(PREC), parameter :: pi = 3.14159265358979324 
  real(PREC), parameter :: ctoa = 7.336E+21 
  real(PREC), parameter :: e = 1.0 
  real(PREC), parameter :: N_A = 6.022E+23  ! molec/mol
  real(PREC), parameter :: rho_ice = 0.92   ! g/cm3
  real(PREC), parameter :: m_h2o = 18.0     ! g/mol
  
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: knat, ksat, ph2o, phcl, sice, vhno3, radnat, tice, delt, c1, c2, &
       theta, c3, vh2o, radice, ch2so4, vsat, radsat, gamhclsathr0, c31, c32 
  real(PREC) :: slip_correction,viscosity
  real(PREC) :: rhonat

  !-----------------------------------------------

  viscosity = 1.72e-4      ! [g/(cm s)]
  rhonat = 1.6252-2.3585e-5*t ! [g/cm^3]
  ! pre-define terminal_velocity (jug, 26.06.2002)
  terminal_velocity=0.0

  ! NOTE: 'IF' STATEMENTS CONTAIN DENSNAT, DENSSAT AND DENSICE, RATHER
  ! THAN KSTATE, SINCE ALL THREE STATES CAN EXIST TOGETHER.
  !
  ph2o = ch2o/ctoa*t                         ! atm 
  phcl = chcl/ctoa*t*1.E+5                   ! pa 
  sice = ph2o/phmh2o(t) 
  !
  ! *--------------------------
  ! *  REACTIONS ON NAT
  ! *--------------------------
  if (densnat > 0.0) then 
     vhno3 = chno3*(1.0 - parthno3)*117.0/(N_A*rhonat)/(densnat + densice) 
     ! (DENSICE+DENSNAT) BECAUSE SOME NAT ARE COVERED IN ICE. WITHOUT THIS,
     ! THE NAT SIZE WOULD CHANGE WITH THE NUMBER OF ICE PARTICLES NUCLEATED!
     radnat = (3.0*vhno3/4.0/pi)**(1.0/3.0) 
     anat = densnat*4.0*pi*radnat**2 

     ! Berechnung der Fallgeschwindigkeit des Teilchens (Seinfeld, S. 466)
     ! Dichte : 1.6 g/cm^3 ???
     ! Genauer 1.6252 -2.3585e-5*T [g/cm^3]
     ! Dynamische Viskositaet : 1.72 g/(cm s)

     slip_correction = 1.+(mean_free_path/radnat)*(1.257+.4&
          &*exp(-1.1*radnat/mean_free_path))
     terminal_velocity(1) = (4.*radnat*radnat*rhonat*981.*slip_correction)&
          &/(m_h2o *viscosity)

     if (param_nat_HR) then 
        ! ==== REACTIONS FROM HANSON AND RAVISHANKARA DATA =====
        ! THE NON-CONSTANT GAMMA VALUES GIVEN HERE ARE GIVEN IN TABLE 1
        ! OF CARSLAW ET AL., GEOPHYS. RES. LETT., 24, 1747-, 1997.
        ! SEE ALSO CARSLAW AND PETER, GEOPHYS. RES. LETT., 24, 1743-, 1997.
        tice = 2668.7/(10.431 - (log(ph2o) + log(760.0))/log(10.0)) 
        delt = t - tice 
        ! H+R PARAMETRISATION IS VALID FOR PHCL = 0.75E-7 TORR ...
        ! corrected jug/tw 03/2009
        c1 = 1.0/(4.34 + 1.42*exp(0.518*delt)*1./(phcl*1.013E+5)**0.6) 
        !...SO INCLUDE WEAK PHCL DEPENDENCE DERIVED FROM THEIR OTHER DATA...
        ! THE FACTOR 1.013E+5 (1/PA) IS 1/(0.75E-7 TORR) - THE PHCL
        ! CORRESPONDING TO GAMHCLNATHR0 USED BY H+R.
        c2 = 1.0/(4.34 + 8403.0*exp((-2.81*sice))) 
        !
     else 
        !
        ! ==== REACTIONS FROM ABBATT AND MOLINA DATA =====
        ! (PARAMETRISED IN CARSLAW+PETER AND CARSLAW ET AL)
        knat = exp(-1.8 + 8.7*sice)          ! 1/pa 
        theta = knat*phcl/(1.0 + knat*phcl) 
        ! REARRANGED TO PREVENT EXPLOSION
        c1 = 8.2*theta/(3.3333*8.2*theta + 1.0) 
        !THE A+M FIT
        c2 = 1.0/(4.34 + 95798.3*exp((-4.97*sice))) 
     endif
     !
     ! HOCL + HCL REACTION SAME FOR HR AND A+M
     knat = exp(-1.8 + 8.7*sice)             ! 1/Pa 
     theta = knat*phcl/(1.0 + knat*phcl) 
     c3 = 5.1*theta/(6.6666*5.1*theta + 1.0) !REARRANGED TO PREVENT EXPLOSION 
     !
     !     C1=0.3       !VALUE GIVEN IN JPL
     !     C2=0.006     !VALUE GIVEN IN JPL
     !     C3=0.1       !VALUE GIVEN IN JPL
     !
     khet1(1) = gamma(1)*c1*4.56E+4*sqrt(t/97.46)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(1)*c1*radnat*press/t)/(max(chcl,cclno3) + e) 
     khet1(2) = gamma(2)*c2*4.56E+4*sqrt(t/97.46)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(2)*c2*radnat*press/t)/ch2o 
     khet1(3) = gamma(3)*c3*4.56E+4*sqrt(t/52.46)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(3)*c3*radnat*press/t)/(max(chocl,chcl) + e) 
     khet1(4) = gamma(4)*4.56E+4*sqrt(t/108.0)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(4)*radnat*press/t)/abs(chcl + e) 
     khet1(5) = gamma(5)*4.56E+4*sqrt(t/108.0)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(5)*radnat*press/t)/ch2o 
     khet1(6) = gamma(6)*4.56E+4*sqrt(t/97.46)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(6)*radnat*press/t)/(max(chbr,cclno3) + e) 
     khet1(7) = gamma(7)*4.56E+4*sqrt(t/142.0)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(7)*radnat*press/t)/(max(cbrno3,chcl) + e) 
     khet1(8) = gamma(8)*4.56E+4*sqrt(t/52.46)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(8)*radnat*press/t)/(max(chbr,chocl) + e) 
     ! NOW NEW REACTIONS THAT WERE NOT IN HETERO1-5...
     khet1(9) = gamma(9)*4.56E+4*sqrt(t/96.91)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(9)*radnat*press/t)/(max(chobr,chcl) + e) 
     khet1(10) = gamma(10)*4.56E+4*sqrt(t/96.91)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(10)*radnat*press/t)/(max(chobr,chbr) + e) 
     khet1(11) = gamma(11)*4.56E+4*sqrt(t/142.0)*radnat**2*densnat/ &
          (1.0 + 3.3E+4*gamma(11)*radnat*press/t)/ch2o 
  else 
     anat = 0.0 
     khet1(:11) = 0. 
  endif 
 
 
  ! *--------------------------
  ! *  REACTIONS ON ICE
  ! *--------------------------
  if (densice > 0.) then 
     vh2o = ch2o*(1.0 - parth2o)*m_h2o/(N_A*rho_ice)/densice 
     radice = (3.0*vh2o/4.0/pi)**(1.0/3.0) 
     aice = densice*4.0*pi*radice**2 
 
     ! Calculation for terminal_velocity and slip correction for ice (s.above)

     slip_correction = 1.+(mean_free_path/radice)*(1.257+.4&
          &*exp(-1.1*radice/mean_free_path))
     terminal_velocity(2) = (4.*radice*radice*rho_ice*981.*slip_correction)&
          &/(m_h2o*viscosity)
     
     khet1(12) = gamma(12)*4.56E+4*sqrt(t/97.46)*radice**2*densice/ &
          (1.0 + 3.3E+4*gamma(12)*radice*press/t)/(max(cclno3,chcl) + e) 
     khet1(13) = gamma(13)*4.56E+4*sqrt(t/97.46)*radice**2*densice/ &
          (1.0 + 3.3E+4*gamma(13)*radice*press/t)/ch2o 
     khet1(14) = gamma(14)*4.56E+4*sqrt(t/52.46)*radice**2*densice/ &
          (1.0 + 3.3E+4*gamma(14)*radice*press/t)/(max(chocl,chcl) + e) 
     khet1(15) = gamma(15)*4.56E+4*sqrt(t/108.0)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(15)*radice*press/t)/abs(chcl + e) 
     khet1(16) = gamma(16)*4.56E+4*sqrt(t/108.0)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(16)*radice*press/t)/ch2o 
     khet1(17) = gamma(17)*4.56E+4*sqrt(t/97.46)*radice**2*densice/ &
          (1.0 + 3.3E+4*gamma(17)*radice*press/t)/(max(chbr,cclno3) + e) 
     khet1(18) = gamma(18)*4.56E+4*sqrt(t/142.0)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(18)*radice*press/t)/(max(cbrno3,chcl) + e) 
     khet1(19) = gamma(19)*4.56E+4*sqrt(t/52.46)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(19)*radice*press/t)/(max(chbr,chocl) + e) 
     ! NOW NEW REACTIONS THAT WERE NOT IN HETERO1-5...
     khet1(20) = gamma(20)*4.56E+4*sqrt(t/96.91)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(20)*radice*press/t)/(max(chobr,chcl) + e) 
     khet1(21) = gamma(21)*4.56E+4*sqrt(t/96.91)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(21)*radice*press/t)/(max(chobr,chbr) + e) 
     khet1(22) = gamma(22)*4.56E+4*sqrt(t/142.0)*radice**2*densice/ &
          (1.0 +3.3E+4*gamma(22)*radice*press/t)/ch2o 
 
  else 
     aice = 0.0 
     khet1(12:22) = 0. 
  endif
 
  ! *--------------------------
  ! *  REACTIONS ON SAT
  ! *--------------------------
  if (denssat > 0) then 
     ch2so4 = aer_h2so4_default*press*7.216E+9/t*denssat/n 
     ! (DENSSAT/N) ENSURES THAT H2SO4 IN FORM OF SAT IS SCALED
     ! IN TERMS OF THE FRACTION OF H2SO4 AEROSOLS THAT ARE SAT
     vsat = ch2so4*170.0/(N_A*1.6)/(denssat + densnat + densice) 
     ! SEE COMMENT FOR VHNO3 ABOVE
     radsat = (3.0*vsat/4.0/pi)**(1.0/3.0) 
     asat = denssat*4.0*pi*radsat**2 
     if (param_nat_HR) then 
        ! ==== REACTIONS FROM HANSON AND RAVISHANKARA DATA =====
        tice = 2668.7/(10.431 - (log(ph2o) + log(760.0))/log(10.0)) 
        delt = t - tice 
        ! H+R PARAMETRISATION IS VALID FOR PHCL = 1.5E-7 TORR (NO PHCL DEPENDENCE)
        gamhclsathr0 = exp(-0.636 - 0.4802*delt) 
        c31 = gamhclsathr0/(1.0 + gamhclsathr0/0.23) 
        ! ABOVE EXPRESSION WAS REARRANGED TO APPROACH 0.23 BELOW TICE WHICH WAS
        ! GIVEN BY H+R.
        c32 = 1.0/(4.34 + 1652.0*exp((-2.94*sice))) 
        !
     else 
        !
        ! ==== REACTIONS FROM ABBATT AND MOLINA DATA =====
        ! (PARAMETRISED IN CARSLAW+PETER)
        ksat = exp(0.59 + 7.3*sqrt(sice))    ! 1/pa 
        theta = ksat*phcl/(1.0 + ksat*phcl) 
        !!rearranged to prevent
        c31 = 1.1*theta/(10.0*1.1*theta + 1.0) 
        c32 = 1.0/(4.34 + 1652.0*exp((-2.94*sice))) 
     endif
     !
     !  c25=0.0
     !  c26=0.0
     !
     khet1(31) = gamma(31)*c31*4.56E+4*sqrt(t/97.46)*radsat**2*denssat/(1.0&
          + 3.3E+4*gamma(31)*c31*radsat*press/t)/(max(cclno3,chcl) + e) 
     khet1(32) = gamma(32)*c32*4.56E+4*sqrt(t/97.46)*radsat**2*denssat/(1.0&
          + 3.3E+4*gamma(32)*c32*radsat*press/t)/ch2o 
     khet1(33) = gamma(33)*4.56E+4*sqrt(t/108.0)*radsat**2*denssat/(1.0 + &
          3.3E+4*gamma(33)*radsat*press/t)/ch2o 
 
     !  write(112,*)t,c31
     !  write(113,*)t,c32
 
  else 
     asat = 0.0 
     khet1(31:33) = 0. 
!
  endif
  return  
end subroutine hetsolid 


! ======================================================================

function phmhno3 (t, ph2o) 

  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none

  real(PREC) :: phmhno3
  !
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: t 
  real(PREC) , intent(in) :: ph2o 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: dum1, dum2 
  !-----------------------------------------------
  !     CALCULATES PHNO3 IN ATM FROM HANSON AND MAUERSBERGER EQUATION
  dum1 = (-2.7836) - 0.00088*t 
  dum2 = 38.9855 - 11397.0/t + 0.009179*t 

  if (ph2o <= 0.) write(*,*) 'warning phmhno3 ',ph2o, t, dum1,dum2

  phmhno3 = 10.0**(dum1*log10(ph2o*760.0) + dum2)/760.0 
  return  
end function phmhno3 


!=======================================================================
function phmh2o (t) 
  ! CALCULATES PH2O IN ATM FROM HM EQUATION
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  real(PREC) :: phmh2o
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC) , INTENT(IN) :: T 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  !-----------------------------------------------
  phmh2o = 10**(10.431 - 2668.7/t)/760.0 
  return  
end function phmh2o 


! ======================================================================
function dhocl (wtsb, t) 

  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  real(PREC) :: dhocl
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC), intent(in)  :: wtsb 
  real(PREC), intent(in)  :: t 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: c, a, t0,  visc 
  !-----------------------------------------------
  !
  !
  !
  !**********************************************************************
  !DIFFUSION COEFFICIENT HOCL in CM2/S
  !TAKEN FROM SHI ET AL, JGR, 106, 24259-, 2001.
  !C:WT % OF H2SO4
  !VISC:VICOSITY IN CP
  !**********************************************************************
  !
  !
  !
  c = wtsb*100.0 
  a = 169.5+5.18*c-0.0825*c**2.+3.27E-3*c**3. 
  t0 = 144.11+0.166*c-0.015*c**2.+2.18E-4*c**3. 
  visc = a*t**(-1.43)*exp(448./(t-t0)) 
  !
  dhocl = 6.4E-8*t/visc
  return  

end function dhocl 


! ====================================================================
function dclono2 (wtsb, t) 

  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  real(PREC) :: dclono2
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  real(PREC), intent(in)  :: wtsb  
  real(PREC), intent(in)  :: t  
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  real(PREC) :: c, a, t0,  visc 
  !-----------------------------------------------
  !
  !
  !
  !**********************************************************************
  !DIFFUSION COEFFICIENT ClONO2 in CM2/S
  !TAKEN FROM SHI ET AL, JGR, 106, 24259-, 2001.
  !C:WT % OF H2SO4
  !VISC:VICOSITY IN CP!
  !**********************************************************************

  c = wtsb*100.0 
  a = 169.5+5.18*c-0.0825*c**2.+3.27E-3*c**3. 
  t0 = 144.11+0.166*c-0.015*c**2.+2.18E-4*c**3.
  visc = a*t**(-1.43)*exp(448./(t-t0)) 
  !
  dclono2 = 5E-8*t/visc
  return  

end function dclono2

! ====================================================================

function density (msb, t) 
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  real(PREC) :: density

  !-----------------------------------------------
  !   d u m m y   a r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: msb  
  real(PREC) , intent(in) :: t 
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s
  !----------------------------------------------- 
  real(PREC) :: z1, z2, z3 
  !-----------------------------------------------
  !
  !
  !*********************************************************************
  !     DENSITY OF BINARY SOLUTION IN G/CM3
  !     TAKEN FROM SHI ET AL, JGR, 106, 24259-, 2001.
  !     MS IS MOLALITY H2SO4 in mol kg-1
  !*********************************************************************
  !
  !
  !
  z1 = 0.12364-5.6E-7*t**2.
  z2 = -0.02954+1.814E-7*t**2.
  z3 = 2.343E-3-1.487E-6*t-1.324E-8*t**2.
  density = 1.+z1*msb+z2*msb**(1.5)+z3*msb**2. 
  !
  return  
end function density 

function density_ter (wtst, wn, t) 
  !...Translated by Pacific-Sierra Research 77to90  4.3E  12:46:19   7/ 2/99  
  !...Switches: -ymhf               
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
  USE messy_clams_global,     ONLY: prec
  implicit none
  real(PREC) :: density_ter

  !-----------------------------------------------
  !   d u m m y   a r g u m e n t s
  !-----------------------------------------------
  real(PREC) , intent(in) :: wtst
  real(PREC) , intent(in) :: wn 
  real(PREC) , intent(in) :: t 
  !-----------------------------------------------
  !   l o c a l   v a r i a b l e s
  !-----------------------------------------------
  real(PREC) , dimension(22) :: x 
  real(PREC) :: w, wh, v1, vs, vn, vmcal 
  !-----------------------------------------------
  !
  !
  !*********************************************************************
  !     DENSITY OF TERNARY SOLUTION IN G/CM3
  !     TAKEN FROM LUO ET AL, GRL, 23, 3707-, 1996.
  !     WS ,WN ARE WT FRACTION,
  !     FITTED TO 0.05<WS+WN<0.70 WT FRACTION, BUT EXTRAPOLATES WELL
  !     185 < T (K)
  !*********************************************************************
  !
  !
  !
  DATA X/ 2.393284E-02, -4.359335E-05, 7.961181E-08, 0.0, -0.198716351, &
       1.39564574E-03, -2.020633E-06, 0.51684706, -3.0539E-03, 4.505475E-06, &
       -0.30119511, 1.840408E-03, -2.7221253742E-06, -0.11331674116, &
       8.47763E-04, -1.22336185E-06, 0.3455282, -2.2111E-03, 3.503768245E-06, &
       -0.2315332, 1.60074E-03,  -2.5827835E-06/  
  w = wtst + wn 
  wh = 1.0 - w 
  v1 = x(1) + x(2)*t + x(3)*t**2 + x(4)*t**3 
  vs = x(5) + x(6)*t + x(7)*t**2 + (x(8)+x(9)*t+x(10)*t**2)*w + &
       (x(11)+x(12)*t+x(13)*t**2)*w*w 
  vn = x(14) + x(15)*t + x(16)*t**2 + (x(17)+x(18)*t+x(19)*t**2)*w + &
       (x(20)+x(21)*t+x(22)*t**2)*w*w 
  vmcal = wh/18.0160*v1 + vs*wtst/98.080 + vn*wn/63.0160 
  density_ter = 1.0/vmcal*0.001 
  !
  return  
end function density_ter 


end Module messy_clamschem_hetero_shi
