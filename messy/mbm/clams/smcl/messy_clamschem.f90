!*******************************************************************************
! MODULE FOR CLaMS CHEMISTRY MODULE 'CHEM'
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2010
! Jens-Uwe Grooss, Daniel S. McKenna, Rolf Mueller, Nicole Thomas,
! Glenn Carver, David Lary, Kenneth S. Carslaw, Tobias Wegner
! Forschungszentrum Juelich GmbH
!
! Last Modified By: N.Thomas
! Last Modified On: Fri Oct  2 10:01:11 2020
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
!
! Commented by DMcK 29.08.96
! 
! Changed version for MESSy interface 19.10.2011 (C.Hoppe) 
!           based on Version from 13.07.1999 (J.-U.Grooss)
!
! This is a box model program capabile of running on a multiple
! trajectory dataset written in netCDF format. It is the main program
! of the CLaMS module chem.
!
! The integration is achieved using the ASAD package provided by
! Glenn Carver at Cambridge.
!
! The photolysis code was provide by Rolf Mueller and is basically
! the David Lary code with a few additions.
!
! the heterogeneous chemisty is based on the code of Ken Carslaw
! and used analytical expressions to determine the condensation of
! sulphate, nitrate and water vapour.
!
! ------------------------------------------------------------------------------
MODULE messy_clamschem

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamschem'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

CONTAINS
  
!**************************************************************************
!
!**************************************************************************
  SUBROUTINE chem (status, CHEMSPECARR, LAT, LON, LEV, TEMPERATURE, PRESSURE,  &
                   LAT_OLD, LON_OLD, LEV_OLD, TEMPERATURE_OLD, PRESSURE_OLD, &
                   BRATES, TRATES, JRATES, HRATES, &
                   BCONST, TCONST, JCONST, HCONST, HETPAR)

  !...Translated by Pacific-Sierra Research 77to90  4.3E  15:40:33   7/23/98  
  !

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! INPUTS
  !
  ! CONTROLS
  !
  ! the main control inputs including initializations come from
  ! namelist read statements which are included in the .exec file
  ! and are written to a temporary file chem.inp
  !
  ! namelist CNTL various parameters to control the trajectory
  ! integration operation etc
  !
  ! ldepw, ldepd
  ! logical flags set .true. to if any deposition losses are
  ! to be calculated.NB this does not work here.
  !
  ! lemit
  ! logical flag set .true. to if any emissions are
  ! to be calculated.NB this does not work here.
  !
  ! lhet
  ! logical flag set .true. to if any heterogeneous reactions
  ! are to be calculated.
  !
  ! lphotol
  ! logical flag set .true. to if any photolysis reactions
  ! are to be calculated.
  !
  ! lvmr
  ! logical flag set .true. to if volume mixing rations and not
  ! number densities are to be stored in output.
  !
  ! method
  ! integer indicating which integration method to employ
  !      method = 1 :  Standard ASAD integrator
  !
  ! ncdt
  ! constant chemical time step in seconds to be
  ! used by ASAD
  !
  ! nit0 constant  used by ASAD to control integration
  !
  ! nitfg constant  used by ASAD to control integration
  !
  ! nrsteps  constant  used by ASAD to control integration
  !
  ! IODUMP  logical if diagnostic dumps are required.
  !
  !
  ! IODUMPo  logical if larger dumps are required.
  !
  ! rates
  ! logical flag set to .true. if the elementary reaction
  ! rates are required to be stored.
  !
  ! const
  ! logical flag set to .true. if the elementary reaction
  ! rate constants are required to be stored.
  !
  ! trajectory data:
  ! This data which should comform to the standard netCDF trajectory
  ! format and should include at least pressure temperature data will 
  ! normally reside in INPDATA.
  ! The solar zenith angle (SZA) SZA will be calculated by the model.
  !
  ! OUTPUTS:
  !
  ! a netCDF file containing an almost total description of the
  ! the trajectory integration 
  !
  !------------------------------------------------------------------------------
  !
  ! Timesteps:
  ! There are two different timesteps in the model:
  ! 1) The chemical input timestep: timestep_chem (seconds). 
  !    This corresponds to the frequency at which the chemistry routine is called. 
  ! 2) The regular internal chemistry timestep: ncdt (seconds). It is given by
  !    the namelist cntl.
  !
  !-------------------------------------------------------------------------------
  !
  !   M o d u l e s 
  !
  !-------------------------------------------------------------------------------

  ! MESSy Main
  USE messy_main_constants_mem, ONLY: DP

  ! ASAD
  use messy_clamschem_asad_mod, ONLY: method, &
                                     asad_mod_init_loop, asad_mod_final_loop, &
                                     speci, advt, family, madvtr, majors, &
                                     moffam, ctype, p, t, converged, positive
  use messy_clamschem_asad_mod_clams, ONLY: jpspec, jpctr, &
                                     lhet, theta_field_size, mype

  USE messy_clamschem_asad_cinit, ONLY: asad_cinit_loop
  USE messy_clamschem_asad_ftoy,  ONLY: ASAD_FTOY
  USE messy_clamschem_asad_cdrive,ONLY: ASAD_CDRIVE

  ! CLaMS
  use messy_clamschem_global,    ONLY: rc_type, ntraj, ipart, npacks, jpnl_count, &
                                       iodump, iodumpo, zangle, &
                                       ftr, missing_index

  USE messy_clamschem_dynamic, ONLY: indynam, dynamic
  USE messy_clamschem_specinit,ONLY: specinit
  USE messy_clamschem_mixinit, ONLY: mixinit
  USE messy_clamschem_mixadd,  ONLY: mixadd
  USE messy_clamschem_inhet,   ONLY: inhet

  USE messy_clams_global,      ONLY: rank, species_type, nchemspec, init_vertcoorname

  implicit none

  INTEGER :: status
  
  ! MESSy species array:
  TYPE(species_type), DIMENSION(:), pointer :: CHEMSPECARR

  TYPE(species_type), DIMENSION(:), POINTER :: HETPAR 

  TYPE(rc_type),      DIMENSION(:), POINTER :: BRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: TRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: JRATES
  TYPE(rc_type),      DIMENSION(:), POINTER :: HRATES

  TYPE(rc_type),      DIMENSION(:), POINTER :: BCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: TCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: JCONST
  TYPE(rc_type),      DIMENSION(:), POINTER :: HCONST

  ! other MESSy arrays:
  REAL(DP) :: LAT(:)
  REAL(DP) :: LON(:) 
  REAL(DP) :: LEV(:)
  REAL(DP) :: TEMPERATURE(:)
  REAL(DP) :: PRESSURE(:)
  REAL(DP) :: LAT_OLD(:)
  REAL(DP) :: LON_OLD(:)
  REAL(DP) :: LEV_OLD(:)
  REAL(DP) :: TEMPERATURE_OLD(:)
  REAL(DP) :: PRESSURE_OLD(:)

  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: i, ihelp
  real(DP) :: prdim
  real(DP), dimension(:,:), allocatable :: cdot 
  integer       :: startclock,endclock,counts_per_sec
  real(DP)      :: seconds

  CALL system_clock (count=startclock)

  status = 0 ! no error

  if(rank==0)THEN

     WRITE (*,*) 
     WRITE (*,*) 'START OF CHEM'
     WRITE (*,*) 

     ! 1=IMPACT; 3=N-R solver; 5=Backward-Euler; 11=SVODE
     SELECT CASE (method)
     CASE (1)
        write (*,*) 'Method of integration:  IMPACT'
     CASE (3)
        write (*,*) 'Method of integration:  N-R solver'
     CASE (5)
        write (*,*) 'Method of integration:  Backward-Euler'
     CASE (11)
        write (*,*) 'Method of integration:  SVODE'
     CASE DEFAULT
        write (*,*) 'Method of integration:  UNKWOWN !!!'
     END SELECT
     WRITE (*,*) 

  ENDIF

!!!!!! => SMIL
!   ! Get chemistry input
!   call clams_chem_init

!!!!!!  => SMIL
!   ! Allocate and initialize ASAD arrays and variables
!   call messy_clamschem_asad_mod_init


!!!!!!  => SMIL
!!!!! Nach messy_clamschem_asad_mod_init belegen:
!!$  ldepd = .FALSE. ! in messy_clamschem_asad_mod.f90 deklariert
!!$  ldepw = .FALSE. ! in messy_clamschem_asad_mod.f90 deklariert
!!$  nit0    = 20    ! in messy_clamschem_asad_mod: mit 20 initialisiert
!!$  nitfg   = 20    ! in messy_clamschem_asad_mod_init: mit 10 initialisiert
!!$  nitnr   = 20    ! in messy_clamschem_asad_mod_init: mit 10 initialisiert
!!$  nrsteps = 50    ! in messy_clamschem_asad_mod_init: mit 45 initialisiert
!!$  cdt     = ncdt  ! in messy_clamschem_asad_mod.f90 deklariert
!!$
!!$  dtime = timestep_chem   ! dtime in messy_clamschem_asad_mod_init initialisiert !
 
  ! Initialize the trajectory calling routine
  call indynam (status, prdim)    
  if (status /= 0) return

    
!!!!! => SMIL
  ! ! Allocate arrays 
  ! call allocate_chem_vars

   ! loop over (potentially parallel) packages of trajectories
   do ipart = 1, npacks 

      if (iodump) write(*,'(4(A,I4),A)') 'rank', rank, ':  part ', ipart, &
                     ' of ',npacks,' trajectory packages'

      ntraj=jpnl_count(ipart)
      !write (*,'(A,I4,A,I12)') 'rank ',rank,': ntraj=',ntraj
     
!!!!! Dimension theta_field_size fuer ASAD !!!
      theta_field_size = ntraj 

     if (method==3) allocate(positive(ntraj),converged(ntraj))


!!!!!
     ! Allocate and initialize ASAD arrays and variables
     if (iodump) write (*,*) 'call messy_clamschem_asad_mod_init_loop'
     call asad_mod_init_loop

!!!!!
     allocate (cdot(ntraj,jpctr))
     allocate (ftr(ntraj,jpctr))
     allocate (zangle(ntraj))
     allocate (missing_index(ntraj))

     call dynamic(status, 0, prdim, &
                  LAT, LON, LEV, TEMPERATURE, PRESSURE,&
                  LAT_OLD, LON_OLD, LEV_OLD, TEMPERATURE_OLD, PRESSURE_OLD)
     if (status /= 0) return

!!!!! => SMIL: clamschem_initialize 
!!$     ! --------------------------------------
!!$     !    Initialise the Chemistry
!!$     ! --------------------------------------
!!$     !
!!$     ! Call ASAD routine CINIT
!!$     !
!!$     call asad_cinit (1)
     if (rank==0) write (*,*) 'call asad_cinit_loop !!!'
     call asad_cinit_loop 
     
!!$     if (iodump .and. rank==0) then 
!!$        write (6, *) 'speci= ', (speci(i),i=1,jpspec) 
!!$        write (6, *) 'advt= ', (advt(i),i=1,jpctr) 
!!$        write (6, *) 'family= ', (family(i),i=1,jpspec) 
!!$        write (6, *) 'madvtr= ', (madvtr(i),i=1,jpspec) 
!!$        write (6, *) 'majors= ', (majors(i),i=1,jpctr) 
!!$        write (6, *) 'moffam= ', (moffam(i),i=1,jpspec) 
!!$        write (6, *) 'ctype= ', (ctype(i),i=1,jpspec) 
!!$     endif
     
     !  Set up output array for the chemical species
     if (ipart == 1) then
        if (iodump .and. rank==0)  write (*,*) 'call specinit'
        call specinit 
     endif

 
     ! Initialize the composition array FTR
     if (iodump .and. rank==0)  write (*,*) 'call mixinit'
     call mixinit (status, CHEMSPECARR)
     if (status /= 0) return
     
     ! Initialise heterogeneous chemistry
     ! (should be called from cinit according to
     ! ASAD standard, but we need the concentrations
     ! of HNO3 etc and thus ftr from mixinit --
     ! cinit must be called before mixinit...)
!!$ LHET
     if (lhet) then
        if (iodump .and. rank==0)  write (*,*) 'call inhet'
        call inhet (CHEMSPECARR) 
     endif


     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Chemistry integration starts here
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     if (iodump .and. rank==0)  write (*,*) 'call dynamic'
     call dynamic (status, 1, prdim, &
                   LAT, LON, LEV, TEMPERATURE, PRESSURE, &
                   LAT_OLD, LON_OLD, LEV_OLD, TEMPERATURE_OLD, PRESSURE_OLD) 
     if (status /= 0) return

     ! call chemistry integration
     ! write(*,*)'chem 1:',ftr(1,1:5)
     call asad_cdrive (cdot, ftr, p, t, 1, ntraj) 
     !write(*,*)'chem 2:',ftr(1,1:5)
     
    
        
     ! make sure that the species array y contains the right numbers after
     ! one timestep in the case of a stiff integrator, since it is used
     ! for the output.        --- J.U. Grooss, 29.05.1998
     ihelp = 0
     if (method >= 10) call asad_ftoy (.FALSE., ihelp, ntraj) 
    
       
     ! Add calculated data points into the output netCDF dataset
     !write(*,*)'chem 3:',ftr(1,1:5)
     call mixadd (CHEMSPECARR, &
                  BRATES, TRATES, JRATES, HRATES, &
                  BCONST, TCONST, JCONST, HCONST, HETPAR)
 
     ! clean up
     if (iodump) write (*,*) 'clean up'
     call clean_up
!!!!! lokal deklariert:
     deallocate (cdot)

     ! Deallocate ASAD arrays
     if (iodump) write (*,*) 'call asad_mod_final_loop',rank
     call asad_mod_final_loop
     if (iodump) write (*,*) 'nach asad_mod_final_loop',rank

     
  enddo  !ipart
  

!!!!! => SMIL
  ! ! Deallocate ASAD arrays
  ! if (iodump) write (*,*) 'call asad_mod_final'
  ! call asad_mod_final

!!!!! Deallocate => SMIL
!  call deallocate_chem_vars

  IF (rank==0) THEN
     WRITE (*,*) 'Normal termination of chem'  
                  
     CALL system_clock (COUNT_RATE=counts_per_sec)
     CALL system_clock (count=endclock)
     IF (endclock > startclock) THEN
        seconds = float(endclock-startclock) / float(counts_per_sec)
     ELSE
        seconds = 0.0
     ENDIF
     
     !WRITE (*,*)
     !WRITE (*,*) 'System clock runs at ', counts_per_sec, 'ticks per second'
     WRITE (*,*)
     WRITE (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.' 
  ENDIF
  
END SUBROUTINE chem

!*********************************************************************

  SUBROUTINE clamschem_read_nml(status, iou)

    !Set default values and read namelist

    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, &
                                        read_nml_close
    USE messy_main_constants_mem, ONLY: DP
    USE messy_clams_global,       ONLY: YEAR, MONTH, DAY, HOUR

    USE messy_clamschem_asad_mod,                 ONLY: nfphot, method
    USE messy_clamschem_asad_mod_clams,           ONLY: lhet, lphotol, nfhet, &
                                        chemdata_type, mype

    USE messy_clams_global,       ONLY: rank, ldiagout
    USE messy_clamschem_global,   ONLY: ncdt, timestep_chem, &
                                        dsn_twodavg, iodump, iodumpo, &
                                        rates, const, emrates, hetparam, emit
    USE messy_clamschem_globalhet,ONLY: liq_sdist_sigma,  &
                                        sat_meltallowed, param_nat_HR, &
                                        transform, saturation_criteria, gamma, &
                                        aer_h2so4_default, densaero_default, &
                                        ciceinit, cnatinit

    IMPLICIT NONE
    
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamschem_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output
    INTEGER              :: i, k
    
    NAMELIST /CTRL/ timestep_chem, lhet, lphotol, method, nfphot, ncdt, &
                    iodump, iodumpo, &
                    dsn_twodavg, rates, const, emrates, hetparam, &
                    chemdata_type, emit

    NAMELIST /CTRL_HETERO/ liq_sdist_sigma, densaero_default, aer_h2so4_default, &
                           ciceinit, cnatinit, sat_meltallowed, param_nat_HR, &
                           transform, saturation_criteria, gamma

   
    status = 1 !ERROR 
    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default values (taken from original CHEM main program):
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    mype = rank  ! rank: module global, mype: messy_clamschem_asad_mod_clams => ASAD

    ! Initialize variables for call CINIT
    lhet    = .TRUE.    ! messy_clamschem_asad_mod_clams.f90 -> wird in cdrive genutzt
    lphotol = .TRUE.    ! messy_clamschem_asad_mod_clams.f90 -> wird in cdrive genutzt
    
    ! -------------------------------------------------------------------
    ! Setting default values for the control parameters

    method = 1      ! method of time integration: 
                    ! (1=IMPACT; 3=N-R solver; 5=Backward-Euler; 11=SVODE)

    timestep_chem = 3600 ! timestep for calling chem (s)

    ncdt = 600      ! internal chemistry timestep (s)

    nfphot = 0      ! call photolysis every chemical timestep but not substep
                    ! bisher: in chcctl.h dekl., in cinit.f genutzt
                    ! jetzt: in messy_clamschem_asad_mod.f90 dekl., in asad_cinit.f90 genutzt

    nfhet  = 0      ! call hetero every chemical timestep but not substep
                    ! bisher: in cmcctl.h dekl. -> messy_clamschem_asad_mod_clams.f90


    ! Default settings NML IO (from dynamic.f90):
    dsn_twodavg='/dat_icg1/twod_mz/avg00T05b.ref.nc'

    ! Default settings NML SIO (from specinit.f90):
    rates=.false.
    const=.false.
    emrates=.false.
    hetparam=.false.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set default values for heterogenous chemistry
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    liq_sdist_sigma = 1.8     ! sigma of log normal size distribution of liq. aerosol
    densaero_default = 10.0   ! default liquid aerosol number density [cm^-3] 
    aer_h2so4_default = 0.2   ! default gasphase H2SO4 equivalent mixing ratio [ppbv]
    ciceinit = 0.003          ! initial ice particle number density [cm^-3] 
    cnatinit = 0.003          ! initial NAT particle number density [cm^-3] 
    sat_meltallowed = .false. ! allow SAT melting (Koop and Carslaw, 1996)
    param_nat_HR = .true.     ! Hanson and Ravishankara param. for NAT/SAT reactions
                              ! instead of Abbatt and Molina

    ! 4x4 transformation matrix for particle phases
    transform(1,:) = (/  0,  0,  0,  0 /)        
    transform(2,:) = (/  1,  0,  0,  1 /)
    transform(3,:) = (/ 99,  1,  0,  1 /)
    transform(4,:) = (/ 99, 99,  1,  0 /)

    ! corresponding saturation criteria for phase transformations
    saturation_criteria(1,:) = (/   0., 1., 100., 1. /)
    saturation_criteria(2,:) = (/   1., 0.,   1., 1. /)
    saturation_criteria(3,:) = (/   1., 1.,   0., 1. /)
    saturation_criteria(4,:) = (/   0., 0.,   1., 0. /)

    ! gamma values for NAT, ICE, liquid and SAT
    !  #### Reactions on NAT ####
    gamma(1) = 1.0      ! ClONO2 + HCl   / set to 1 or 0 as a switch
    gamma(2) = 1.0      ! ClONO2 + H2O   / set to 1 or 0 as a switch
    gamma(3) = 1.0       ! HOCl + HCl     / set to 1 or 0 as a switch
    gamma(4) = 0.003     ! N2O5 + HCl
    gamma(5)= 0.0003     ! N2O5 + H2O
    gamma(6) = 0.3       ! ClONO2 + HBr
    gamma(7) = 0.3       ! BrNO3 + HCl
    gamma(8) = 0.3       ! HBr + HOCl
    gamma(9) = 0.1       ! HOBr + HCl
    gamma(10) = 0.1      ! HOBr + HBr
    gamma(11) = 0.001    ! BrONO2 + H2O
    !  #### Reactions on ICE ####
    gamma(12) = 0.3     ! ClONO2 + HCl
    gamma(13) = 0.3      ! ClONO2 + H2O
    gamma(14) = 0.3      ! HOCl + HCl
    gamma(15) = 0.03     ! N2O5 + HCl
    gamma(16) = 0.01     ! N2O5 + H2O
    gamma(17) = 0.3      ! ClONO2 + HBr
    gamma(18) = 0.3      ! BrNO3 + HCl
    gamma(19) = 0.3      ! HBr + HOCl
    gamma(20) = 0.3      ! HOBr + HCl
    gamma(21) = 0.1      ! HOBr + HBr
    gamma(22) = 0.3      ! BrONO2 + H2O
    !  #### Reactions on liquid aerosol ####
    gamma(23) = 1.0     ! HOCl + HCl    / set to 1 or 0 as a switch
    gamma(24) = 1.0      ! ClONO2 + HCl   / set to 1 or 0 as a switch
    gamma(25) = 1.0      ! ClONO2 + H2O   / set to 1 or 0 as a switch
    gamma(26) = 1.0      ! N2O5 + H2O     / set to 1 or 0 as a switch
    gamma(27) = 1.0      ! HOBr + HCl     / set to 1 or 0 as a switch
    gamma(28) = 1.0      ! HBr + HOBr     / set to 1 or 0 as a switch
    gamma(29) = 1.0      ! HBr + HOCl     / set to 1 or 0 as a switch
    gamma(30) = 1.0      ! BrONO2 + H2O   / set to 1 or 0 as a switch
    !  #### Reaction on SAT ####
    gamma(31) = 1.0     ! ClONO2 + HCl   / set to 1 or 0 as a switch
    gamma(32) = 1.0      ! ClONO2 + H2O   / set to 1 or 0 as a switch
    gamma(33) = 0.006    ! N2O5 + H2O
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read namelist variables:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! read namelist CTRL_HETERO if heterogenous chemistry is switch on
    IF (lhet) THEN

       if (rank==0 .and. iodump) then
          write (*,*) 'Matrizen nach Initialisierung:'
          write (*,*) 'transform:'
          do i = 1,4
             write (*,'(4F9.3)') (transform(i,k),k=1,4)
          enddo
          write (*,*) 'saturation:'
          do i = 1,4
             write (*,'(4F9.3)') (saturation_criteria(i,k),k=1,4)
          enddo
       endif

       ! transpose matrices
       ! Falls die CTRL_HETERO vorhanden ist, aber die Matrizen dort nicht gesetzt werden,
       ! werden sie weiter oben zunaechst richtig (2. Dimension zuerst laufend) initialisiert.
       ! Nach dem Einlesen von CTRL_HETERO (ohne Matrizen) wird auch in diesem Fall dann das 
       ! "transpose" ausgefuehrt. Dann sind sie aber (falsch!) transponiert abgelegt.
       ! Daher wird vor dem read_nml bereits ein transpose ausgefuehrt, damit sie in diesem
       ! Fall wieder richtig gespeichert sind.
       transform = transpose(transform)
       saturation_criteria = transpose(saturation_criteria)

       READ(iou, NML=CTRL_HETERO, IOSTAT=fstat)
       CALL read_nml_check(fstat, substr, iou, 'CTRL_HETERO', modstr, l_print)
       IF (fstat /= 0) RETURN  ! error while reading namelist

       ! transpose matrices 
       ! Die in der CTRL_HETERO angegebenen Matrizen muessen dort - wie oben in der Initialisierung -
       ! zeilenweise angegeben werden. Fortran liest sie dann spaltenweise! (1. Dim. zuerst laufend)
       ! ein. Daher muss in diesem Fall ein transpose ausgefuehrt werden.
       transform = transpose(transform)
       saturation_criteria = transpose(saturation_criteria)

       if (rank==0 .and. iodump) then
          write (*,*) 'Matrizen nach Einlesen der Namelist:'
          write (*,*) 'transform:'
          do i = 1,4
             write (*,'(4F9.3)') (transform(i,k),k=1,4)
          enddo
          write (*,*) 'saturation:'
          do i = 1,4
             write (*,'(4F9.3)') (saturation_criteria(i,k),k=1,4)
          enddo
       endif

    ENDIF

    CALL read_nml_close(substr, iou, modstr, l_print)


    status = 0 !NO ERROR
    
  END SUBROUTINE clamschem_read_nml
  
!*********************************************************************
subroutine allocate_chem_vars

  use messy_clamschem_asad_mod_clams,          only: jpctr, jpnr, &
                                     jpbk, jptk, jppj, jphkp1
  use messy_clamschem_asad_mod,                only: jpspb, jpspt, jpspj, jpsph

  use messy_clamschem_global,  only: ftrindex, fnames, &
                                     therm_flag, cindex, &
                                     slenb, slent, slenj, slenh

  implicit none

  allocate (ftrindex(jpctr))
  allocate (fnames(jpctr))
  allocate (therm_flag(jptk))
  allocate (slenb(jpspb),slent(jpspt),slenj(jpspj),slenh(jpsph))
  allocate (cindex(jppj))

end subroutine allocate_chem_vars

!*********************************************************************
subroutine deallocate_chem_vars

  use messy_clamschem_global,  only: ftrindex, fnames, &
                                     therm_flag, cindex, &
                                     slenb, slent, slenj, slenh

  implicit none

  deallocate (ftrindex)
  deallocate (fnames)
  deallocate (therm_flag)
  deallocate (slenb,slent,slenj,slenh)
  deallocate (cindex)

end subroutine deallocate_chem_vars

!*********************************************************************
subroutine clean_up

  use messy_clamschem_asad_mod_clams, only: failed, lhet
  use messy_clamschem_asad_mod,     only: converged, positive

  use messy_clamschem_globalhet, only: teold, prold, densnat_old, densice_old,  &
                                       denssat_old, cnat_old,cice_old,vliq_save, &
                                       astate,lstate, log_state
  use messy_clamschem_global,  only: missing_index, ftr, ftr_ini, zangle, &
                                     con, chindex, shindex, &
                                     densaero, aerh2so4, parth, wt, ar, &
                                     sedinucl, &
                                     lats,lons,theta,press,temps, sza, slt

  implicit none

  integer :: i

     deallocate (ftr)
!!!!! ???
     if (allocated(ftr_ini)) deallocate(ftr_ini)
     if (allocated(converged)) deallocate(converged)
     if (allocated(positive)) deallocate(positive)

     deallocate (zangle)

     deallocate (lats,lons,theta,press,temps)
     deallocate (sza,slt)
     deallocate (failed)

     deallocate(missing_index)
     
     if (lhet) then
        deallocate (con, shindex, chindex)
        deallocate (densaero, aerh2so4, parth, wt, ar)
        !deallocate (sedinucl)
        deallocate (teold, prold, densnat_old, densice_old,  &
                    denssat_old, cnat_old,cice_old,vliq_save, &
                    astate,lstate, log_state)
    endif

     

end subroutine clean_up

END MODULE messy_clamschem
