!        1         2         3         4         5         6         7         8
!2345678901234567890123456789012345678901234567890123456789012345678901234567890

MODULE messy_made3_box

   ! DESCRIPTION
   ! Interface layer for MADE3 box model
   !
   ! AUTHOR
   ! Axel Lauer, DLR Oberfpaffenhofen, Germany
   ! questions/suggestions: axel.lauer@dlr.de
   !
   ! LAST CHANGES
   ! 2008 by Valentina Aquila, DLR Oberpfaffenhofen (valentina.aquila@dlr.de)
   !      - introduction of akns, accs, sooti and sootj modes
   ! 2012, 2013 by Christopher Kaiser, DLR Oberpfaffenhofen
   !         (christopher.kaiser@dlr.de)
   !      - renamed from made to madein
   !      - adapted to MESSy2
   !      - modifications to enable initialization via namelist
   !      - introduction of several tests that can be switched on/off in
   !        namelist
   !      - renamed from madein to made3
   !      - cleaned up a little bit and adapted to include new modes
   !      - included Cl-/HCl
   !      - included HNO3 formation rate and BC emissions
   !
   ! TODO: Output diameters and species masses vs. time in nc files so that
   !       everything can be diagnosed offline.

   ! MESSy
   USE messy_main_constants_mem, ONLY: dp

   ! MADE3
   USE messy_made3, ONLY: modstr, dim1_cblk, dim2_cblk, nmod, nspec

   IMPLICIT NONE
   PRIVATE
   SAVE

   ! *** Constants ***

   INTEGER, PARAMETER, PUBLIC :: blksize = 1   ! array size
   INTEGER, PARAMETER, PUBLIC :: numcells = 1  ! number of cells (boxes)

   REAL(dp), PARAMETER, PUBLIC :: conmin = 1.e-30_dp  ! conc. lower lim. [ug/m3]
   REAL(dp), PARAMETER, PUBLIC :: dgmin  = 1.e-9_dp   ! lowest median diam. [m]
   REAL(dp), PARAMETER, PUBLIC :: one3   = 1._dp / 3._dp

   CHARACTER(LEN=10), PUBLIC   :: specname(nspec)     ! chemical species names
   CHARACTER(LEN=5),  PUBLIC   :: modname(nmod)       ! mode names

   ! *** Namelist-set parameters ***

   INTEGER,  PUBLIC :: timesteps              ! number of timesteps
   REAL(dp), PUBLIC :: tmst                   ! time step duration [s]

   REAL(dp), PUBLIC :: pressure(blksize)      ! pressure [Pa]
   REAL(dp), PUBLIC :: temperature(blksize)   ! temperature [K]
   REAL(dp), PUBLIC :: relhum(blksize)        ! relative humidity [0-1]
   REAL(dp), PUBLIC :: cloudcover(blksize)    ! fract. cloud cover [0-1]
   REAL(dp), PUBLIC :: rh_hist(nmod,blksize)  ! rel. hum. history
                                              ! (---> hysteresis)
   REAL(dp), PUBLIC :: so4rat(blksize)        ! H2SO4(g), rate of change
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: no3rat(blksize)        ! HNO3(g), rate of change
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: soa(blksize)           ! SOA (gas phase) "emissions"
                                              ! [ug/(m3s)]
   REAL(dp), PUBLIC :: em_bc_mass(3,blksize)  ! BC mass emission rate [ug m-3 s-1]
   REAL(dp), PUBLIC :: em_bc_num(3,blksize)   ! BC number emission rate [# m-3 s-1]

   LOGICAL, PUBLIC  :: ltest_mcon             ! test mass conservation
   LOGICAL, PUBLIC  :: ltest_ncon             ! test number consistency
   LOGICAL, PUBLIC  :: ltest_so4              ! test SO4 production
! op_ck_20120411+
   LOGICAL, PUBLIC  :: ltest_adapdt           ! test if adap. time step used
! op_ck_20120411-

   ! *** Global variables ***

   REAL(dp), PUBLIC :: tracer(dim1_cblk,dim2_cblk,blksize)! array passed to core
   REAL(dp), PUBLIC :: conv(dim1_cblk+3)                  ! conversion factors
   REAL(dp), PUBLIC :: dg(nmod,blksize)                   ! median diameters
   REAL(dp), PUBLIC :: h2so4, nh3, hno3, hcl

   ! *** Subroutines ***

   PUBLIC :: made3_init, made3_box_init, made3_log_init, made3_finalize


   CONTAINS

   ! --------------------------------------------------------------------------

   SUBROUTINE made3_init

      ! MADE3 parameters
      USE messy_made3,        ONLY: i_so4, i_nh4, i_no3, i_cl, i_ss, i_pom, &
                                    i_bc, i_bctag, i_du, i_h2o,             &
                                    akn, akns, sooti, acc, accs, sootj,     &
                                    cor, cors, sootc

      ! Subroutines
      USE messy_made3,        ONLY: made3_read_nml, made3_initialize_core
      USE messy_main_blather, ONLY: start_message, end_message

      IMPLICIT NONE
      INTRINSIC TRIM

      CHARACTER(LEN=*), PARAMETER :: substr = 'made3_init'
      INTEGER :: status ! error status

      CALL start_message(TRIM(modstr),'read namelist and initialize made3', &
                         substr)

      specname(i_so4)   = 'SO4'
      specname(i_nh4)   = 'NH4'
      specname(i_no3)   = 'NO3'
      specname(i_ss)    = 'SS'
      specname(i_cl)    = 'Cl'
      specname(i_pom)   = 'POM'
      specname(i_bc)    = 'BC'
      specname(i_bctag) = 'BCtag'
      specname(i_du)    = 'DU'
      specname(i_h2o)   = 'H2O'

      modname(akn)      = 'akn'
      modname(akns)     = 'akns'
      modname(sooti)    = 'sooti'
      modname(acc)      = 'acc'
      modname(accs)     = 'accs'
      modname(sootj)    = 'sootj'
      modname(cor)      = 'cor'
      modname(cors)     = 'cors'
      modname(sootc)    = 'sootc'

      ! *** ATTENTION ***
      ! Namelist MUST be read BEFORE calling 'made3_initialize_core'!!!

      ! read CTRL namelist
      CALL made3_read_nml(status, 99)
      IF (status /= 0) STOP

      WRITE(*,*) ' '
      WRITE(*,*) ' input data successfully read in'
      WRITE(*,*) ' '

      CALL made3_initialize_core
      CALL end_message(TRIM(modstr),'read namelist and initialize made3', &
           substr)

   END SUBROUTINE made3_init

   ! --------------------------------------------------------------------------

   SUBROUTINE made3_box_init   ! ATTENTION: Must be called AFTER made3_init.

     ! MESSy
     USE messy_main_constants_mem, ONLY: R_gas
     USE messy_main_blather, ONLY: start_message, end_message
     USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

     ! MADE3
     USE messy_made3, ONLY: akn, acc, cor, akns, accs, sooti, sootj, gas,     &
           i_num, i_mom3, nspec, i_h2so4, i_nh3, i_hno3, i_hcl,               &
           MW, i_mwaero, i_mwgas, rho, f6dpim9, sigma, es36

      IMPLICIT NONE
      INTRINSIC EXP, LOG, MAX, TRIM

      CHARACTER(LEN=*), PARAMETER :: substr = 'made3_box_init'
      INTEGER,          PARAMETER :: iou = 99     ! logical I/O unit
      LOGICAL                     :: lex          ! file exists ?
      INTEGER                     :: fstat        ! file status
      INTEGER                     :: i, j, k      ! loop indices

      NAMELIST /BOXINIT/ TIMESTEPS, TMST,                       &
           PRESSURE, TEMPERATURE, RELHUM, CLOUDCOVER,           &
           RH_HIST, SO4RAT, NO3RAT, SOA, EM_BC_MASS, EM_BC_NUM, &
           TRACER, LTEST_MCON, LTEST_NCON, LTEST_SO4, LTEST_ADAPDT

      ! Array initializations
      rh_hist    = -1.0_dp
      dg         = -1.0_dp
      no3rat     = 0.0_dp
      em_bc_mass = 0.0_dp
      em_bc_num  = 0.0_dp
      tracer     = 0.0_dp

      ! read BOXINIT namelist

      CALL start_message(TRIM(modstr),'read namelist and initialize ' &
           // 'parameters and variables     for made3 box model', substr)

      CALL read_nml_open(lex, substr, iou, 'BOXINIT', modstr)
      IF (.not.lex) STOP 1   ! <modstr>.nml does not exist

      READ(iou, NML=BOXINIT, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'BOXINIT', modstr)
      IF (fstat /= 0) STOP 1  ! error while reading namelist

      CALL read_nml_close(substr, iou, modstr)

      WRITE(*,*) ' '
      WRITE(*,*) ' initialization data successfully read in'
      WRITE(*,*) ' '

      CALL end_message(TRIM(modstr),'read namelist and initialize ' &
           // 'variables     for made3 box model', substr)

      ! set conversion factors
      ! (from 1/mol(air) to #/m3 and mol/mol(air) to ug/m3, respectively)

      conv(i_num)  = pressure(1) / (R_gas * temperature(1))

      DO i = 1, nspec
         conv(i) = 1.e6_dp * MW(i,i_mwaero) * conv(i_num)
      END DO

      do i = 1, blksize

         ! convert tracers to #/m3 and ug/m3, respectively,
         ! and calculate 3rd mom. conc. [m3 m-3]
         do j = 1, nmod
            do k = 1, nspec
               tracer(k,j,i)      = conv(k) * tracer(k,j,i)
               tracer(i_mom3,j,i) = tracer(i_mom3,j,i) &
                    + f6dpim9 * tracer(k,j,i) / rho(k)
            end do
            tracer(i_num,j,i)  = conv(i_num) * tracer(i_num,j,i)
            tracer(i_mom3,j,i) = MAX(conmin, tracer(i_mom3,j,i))
            ! Standard deviation fixed in all modes, so diagnose diameter
            ! from 3rd moment and number concentrations:
            dg(j,i) = MAX(dgmin, (tracer(i_mom3,j,i) / (tracer(i_num,j,i) &
                 * es36(j)))**one3)
         end do
         tracer(i_h2so4,gas,1) = 1.e6 * MW(i_h2so4,i_mwgas) * conv(i_num) &
              * tracer(i_h2so4,gas,1)
         h2so4 = tracer(i_h2so4,gas,1)
         tracer(i_nh3,gas,1)   = 1.e6 * MW(i_nh3,i_mwgas)   * conv(i_num) &
              * tracer(i_nh3,gas,1)
         nh3 = tracer(i_nh3,gas,1)
         tracer(i_hno3,gas,1)  = 1.e6 * MW(i_hno3,i_mwgas)  * conv(i_num) &
              * tracer(i_hno3,gas,1)
         hno3 = tracer(i_hno3,gas,1)
         tracer(i_hcl,gas,1)  = 1.e6 * MW(i_hcl,i_mwgas)  * conv(i_num) &
              * tracer(i_hcl,gas,1)
         hcl = tracer(i_hcl,gas,1)

      end do

   END SUBROUTINE made3_box_init

   ! --------------------------------------------------------------------------

   SUBROUTINE made3_finalize

      !**** *SUBROUTINE* *MADE3_FINALIZE*  Finalize MADE3 after last call.
      !
      !         A. Lauer   DLR Oberpfaffenhofen  11/2004
      !
      !     PURPOSE.
      !     --------
      !
      !         1. Finalize MADE3 after last call.
      !         2. Clean up.
      !
      !**   INTERFACE.
      !     ----------
      !
      !       *MADE3_FINALIZE* IS CALLED FROM PROGRAM *MADE3*
      !
      !       *PARAMETERS*: none
      !
      !
      !     METHODS.
      !     --------
      !
      !         -
      !
      !     EXTERNALS.
      !     ----------
      !
      !         None.
      !
      !     REFERENCE.
      !     ----------
      !
      !         None.
      !
      !     ------------------------------------------------------------

      IMPLICIT NONE

      !*        1. NOTHING TO DO YET
      !            ------- -- -- ---      

   END SUBROUTINE made3_finalize

   ! --------------------------------------------------------------------------

   SUBROUTINE made3_log_init

      !**** *SUBROUTINE* *MADE3_LOG_INIT*  Write initial setup to file
      !
      !         C. Kaiser   DLR Oberpfaffenhofen  04/2012
      !
      !     PURPOSE.
      !     --------
      !
      !         1. Integrate initial tracer concentrations over size ranges
      !            (Aitken, accumulation, and coarse size range).
      !         2. Write these concentrations [ug/m3] to file.
      !
      !**   INTERFACE.
      !     ----------
      !
      !       *MADE3_LOG_INIT* IS CALLED FROM PROGRAM *MADE3*
      !
      !       *PARAMETERS*: none
      !
      !
      !     METHODS.
      !     --------
      !
      !         -
      !
      !     EXTERNALS.
      !     ----------
      !
      !         None.
      !
      !     REFERENCE.
      !     ----------
      !
      !         None.
      !
      !     ------------------------------------------------------------

      USE messy_made3, ONLY: i_num, nspec,                                    &
                             akn, akns, sooti, acc, accs, sootj, cor, cors,   &
                             sootc

      IMPLICIT NONE

      INTEGER :: initout = 105   ! unit for output file for initialization
                                 ! values
      INTEGER :: k               ! loop index

      open (initout, action='WRITE', file='made3_ini_agg.dat', &
           form='FORMATTED')

      write (initout, '(A)') '# This file contains the initialization '  &
           // 'values for number and mass concentrations'
      write (initout, '(A)') '# in the aggregated Aitken, accumulation ' &
           // 'and coarse modes as well as the'
      write (initout, '(A)') '# mass contributions of the individual '   &
           // 'chemical components. It also lists the'
      write (initout, '(A,/)') '# initial values of the meteorological ' &
           // 'parameters.'

      write (initout,'(A6,3(2X,A12))') ' ', 'akn_tot', 'acc_tot', 'cor_tot'

      write (initout,'(A6,3(2X,es12.6))') 'N',      &
           SUM(tracer(i_num,(/akn,akns,sooti/),1)), &
           SUM(tracer(i_num,(/acc,accs,sootj/),1)), &
           SUM(tracer(i_num,(/cor,cors,sootc/),1))
      write (initout,'(A6,3(2X,es12.6))') 'm',        &
           SUM(tracer(1:nspec,(/akn,akns,sooti/),1)), &
           SUM(tracer(1:nspec,(/acc,accs,sootj/),1)), &
           SUM(tracer(1:nspec,(/cor,cors,sootc/),1))
      DO k = 1, nspec
         write (initout,'(A6,3(2X,es12.6))') specname(k), &
              SUM(tracer(k,(/akn,akns,sooti/),1)),        &
              SUM(tracer(k,(/acc,accs,sootj/),1)),        &
              SUM(tracer(k,(/cor,cors,sootc/),1))
      END DO

      write (initout,'(/,A2,2X,es12.6)') 'p', pressure(1)
      write (initout,'(A2,2X,es12.6)')   'T', temperature(1)
      write (initout,'(A2,2X,es12.6)')  'RH', relhum(1)

      close (initout)
      
    END SUBROUTINE made3_log_init
    
   ! --------------------------------------------------------------------------

END MODULE messy_made3_box

!*****************************************************************************

PROGRAM made3

   USE messy_main_constants_mem, ONLY: dp

   USE messy_made3, ONLY: i_so4, i_nh4, i_no3, i_cl, i_bc, i_bctag, i_h2o, &
        i_num, i_mom3, i_h2so4, i_nh3, i_hno3, i_hcl, nspec,               &
        akn, akns, sooti, acc, accs, sootj, cor, cors, sootc, gas, nmod,   &
        MW, i_mwaero, i_mwgas, sigma, es36, rho, f6dpim9

   USE messy_made3_box, ONLY: blksize, numcells, timesteps, tmst,  &
        pressure, temperature, relhum, cloudcover,                 &
        rh_hist, so4rat, no3rat, soa, em_bc_mass, em_bc_num,       &
        tracer, dg, one3, dgmin, ltest_mcon, ltest_ncon, ltest_so4,&
        ltest_adapdt, specname, modname, h2so4, nh3, hno3, hcl

   ! SUBROUTINES

   USE messy_made3_box, ONLY: made3_init, made3_box_init,    &
                              made3_log_init, made3_finalize
   USE messy_made3,     ONLY: made3_main

   IMPLICIT NONE
   INTRINSIC :: MAX, SUM

   integer :: status                   ! error status
   integer :: i, j                     ! loop indices
   integer :: outunit = 104            ! unit for output file for data of
                                       ! aggregated modes
   integer :: testout = 106            ! unit for optional output of tests
   integer :: o_species = 107          ! unit for output of species mass
                                       ! time series
   integer :: o_dry = 108              ! unit for output of dry diameters/masses
   integer :: o_spec_mod(nmod)         ! units for output of species mass time
                                       ! series per mode

   ! modal diameters (diagnostic output for box-version)
   REAL(dp) :: DGDRY(nmod,BLKSIZE) ! modal geometric mean diameters (dry) [m]
   REAL(dp) :: DENS(nmod,BLKSIZE)  ! average modal densities [kg m-3]

   ! for aggregated output
   REAL(dp) :: sumakn0, sumacc0,   & ! sum number conc. over all Ait./acc. modes
               sumcor0,            & ! [#/m3]
               sumakn3, sumacc3,   & ! sum 3rd mom. conc. [m3/m3]
               sumcor3,            &
               dgnuctot, dgacctot, & ! median diam. of aggregated modes [m]
               dgcortot
   INTEGER  :: k                            ! loop index
   REAL(dp) :: oldmass, newmass             ! for mass conservation test
   REAL(dp) :: oldnum, newnum               ! for number consistency test
   REAL(dp) :: oldso4, newso4               ! for SO4 production test
! op_ck_20120411+
   REAL(dp) :: oldtmst                      ! for adaptive time step test
! op_ck_20120411-

   ! Array initializations
   DGDRY = -1.0_dp
   DENS  = -1.0_dp

   ! ----------------------------------------------------------------------
   !     Set initial values - box version only.
   ! ----------------------------------------------------------------------

   call made3_init       ! Initialization of MADE3 before 1st call
   call made3_box_init   ! Read initial parameters and variables for box
                         ! model run from namelist file

   ! ----------------------------------------------------------------------

   call made3_log_init

   ! Output of aggregated modes for comparison with MADE
   open (outunit, action='WRITE', file='made3_agg_modes.dat', form='FORMATTED')
   write (outunit, '(A)') '# This file contains the time evolution of the ' &
        // 'number concentrations and median'
   write (outunit, '(A)') '# diameters of the aggregated Aitken, '          &
        // 'accumulation, and coarse modes.'
   write (outunit, '(A)') '# The aggregated Aitken mode (akn_tot) '         &
        // 'comprises modes akn, akns, and sooti,'
   write (outunit, '(A)') '# the aggregated accumulation mode (acc_tot) '   &
        // 'modes acc, accs, and sootj,'
   write (outunit, '(A)') '# and the aggregated coarse mode (cor_tot) '     &
        // 'modes cor, cors, and sootc.'
   write (outunit, '(A)') '# Particle number concentrations are summed '    &
        // 'over the respective modes, median'
   write (outunit, '(A)') '# diameters are calculated from the sum of the ' &
        // 'respective 0th and 3rd moments.'
   write (outunit, *) ''
   write (outunit,'(A5,2(3(2X,A12)))') 'step',    &
        'num, akn_tot', 'num, acc_tot', 'num, cor_tot', &
        'dg, akn_tot', 'dg, acc_tot', 'dg, cor_tot'

   ! print some diagnostics of aerosol components...
   write (*,*)
   write (*,'(A10,9(2X,A9))') 'component','akn','akns','sooti','acc','accs',  &
        'sootj','cor','cors','sootc'
   DO k = 1, nspec
      write (*,'(A10,9(2X,es9.3))') specname(k), &
           tracer(k,akn,1), tracer(k,akns,1), tracer(k,sooti,1),   &
           tracer(k,acc,1), tracer(k,accs,1), tracer(k,sootj,1),   &
           tracer(k,cor,1), tracer(k,cors,1), tracer(k,sootc,1)
   END DO

   ! ...and gas-phase species
   write (*,*)
   write (*,'(A10,2X,es9.3)') 'NH3(g)',   tracer(i_nh3,gas,1)
   write (*,'(A10,2X,es9.3)') 'H2SO4(g)', tracer(i_h2so4,gas,1)
   write (*,'(A10,2X,es9.3)') 'HNO3(g)',  tracer(i_hno3,gas,1)
   write (*,'(A10,2X,es9.3)') 'HCl(g)',   tracer(i_hcl,gas,1)

   ! Initialize time series output file
   open (o_species, action='WRITE', file='made3_species.dat', form='FORMATTED')
   write (o_species, '(A)') '# This file contains the time [s] evolution of the ' &
        // 'individual species mass concentrations [ug m-3]'
   write (o_species, *) ''
   write (o_species, '(A6)', advance='no') 'time'
   do k = 1, nspec
      write (o_species, '(1X,A9)', advance='no') specname(k)
   end do
   write (o_species, *) ''
   write (o_species, '(i6)', advance='no') 0
   do k = 1, nspec
      write (o_species, '(1X,es9.3)', advance='no') sum(tracer(k,1:nmod,1))
   end do
   write (o_species, *) ''

   ! Initialize dry diameter/mass output file
   open (o_dry, action='WRITE', file='made3_dry.dat', form='FORMATTED')
   write (o_dry, '(A)') '# This file contains the time [s] evolution of the ' &
        // 'dry diameters [m] and modal dry mass concentrations [ug m-3]'
   write (o_dry, *) ''
   write (o_dry, '(A6)', advance='no') 'time'
   do j = 1, nmod
      write (o_dry, '(1X,A9)', advance='no') 'dg, '//modname(j)
   end do
   do j = 1, nmod
      write (o_dry, '(1X,A9)', advance='no') 'm, '//modname(j)
   end do
   write (o_dry, *) ''
   write (o_dry, '(i6)', advance='no') 0
   do j = 1, nmod
      write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
           ((tracer(i_mom3,j,1) - tracer(i_h2o,j,1) * f6dpim9 / rho(i_h2o)) &
           / (tracer(i_num,j,1) * es36(j)))**one3)
   end do
   do j = 1, nmod
      write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1:nspec,j,1)) &
           - tracer(i_h2o,j,1)
   end do
   write (o_dry, *) ''

   ! Initialize time series files for modal species mass concentrations
   do j = 1, nmod
      o_spec_mod(j) = o_dry + j
      open (o_spec_mod(j), action='WRITE', &
           file='made3_species_'//trim(modname(j))//'.dat', form='FORMATTED')
      write (o_spec_mod(j), '(A)') '# This file contains the time' &
           // ' [s] evolution of the wet volume' &
           // ' median diameter [m] and corresponding species' &
           // ' mass concentrations [ug m-3] in mode '//modname(j)
      write (o_spec_mod(j), *) ''
      write (o_spec_mod(j), '(A6)', advance='no') 'time'
      write (o_spec_mod(j), '(1X,A9)', advance='no') 'dgv'
      do k = 1, nspec
         if (k .eq. i_bctag) cycle
         write (o_spec_mod(j), '(1X,A9)', advance='no') 'm, '//specname(k)
      end do
      write (o_spec_mod(j), *) ''
      write (o_spec_mod(j), '(i6)', advance='no') 0
      write (o_spec_mod(j), '(1X,es9.3)', advance='no') &
           dg(j,1) * exp(3._dp * log(sigma(j)) * log(sigma(j)))
      do k = 1, nspec
         if (k .eq. i_bc) then 
            write (o_spec_mod(j), '(1X,es9.3)', advance='no') tracer(k,j,1) &
                 + tracer(k+1,j,1)
         elseif (k .eq. i_bctag) then
            cycle
         else
            write (o_spec_mod(j), '(1X,es9.3)', advance='no') tracer(k,j,1)
         end if
      end do
      write (o_spec_mod(j), *) ''
   end do

   write (*,*)
   write (*,'(A5,18(1X,A10))') 'step', &
        'num, akn', 'num, akns', 'num, sooti', &
        'num, acc', 'num, accs', 'num, sootj', &
        'num, cor', 'num, cors', 'num, sootc', &
        'dg, akn', 'dg, akns', 'dg, sooti', &
        'dg, acc', 'dg, accs', 'dg, sootj', &
        'dg, cor', 'dg, cors', 'dg, sootc'

   ! Write initial number conc.s and diameters to standard output
   write (*,'(i5,18(1X,es9.3))') 0,                                       &
        tracer(i_num,akn,1), tracer(i_num,akns,1), tracer(i_num,sooti,1), &
        tracer(i_num,acc,1), tracer(i_num,accs,1), tracer(i_num,sootj,1), &
        tracer(i_num,cor,1), tracer(i_num,cors,1), tracer(i_num,sootc,1), &
        dg(akn,1), dg(akns,1), dg(sooti,1),                               &
        dg(acc,1), dg(accs,1), dg(sootj,1),                               &
        dg(cor,1), dg(cors,1), dg(sootc,1)

   ! Write aggregated initial values to outunit
   sumakn0 = SUM(tracer(i_num,(/akn,akns,sooti/),1))
   sumacc0 = SUM(tracer(i_num,(/acc,accs,sootj/),1))
   sumcor0 = SUM(tracer(i_num,(/cor,cors,sootc/),1))
   sumakn3 = SUM(tracer(i_mom3,(/akn,akns,sooti/),1))
   sumacc3 = SUM(tracer(i_mom3,(/acc,accs,sootj/),1))
   sumcor3 = SUM(tracer(i_mom3,(/cor,cors,sootc/),1))
   ! NOTE: The following lines are correct only if all Aitken modes have the
   ! same width, all accumulation modes have the same width, and all coarse
   ! modes have the same width.
   dgnuctot = MAX(dgmin, (sumakn3 / (sumakn0 * es36(akn)))**one3)
   dgacctot = MAX(dgmin, (sumacc3 / (sumacc0 * es36(acc)))**one3)
   dgcortot = MAX(dgmin, (sumcor3 / (sumcor0 * es36(cor)))**one3)
   write (outunit,'(i5,2(2(2X,es12.6),2X,es12.6))') 0, &
        sumakn0, sumacc0, sumcor0, &
        dgnuctot, dgacctot, dgcortot

   ! File for test output
   if (ltest_mcon .or. ltest_ncon .or. ltest_so4 .or. ltest_adapdt) then
      open(testout, action='WRITE', file='made3_test.log', form='FORMATTED')
   end if

   ! *** For mass conservation test
   if (ltest_mcon) then
      oldmass = 0.
      do j = 1, nmod
         do k = 2, 9
            oldmass = oldmass + tracer(k,j,1)
         end do
      end do
      oldmass = oldmass &
           + tracer(i_nh3,gas,1)  * MW(i_nh4,i_mwaero) / MW(i_nh3,i_mwgas)  &
           + tracer(i_hno3,gas,1) * MW(i_no3,i_mwaero) / MW(i_hno3,i_mwgas) &
           + tracer(i_hcl,gas,1)  * MW(i_cl,i_mwaero)  / MW(i_hcl,i_mwgas)
      write (testout,'(A40,2X,es13.6)') 'total initial dry mass conc. ' &
           // '(incl. gases)', oldmass
   end if

   ! *** For number consistency test
   if (ltest_ncon) then
      oldnum = SUM(tracer(i_num,1:nmod,1))
      write (testout,'(A40,2X,es13.6)') 'total initial number conc.', &
           oldnum
   end if

   ! *** For SO4 production test
   if (ltest_so4) then
      oldso4 = tracer(i_h2so4,gas,1) * MW(i_so4,i_mwaero) / MW(i_h2so4,i_mwgas)&
           + SUM(tracer(i_so4,1:nmod,1))
   end if

! op_ck_20120411+
   ! *** For adaptive time step test
   oldtmst = tmst
   if (ltest_adapdt) then
      write (testout,'(A40,2X,es13.6,A)') 'time step initially set to', &
           tmst, 's'
   end if
! op_ck_20120411-

   do i = 1, timesteps
!!$      tracer(i_h2so4,gas,1) = h2so4 + so4rat(1) * tmst
!!$      tracer(i_nh4,gas,1)   = nh3
!!$      tracer(i_hno3,gas,1)  = hno3
!!$      tracer(i_hcl,gas,1)   = hcl
      ! *** Emissions ***
      tracer(i_h2so4,gas,1) = tracer(i_h2so4,gas,1) + so4rat(1) * tmst
      tracer(i_hno3,gas,1)  = tracer(i_hno3,gas,1)  + no3rat(1) * tmst
      tracer(i_bc,sooti,1)  = tracer(i_bc,sooti,1)  + em_bc_mass(1,1) * tmst
      tracer(i_bc,sootj,1)  = tracer(i_bc,sootj,1)  + em_bc_mass(2,1) * tmst
      tracer(i_bc,sootc,1)  = tracer(i_bc,sootc,1)  + em_bc_mass(3,1) * tmst
      tracer(i_num,sooti,1) = tracer(i_num,sooti,1) + em_bc_num(1,1)  * tmst
      tracer(i_num,sootj,1) = tracer(i_num,sootj,1) + em_bc_num(2,1)  * tmst
      tracer(i_num,sootc,1) = tracer(i_num,sootc,1) + em_bc_num(3,1)  * tmst
      call made3_main(status, BLKSIZE, NUMCELLS, PRESSURE, TEMPERATURE, &
                      RELHUM, TMST, SO4RAT, SOA, CLOUDCOVER, TRACER,     &
                      RH_HIST, DG, DGDRY, DENS)
      if (status /= 0) stop

      write (*,'(i5,18(1X,es9.3))') i,                                       &
           tracer(i_num,akn,1), tracer(i_num,akns,1), tracer(i_num,sooti,1), &
           tracer(i_num,acc,1), tracer(i_num,accs,1), tracer(i_num,sootj,1), &
           tracer(i_num,cor,1), tracer(i_num,cors,1), tracer(i_num,sootc,1), &
           dg(akn,1), dg(akns,1), dg(sooti,1),                               &
           dg(acc,1), dg(accs,1), dg(sootj,1),                               &
           dg(cor,1), dg(cors,1), dg(sootc,1)

      ! Write time series of species mass concentrations
      write (o_species, '(i6)', advance='no') int(i * oldtmst)
      do k = 1, nspec
         write (o_species, '(1X,es9.3)', advance='no') sum(tracer(k,1:nmod,1))
      end do
      write (o_species, *) ''

      ! Write time series of dry diameters/masses
      write (o_dry, '(i6)', advance='no') int(i * oldtmst)
      do j = 1, nmod
         write (o_dry, '(1X,es9.3)', advance='no') MAX(dgmin, &
              ((tracer(i_mom3,j,1) - tracer(i_h2o,j,1) * f6dpim9 / rho(i_h2o)) &
              / (tracer(i_num,j,1) * es36(j)))**one3)
      end do
      do j = 1, nmod
         write (o_dry, '(1X,es9.3)', advance='no') SUM(tracer(1:nspec,j,1)) &
              - tracer(i_h2o,j,1)
      end do
      write (o_dry, *) ''

      ! Write modal time series of species mass concentrations
      do j = 1, nmod
         write (o_spec_mod(j), '(i6)', advance='no') int(i * oldtmst)
         write (o_spec_mod(j), '(1X,es9.3)', advance='no') &
              dg(j,1) * exp(3._dp * log(sigma(j)) * log(sigma(j)))
         do k = 1, nspec
            if (k .eq. i_bc) then 
               write (o_spec_mod(j), '(1X,es9.3)', advance='no') tracer(k,j,1) &
                    + tracer(k+1,j,1)
            elseif (k .eq. i_bctag) then
               cycle
            else
               write (o_spec_mod(j), '(1X,es9.3)', advance='no') tracer(k,j,1)
            end if
         end do
         write (o_spec_mod(j), *) ''
      end do

      ! For aggregated output
      sumakn0 = SUM(tracer(i_num,(/akn,akns,sooti/),1))
      sumacc0 = SUM(tracer(i_num,(/acc,accs,sootj/),1))
      sumcor0 = SUM(tracer(i_num,(/cor,cors,sootc/),1))
      sumakn3 = SUM(tracer(i_mom3,(/akn,akns,sooti/),1))
      sumacc3 = SUM(tracer(i_mom3,(/acc,accs,sootj/),1))
      sumcor3 = SUM(tracer(i_mom3,(/cor,cors,sootc/),1))
      ! NOTE: The following lines are correct only if all Aitken modes have the
      ! same width, all accumulation modes have the same width, and all coarse
      ! modes have the same width.
      dgnuctot = MAX(dgmin, (sumakn3 / (sumakn0 * es36(akn)))**one3)
      dgacctot = MAX(dgmin, (sumacc3 / (sumacc0 * es36(acc)))**one3)
      dgcortot = MAX(dgmin, (sumcor3 / (sumcor0 * es36(cor)))**one3)
      write (outunit,'(i5,2(2(2X,es12.6),2X,es12.6))') i, &
           sumakn0, sumacc0, sumcor0, &
           dgnuctot, dgacctot, dgcortot

      ! *** Mass conservation test
      if (ltest_mcon) then
         newmass = 0.
         do j = 1, nmod
            do k = 2, 9
               newmass = newmass + tracer(k,j,1)
            end do
         end do
         newmass = newmass &
              + tracer(i_nh3,gas,1)  * MW(i_nh4,i_mwaero) / MW(i_nh3,i_mwgas)  &
              + tracer(i_hno3,gas,1) * MW(i_no3,i_mwaero) / MW(i_hno3,i_mwgas) &
              + tracer(i_hcl,gas,1)  * MW(i_cl,i_mwaero)  / MW(i_hcl,i_mwgas)
         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative mass conc. ' &
              // 'change step', i, (newmass - oldmass) * 100. / oldmass, '%'
         write (testout,'(A36,1X,i3,2X,es13.6)') 'total mass conc. step', &
              i, newmass
         oldmass = newmass
      end if

      ! *** Number consistency test
      if (ltest_ncon) then
         newnum = SUM(tracer(i_num,1:nmod,1))
         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative number conc. ' &
              // 'change step', i, (newnum - oldnum) * 100. / oldnum, '%'
         write (testout,'(A36,1X,i3,2X,es13.6)') 'total number conc. step', &
              i, newnum
         oldnum = newnum
      end if

      ! *** SO4 production test
      if (ltest_so4) then
         newso4 = tracer(i_h2so4,gas,1)*MW(i_so4,i_mwaero)/MW(i_h2so4,i_mwgas) &
              + SUM(tracer(i_so4,1:nmod,1))
         write (testout,'(A36,1X,i3,2X,f13.9,A)') 'relative SO4 prod. ' &
              // 'deviation step', i, ((newso4 - oldso4) * 100._dp &
              / (so4rat * oldtmst * MW(i_so4,i_mwaero) / MW(i_h2so4,i_mwgas))) &
              - 100._dp, '%'
         write (testout,'(A36,1X,i3,2X,es13.6)') 'absolute SO4 prod. step', &
              i, (newso4 - oldso4) / oldtmst
         oldso4 = newso4
      end if

      ! *** Adaptive time step test
      if (ltest_adapdt) then
         write(testout,'(A36,1X,i3,2X,es13.6,A)') 'actual step duration ' &
              // 'during step', i, tmst, 's'
      end if
      tmst = oldtmst

   end do

   ! print some diagnostics of aerosol components...

   write (*,*)
   write (*,'(A10,9(2X,A9))') 'component','akn','akns','sooti','acc','accs',  &
        'sootj','cor','cors','sootc'
   DO k = 1, nspec
      write (*,'(A10,9(2X,es9.3))') specname(k), &
           tracer(k,akn,1), tracer(k,akns,1), tracer(k,sooti,1),   &
           tracer(k,acc,1), tracer(k,accs,1), tracer(k,sootj,1),   &
           tracer(k,cor,1), tracer(k,cors,1), tracer(k,sootc,1)
   END DO

   ! ...and gas-phase species

   write (*,*)
   write (*,'(A10,2X,es9.3)') 'NH3(g)',   tracer(i_nh3,gas,1)
   write (*,'(A10,2X,es9.3)') 'H2SO4(g)', tracer(i_h2so4,gas,1)
   write (*,'(A10,2X,es9.3)') 'HNO3(g)',  tracer(i_hno3,gas,1)
   write (*,'(A10,2X,es9.3)') 'HCl(g)',   tracer(i_hcl,gas,1)

   call made3_finalize   ! Finalization of MADE after last call

   do j = 1, nmod
      close (o_spec_mod(j)) ! Files for modal mass conc.s' time series
   end do
   close (o_dry)     ! File for dry diameter/mass time series
   close (o_species) ! File for mass concentrations' time series
   close (outunit)   ! File for aggregated output
   if (ltest_mcon .or. ltest_ncon .or. ltest_so4 .or. ltest_adapdt) &
        close(testout)   ! File for test output

   write (*,*)
   write (*,*) 'BoxMADE3 run finished'

END PROGRAM made3

!*****************************************************************************
