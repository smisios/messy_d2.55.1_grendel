 
!**********************************************************************
MODULE messy_clamschem_si
!**********************************************************************
!  Submodel interface for clamschem 
!**********************************************************************
! SMCL
  USE messy_clams_global,     ONLY: PREC, DP, species_type
  USE messy_clamschem
  USE messy_clamschem_global, ONLY: rc_type, rate_type
  USE messy_main_timer_event, ONLY: time_event, io_time_event, event_is_active 

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  PUBLIC :: clamschem_initialize
  PUBLIC :: clamschem_init_memory
  PUBLIC :: clamschem_init_coupling
  PUBLIC :: clamschem_global_end
  PUBLIC :: clamschem_free_memory

  ! MODULE VARIABLES

  REAL(PREC), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV        => NULL()
  !REAL(DP),   DIMENSION(:), POINTER :: JULSEC     => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: TEMP       => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: PRESS      => NULL()
  REAL(DP),   POINTER               :: JULTIME

  REAL(PREC), DIMENSION(:), POINTER :: LAT_OLD    => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_OLD    => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_OLD    => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: TEMP_OLD   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: PRESS_OLD  => NULL()
  REAL(DP),   POINTER               :: JULTIME_OLD


  TYPE(rc_type), DIMENSION(:), POINTER :: BRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: TRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: JRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: HRATES => NULL()

  TYPE(rc_type), DIMENSION(:), POINTER :: BCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: TCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: JCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: HCONST => NULL()

  TYPE(species_type), DIMENSION(:), POINTER :: CHEMSPECARR => NULL()

  TYPE(species_type), DIMENSION(:), POINTER :: HETPAR => NULL()

! op_pj_20160606+
!!$  TYPE(time_event), PUBLIC :: chemevent
!!$  TYPE(io_time_event):: io_chemevent
  TYPE(time_event), PUBLIC, SAVE :: chemevent
  TYPE(io_time_event), SAVE :: io_chemevent
! op_pj_20160606-

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE clamschem_initialize

    ! BMIL
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_timer_bi,      ONLY: timer_event_init

    ! SMCL
    USE messy_clams_global,          ONLY: nspec, nchemspec, maxspec, SPECARR, &
                                           rank, dates30, ldiagout, specnames
! op_pj_20170110+
!!$ USE messy_clamschem_global,      ONLY: timestep_chem, ncdt, asad_gfirst, &
    USE messy_clamschem_global,      ONLY: timestep_chem, ncdt, &
                                           ip_messy, &
                                           hetvar, nhetspec, nhetpar
    USE messy_clams_global,          ONLY: asad_gfirst
! op_pj_20170110-
    USE messy_clamschem_defs_mod,    ONLY: chch_defs
    USE messy_clamschem_data,        ONLY: clams_chem_init
    USE messy_clamschem_data_hetero, ONLY: clams_chem_init_hetero
    USE messy_clamschem,             ONLY: allocate_chem_vars
    USE messy_clams_tools_utils,     ONLY: uppercase

    USE messy_clamschem_asad_mod,    ONLY: asad_mod_init, &
                                           speci, advt, family, madvtr, majors, &
                                           moffam, ctype, &
                                           nit0, nitfg, nitnr, nrsteps, &
                                           ldepd, ldepw, dtime, cdt
    USE messy_clamschem_asad_mod_clams,  ONLY: jpspec, jpctr, lhet
    USE messy_clamschem_asad_cinit,  ONLY: asad_cinit

    USE messy_main_timer,            ONLY: delta_time
    USE messy_main_tools,            ONLY: find_next_free_unit

    USE messy_cmn_photol_mem 


    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamschem_initialize'

    INTEGER :: status, iou, ios, n, nodd, i
    !CHARACTER(2)  :: ctype
    CHARACTER(20) :: specname

    IF (rank==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
    ENDIF
    
    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    CALL clamschem_read_nml(status, iou)
    IF (status /= 0) CALL error_bi('Error in clamschem_read_nml !',substr)

    ! Define CHEM event / check chemistry timestep
    if (timestep_chem < delta_time) then
       if (rank==0) then
          write (*,*)
          write (*,*) 'timestep_chem < delta_time'
          write (*,*) 'timestep_chem:', timestep_chem
          write (*,*) 'delta_time:', delta_time
          write (*,*) 'timestep_chem is set to:', delta_time
          write (*,*)
       endif
       timestep_chem = delta_time
    elseif (mod(timestep_chem,int(delta_time)) /= 0) then
       call error_bi ("Wrong chemistry timestep !!!",substr)
    endif
    if (timestep_chem > delta_time) then
       io_chemevent%counter = timestep_chem    
       io_chemevent%unit = 'seconds'
       io_chemevent%adjustment = 'exact'
       io_chemevent%offset = -delta_time
       CALL timer_event_init (chemevent, io_chemevent, 'CHEM_Event', 'present')
    endif
    if (mod(timestep_chem,ncdt) /= 0) then
       call error_bi &
            ("Wrong timesteps: chemistry timestep must be equal to or multiple of internal chemistry timestep (NCDT) !!!",substr)
    endif


    ! For calculation using  30-days-month dates30 has to be set to "true"
    if (dates30) then
       write (*,*)
       write (*,*) 'ATTENTION: 30-days-months'
       write (*,*)
    endif

!    IF (emit) write(*,*) 'emissions of NO through cosmic rays calculated'

    ! Get chemistry input
    IF (rank==0 .and. ldiagout) write (*,*) 'call clams_chem_init'
    call clams_chem_init
    
    ! Allocate and initialize ASAD arrays and variables
    IF (rank==0 .and. ldiagout) write (*,*) 'call asad_mod_init'
    call asad_mod_init
    
!!!!! Nach asad_mod_init belegen:
    ldepd = .FALSE. ! in asad_mod.f90 deklariert
    ldepw = .FALSE. ! in asad_mod.f90 deklariert
    nit0    = 20    ! in asad_mod: mit 20 initialisiert
    nitfg   = 20    ! in asad_mod_init: mit 10 initialisiert
    nitnr   = 20    ! in asad_mod_init: mit 10 initialisiert
    nrsteps = 50    ! in asad_mod_init: mit 45 initialisiert
    cdt     = ncdt  ! in asad_mod.f90 deklariert

    dtime = timestep_chem   ! dtime in asad_mod_init initialisiert !

    asad_gfirst = .true.

    ! Allocate arrays 
    call allocate_chem_vars

!!!!!
     ! --------------------------------------
     !    Initialise the Chemistry
     ! --------------------------------------
     !
     ! Call ASAD routine CINIT
     !
    if (rank==0 .and. ldiagout) write (*,*) 'call asad_cinit !!!'
     call asad_cinit (1)

     if (rank==0) then 
        write (6, *) 'speci= ', (speci(i),i=1,jpspec) 
        write (6, *) 'advt= ', (advt(i),i=1,jpctr) 
        write (6, *) 'family= ', (family(i),i=1,jpspec) 
        write (6, *) 'madvtr= ', (madvtr(i),i=1,jpspec) 
        write (6, *) 'majors= ', (majors(i),i=1,jpctr) 
        write (6, *) 'moffam= ', (moffam(i),i=1,jpspec) 
        write (6, *) 'ctype= ', (ctype(i),i=1,jpspec) 
     endif


    ! Read chemical species from structure chch_defs
    nchemspec = jpspec
    if (nchemspec > maxspec) &
         call error_bi (substr,"To many species: increase MAXSPEC !!!")
    DO i = 1, jpspec
       specnames(i)        = chch_defs(i)%speci
       SPECARR(i)%name     = chch_defs(i)%speci
       SPECARR(i)%ctype    = chch_defs(i)%ctype
       SPECARR(i)%longname = chch_defs(i)%speci
       SPECARR(i)%units    = 'm^3/m^3'
    ENDDO

    ! for heterogeneous chemistry: add species
    if (lhet) then

       IF (rank==0 .and. ldiagout) write (*,*) 'call clams_chem_init_hetero'
       call clams_chem_init_hetero
       DO i = 1, nhetspec
          specnames(nchemspec+i)        = hetvar(i)%name
          SPECARR(nchemspec+i)%name     = hetvar(i)%name
          SPECARR(nchemspec+i)%ctype    = ' '
          SPECARR(nchemspec+i)%longname = hetvar(i)%longname
          SPECARR(nchemspec+i)%units    = hetvar(i)%units
       ENDDO
       nchemspec = nchemspec + nhetspec

       allocate (HETPAR(nhetpar))
       do i = 1, nhetpar
          HETPAR(i)%name     = hetvar(nhetspec+i)%name
          HETPAR(i)%longname = hetvar(nhetspec+i)%longname
          HETPAR(i)%units    = hetvar(nhetspec+i)%units
       enddo

    endif

    nspec = nchemspec
    if (rank==0) then
       write (*,*) 'nchemspec=',nchemspec
       do i = 1, nchemspec
          write(*,*) "i,speciesname(i),type ",i,SPECARR(i)%name,SPECARR(i)%ctype
       enddo
       if (lhet) then
          write (*,*) 'nhetpar=',nhetpar
          do i = 1, nhetpar
             write(*,*) "i,hetpar(i) ",i,HETPAR(i)%name
          enddo
       endif
    endif
    CLOSE(iou)

    ALLOCATE (CHEMSPECARR(nchemspec))


  ! corresponding MESSy photolysis indices from messy_cmn_photol_mem
  !
  ! *** must be in the same order than the above list ***
  ip_messy = (/ip_O2, ip_O3P, ip_O1D, ip_H2O2, ip_Cl2, ip_Cl2O2, &
       ip_HOCl, ip_ClNO2, ip_ClNO3, ip_HNO3, ip_NO2, ip_N2O5, &
       ip_HO2NO2, ip_NO2O, ip_NOO2, ip_CHOH, ip_COH2, ip_CH3OOH, &
       ip_BrONO2, ip_BrCl, ip_OClO, ip_H2O, ip_HCl, ip_NO, &
       ip_N2O, ip_CH3OCl, ip_HOBr, ip_Br2, ip_MEO2NO2, ip_BrO, &
       ip_ClONO2, ip_BrNO3, ip_OHNO3, ip_CFCl3, ip_CF2CL2, ip_CHF2Cl, &
       ip_F113, ip_CH3Cl, ip_CCl4, ip_CH3Br, ip_CF2ClBr, ip_CF3Br,&
       ip_CHClBr2, ip_CH2ClBr, ip_CHCl2Br, ip_CH2Br2, ip_CBrF2CBrF2, ip_CHBr3/)



  END SUBROUTINE clamschem_initialize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamschem_init_memory

    USE messy_main_blather_bi,          ONLY: error_bi
    USE messy_main_channel_error_bi,    ONLY: channel_halt
    USE messy_main_channel_bi,          ONLY: REPR_LG_CLAMS
    USE messy_main_mpi_bi,              ONLY: p_pe

    USE messy_main_channel,             ONLY: new_channel, new_attribute, &
                                              new_channel_object
    USE messy_main_switch,              ONLY: USE_CLAMSMIX
    USE messy_clams_tools_utils,        ONLY: uppercase
    USE messy_clamschem_global,         ONLY: therm_flag, rates, const, &
                                              nhetpar, hetparam
    USE messy_clams_global,             ONLY: username
    USE messy_clamschem_asad_mod_clams, ONLY: jpbk, jptk, jppj, jphk, lhet
    USE messy_clamschem_reacpro,        ONLY: breac, treac, jreac, hreac, &
                                              check_dup_reac

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, TRIM

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamschem_init_memory'

    CHARACTER(30), dimension(:), pointer :: allnames
    INTEGER       :: lpro, lrec
    INTEGER       :: status
    INTEGER       :: ir, i
    CHARACTER(50) :: channelname
    CHARACTER(30) :: name
    CHARACTER(30) :: reactants, products

    IF (p_pe==0) WRITE(*,*) uppercase(substr)
     
    !-----------------------------------------------------------------
    ! Define channel CHEM_RATES
    !----------------------------------------------------------------
    if (rates) then

       channelname = modstr//'_RATES'

       ! Create new channel

       CALL new_channel(status, channelname, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, channelname, modstr//'_version', c = modver)
       CALL channel_halt(substr, status)

       !-----------------------------------------------------------------
       ! bimolecular reaction rates
       !-----------------------------------------------------------------

       ALLOCATE (BRATES(jpbk))
       allocate (allnames(jpbk))
       allnames = ''
       do ir = 1, jpbk 
          ! get name 
          call breac (ir, 'B', name, reactants, products, lpro, lrec) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          BRATES(ir)%longname='B('//trim(reactants)//' ->'//trim(products)//')'
          BRATES(ir)%units='cm^-3 s^-1' 
          BRATES(ir)%descr='Bimolecular reaction rate: ' &
                            //trim(reactants)//' ->'//trim(products)
       enddo
       do ir = 1, jpbk
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=BRATES(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(BRATES(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(BRATES(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(BRATES(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)

       !-----------------------------------------------------------------
       ! termolecular reaction rates
       !-----------------------------------------------------------------

       ALLOCATE (TRATES(jptk))
       allocate (allnames(jptk))
       allnames = ''
       do ir = 1, jptk 
          call treac (ir, 'T', name, reactants, products, lpro, lrec, therm_flag(ir)) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          TRATES(ir)%longname='T('//trim(reactants)//' ->'//trim(products)//')'
          TRATES(ir)%units='cm^-3 s^-1' 
          if (therm_flag(ir) == 0) then
             TRATES(ir)%descr='Termolecular reaction rate: ' &
                            //trim(reactants)//' ->'//trim(products)
          else
             TRATES(ir)%descr='Thermal decay rate: ' &
                            //trim(reactants)//' ->'//trim(products)
          endif
       enddo
       do ir = 1, jptk
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=TRATES(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL channel_halt (substr, status)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(TRATES(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(TRATES(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(TRATES(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)

       !-----------------------------------------------------------------
       ! photolysis rates
       !-----------------------------------------------------------------

       ALLOCATE (JRATES(jppj))
       allocate (allnames(jppj))
       allnames = ''
       do ir = 1, jppj 
          call jreac (ir, 'J', name, reactants, products, lpro, lrec) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          JRATES(ir)%longname='J('//trim(reactants)//' ->'//trim(products)//')'
          JRATES(ir)%units='cm^-3 s^-1' 
          JRATES(ir)%descr='Photolysis rate: ' &
                            //trim(reactants)//' ->'//trim(products)
       enddo
       do ir = 1, jppj
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=JRATES(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL channel_halt (substr, status)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(JRATES(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(JRATES(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(JRATES(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)

       !-----------------------------------------------------------------
       ! heterogeneous rates
       !-----------------------------------------------------------------

       if (lhet) then    
          ALLOCATE (HRATES(jphk))
          allocate (allnames(jphk))
          allnames = ''
          do ir = 1, jphk 
             call hreac (ir, 'H', name, reactants, products, lpro, lrec) 
             ! check duplicate reactants and add name to allnames
             call check_dup_reac (status,allnames,name,ir)
             IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
             ! set attributes
             HRATES(ir)%longname='H('//trim(reactants)//' ->'//trim(products)//')'
             HRATES(ir)%units='cm^-3 s^-1' 
             HRATES(ir)%descr='Heterogeneous reaction rate: ' &
                               //trim(reactants)//' ->'//trim(products)
          enddo
          do ir = 1, jphk
             if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
             CALL new_channel_object(status, channelname, allnames(ir), &
                                     p1=HRATES(ir)%values, &
                                     reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
             CALL channel_halt (substr, status)
             CALL new_attribute(status, channelname, allnames(ir), &
                  'creator_of_parameter', c = TRIM(username))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'longname', c = TRIM(HRATES(ir)%longname))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'units', c = TRIM(HRATES(ir)%units))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'description', c = TRIM(HRATES(ir)%descr))
             CALL channel_halt (substr, status)
          enddo
          deallocate (allnames)
       endif

    endif

    !-----------------------------------------------------------------
    ! Define channel CHEM_CONST
    !-----------------------------------------------------------------
    if (const) then

       !
       ! THE RATE CONSTANTS:
       !
       ! in (cm^3/molec) sec^-1 for bimol and termol rections
       ! in sec^-1 for photolysis and heterogeneous reactions
       ! (denoted by lower case letters)
       !

       channelname = modstr//'_CONST'
      
       ! Create new channel
       CALL new_channel(status, channelname, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, channelname, modstr//'_version', c = modver)
       CALL channel_halt(substr, status)

       !-----------------------------------------------------------------
       ! bimolecular reaction rate constants
       !-----------------------------------------------------------------

       ALLOCATE (BCONST(jpbk))
       allocate (allnames(jpbk))
       allnames = ''
       do ir = 1, jpbk 
          ! get name 
          call breac (ir, 'b', name, reactants, products, lpro, lrec) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          BCONST(ir)%longname='b('//trim(reactants)//' ->'//trim(products)//')'
          BCONST(ir)%units='cm^-3 s^-1' 
          BcONST(ir)%descr='Bimolecular reaction rate constant: ' &
                            //trim(reactants)//' ->'//trim(products)
       enddo
       do ir = 1, jpbk
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=BCONST(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(BCONST(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(BCONST(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(BCONST(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)


       !-----------------------------------------------------------------
       ! termolecular reaction rate constants
       !-----------------------------------------------------------------

       ALLOCATE (TCONST(jptk))
       allocate (allnames(jptk))
       allnames = ''
       do ir = 1, jptk 
          call treac (ir, 't', name, reactants, products, lpro, lrec, therm_flag(ir)) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          TCONST(ir)%longname='t('//trim(reactants)//' ->'//trim(products)//')'
          if (therm_flag(ir) == 0) then
             TCONST(ir)%units='cm^-3 s^-1' 
             TCONST(ir)%descr='Termolecular reaction rate constant: ' &
                            //trim(reactants)//' ->'//trim(products)
          else
             TCONST(ir)%units='s^-1' 
             TCONST(ir)%descr='Thermal decay rate constant: ' &
                            //trim(reactants)//' ->'//trim(products)
          endif
       enddo
       do ir = 1, jptk
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=TCONST(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL channel_halt (substr, status)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(TCONST(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(TCONST(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(TCONST(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)

       !-----------------------------------------------------------------
       ! photolysis rate constants
       !-----------------------------------------------------------------

       ALLOCATE (JCONST(jppj))
       allocate (allnames(jppj))
       allnames = ''
       do ir = 1, jppj 
          call jreac (ir, 'j', name, reactants, products, lpro, lrec) 
          ! check duplicate reactants and add name to allnames
          call check_dup_reac (status,allnames,name,ir)
          IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
          ! set attributes
          JCONST(ir)%longname='j('//trim(reactants)//' ->'//trim(products)//')'
          JCONST(ir)%units='s^-1' 
          JCONST(ir)%descr='Photolysis rate constant: ' &
                            //trim(reactants)//' ->'//trim(products)
       enddo
       do ir = 1, jppj
          if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
          CALL new_channel_object(status, channelname, allnames(ir), &
                                  p1=JCONST(ir)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL channel_halt (substr, status)
          CALL new_attribute(status, channelname, allnames(ir), &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, allnames(ir), &
               'longname', c = TRIM(JCONST(ir)%longname))
          CALL new_attribute(status, channelname, allnames(ir), &
               'units', c = TRIM(JCONST(ir)%units))
          CALL new_attribute(status, channelname, allnames(ir), &
               'description', c = TRIM(JCONST(ir)%descr))
          CALL channel_halt (substr, status)
       enddo
       deallocate (allnames)

       !-----------------------------------------------------------------
       ! heterogeneous rate constants
       !-----------------------------------------------------------------

       if (lhet) then    
          ALLOCATE (HCONST(jphk))
          allocate (allnames(jphk))
          allnames = ''
          do ir = 1, jphk 
             call hreac (ir, 'h', name, reactants, products, lpro, lrec) 
             ! check duplicate reactants and add name to allnames
             call check_dup_reac (status,allnames,name,ir)
             IF (status/=0) call error_bi ("Error check_dup_reac !!!",substr)
             ! set attributes
             HCONST(ir)%longname='h('//trim(reactants)//' ->'//trim(products)//')'
             HCONST(ir)%units='s^-1' 
             HCONST(ir)%descr='Heterogeneous reaction rate constant: ' &
                               //trim(reactants)//' ->'//trim(products)
          enddo
          do ir = 1, jphk
             if (p_pe==0) write (*,*) 'create channel object:', allnames(ir)
             CALL new_channel_object(status, channelname, allnames(ir), &
                                     p1=HCONST(ir)%values, &
                                     reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
             CALL channel_halt (substr, status)
             CALL new_attribute(status, channelname, allnames(ir), &
                  'creator_of_parameter', c = TRIM(username))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'longname', c = TRIM(HCONST(ir)%longname))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'units', c = TRIM(HCONST(ir)%units))
             CALL new_attribute(status, channelname, allnames(ir), &
                  'description', c = TRIM(HCONST(ir)%descr))
             CALL channel_halt (substr, status)
          enddo
          deallocate (allnames)
       endif

    endif


    !-----------------------------------------------------------------
    ! Define channel CHEM_HETPAR
    !-----------------------------------------------------------------
    if (hetparam) then 
       if (lhet) then

       channelname = modstr//'_HETPAR'

       ! Create new channel
       CALL new_channel(status, channelname, lrestreq=.FALSE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, channelname, modstr//'_version', c = modver)
       CALL channel_halt(substr, status)

       do i = 1, nhetpar

          if (p_pe==0) write (*,*) 'create channel object ',trim(HETPAR(i)%name)
          CALL new_channel_object(status, channelname, HETPAR(i)%name, &
                                  p1=HETPAR(i)%values, &
                                  reprid=REPR_LG_CLAMS, lrestreq=.FALSE.)
          CALL channel_halt (substr, status)
          CALL new_attribute(status, channelname, HETPAR(i)%name, &
               'creator_of_parameter', c = TRIM(username))
          CALL new_attribute(status, channelname, HETPAR(i)%name, &
               'longname', c = TRIM(HETPAR(i)%longname))
          CALL new_attribute(status, channelname, HETPAR(i)%name, &
               'units', c = TRIM(HETPAR(i)%units))
          CALL channel_halt (substr, status)

       enddo

       endif
    endif

  END SUBROUTINE clamschem_init_memory
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamschem_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_pe
    ! SMCL: CLAMS
    USE messy_clams_global,          ONLY: nchemspec, SPECARR, init_vertcoorname
    USE messy_clamschem_global,      ONLY: jpdim, ip_messy, dissoc_rate
    USE messy_clams_tools_utils,     ONLY: uppercase
    ! SMCL: MESSY
    USE messy_main_channel,          ONLY: get_channel_object, &
                                           set_channel_output, &
                                           STRLEN_OBJECT
    USE messy_cmn_photol_mem,        ONLY: jname

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamschem_init_coupling'
    CHARACTER(LEN=STRLEN_OBJECT) :: objname = ''
    INTEGER  :: status, i

    IF (p_pe==0) WRITE(*,*) uppercase(substr)

    ! Get arrays from CLAMS submodel:
    CALL get_channel_object(status, 'clams', 'TEMP', p1=TEMP)
    CALL channel_halt(substr, status)
    IF (TRIM(init_vertcoorname) == 'press') THEN 
       CALL get_channel_object(status, 'clams', 'LEV', p1=PRESS)
       CALL channel_halt(substr, status)
    ELSE
       CALL get_channel_object(status, 'clams', 'PRESS', p1=PRESS)
       CALL channel_halt(substr, status)
    END IF
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)
    !CALL get_channel_object(status, 'clams', 'JULSEC', p1=JULSEC)
    !CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME', p0=JULTIME)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'clams', 'LAT_OLD', p1=LAT_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON_OLD', p1=LON_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV_OLD', p1=LEV_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULTIME_OLD', p0=JULTIME_OLD)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'TEMP_OLD', p1=TEMP_OLD)
    CALL channel_halt(substr, status)
    IF (TRIM(init_vertcoorname) == 'press') THEN 
       CALL get_channel_object(status, 'clams', 'LEV_OLD', p1=PRESS_OLD)
       CALL channel_halt(substr, status)
    ELSE
       CALL get_channel_object(status, 'clams', 'PRESS_OLD', p1=PRESS_OLD)
       CALL channel_halt(substr, status)
    END IF


    ! couple photolysis rates
    allocate (DISSOC_RATE(jpdim))
    DO i = 1, jpdim
       objname = 'J'//jname(ip_messy(i))
       CALL get_channel_object(status, 'dissoc', objname, p1=DISSOC_RATE(i)%values)
       if (p_pe==0) write (*,*) 'in clamschem_init_coupling ',i, ' ', trim(objname)
    ENDDO



    ! Get chemical species (from CLAMS):
    DO i = 1, nchemspec
       CALL get_channel_object(status, 'clams', &
            trim(SPECARR(i)%name), p1=CHEMSPECARR(i)%values)
       CALL channel_halt(substr, status)
       CHEMSPECARR(i)%name     = SPECARR(i)%name
       CHEMSPECARR(i)%longname = SPECARR(i)%longname
       CHEMSPECARR(i)%units    = SPECARR(i)%units
       CHEMSPECARR(i)%ctype    = SPECARR(i)%ctype
    ENDDO


  END SUBROUTINE clamschem_init_coupling
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamschem_global_end

    ! BMIL
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt

    ! SMCL
    USE messy_main_channel,    ONLY: set_channel_output
    USE messy_main_switch,     ONLY: USE_CLAMSMIX
    USE messy_main_timer,      ONLY: YEAR, MONTH

    USE messy_clams_global,   ONLY: rank, dnparts_max, dnparts, ldiagout, &
                                    nchemspec, lchemevent,  &
                                    mdi, eps, dates30

    USE messy_clamschem,            ONLY: chem
    USE messy_clamschem_global,     ONLY: rates, const, hetparam, &
                                          js_monthmean
    USE messy_clams_tools_dateconv, ONLY: ymds2js_interface

!    USE messy_dissoc,          ONLY: 

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamschem_global_end'
    INTEGER :: status
    INTEGER :: i

    IF (lchemevent) THEN
        IF (rank == 0 .and. ldiagout) THEN
          WRITE(*,*) 'Active CHEM event'
          DO i = 1, nchemspec
             WRITE (*,*) 'i,chemspecarr(i)%name=',i,chemspecarr(i)%name
          END DO
       END IF

       ! js_monthmean (15th of current month) for dissoc:  
       js_monthmean = ymds2js_interface(YEAR, MONTH, 15, 43200, dates30) 

!!!!!! ???
       if (dnparts>0) then

          ! ju_jug_20150911 This filtering out of air parcels with negative
          !                 concentrations may be not ideal
          DO i=1, nchemspec
             IF (chemspecarr(i)%ctype=='TR') THEN
                WHERE (chemspecarr(i)%values .LT. 0.)
                   LAT = mdi
                   LON = mdi
                   LEV = mdi
                END WHERE
             END IF
          END DO


          CALL chem (status, CHEMSPECARR, LAT, LON, LEV, TEMP, PRESS, &
                     LAT_OLD, LON_OLD, LEV_OLD, TEMP_OLD, PRESS_OLD, &
                     BRATES, TRATES, JRATES, HRATES, &
                     BCONST, TCONST, JCONST, HCONST, HETPAR)
          IF (status/=0) call error_bi ("Error in CHEM !!!",substr)
 

       endif

       ! Set _OLD arrays for next call of CHEM:
       if (rank==0) write (*,*) 'dnparts auf rank 0 nach CHEM:',dnparts
       if (rank==0) write (*,*) 'Save OLD positions !'
       LAT_OLD = LAT
       LON_OLD = LON
       LEV_OLD = LEV
       JULTIME_OLD = JULTIME
       TEMP_OLD = TEMP
       PRESS_OLD = PRESS

!!$       if (rates) then
!!$          CALL set_channel_output(status, 'clamschem_RATES', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       endif
!!$       if (const) then
!!$          CALL set_channel_output(status, 'clamschem_CONST', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       endif
!!$       if (hetparam) then
!!$          CALL set_channel_output(status, 'clamschem_HETPAR', .TRUE.)
!!$          CALL channel_halt(substr, status)
!!$       endif

!!!!!! => Output in clams_main.f90   
!!$       if (rates .or. const .or. hetparam) then
!!$          CALL set_channel_output(status, 'clams', .FALSE.)
!!$          CALL channel_halt(substr, status)
!!$          CALL messy_write_output
!!$          if (lclamsoutevent) then
!!$             CALL set_channel_output(status, 'clams', .TRUE.)
!!$             CALL channel_halt(substr, status)
!!$          endif
!!$       endif
       
    ENDIF

  END SUBROUTINE clamschem_global_end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamschem_free_memory

    use messy_clamschem_asad_mod,        ONLY: asad_mod_final
    USE messy_clamschem_asad_mod_clams,  ONLY: lhet
    USE messy_clams_global,      ONLY: rank, ldiagout
    USE messy_clamschem,         ONLY: deallocate_chem_vars
    USE messy_clamschem_global,  ONLY: rates, const

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamschem_free_memory'

    ! Deallocate ASAD arrays
    if (rank==0 .and. ldiagout) write (*,*) 'call asad_mod_final'
    call asad_mod_final

    ! Deallocate
    call deallocate_chem_vars
 
    DEALLOCATE(CHEMSPECARR)

    if (rates) then
       DEALLOCATE (BRATES, TRATES, JRATES)
       if (lhet) DEALLOCATE (HRATES)
    endif
    if (const) then
       DEALLOCATE (BCONST, TCONST, JCONST)
       if (lhet) DEALLOCATE (HCONST)
    endif

  END SUBROUTINE clamschem_free_memory
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
END MODULE messy_clamschem_si
