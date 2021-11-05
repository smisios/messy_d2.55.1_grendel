!**********************************************************************
MODULE messy_clamscheme5_si

#if defined(ECHAM5)
!**********************************************************************
!  Submodel interface for clamschem 
!**********************************************************************

  USE messy_clamscheme5
  USE messy_main_blather_bi,  ONLY: start_message_bi, end_message_bi, error_bi
  USE messy_clams_global,     ONLY: species_type, species_type_3d
  USE messy_clamschem_global, ONLY: rc_type

  IMPLICIT NONE

  ! Module variables:
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LAT3D_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_LON3D_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_ZETA_D => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_TEMP => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: E5_PRESS => NULL()

  TYPE(species_type), DIMENSION(:), POINTER :: E5CHEMSPECARR ! ECHAMTRACER
  TYPE(species_type_3d), DIMENSION(:), POINTER :: E5CHEMSPECARR_tte ! ECHAMTRACER Tendency

  TYPE(rc_type), DIMENSION(:), POINTER :: BRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: TRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: JRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: HRATES => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: BCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: TCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: JCONST => NULL()
  TYPE(rc_type), DIMENSION(:), POINTER :: HCONST => NULL()
  TYPE(species_type), DIMENSION(:), POINTER :: HETPAR => NULL()


  PUBLIC :: clamscheme5_initialize
  PUBLIC :: clamscheme5_init_memory
  PUBLIC :: clamscheme5_init_coupling
  PUBLIC :: clamscheme5_new_tracer
  PUBLIC :: clamscheme5_init_tracer
  PUBLIC :: clamscheme5_global_end
  PUBLIC :: clamscheme5_free_memory

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

  SUBROUTINE clamscheme5_initialize

    ! SMCL
    USE messy_main_timer,        ONLY: delta_time
    USE messy_main_tools,        ONLY: find_next_free_unit

    USE messy_clams_global,          ONLY: nspec, nchemspec, maxspec, SPECARR, &
                                           rank, initfile, dnparts_max, specnames
    USE messy_clamsrdfrc_tools,      ONLY: nc_read_ap_s_info
    USE messy_clamsrdfrc_global,     ONLY: lev_window 

    ! USE messy_clamschem_global,      ONLY: timestep_chem, ncdt, asad_gfirst
    ! USE messy_clamschem_defs_mod,    ONLY: chch_defs
    ! USE messy_clamschem_data,        ONLY: clams_chem_init
    ! USE messy_clamschem_data_hetero, ONLY: clams_chem_init_hetero
    ! USE messy_clamschem,             ONLY: allocate_chem_vars
    ! USE messy_clamschem_usrhet,      ONLY: hetvar, nhetspec, nhetpar
    ! USE asad_mod,                    ONLY: asad_mod_init, &
    !                                        speci, advt, family, madvtr, majors, &
    !                                        moffam, ctype, &
    !                                        nit0, nitfg, nitnr, nrsteps, &
    !                                        ldepd, ldepw, dtime, cdt
    !USE asad_mod_clams,              ONLY: jpspec, jpctr, lhet


    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_initialize'

    INTEGER :: status, iou, ichem
!!$     INTEGER :: status, iou, ios, n, nodd, i, ichem
!!$    CHARACTER(2)  :: ctype
!!$    CHARACTER(20) :: specname

    !WRITE(*,*) substr
    
    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    CALL clamscheme5_read_nml(status, iou)
    CALL nc_read_ap_s_info(initfile, lev_window)


!!$ ASAD-F77: START
!!$ 
!!$    ALLOCATE (E5CHEMSPECARR(nchemspec))
!!$    ALLOCATE (E5CHEMSPECARR_tte(nchemspec))
!!$
!!$    ! Read chemical species from file data/chch.d
!!$    iou = find_next_free_unit(100,200)
!!$    OPEN(iou,file='data/chch.d',iostat=ios)
!!$    IF (ios/=0) THEN
!!$       CALL error_bi ("File data/chch.d could not be opened !!!",substr)
!!$    ENDIF
!!$    ichem = 0
!!$    DO
!!$       READ (iou,*,iostat=ios) n, specname, nodd, ctype
!!$       IF (ios/=0) EXIT
!!$       IF (specname/='') THEN
!!$          ichem = ichem + 1
!!$          IF (nchemspec > maxspec) &
!!$               CALL error_bi (substr,"To many species: increase MAXSPEC !!!")
!!$          E5CHEMSPECARR(ichem)%name = TRIM(specname)
!!$          E5CHEMSPECARR(ichem)%ctype = ctype
!!$          E5CHEMSPECARR(ichem)%longname = TRIM(specname)
!!$          E5CHEMSPECARR(ichem)%units = 'm^3/m^3'
!!$          E5CHEMSPECARR_tte(ichem)%name = TRIM(specname)
!!$          E5CHEMSPECARR_tte(ichem)%ctype = ctype
!!$          E5CHEMSPECARR_tte(ichem)%longname = TRIM(specname)
!!$       ENDIF
!!$    ENDDO
!!$
!!$ ASAD-F77: END


!!!!!!!!!!!!!!!!!!!! ASAD-F90: START

!!!!! Wird bereits in clamschem_initialize gesetzt ?!?
!
!     ! Get chemistry input
!     IF (rank==0) write (*,*) 'call clams_chem_init'
!     call clams_chem_init
    
!     ! Allocate and initialize ASAD arrays and variables
!     IF (rank==0) write (*,*) 'call asad_mod_init'
!     call asad_mod_init
    
! !!!!! Nach asad_mod_init belegen:
!     ldepd = .FALSE. ! in asad_mod.f90 deklariert
!     ldepw = .FALSE. ! in asad_mod.f90 deklariert
!     nit0    = 20    ! in asad_mod: mit 20 initialisiert
!     nitfg   = 20    ! in asad_mod_init: mit 10 initialisiert
!     nitnr   = 20    ! in asad_mod_init: mit 10 initialisiert
!     nrsteps = 50    ! in asad_mod_init: mit 45 initialisiert
!     cdt     = ncdt  ! in asad_mod.f90 deklariert

!     dtime = timestep_chem   ! dtime in asad_mod_init initialisiert !

!     asad_gfirst = .true.

!     ! Allocate arrays 
!     call allocate_chem_vars
    
!     ! --------------------------------------
!     !    Initialise the Chemistry
!     ! --------------------------------------
!     !
!     ! Call ASAD routine CINIT
!     !
!     if (rank==0) write (*,*) 'call asad_cinit !!!'
!     call asad_cinit (1)

!     if (rank==0) then 
!        write (6, *) 'speci= ', (speci(i),i=1,jpspec) 
!        write (6, *) 'advt= ', (advt(i),i=1,jpctr) 
!        write (6, *) 'family= ', (family(i),i=1,jpspec) 
!        write (6, *) 'madvtr= ', (madvtr(i),i=1,jpspec) 
!        write (6, *) 'majors= ', (majors(i),i=1,jpctr) 
!        write (6, *) 'moffam= ', (moffam(i),i=1,jpspec) 
!        write (6, *) 'ctype= ', (ctype(i),i=1,jpspec) 
!     endif

    ! ! Read chemical species from structure chch_defs
    ! nchemspec = jpspec
    ! if (nchemspec > maxspec) &
    !      call error_bi (substr,"To many species: increase MAXSPEC !!!")
    ! DO i = 1, jpspec
    !    specnames(i)        = chch_defs(i)%speci
    !    SPECARR(i)%name     = chch_defs(i)%speci
    !    SPECARR(i)%ctype    = chch_defs(i)%ctype
    !    SPECARR(i)%longname = chch_defs(i)%speci
    !    SPECARR(i)%units    = 'm^3/m^3'
    ! ENDDO

    ! ! for heterogeneous chemistry: add species
    ! if (lhet) then

    !    IF (rank==0) write (*,*) 'call clams_chem_init_hetero'
    !    call clams_chem_init_hetero
    !    DO i = 1, nhetspec
    !       specnames(nchemspec+i)        = hetvar(i)%name
    !       SPECARR(nchemspec+i)%name     = hetvar(i)%name
    !       SPECARR(nchemspec+i)%ctype    = ' '
    !       SPECARR(nchemspec+i)%longname = hetvar(i)%longname
    !       SPECARR(nchemspec+i)%units    = hetvar(i)%units
    !    ENDDO
    !    nchemspec = nchemspec + nhetspec

    !    allocate (HETPAR(nhetpar))
    !    do i = 1, nhetpar
    !       HETPAR(i)%name     = hetvar(nhetspec+i)%name
    !       HETPAR(i)%longname = hetvar(nhetspec+i)%longname
    !       HETPAR(i)%units    = hetvar(nhetspec+i)%units
    !    enddo

    ! endif

    ! nspec = nchemspec
    ! if (rank==0) then
    !    write (*,*) 'nchemspec=',nchemspec
    !    do i = 1, nchemspec
    !       write(*,*) "i,speciesname(i),type ",i,SPECARR(i)%name,SPECARR(i)%ctype
    !    enddo
    !    if (lhet) then
    !       write (*,*) 'nhetpar=',nhetpar
    !       do i = 1, nhetpar
    !          write(*,*) "i,hetpar(i) ",i,HETPAR(i)%name
    !       enddo
    !    endif
    ! endif
    ! CLOSE(iou)

    ALLOCATE (E5CHEMSPECARR(nchemspec))
    ALLOCATE (E5CHEMSPECARR_tte(nchemspec))

    DO ichem = 1, nchemspec
   
       E5CHEMSPECARR(ichem)%name     = SPECARR(ichem)%name
       E5CHEMSPECARR(ichem)%ctype    = SPECARR(ichem)%ctype
       E5CHEMSPECARR(ichem)%longname = SPECARR(ichem)%longname
       E5CHEMSPECARR(ichem)%units    = SPECARR(ichem)%units
       E5CHEMSPECARR_tte(ichem)%name     = SPECARR(ichem)%name
       E5CHEMSPECARR_tte(ichem)%ctype    = SPECARR(ichem)%ctype
       E5CHEMSPECARR_tte(ichem)%longname = SPECARR(ichem)%longname
       
    ENDDO

!!!!!!!!!!!!!!!!!!! ASAD-F90: END

    IF (rank==0) THEN
       DO ichem = 1, nchemspec
          WRITE(*,*) " E5: ichem,speciesname(ichem),type ",ichem,E5CHEMSPECARR(ichem)%name,&
               E5CHEMSPECARR(ichem)%ctype
       ENDDO
    ENDIF
    CLOSE(iou)

  END SUBROUTINE clamscheme5_initialize
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, REPR_LG_CLAMS
    USE messy_main_channel,          ONLY: new_channel, new_attribute &
                                         , new_channel_object
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_clams_global,          ONLY: nchemspec, username
    USE messy_clamschem_global,      ONLY: rates, const, hetparam &
                                        , therm_flag, nhetpar
    USE messy_clamschem_asad_mod_clams, ONLY: jpbk, jptk, jppj, jphk, lhet
    USE messy_clamschem_reacpro,     ONLY: breac, treac, jreac, hreac, &
                                           check_dup_reac

    IMPLICIT NONE

    ! LOCAL 
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_init_memory'
    CHARACTER(30), dimension(:), pointer :: allnames
    CHARACTER(50) :: channelname
    CHARACTER(30) :: name
    CHARACTER(30) :: reactants, products
    INTEGER       :: lpro, lrec
    INTEGER       :: status, i, ir
 
    ! Define channel CHEM
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

    DO i=1, nchemspec
       IF (E5CHEMSPECARR_tte(i)%ctype == 'TR') THEN
          CALL new_channel_object(status, modstr, TRIM(E5CHEMSPECARR_tte(i)%name)//'_tte',&
               p3=E5CHEMSPECARR_tte(i)%values, reprid=GP_3D_MID, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
       END IF
    END DO

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
!!!!!   ???
    if (hetparam) then 
!!$    if (lhet) then

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

  END SUBROUTINE clamscheme5_init_memory
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_init_coupling'
    INTEGER :: status

    CALL start_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)

    CALL get_channel_object(status, 'ECHAM5', 'tm1', p3=E5_TEMP)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'ECHAM5', 'press', p3=E5_PRESS)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_ZETA', p3=E5_ZETA_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_LAT3D', p3=E5_LAT3D_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_LON3D', p3=E5_LON3D_D)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'COUPLING TO DRIVER FIELDS', substr)

  END SUBROUTINE clamscheme5_init_coupling
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_new_tracer

    USE messy_main_tracer,        ONLY: new_tracer, tracer_error_str
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG
    USE messy_clams_global,       ONLY: nchemspec
    USE messy_main_mpi_bi,        ONLY: p_pe

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_new_tracer'

    INTEGER :: status, iou, i
    CHARACTER(LEN=STRLEN_VLONG):: msg

    IF (p_pe==0) WRITE(*,*) substr

    DO i=1, nchemspec
       IF(TRIM((E5CHEMSPECARR(i)%name))=='H2O_100') E5CHEMSPECARR(i)%name='H2O100'
       IF(TRIM((E5CHEMSPECARR(i)%name))=='IWC_100') E5CHEMSPECARR(i)%name='IWC100'

       CALL new_tracer(status, 'gp', E5CHEMSPECARR(i)%name, modstr, &
            subname=GPTRSTR)
       IF (status .NE. 0) THEN
          msg=tracer_error_str(status)
          WRITE(*,*) msg
          STOP
       END IF
    END DO

  END SUBROUTINE clamscheme5_new_tracer
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_init_tracer

    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl
    USE netcdf
    USE messy_main_mpi_bi,        ONLY: p_pe, p_barrier
    USE messy_clams_tools_ncutils,ONLY: nc_check_error
    USE messy_clams_global,       ONLY: nchemspec, mdi
    USE messy_main_timer,         ONLY: lstart

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_init_tracer'
    INTEGER :: rcode, ncid, status, varid
    INTEGER :: i
    REAL(DP), DIMENSION(:,:,:), POINTER:: E5_TRACER_TEMP_G
    REAL(DP), DIMENSION(:,:,:), POINTER:: E5_TRACER_TEMP_D
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE:: E5_INIT_TEMP
    INTEGER,DIMENSION(NF90_MAX_VAR_DIMS) :: dim_array
    INTEGER :: dim_len1, dim_len2, dim_len3
    INTEGER :: ilon, ilat, ilev

    IF (lstart) THEN
       !Read tracer from initfile:
       rcode= nf90_open(TRIM(e5chem_initfile), nf90_nowrite, ncid)
       call nc_check_error (rcode,'Cannot open file',abort=.true.)
       
       DO i=1,nchemspec
          IF (p_pe==0) WRITE(*,*) 'Initialize EMAC Tracer ',E5CHEMSPECARR(i)%name
          ! Set pointer to EMAC Tracer Array:
          CALL get_tracer(status, 'gp', E5CHEMSPECARR(i)%name, subname=GPTRSTR,&
               & pxt=E5_TRACER_TEMP_D)
          CALL trp_gpdc_gpgl (1,E5_TRACER_TEMP_D, E5_TRACER_TEMP_G)

          IF (E5CHEMSPECARR(i)%ctype=='TR') THEN
             ! Read tracer from e5chem_initfile:
             rcode = nf90_inq_varid(ncid,TRIM(E5CHEMSPECARR(i)%name)//'_gp',varid)
             call nc_check_error (rcode,'Cannot find variable'//E5CHEMSPECARR(i)%name,abort=.true.)
             rcode=nf90_inquire_variable(ncid,varid,dimids=dim_array)
             rcode=nf90_inquire_dimension(ncid,dim_array(1),len=dim_len1)
             rcode=nf90_inquire_dimension(ncid,dim_array(2),len=dim_len2)
             rcode=nf90_inquire_dimension(ncid,dim_array(3),len=dim_len3)
             ALLOCATE(E5_INIT_TEMP(dim_len1, dim_len2, dim_len3))
             rcode = nf90_get_var(ncid,varid,E5_INIT_TEMP)
             call nc_check_error (rcode,'Cannot read variable'//E5CHEMSPECARR(i)%name,abort=.true.)
          END IF

          E5_TRACER_TEMP_G(:,:,:) = 0.

          IF (E5CHEMSPECARR(i)%ctype=='TR') THEN
             DO ilon=1,dim_len1
                DO ilat=1, dim_len2
                   DO ilev =1, dim_len3
                      E5_TRACER_TEMP_G(ilon,ilev,ilat) = E5_INIT_TEMP(ilon,ilat,ilev)
                      IF (E5_INIT_TEMP(ilon,ilat,ilev) .LT. 0. ) THEN
                         IF (ilat==1) THEN
                            E5_TRACER_TEMP_G(ilon,ilev,ilat) = 0.
                         ELSE
                            E5_TRACER_TEMP_G(ilon,ilev,ilat) = E5_TRACER_TEMP_G(ilon,ilev,ilat-1)
                         END IF
                      END IF
                   END DO
                END DO
             END DO
             DEALLOCATE (E5_INIT_TEMP)
          END IF

          CALL trp_gpdc_gpgl (-1,E5_TRACER_TEMP_D, E5_TRACER_TEMP_G)
       END DO

       ! Close INIT-File
       rcode=nf90_close(ncid)
    END IF

  END SUBROUTINE clamscheme5_init_tracer
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_global_end

    USE messy_main_timer,    ONLY: time_step_len, MONTH, DAY, HOUR, MINUTE
    USE messy_clams_global,  ONLY: nchemspec, nparts_max, mdi, eps, &
                                   dnparts_max, lperpetuum
    USE messy_main_mpi_bi,   ONLY: p_pe
    USE messy_clams_global,  ONLY: mdi, lchemevent
    USE messy_main_tracer_mem_bi,   ONLY: xtte, xtm1
    USE messy_clamschem,            ONLY: chem
    USE messy_clamschem_global,     ONLY: ncdt
    USE messy_main_grid_def_mem_bi, ONLY: nlev, ngpblks, nproma, npromz

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamscheme5_global_end'

    INTEGER :: i, ilat, ilev, icnt, counter, ncnt, status
    REAL(DP), DIMENSION(:), POINTER :: E5_TEMP_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_PRESS_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_LAT_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_LON_1D
    REAL(DP), DIMENSION(:), POINTER :: E5_ZETA_1D

    IF (lchemevent) THEN

       IF (p_pe==0) WRITE(*,*) 'Active E5CHEM Event'
       
       ! 1D-Arrays of positions, temperature, pressure
       ALLOCATE (E5_TEMP_1D(dnparts_max))
       ALLOCATE (E5_PRESS_1D(dnparts_max))
       ALLOCATE (E5_LAT_1D(dnparts_max))
       ALLOCATE (E5_LON_1D(dnparts_max))
       ALLOCATE (E5_ZETA_1D(dnparts_max))
       E5_TEMP_1D = mdi
       E5_PRESS_1D = mdi
       E5_LAT_1D = mdi
       E5_LON_1D = mdi
       E5_ZETA_1D = mdi
       counter = 0
       DO ilev = 0, nlev-1
          DO ilat = 0, ngpblks-1
             IF (ilat==ngpblks-1) THEN
                ncnt = npromz
             ELSE
                ncnt = nproma
             END IF
             DO icnt = 0, ncnt-1
                E5_TEMP_1D (counter+1) = E5_TEMP   (icnt+1,ilev+1,ilat+1)
                E5_PRESS_1D(counter+1) = E5_PRESS  (icnt+1,ilev+1,ilat+1)/100.
                E5_LAT_1D  (counter+1) = E5_LAT3D_D(icnt+1,ilev+1,ilat+1)
                E5_LON_1D  (counter+1) = E5_LON3D_D(icnt+1,ilev+1,ilat+1)
                E5_ZETA_1D (counter+1) = E5_ZETA_D (icnt+1,ilev+1,ilat+1)
                counter = counter + 1
             END DO
          END DO
       END DO

       ! Set values in E5CHEMSPECARR
       DO i=1, nchemspec
          ALLOCATE (E5CHEMSPECARR(i)%values(dnparts_max))
          E5CHEMSPECARR(i)%values=mdi
          ! Value of last timestep + tendency * timestep-length
          counter = 0
          DO ilev = 0, nlev-1
             DO ilat = 0, ngpblks-1
                IF (ilat==ngpblks-1) THEN
                   ncnt = npromz
                ELSE
                   ncnt = nproma
                END IF
                DO icnt = 0, ncnt-1
                   E5CHEMSPECARR(i)%values(counter+1) = &
                        xtm1(icnt+1,ilev+1,i,ilat+1) + xtte(icnt+1,ilev+1,i,ilat+1) * time_step_len
                   IF (E5CHEMSPECARR(i)%values(counter+1) .LT.0.) THEN
                      E5CHEMSPECARR(i)%values(counter+1) = 0.
                   END IF
                   counter = counter + 1
                END DO
             END DO
          END DO
       END DO
      
       CALL chem (status, E5CHEMSPECARR, E5_LAT_1D, E5_LON_1D, E5_ZETA_1D, E5_TEMP_1D, &
             E5_PRESS_1D, E5_LAT_1D, E5_LON_1D, E5_ZETA_1D, E5_TEMP_1D, E5_PRESS_1D, &
             BRATES, TRATES, JRATES, HRATES, &
             BCONST, TCONST, JCONST, HCONST, HETPAR)
             

       ! Reset mean age tracer in perpetuum run:
       IF (lperpetuum) THEN
          IF (MONTH==1 .AND. DAY==1) THEN
             DO i=1,nchemspec
                IF (trim(E5CHEMSPECARR(i)%name)/='BA') THEN
                   CYCLE
                ELSE
                   counter = 0
                   DO ilev = 0, nlev-1
                      DO ilat = 0, ngpblks-1
                         IF (ilat==ngpblks-1) THEN
                            ncnt = npromz
                         ELSE
                            ncnt = nproma
                         END IF
                         DO icnt = 0, ncnt-1
                            E5CHEMSPECARR(i)%values(counter+1) = &
                                 E5CHEMSPECARR(i)%values(counter+1)-0.365
                            counter = counter + 1
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END IF
       END IF
       
       DO i=1, nchemspec
          IF (E5CHEMSPECARR(i)%ctype == 'TR') THEN
             counter = 0
             DO ilev = 0, nlev-1
                DO ilat = 0, ngpblks-1
                   IF (ilat==ngpblks-1) THEN
                      ncnt = npromz
                   ELSE
                      ncnt = nproma
                   END IF
                   DO icnt = 0, ncnt-1
                      E5CHEMSPECARR_tte(i)%values(icnt+1,ilev+1,ilat+1) = (((E5CHEMSPECARR(i)%values(counter+1)&
                           -xtm1(icnt+1,ilev+1,i,ilat+1))/ time_step_len) - xtte(icnt+1,ilev+1,i,ilat+1)) / (ncdt/time_step_len)
                      counter = counter + 1
                   END DO
                END DO
             END DO
          END IF
          DEALLOCATE (E5CHEMSPECARR(i)%values)
       END DO
  
       DEALLOCATE(E5_TEMP_1D)
       DEALLOCATE(E5_PRESS_1D)
       DEALLOCATE(E5_LAT_1D)
       DEALLOCATE(E5_LON_1D)
       DEALLOCATE(E5_ZETA_1D)
    END IF !lchemevent

    ! Add tendencies due to chemistry:
    DO i=1, nchemspec
       IF (E5CHEMSPECARR(i)%ctype == 'CT') THEN
          xtte(:,:,i,:) = 0. 
       ELSE
          xtte(:,:,i,:) = xtte(:,:,i,:) + E5CHEMSPECARR_tte(i)%values(:,:,:) 
          WHERE (xtm1(:,:,i,:)+xtte(:,:,i,:)*time_step_len .LT. 0.)
             xtte(:,:,i,:) = - xtm1(:,:,i,:) / time_step_len   
          END WHERE
       END IF
    END DO


  END SUBROUTINE clamscheme5_global_end
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  SUBROUTINE clamscheme5_free_memory

    USE messy_clamschem_asad_mod_clams, ONLY: lhet
    USE messy_clamschem_global,         ONLY: rates, const

    IMPLICIT NONE

    DEALLOCATE(E5CHEMSPECARR)
    DEALLOCATE (E5CHEMSPECARR_tte)

    if (rates) then
       DEALLOCATE (BRATES, TRATES, JRATES)
       if (lhet) DEALLOCATE (HRATES)
    endif
    if (const) then
       DEALLOCATE (BCONST, TCONST, JCONST)
       if (lhet) DEALLOCATE (HCONST)
    endif

  END SUBROUTINE clamscheme5_free_memory
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
#endif
END MODULE messy_clamscheme5_si
