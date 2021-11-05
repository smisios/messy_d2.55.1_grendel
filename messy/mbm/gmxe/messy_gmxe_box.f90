PROGRAM GMXE_BOX

  USE messy_main_constants_mem,  ONLY: dp
  USE messy_gmxe_mem
  USE messy_gmxe

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: nmaxtrac=1000

  TYPE trac_struct
    REAL(dp)          :: field = 0._dp
    CHARACTER(LEN=32) :: name=''
    CHARACTER(LEN=2)  :: subname=''
    CHARACTER(LEN=32) :: S_aerosol_model='gmxe'
    CHARACTER(LEN=32) :: unit=''
    INTEGER           :: I_aerosol_mode
    INTEGER           :: I_aerosol_sol
    INTEGER           :: I_aerosol_method
!    REAL              :: R_aerosol_density
    REAL              :: R_molarmass
!    REAL              :: R_henry
!    REAL              :: R_dryreac_sf
    INTEGER           :: medium
    INTEGER           :: index
    INTEGER           :: quantity
  END type trac_struct

  TYPE (trac_struct), DIMENSION(nmaxtrac) :: xt
  INTEGER :: tcount

  REAL(dp), DIMENSION(:),   POINTER  :: rdryaer, rwetaer, ddryaer, aernumb
  REAL(dp), DIMENSION(:,:), POINTER  :: diagaer

  REAL(dp) :: temp_in, press_in, rh_in
  LOGICAL,  DIMENSION(:),     ALLOCATABLE :: lrhoa            ! true if air density is needed for unit conversion
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zconvert,      & ! factor to convert unit including air density
                                             xconvert         ! factor to convert unit excluding air density

  ! for output

  INTEGER, DIMENSION(:), POINTER, SAVE :: idvar
  INTEGER, SAVE :: it, idfile


  
  CALL gmxe_initialize
  print*, "Inititalisation complete!"

  CALL define_pseudo_tracers

  print*, "Pseudo-Tracer definition complete!"

  print*, "The next SR reads the boundary conditions, i.e. temperature, "
  print*, "pressure and RH and the input concentrations of all given species "
  print*, "from a data file: gmxe_input.txt"
  print*, "In the first row of this value T, p, and RH must be written!"
  print*, "All following rows contain Speciesname, subname and mixing ratio!"
  
  CALL init_fields
  print*, "Initialisation complete!"

  call gmxe_box_init_mem
  CAll gmxe_box_init_cpl

  print*, "calling driver"
  call gmxe_box_driver
  print*, "done!"
  
CONTAINS
  !==============================================================================

  SUBROUTINE gmxe_initialize
  !
    USE messy_main_tools,        ONLY: find_next_free_unit
!!$    USE messy_ncregrid_tools_bi, ONLY: rgtevent_init_nml ! op_pj_20160815
    USE messy_gmxe_aerchem,      ONLY: init_aerchem
    USE messy_gmxe_soa,          ONLY: set_soa_params, l_gasprec, l_oxi, &
                                       ngas_max_SOA => ngas_max, noxi_MAX, &
                                       excl_str_soa, exclude_max_soa,    &
                                       noxi, ngas_soa => ngas, nsoa
    USE messy_gmxe_isorropia2,   ONLY: init_species_isorropia
    USE messy_gmxe_eqsam4clim,   ONLY: init_species_eqsam4clim


    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'gmxe_initialize' ! name of subroutine
    INTEGER                         :: iou    ! I/O unit
    INTEGER                         :: status ! error status
    INTEGER                         :: jm, jc   


    ! Initialise namelist variables
    L_GASPREC(:,:) = .FALSE.
    L_OXI(:)       = .FALSE.
    
    !--- Read namelist and control variables:

    ! INITIALIZE MAIN-CTRL
    iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
    CALL gmxe_read_nml_ctrl(status, iou)
    IF (status /= 0)  print*, 'error in gmxe_read_nml_ctrl'
    
    l_io=.FALSE.
    IF (loutput) l_io=.TRUE.
!--- Initialize core:
    CALL gmxe_initialize_core(status)
    IF (status /= 0) print*, "error in core initialisation"

    IF (l_aerchem) THEN
      DO jm=1,nmod
        IF ( (JM >= lmode) .AND. (JM <= umode) ) AERCHEM(JM) = .TRUE.
      END DO
      CALL init_aerchem
    END IF

    IF (L_soa) THEN
      IF (L_GASPREC(7,1) .OR. L_GASPREC(8,1)) L_GASPREC(6,1) =.TRUE.
      DO jm=1,nmod
        IF ( (JM >= lmode_soa) .AND. (JM <= umode_soa) ) LSOA(JM) = .TRUE.
      END DO
      call set_soa_params
    ENDIF

    CALL gmxe_initialize_species(0,0)

    SELECT CASE (neqm)
    CASE(1)
       CAll init_species_eqsam4clim(1)    
    CASE(2)
       CAll init_species_isorropia(1)
    END SELECT
   
  END SUBROUTINE gmxe_initialize

  !============================================================================
  SUBROUTINE define_pseudo_tracers

    USE messy_main_tools,         ONLY: strcrack
    USE MESSY_GMXE_OC_AGING,      ONLY: num_wsoc
    USE MESSY_GMXE_SOA,           ONLY: ngas_soa => ngas, nsoa, noxi, mw,     &
                                        names
    USE MESSY_GMXE_AERCHEM_LIQ
    USE MESSY_GMXE_AERCHEM_KPP,   ONLY: str_field_kpp_l => spc_names

    INTEGER :: jm, jc, jt, idx1, dummy
    CHARACTER(LEN=10)           :: cunit, cname
    CHARACTER(LEN=26), POINTER  :: strname(:) => NULL()
    CHARACTER(LEN=26)           :: c_name
    CHARACTER(LEN=2)            :: str_wsoc, str_pa, csubname

    INTEGER                     :: isol, nmedium, nquantity

    nquantity = 1
    cunit='mol/mol'
    nmedium=1

    tcount = 1
    
    DO jm=0,nmod
      nmedium=1
      csubname = ''
      IF (jm > 0) THEN
        csubname=TRIM(cmodes(jm))
        nmedium=2
      END IF
      isol = 1
      IF (JM > nsoluble) isol = 0

      DO jc=1,ngas
        IF (.NOT. td%gas(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%gas(jc)%name)

        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = nmedium
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = jm
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = td%gas(jc)%density*1.e3_dp
        xt(tcount)%R_molarmass       = td%gas(jc)%molmass
 !       IF (jm == 0) THEN
 !         xt(tcount)%R_henry         = td%gas(jc)%henry
 !         xt(tcount)%R_dryreac_sf    = td%gas(jc)%dryreac
 !       ENDIF

        tcount = tcount + 1
      END DO

      DO jc=1,nanions
        IF (jm == 0) CYCLE
        IF (.NOT. td%anion(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%anion(jc)%name)
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = nmedium
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = jm
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = td%anion(jc)%density*1.e3_dp
        xt(tcount)%R_molarmass       = td%anion(jc)%molmass

        tcount = tcount + 1
      END DO

      DO jc=1,ncations
        IF (jm == 0) CYCLE
        IF (.NOT. td%cation(jc)%ltreat(jm)) CYCLE
        IF (jm > nsoluble) CYCLE
        cname = TRIM(td%cation(jc)%name)
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = nmedium
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = jm
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = td%cation(jc)%density*1.e3_dp
        xt(tcount)%R_molarmass       = td%cation(jc)%molmass

        tcount = tcount + 1
      END DO

      DO jc=1,nsolutes
        IF (jm > nsoluble) CYCLE
        IF (.NOT. td%solute(jc)%ltreat(jm)) CYCLE
        cname = TRIM(td%solute(jc)%name)
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = nmedium
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = jm
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = td%solute(jc)%density*1.e3_dp
        xt(tcount)%R_molarmass       = td%solute(jc)%molmass

        tcount = tcount + 1
      END DO

      DO jc=1,nbulk
        IF (jm == 0) CYCLE
        IF (.NOT. bulk(jc)%ltreat(jm)) CYCLE
        cname = TRIM(bulk(jc)%name)
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = nmedium
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = jm
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = bulk(jc)%density*1.e3_dp
        xt(tcount)%R_molarmass       = bulk(jc)%molmass

        tcount = tcount + 1
      END DO

    END DO

    DO jm=1,nmod
      csubname=TRIM(cmodes(jm))
      nmedium=2
      isol = 1
      IF (JM > nsoluble) isol = 0
      nquantity = 2
      cname = "N"
      cunit='1/mol'
      xt(tcount)%name              = TRIM(cname)
      xt(tcount)%subname           = TRIM(csubname)
      xt(tcount)%unit              = TRIM(cunit)
      xt(tcount)%medium            = nmedium
      xt(tcount)%quantity          = nquantity
      xt(tcount)%index             = tcount
      xt(tcount)%S_aerosol_model   = 'gmxe'
      xt(tcount)%I_aerosol_mode    = jm
      xt(tcount)%I_aerosol_sol     = isol
      xt(tcount)%I_aerosol_method  = 2
!      xt(tcount)%R_aerosol_density = 1._dp
      xt(tcount)%R_molarmass       = 1._dp

      tcount = tcount + 1
    END DO

    cunit='mol/mol'
    nquantity = 1
    IF (L_aerchem) THEN

      nspec_aer = 0            ! counter for additional aerosol phase compounds 
                               ! in paerml array
      do jm = 1, nsoluble
        IF (.NOT.AERCHEM(jm) ) CYCLE
        ! use all modes in which aerchem shall be applied
        ! set subname according to the actual modes
        csubname=TRIM(cmodes(jm))
        do jt = 1,lspec_liq
          idx1 = kpp_l_idx%liq_spec(jt, liq_idx)
          ! check if a tracer with the same name already exists for each mode 
          ! treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          IF (TRIM(strname(1)) == "Prod") CYCLE
          c_name=TRIM(strname(1))
          xt(tcount)%name              = TRIM(cname)
          xt(tcount)%subname           = TRIM(csubname)
          xt(tcount)%unit              = TRIM(cunit)
          xt(tcount)%medium            = nmedium
          xt(tcount)%quantity          = nquantity
          xt(tcount)%index             = tcount
          xt(tcount)%S_aerosol_model   = 'gmxe'
          xt(tcount)%I_aerosol_mode    = jm
          xt(tcount)%I_aerosol_sol     = isol
          xt(tcount)%I_aerosol_method  = 2
!          xt(tcount)%R_aerosol_density = kpp_l_idx%liq_attr(jt,liq_dens)*1.e3_dp
          xt(tcount)%R_molarmass       = kpp_l_idx%liq_attr(jt,liq_mw)
          
          nspec_aer = nspec_aer + 1
          tcount = tcount + 1

        END do
      END do
      ! check for gas phase compounds necessary to guarantee mass conservation
      ! for kpp chemistry
      nspec_gas = 0               ! counter for additional gas phase compounds 
      ! in paerml array
      do jt = 1, lspec_gas  
        !      print*, "new tracer gas: ", str_field_kpp_l
        IDX1 = kpp_l_idx%gas_spec(jt,gas_idx)
        cname=TRIM(str_field_kpp_l(idx1))

        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = 1
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = 0
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
        xt(tcount)%R_molarmass       = KPP_L_IDX%GAS_ATTR(JT,GAS_MW)
          
        nspec_gas = nspec_gas + 1
        tcount = tcount + 1
        
      end do
    END IF

    passive:IF (L_PASSIVE_AER) THEN
      mode_loop: DO jm=pamode1,pamode2
        csubname=TRIM(cmodes(jm)) ! mode name, e.g. ki, etc.
        isol = 1
        IF ( jm > nsoluble ) isol = 0 ! if insoluble tracer
        DO jc = 1,num_pa
          IF (jc < 10) then
            write (str_pa,'(A,I1)') "0",jc
          ELSE
            write (str_pa,'(I2)') jc
          END IF
          cname="PASSAER"//str_pa
          xt(tcount)%name              = TRIM(cname)
          xt(tcount)%subname           = TRIM(csubname)
          xt(tcount)%unit              = TRIM(cunit)
          xt(tcount)%medium            = nmedium
          xt(tcount)%quantity          = nquantity
          xt(tcount)%index             = tcount
          xt(tcount)%S_aerosol_model   = 'gmxe'
          xt(tcount)%I_aerosol_mode    = jm
          xt(tcount)%I_aerosol_sol     = isol
          xt(tcount)%I_aerosol_method  = 2
!          xt(tcount)%R_aerosol_density = 1000._dp
          xt(tcount)%R_molarmass       = 1._dp
          
          tcount = tcount + 1

        END DO
      ENDDO mode_loop
    END IF passive

    IF (L_SOA) THEN
     ! gaseous species first
      DO jt=1,NSOA ! gaseous SOA compounds
        cname=TRIM(names(jt))
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = 1
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = 0
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = 0._dp
        xt(tcount)%R_molarmass       = MW(jt)

        tcount = tcount + 1
      END DO
      DO jt=1,NGAS_SOA ! gaseous SOA precursors
        cname=TRIM(names(NSOA + NOXI + jt))
        xt(tcount)%name              = TRIM(cname)
        xt(tcount)%subname           = TRIM(csubname)
        xt(tcount)%unit              = TRIM(cunit)
        xt(tcount)%medium            = 1
        xt(tcount)%quantity          = nquantity
        xt(tcount)%index             = tcount
        xt(tcount)%S_aerosol_model   = 'gmxe'
        xt(tcount)%I_aerosol_mode    = 0
        xt(tcount)%I_aerosol_sol     = isol
        xt(tcount)%I_aerosol_method  = 2
!        xt(tcount)%R_aerosol_density = 0._dp
        xt(tcount)%R_molarmass       = MW(jt)

        tcount = tcount + 1
      END DO
    END IF
    tcount = tcount -1
    
    print*, "number of pseudo tracers: ", tcount
      
  END SUBROUTINE define_pseudo_tracers

!============================================================================  

  SUBROUTINE init_fields

    USE messy_main_tools,        ONLY: find_next_free_unit

    INTEGER           :: iou, jc, check
    REAL(dp)          :: temp, press, rh
    CHARACTER(LEN=14) :: tname
    CHARACTER(LEN=4)  :: tsubname
    REAL(dp)          :: xtval
    

    do jc=1,tcount
      xt(jc)%field=0._dp
    END do
    iou = find_next_free_unit(100,200)
    open(iou, file='gmxe_input.txt')
    read(iou,*,IOSTAT=check) temp_in, press_in, rh_in
    IF (check<0) print*, "No temperature data - file not found!"
    DO
      tname=''
      tsubname=''
      xtval = 0._dp
      read(iou,fmt='(A14,A4,E14.4)',IOSTAT=check) tname, tsubname, xtval
      IF (check<0) EXIT
      DO jc=1,tcount
        IF (TRIM(tname) == TRIM(xt(jc)%name)) THEN
          IF (TRIM(ADJUSTL(tsubname)) == TRIM(xt(jc)%subname)) THEN
            xt(jc)%field = xtval
          ENDIF
        ENDIF
      END DO
    ENDDO


    print*, "Initialisation final status!"
    DO jc=1,tcount
      IF (TRIM(xt(jc)%subname) == '') THEN
        print*, TRIM(xt(jc)%name), xt(jc)%field
      ELSE
        print*, TRIM(xt(jc)%name),"_", TRIM(xt(jc)%subname), xt(jc)%field
      END IF
    ENDDO
    
  END SUBROUTINE init_fields

!============================================================================  
  SUBROUTINE gmxe_box_init_mem

    ! This subroutine sets up memory structures for in/output variables

!    ALLOCATE(temperature(1,1))
!    ALLOCATE(pressure(1,1))
!    ALLOCATE(rhumidity(1,1))

    ALLOCATE(rdryaer(nmod))
    ALLOCATE(rwetaer(nmod))
    ALLOCATE(ddryaer(nmod))
    ALLOCATE(aernumb(nmod))

    ALLOCATE(diagaer(0:naerodiag,nlowermode:nuppermode))

!================================================================  
! Design a netcdf output-file with all variables with time dependence

    
    
    
  END SUBROUTINE gmxe_box_init_mem
  
!============================================================================  
  SUBROUTINE gmxe_box_init_cpl

    USE messy_gmxe_aerchem_liq,   ONLY: jval_h2o2, jval_no3, jval_o3
    USE messy_gmxe_kappa,         ONLY: initialise_salts
    USE messy_gmxe_isorropia2,    ONLY: init_species_isorropia
    USE messy_gmxe_eqsam4clim,    ONLY: init_species_eqsam4clim

    IMPLICIT NONE
    INTEGER :: jt, idx1, idx2, jm
    LOGICAL :: lfound

    lfound = .false.
    do jt=1,tcount
      IF (TRIM(xt(jt)%name) == "SO4mm") THEN
        IF (TRIM(xt(jt)%subname) == cmodes(1)) THEN
          lfound = .true.
          EXIT
        END IF
      END IF
    end do

    If (.not.lfound) THEN
      print*, "No nucleation mode sulphate tracer available!"
      print*, "Therefore, nucleation is switched off - LNUCL = .FALSE."
      lnucl = .FALSE.
    END If

    CALL INIT_CPL_SPECIES_E5
    CALL initialise_salts(1,1,nmod)
    SELECT CASE (neqm)
    CASE(1)
       CALL init_species_eqsam4clim(2)
    CASE(2)
       CALL init_species_isorropia(2)
    END SELECT
    print*, 'Species structure finished.'

    ALLOCATE(zconvert (1,1,0:naertot)); zconvert(:,:,:) = FACp0
    ALLOCATE(xconvert (1,1,0:naertot)); xconvert(:,:,:) = FACp0
    ALLOCATE(lrhoa    (0:naertot))            ; lrhoa       (:) = .TRUE.

    DO jt=1,spec_number
      DO jm = 0,nmod
        
        idx1 = species(jt)%tracidx(jm)
        IF (idx1 == 0) cycle
        idx2 = species(jt)%aermlidx(jm)
        lrhoa(idx2) =.true.
        
        !--- Assign conversion factors depending on the unit of the tracer:
        
        
        SELECT CASE (TRIM(xt(idx1)%unit))
          
          !--- Particle mass:
        CASE ('kg/kg')
          ! [kg kg^-1(air)]   => [umol m^-3 (air)]
          xconvert(1,1,idx2) = FACp9    & 
            / xt(idx1)%R_molarmass
          
        CASE ('mol/mol')
          !  [mol mol^-1 (air)] => [umol g^-1 (air)]
          xconvert(1,1,idx2) = FACp9 / M_air 
          
          !--- Particle number:
          
        CASE ('1/kg')
          !  [# kg^-1 (air)]  => [# mg^- (air)]
          xconvert(1,1,idx2) = FACm6 
          
        CASE ('1/mol')
          !  [# mol^-1 (air)] => [# cm^-3 (air)]
          xconvert(1,1,idx2) = FACm3 / M_air
          
        CASE ('1/cm3')
          !  [# cm^-3 (air)] => [# cm^-3 (air)]
          lrhoa(idx2)    = .false.
          xconvert(1,1,idx2) = FACp0
          
        END SELECT
      END DO
    END DO
    

    
    
  END SUBROUTINE gmxe_box_init_cpl
!============================================================================  
  SUBROUTINE GMXE_box_driver

    USE messy_main_constants_mem,   ONLY: rd, avo=>N_a
    
    INTEGER, PARAMETER  :: kproma = 1
    INTEGER, PARAMETER  :: nlev   = 1
    REAL(dp), PARAMETER :: zero = 0._dp

    REAL(dp) :: time_step_len = 100.  ! in seconds
    
    REAL(dp) :: ztemp(1,1), zpress(1,1), zrhum(1,1)

! dummy for pseudo cloud coupling (by Swen Metzger)
    
    REAL(dp) :: zsat(1,1), zaclc(1,1)
    REAL(dp) :: zaopt(1,1)
    REAL(dp) :: zcih2o(1,1), zclh2o(1,1)
    REAL(dp) :: zcdnc(1,1), zicnc(1,1), zcrain(1,1)
    REAL(dp) :: zcph(1,1,nmod), zccn(1,1,nmod)

! concentrations    
    REAL(dp) :: zaerml(1,1,0:naertot)
    REAL(dp) :: zxtm(1,1,0:naertot)
    REAL(dp) :: zaernl(1,1,nmod)
    REAL(dp) :: zaerosol(1,1,0:nmod,0:naerodiag)

! additional parameters    
    REAL(dp) :: zconv
    REAL(dp) :: zvap(1,1)
    REAL(dp) :: zvap2(1,1)
    REAL(dp) :: zprod(1,1)
    REAL(dp) :: zminc(1,1)
    REAL(dp) :: zrhoa(1,1)
    REAL(dp) :: zpressi(1,2)
    REAL(dp) :: vervel(1,1)

    REAL(dp) :: zrdry(1,1,nmod)    
    REAL(dp) :: zrwet(1,1,nmod)    
    REAL(dp) :: zddry(1,1,nmod)

    REAL(dp) :: time, end_time

    INTEGER  :: jm, jt, idx1, idx2, idt

    REAL(dp) :: kmod(nmod)
    LOGICAL  :: end_reached
    ! time loop

    TIME = 0._dp
    do jt=1,nmod
      kmod(jt)=REAL(jt,dp)
    enddo
    CALL create_outputfile_netcdf(nmod,kmod,tcount)
    
    END_TIME = 10._dp * 86400._dp

    ztemp  = temp_in
    zpress = press_in
    zrhoa  = zpress / (Rd * ztemp)
    zrhum  = rh_in

    zpressi(1,1) = zpress(1,1) - 25000._dp
    zpressi(1,2) = zpress(1,1) + 25000._dp

!    zdpress     (:,:) = zero
!    zshum       (:,:) = zero
    zsat        (:,:) = zero
!    zgh2o       (:,:) = zero
    zclh2o      (:,:) = zero
    zcih2o      (:,:) = zero
    zaclc       (:,:) = zero
    zaopt       (:,:) = zero
    zcrain      (:,:) = zero
    zcdnc       (:,:) = zero
    zicnc       (:,:) = zero
    zcph      (:,:,:) = NaN
    zccn      (:,:,:) = zero
    zaerml    (:,:,:) = zero
    zaernl    (:,:,:) = zero
    zxtm      (:,:,:) = zero
    zaerosol(:,:,:,:) = ZERO

    zrdry = 0._dp
    zrwet = 0._dp
    zddry = 0._dp

    zconv    = 1._dp / (avo * 1.e-6_dp * 1.e-6_dp)
    
    do jt=1,spec_number
      do jm=0,nmod
        idx1 = species(jt)%tracidx(jm)
        idx2 = species(jt)%aermlidx(jm)
        IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
        zconvert(1,1,idx2) = xconvert (1,1,idx2)
        IF(lrhoa(idx2)) &
          zconvert(1,1,idx2) = xconvert (1,1,idx2) * zrhoa(1,1)
        idt = idx1
        zxtm(1:kproma,1:nlev,idx2) = MAX(0._dp, xt(idt)%field)

        
        zaerml(1:kproma,1:nlev,idx2) = zxtm(1:kproma,1:nlev,idx2)       &  
                                     ! tracer mass       [mol  mol-1]
                                     * zconvert(1:kproma,1:nlev,idx2)  
                                     ! conversion factor [umol m-3 (air)]

    
        IF(.NOT. lnumber) THEN
          idx2 = species(jt)%aernlidx(jm)
          IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
          zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2)
          IF(lrhoa(idx2)) &
            zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2) &
                                           * zrhoa(1:kproma,1:nlev)
          idt = idx1
          
          zxtm(1:kproma,1:nlev,idx2) = MAX(0._dp, xt(idt)%field)

          zaernl(1:kproma,1:nlev,jm) = zxtm(1:kproma,1:nlev,idx2)        &
                                     ! tracer number       [#  mol-1]
                                     * zconvert(1:kproma,1:nlev,idx2)    
          zaerml(1:kproma,1:nlev,idx2) = 0._dp
        END IF
      END DO
    END DO
    
    
    DO
      end_reached = time > end_time
      IF (end_reached) EXIT

      CALL gmxe_main(1,          1,           time_step_len,                &
                 zrhum,          ztemp,          zpress,          zsat,     &
                 zaclc,          zaopt,          zcih2o,          zclh2o,   &
                 zcdnc,          zicnc,          zcph,            zcrain,   &
                 zccn,           zaerml,         zaernl,          zaerosol, &
                 zconv,          zvap,           zvap2,           zprod,    &
                 zminc,          zrhoa*FACp3,    zpressi,         vervel,   &
                 zrdry,          zrwet,          zddry )


      do jt=1,spec_number
        do jm=0,nmod
          idx1 = species(jt)%tracidx(jm)
          idx2 = species(jt)%aermlidx(jm)
          IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
          zconvert(1,1,idx2) = xconvert (1,1,idx2)
          IF(lrhoa(idx2)) &
            zconvert(1,1,idx2) = xconvert (1,1,idx2) * zrhoa(1,1)
          idt = idx1
        
          zxtm(1:kproma,1:nlev,idx2) = zaerml(1:kproma,1:nlev,idx2) &  
                                     ! tracer mass       [mol  mol-1]
                                     / zconvert(1:kproma,1:nlev,idx2)  
                                     ! conversion factor [umol m-3 (air)]
    
          xt(idt)%field = MAX(0._dp,zxtm(1,1,idx2))
          
          IF(.NOT. lnumber) THEN
            idx2 = species(jt)%aernlidx(jm)
            IF ( (idx1 == 0) .or. (idx2 == 0) ) CYCLE
            zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2)
            IF(lrhoa(idx2)) &
              zconvert(1:kproma,1:nlev,idx2) = xconvert (1:kproma,1:nlev,idx2) &
              * zrhoa(1:kproma,1:nlev)
            idt = idx1
            
            zxtm(1:kproma,1:nlev,idx2) = zaernl(1:kproma,1:nlev,jm) &  
              ! tracer mass       [mol  mol-1]
              / zconvert(1:kproma,1:nlev,idx2)  
            ! conversion factor [umol m-3 (air)]
            
            xt(idt)%field = MAX(0._dp,zxtm(1,1,idx2))
            
          END IF
        END DO
      END DO
      
      call write_output_netcdf(nmod,time,zrdry(1,1,:),zrwet(1,1,:),zddry(1,1,:), &
        ztemp(1,1), zpress(1,1),zrhum(1,1),END_reached)
      
      
      TIME = time + time_step_len
      
      
    ENDDO
    
  END SUBROUTINE GMXE_box_driver
!============================================================================  
  SUBROUTINE INIT_CPL_SPECIES_E5

    USE MESSY_MAIN_TOOLS,              ONLY: strcrack
    USE MESSY_GMXE_AERCHEM_LIQ
    USE MESSY_GMXE_OC_AGING,           ONLY: num_wsoc, kappa, init_oc_aging, &
                                             spec_idx_wsoc
    USE MESSY_GMXE_SOA,                ONLY: nsoa, noxi, ngas_soa => ngas, names
    USE MESSY_GMXE_AERCHEM_KPP,        ONLY: str_field_kpp_l => spc_names

    INTEGER :: jm ,jc, jt, jk !mz_sb_20161020: jk added
    INTEGER :: counter
    LOGICAL :: FOUND=.FALSE.

    CHARACTER(LEN=STRLEN_MEDIUM) :: names_array(NMAX_SPEC)
    CHARACTER(LEN=STRLEN_MEDIUM) :: cname

    INTEGER                      :: dummy, idx1
    CHARACTER(LEN=26), POINTER   :: strname(:) => NULL()
    CHARACTER(LEN=2)             :: str_wsoc
    !mz_sb_20161020+
    CHARACTER(LEN=2)             :: str_soa
    CHARACTER(LEN=2)             :: str_mom
    !mz_sb_20161020-

    CHARACTER(LEN=2)             :: str_pa ! mz_dk_20120120
    INTEGER                      :: spec_idx_pa(1:num_pa) ! mz_dk_20120120
!mz_sb_20160211+
!    INTEGER                      :: spec_idx_fPOA(NfPOA),     spec_idx_bbPOA(NbbPOA),      &
!                                    spec_idx_fSOAsv(NfSOAsv), spec_idx_bbSOAsv(NbbSOAsv),  &
!                                    spec_idx_fSOAiv(NfSOAiv), spec_idx_bbSOAiv(NbbSOAiv),  &
!                                    spec_idx_SOAv(NSOAv,maxMOM) 
!mz_sb_20160211-

    cname          = ""
    names_array(:) = ""

    counter = 0
      ! thermodynamic / inorganic partitioning compounds
    DO jc = 1, ngas
      cname = TRIM(td%gas(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,nanions
      cname = TRIM(td%anion(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,ncations
      cname = TRIM(td%cation(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO

    DO jc=1,nsolutes
      cname = TRIM(td%solute(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
    DO jc=1,nbulk
      cname = TRIM(bulk(jc)%name)
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
          found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
    END DO
    
! This one is for numbers per mode
    IF (.NOT. lnumber) THEN
      DO jm=1,nmod
        cname = 'N'
        DO jt=1,counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) &
            found=.true.
        END DO
        IF (.NOT.FOUND) THEN
          counter = counter + 1
          names_array(counter) = cname
        END IF
      END DO
    END IF
      
    IF (l_aerchem) THEN
      ! additional species from aerchem (aerosol phase)
      DO jm=1,nsoluble
        IF (.NOT.AERCHEM(jm)) CYCLE
        
        DO jc = 1, lspec_liq
          idx1 = kpp_l_idx%liq_spec(jc, liq_idx)
          ! check if a tracer with the same name already exists for each mode 
          ! treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          IF (TRIM(strname(1)) == 'Prod') CYCLE
          cname=TRIM(strname(1))
          found=.false.
          DO jt=1,counter
            IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
          END DO
          IF (.NOT.FOUND) THEN
            counter = counter + 1
            names_array(counter) = cname
          END IF
        END DO
      END DO
      
      ! additional species from aerchem (gas phase)
      jm = 0
      
      DO jc = 1, lspec_gas
        idx1 = kpp_l_idx%gas_spec(jc, gas_idx)
        ! check if a tracer with the same name already exists for each mode 
        ! treated by aerchem
        if (associated(strname)) DEALLOCATE (strname)
        NULLIFY(strname) 
        CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
        cname=TRIM(strname(1))
        found=.false.
        DO jt=1,counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.FOUND) THEN
          counter = counter + 1
          names_array(counter) = cname
        END IF
      END DO
      
    END IF ! aerchem

    IF (L_OC_AGING) THEN
      ! get gas phase OH concentration into the species array
      found = .false.
      cname = "OH"
      DO jt=1,counter
        IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
      END DO
      IF (.NOT.FOUND) THEN
        counter = counter + 1
        names_array(counter) = cname
      END IF
       
    ENDIF ! oc_aging
    
    ! mz_dk_20120119+
    ! PASSIVE AEROSOLS
    ! set names_array and counter
    IF (L_PASSIVE_AER) THEN
      DO jm = 1, nmod
        DO jc = 1, num_pa
          IF (jc < 10) then
            write (str_pa,'(A,I1)') "0",jc
          ELSE
            write (str_pa,'(I2)') jc
          END IF
          found=.false.
          cname = "PASSAER"//str_pa
          DO jt = 1, counter
            IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
          END DO
          IF (.NOT.found) THEN
            counter = counter+1
            names_array(counter) = cname
          END IF
        END DO
      END DO
    END IF  ! pass_aer
    ! mz_dk_20120119-

    IF (L_SOA) THEN
      DO jc=1,ngas_soa
        found=.false.
        cname = TRIM(names(NSOA+NOXI+jc))
        DO jt = 1, counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.found) THEN
          counter = counter+1
          names_array(counter) = cname
        END IF
      END DO
      DO jc=1,noxi
        found=.false.
        cname = TRIM(names(NSOA+jc))
        DO jt = 1, counter
          IF ( TRIM(names_array(jt)) == TRIM(cname) ) found=.true.
        END DO
        IF (.NOT.found) THEN
          counter = counter+1
          names_array(counter) = cname
        END IF
      END DO
       ! aerosol species should be already included in the bulk structure
!!$       DO jm=1,nmod
!!$          IF (.NOT.LSOA(jm)) CYCLE
!!$          DO jc=1,nsoa
!!$             cname = TRIM(names(jc))
!!$             found=.false.
!!$             DO jt = 1, counter
!!$                IF ( TRIM(names_array(jt)) == TRIM(cname)) found=.true.
!!$             END DO
!!$             IF (.NOT.found) THEN
!!$                counter = counter+1
!!$                names_array(counter) = cname
!!$             END IF
!!$          END DO
!!$       END DO
    END IF ! SOA

    ! allocate the species index field
    spec_number = counter
    ALLOCATE(SPECIES(spec_number))
    DO jt=1,spec_number
      ALLOCATE(SPECIES(jt)%tracidx(0:nmod))
      ALLOCATE(SPECIES(jt)%aermlidx(0:nmod))
      ALLOCATE(SPECIES(jt)%aernlidx(0:nmod))
      ALLOCATE(SPECIES(jt)%kppidx(0:nmod))
      ALLOCATE(SPECIES(jt)%soaidx(0:nmod))
      ALLOCATE(SPECIES(jt)%zcoup(0:nmod))
      ALLOCATE(SPECIES(jt)%npassive(0:nmod)) ! mz_dk_20120120 PASSIVE TRACER or NOT?
      SPECIES(jt)%name = names_array(jt)
      DO jm=0,nmod
        SPECIES(jt)%tracidx(jm)  = 0
        SPECIES(jt)%aermlidx(jm) = 0
        SPECIES(jt)%aernlidx(jm) = 0
        SPECIES(jt)%kppidx(jm)   = 0
        SPECIES(jt)%soaidx(jm)   = 0
        SPECIES(jt)%zcoup(jm)    = 1.0_dp
        SPECIES(jt)%molmass      = 0.0_dp
!        SPECIES(jt)%density      = 0.0_dp
        SPECIES(jt)%kappa        = 0.0_dp
        SPECIES(jt)%charge       = 0.0_dp
        SPECIES(jt)%npassive(jm) = .false. ! mz_dk_20120120
      END DO
    END DO


    ! determine tracer indices
    DO jt=1,tcount
      IF (xt(jt)%medium == 2) THEN
        IF (TRIM(xt(jt)%S_aerosol_model) == 'gmxe') THEN
          jm = xt(jt)%I_aerosol_mode
        ELSE
          CYCLE
        ENDIF
      ENDIF

      IF (xt(jt)%medium == 1) &
        jm = 0

      IF (xt(jt)%medium == 3) CYCLE
        
      DO jc=1,spec_number
        IF ( TRIM(xt(jt)%name) == TRIM(species(jc)%name) ) THEN
          SPECIES(jc)%tracidx(jm) = jt
          SPECIES(jc)%molmass     = xt(jt)%R_molarmass
!          IF (jm > 0) &
!            SPECIES(jc)%density   = xt(jt)%R_aerosol_density / 1.e3_dp
        END IF
      ENDDO
    ENDDO

    ! determine paerml indices and coupling switches

    counter = 0
    DO jm = 0, nmod
      DO jc=1, ngas
        IF (td%gas(jc)%ltreat(jm)) THEN
          cname = TRIM(td%gas(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == TRIM(species(jt)%name)) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%gas(jc)%aerml_idx(jm) = counter
              IF (lgas) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, nanions
        IF (jm == 0) CYCLE
        IF (td%anion(jc)%ltreat(jm)) THEN
          cname = TRIM(td%anion(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%anion(jc)%aerml_idx(jm) = counter
              species(jt)%charge = td%anion(jc)%charge
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, ncations
        IF (jm == 0) CYCLE
        IF (td%cation(jc)%ltreat(jm)) THEN
          cname = TRIM(td%cation(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%cation(jc)%aerml_idx(jm) = counter
              species(jt)%charge = td%cation(jc)%charge
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
       DO jc=1, nsolutes
        IF (td%solute(jc)%ltreat(jm)) THEN
          cname = TRIM(td%solute(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              td%solute(jc)%aerml_idx(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END DO
      DO jc=1, nbulk
        IF (jm == 0) CYCLE
        IF (bulk(jc)%ltreat(jm)) THEN
          cname = TRIM(bulk(jc)%name)
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              bulk(jc)%aerml_idx(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
              species(jt)%kappa = bulk(jc)%kappa
            END IF
          END DO
        END IF
      END DO
      IF (.NOT. lnumber) THEN
        IF (jm /= 0) THEN
          cname='N'
          DO jt = 1,spec_number
            IF (TRIM(cname) == species(jt)%name) THEN
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              species(jt)%aermlidx(jm) = counter
              species(jt)%aernlidx(jm) = counter
              nnum(jm) = counter
              IF (laerosol) &
                species(jt)%zcoup(jm) = 1._dp
            END IF
          END DO
        END IF
      END IF

    END DO
    
    counter = 0
    do jt=1,spec_number
      DO jm=0,nmod
        counter = MAX(counter, SPECIES(jt)%aermlidx(jm))
      ENDDO
    ENDDO
    
    IF ( l_aerchem ) THEN

      ! aerchem liquid compounds
      
      DO jm= 1,nsoluble
        IF (.NOT. AERCHEM(jm)) CYCLE
        DO jc = 1,lspec_liq
          idx1 = kpp_l_idx%liq_spec(jc, liq_idx)
        ! check if a compound with the same name already exists for each mode 
        ! within aerml array treated by aerchem
          if (associated(strname)) DEALLOCATE (strname)
          NULLIFY(strname) 
          CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
          cname=TRIM(strname(1))
          DO jt=1,spec_number
            IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
              SPECIES(jt)%kppidx(jm)   = idx1
              IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
              counter = counter + 1
              SPECIES(jt)%aermlidx(jm) = counter
            ENDIF
          ENDDO
        ENDDO

      ENDDO   ! aerchem modes

        ! aerchem gas compounds

      jm = 0
      DO jc = 1,lspec_gas
        idx1 = kpp_l_idx%gas_spec(jc, gas_idx)
        ! check if a compound with the same name already exists for each mode 
        ! within aerml array treated by aerchem
        if (associated(strname)) DEALLOCATE (strname)
        NULLIFY(strname) 
        CALL strcrack(STR_FIELD_KPP_L(idx1),'_', strname, dummy)
        cname=TRIM(strname(1))
        DO jt=1,spec_number
          IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
            SPECIES(jt)%kppidx(jm)   = idx1
            IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
            counter = counter + 1
            SPECIES(jt)%aermlidx(jm) = counter
          ENDIF
        ENDDO
      ENDDO
      
    ENDIF

    IF (L_OC_AGING) THEN
       ! WSOC compounds treated already in the bulk structure
       ! gas phase OH
       jm = 0
       cname="OH"
       DO jt=1,spec_number
          IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
             IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
             IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
             counter = counter + 1
             SPECIES(jt)%aermlidx(jm) = counter
             SPECIES(JT)%zcoup(jm) = 0._dp
          ENDIF
       END DO
    END IF
    
    ! mz_dk_20120119+
    IF (L_PASSIVE_AER) THEN
       ! passive aerosol tracers
       ! set aerosol mass index
       DO jm = 1, nmod
          DO jc = 1, num_pa
             IF (jc < 10) then
                write (str_pa,'(A,I1)') "0",jc
             ELSE
                write (str_pa,'(I2)') jc
             END IF
             cname = "PASSAER"//str_pa
             DO jt = 1,spec_number
                IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                   IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                   IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                   counter = counter + 1
                   SPECIES(jt)%aermlidx(jm) = counter
                END IF
             END DO
          END DO
       END DO
    END IF
    ! mz_dk_20120119-

    IF (L_SOA) THEN
       ! gas phase first
       jm = 0
       DO jc=1,NSOA+NOXI+NGAS_SOA
          cname = TRIM(names(jc))
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%npassive(jm) = .false.
                IF (SPECIES(jt)%tracidx(jm) == 0) CYCLE
                IF (SPECIES(jt)%aermlidx(jm) /= 0) CYCLE
                counter = counter + 1
                SPECIES(jt)%aermlidx(jm) = counter
                species(jt)%SOAIDX(jm) = jc
             END IF
          END DO
       END DO
       ! aerosols treated already in the bulk structure
       ! aermlidx is therefore already set, but SOAIDX 
       ! must still be determined
       DO jm = 1,nmod
         IF (.NOT.LSOA(jm)) CYCLE
         DO jc=1,NSOA
           cname = TRIM(names(jc))
           DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) &
               species(jt)%SOAIDX(jm) = jc
           END DO
         END DO
       END DO
         
!!$       DO jt=1,spec_number
!!$          DO jm=0,nmod
!!$             IF (species(jt)%SOAIDX(jm) == 0 ) CYCLE
!!$             IF (p_pe /= p_io) CYCLE
!!$             print*, "names, soa", species(jt)%name, jt, species(jt)%aermlidx(jm), &
!!$                  species(jt)%SOAIDX(jm), jm, species(jt)%tracidx(jm), &
!!$                  species(jt)%npassive(jm)
!!$          END DO
!!$       END DO
    ENDIF ! L_SOA

    DO jt=1,spec_number
      SELECT CASE (TRIM(species(jt)%name))
      CASE("SS")
!        spec_idx_ss = jt
      CASE("OC")
        spec_idx_oc = jt
      CASE("DU")
!        spec_idx_du = jt     
      CASE("BC")
!        spec_idx_bc = jt       
      CASE("Hp")   
        spec_idx_hp = jt
      CASE("OHm")   
        spec_idx_ohm = jt
      CASE("HSO4m")  
        spec_idx_hso4m = jt
      CASE("H2SO4")
        spec_idx_h2so4 = jt
      CASE("SO4mm")
        spec_idx_so4mm = jt
      CASE("NH3")
        spec_idx_nh3 = jt
      CASE("NH4p")
        spec_idx_nh4p = jt
        
      CASE("H2O")
        spec_idx_h2o = jt
        species(jt)%molmass = mwh2o
      CASE("OH")
        spec_idx_oh = jt
      CASE("Number", "N")
        number_idx = jt
      END SELECT
    END DO

    IF (L_OC_AGING) THEN
       jt = spec_idx_oc
       species(jt)%kappa = 0.01_dp
       DO jc=1,nbulk
         IF ( TRIM(species(jt)%name) == TRIM(bulk(jc)%name) ) &
           bulk(jc)%kappa = 0.01_dp
       END DO
       spec_idx_wsoc(0) = spec_idx_oc
       CALL init_oc_aging

       DO jc=1,num_wsoc
         IF (jc < 10) then
           write (str_wsoc,'(A,I1)') "0",jc
         ELSE
           write (str_wsoc,'(I2)') jc
         END IF
         cname = "WSOC"//str_wsoc
         DO jt=1,spec_number
           IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
             species(jt)%kappa = kappa(jc)
             spec_idx_wsoc(jc) = jt
           END IF
         END DO
         DO jt=1,nbulk
           IF (TRIM(bulk(jt)%name) == TRIM(cname) ) &
             bulk(jt)%kappa = kappa(jc)
         END DO
       END DO
    END IF
    
    ! mz_dk_2012019+
    IF (L_PASSIVE_AER) THEN
       ! set kappa, index, and npassive flag
       DO jc=1,num_pa
          IF (jc < 10) then
             write (str_pa,'(A,I1)') "0",jc
          ELSE
             write (str_pa,'(I2)') jc
          END IF
          cname = "PASSAER"//str_pa
          DO jt=1,spec_number
             IF ( TRIM(cname) == TRIM(species(jt)%name) ) THEN
                species(jt)%kappa = 0._dp
                spec_idx_pa(jc) = jt
                species(jt)%npassive(1:nmod)= .true.
                !print*,species(jt)%name,jt,species(jt)%npassive(1:nmod)
             END IF
          END DO
       END DO
    END IF
    ! mz_dk_2012019-
    naertot = counter


    DO jc=1,spec_number
      DO jm=0,nmod
        jt = species(jc)%aermlidx(jm)
        SELECT CASE (TRIM(species(jc)%name))
        CASE ("H2O")
          nwh2o(jm) = jt
        CASE ("H2SO4")
          nh2so4(jm) = jt
        CASE ("HSO4m")
          nhso4m(jm) = jt
        CASE ("SO4mm")  
          nso42m(jm) = jt
        CASE ("Hp")  
          nhp(jm) = jt
        CASE ("OHm")  
          nohm(jm) = jt
        END SELECT
      END DO
    END DO
    
!    counter = 0
!    DO jt=1,spec_number
!      DO jm=0,nmod
!        IF (SPECIES(jt)%tracidx(jm) /= 0) counter = counter + 1
!      END DO
!    END DO
!
!    print*, "end counter: ",counter, naertot

    do jt=1,spec_number
      print*, "name: ",TRIM(species(jt)%name)
      DO jm=0,nmod
        print*, "mode: ", jm, "tracer_index: ", SPECIES(jt)%tracidx(jm),  &
          "aerml_index: ",  SPECIES(jt)%aermlidx(jm), &
          "kpp_index: ",  SPECIES(jt)%kppidx(jm),     &
!          "density: ",  SPECIES(jt)%density,          &
          "molarmass: ",  SPECIES(jt)%molmass,        &
          "kappa: ",  SPECIES(jt)%kappa,              &
          "npassive: ", species(jt)%npassive(jm)
      ENDDO
    ENDDO


  END SUBROUTINE INIT_CPL_SPECIES_E5

!==============================================================================
  SUBROUTINE create_outputfile_netcdf(ndim,kmod,tcount)! nstep,nstep_fin,t,dt)
        
    use netcdf

    IMPLICIT NONE
    INTEGER :: ndim, tcount
    
    REAL(dp), dimension(ndim)             :: kmod
    
    CHARACTER(LEN=35)  :: vname


    integer :: i,j,imin,imax,jmin,jmax,kmin,kmax,istep_rep
    integer :: nx,ny,nz,istep,isto,nstep,nstep_fin
    real*8 :: t,dt
    !        real*8, dimension(1) :: tt
   
    
    INTEGER :: id_ndim, id_rec, jt
    
    INTEGER, DIMENSION(:), POINTER, SAVE :: jddim
    
    ! open output file
    call check(nf90_create("GMXE_box_out.nc",nf90_clobber,idfile))
    ! set a global attribute
    call check(nf90_put_att(idfile,nf90_global,'title','GMXE-Box-Model'))
    
    ! define dimensions
    CALL check(nf90_def_dim(idfile,'time',nf90_unlimited, id_rec))
    
    CALL check(nf90_def_dim(idfile,'ndim',ndim, id_ndim))
    
    ! define variables for axes
    NULLIFY(jddim)
    ALLOCATE(jddim(1))
    ALLOCATE(idvar(8+tcount))
    jddim(1) = id_ndim
    CALL check(nf90_def_var(idfile,'ndim',nf90_float,jddim,idvar(1)))
    CALL check(nf90_put_att(idfile, idvar(1), "longname", "number of modes") )
    CALL check(nf90_put_att(idfile, idvar(1), "units", "-") )
    jddim(1) = id_rec
    CALL check(nf90_def_var(idfile,'time',nf90_float,jddim,idvar(2)))
    CALL check(nf90_put_att(idfile, idvar(2), "longname", &
      "Time") )
    CALL check(nf90_put_att(idfile, idvar(2), "units", "seconds" ) )
    ! define variable for output
    DEALLOCATE(jddim)
    NULLIFY(jddim)
    
    ALLOCATE(jddim(2))
    jddim(1) = id_ndim
    jddim(2) = id_rec
    
    CALL check(nf90_def_var(idfile, "dryrad", nf90_float, jddim, idvar(3)))
    CALL check(nf90_put_att(idfile, idvar(3), "longname", &
      "aerosol dry radius") )
    CALL check(nf90_put_att(idfile, idvar(3), "units", "m") )
    
    CALL check(nf90_def_var(idfile, "wetrad", nf90_float, jddim, idvar(4)))
    CALL check(nf90_put_att(idfile, idvar(4), "longname", &
      "aerosl ambient (wet) radius") )
    CALL check(nf90_put_att(idfile, idvar(4), "units", "m") )
    
    CALL check(nf90_def_var(idfile, "drydens", nf90_float, jddim, idvar(5)))
    CALL check(nf90_put_att(idfile, idvar(5), "longname", &
      "aerosol dry density") )
    CALL check(nf90_put_att(idfile, idvar(5), "units", "kg/m^⁻3") )
    
    DEALLOCATE(jddim)
    NULLIFY(jddim)
    
    ALLOCATE(jddim(1))
    jddim(1) = id_rec
    
    CALL check(nf90_def_var(idfile, "temperature", nf90_float, jddim, idvar(6)))
    CALL check(nf90_put_att(idfile, idvar(6), "longname", &
      "temperature") )
    CALL check(nf90_put_att(idfile, idvar(6), "units", "K") )
    
    CALL check(nf90_def_var(idfile, "pressure", nf90_float, jddim, idvar(7)))
    CALL check(nf90_put_att(idfile, idvar(7), "longname", &
      "pressure") )
    CALL check(nf90_put_att(idfile, idvar(7), "units", "Pa") )
    
    CALL check(nf90_def_var(idfile, "rhum", nf90_float, jddim, idvar(8)))
    CALL check(nf90_put_att(idfile, idvar(8), "longname", &
      "relative humidity") )
    CALL check(nf90_put_att(idfile, idvar(8), "units", "1/s") )
    
    print*, "creating tracer output list!"
    DO jt=1,tcount

      IF (xt(jt)%medium == 1) THEN
        vname = TRIM(xt(jt)%name)
      ELSE
        vname = TRIM(xt(jt)%name)//"_"//TRIM(xt(jt)%subname)
      ENDIF
      CALL check(nf90_def_var(idfile, TRIM(vname), nf90_float, jddim, idvar(8+jt)))
      CALL check(nf90_put_att(idfile, idvar(8+jt), "longname", &
        TRIM(vname)//" mixing ratio") )
      CALL check(nf90_put_att(idfile, idvar(8+jt), "units", TRIM(xt(jt)%unit)) )
      
    END DO
    
    DEALLOCATE(jddim)
    NULLIFY(jddim)

    ALLOCATE(jddim(1))
    jddim(1) = id_ndim

    CALL check(nf90_def_var(idfile,"crdiv_mid", nf90_float,jddim,idvar(8+tcount+1)))
    CALL check(nf90_put_att(idfile, idvar(8+tcount+1), "longname", "size of bin"))
    CALL check(nf90_put_att(idfile, idvar(8+tcount+1), "units", "m"))
    DEALLOCATE(jddim)
    NULLIFY(jddim)
    ! end define mode
    CALL check(NF90_ENDDEF(idfile))
    
    ! write axis variables
    CALL check(nf90_put_var(idfile, idvar(1), kmod))
    call check(nf90_put_var(idfile,idvar(9+tcount), 1.e-2_dp * crdiv_mid))
    it = 0
    
    print*, "successfully created outputfile variables"
  END SUBROUTINE create_outputfile_netcdf

      !--------------------------------------------------------------------------
  SUBROUTINE write_output_netcdf(ndim, time, rdryrad,rwetrad,ddryden, &
                                     temp, pres, rh, end_reached)

    use netcdf
    INTEGER  :: ndim, jt
    REAL(dp) :: time
    REAL(dp), DIMENSION(ndim) :: rdryrad, rwetrad, ddryden
    REAL(dp) :: temp, pres, rh
    logical  :: end_reached
    INTEGER :: nstt(2), ncnt(2)
    it = it + 1
    CALL check(nf90_put_var(idfile, idvar(2), time, start=(/it/)))
    
    ! write ouput variable
    
    nstt = (/ 1,    it/)
    ncnt = (/ ndim, 1 /)
    
    CALL check(nf90_put_var(idfile, idvar(3), rdryrad*1.e-2_dp, nstt, ncnt))
    CALL check(nf90_put_var(idfile, idvar(4), rwetrad*1.e-2_dp, nstt, ncnt))
    CALL check(nf90_put_var(idfile, idvar(5), ddryden, nstt, ncnt))
    
    CALL check(nf90_put_var(idfile, idvar(6), temp,start=(/it/)))
    CALL check(nf90_put_var(idfile, idvar(7), pres,start=(/it/)))
    CALL check(nf90_put_var(idfile, idvar(8), rh, start=(/it/)))
    
    DO jt=1,tcount
      CALL check(nf90_put_var(idfile, idvar(8+jt), xt(jt)%field,start=(/it/)))
    END DO
    
    CALL check(NF90_SYNC(idfile))
    IF (end_reached) CALL check(nf90_close(idfile))
    
    
  END SUBROUTINE write_output_netcdf

  !----------------------------------------------------------------------
  
  subroutine check(status)

    use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
!==============================================================================
  
END PROGRAM GMXE_BOX
