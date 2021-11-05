#include "messy_main_ppd_bi.inc"
! **********************************************************************
MODULE messy_mxl_si
! **********************************************************************
! This submodel is based on the submodel of the CHANNEL box model
! Authors: Patrick Joeckel, DLR, Oberpfaffenhofen (original code)
!          Astrid Kerkweg, UNI-MZ, Mainz (adapted to example submodel)
!    Ruud Janssen, Andrea Pozzer, MPIC, Mainz (adapted to MXL submodel)
! **********************************************************************

  ! SMCL
  USE messy_mxl
  USE messy_mxl_mem
  USE messy_main_data_bi,    ONLY: cvs, cvw,  alake                 &
                                 , slm, slf, seaice, vgrat          &
                                 , cdnl, cfncl,cfml                 &
                                 , um1, vm1, tvir, tvl, ril         & 
                                 , cdnw, cfmw, cfncw, riw, tvw      &
                                 , cdni, cfmi,cfnci, rii, tvi       &
                                 , tslm1, rh_2m                     &
                                 , rco_leaf, fws                    &
                                 , u10, v10, az0                    & 
                                 , press_3d, pressi_3d, geopoti_3d, geopot_3d &
                                 , tm1, qm1                         &
                                 , tsw, rho_air_dry_3d              &
                                 , wind10_2d, srfl, prc, prl, wsoil &
                                 , tsoil, glac, wsmx                &
                                 , zust_2d, tm1_3d, qm1_3d          &
                                 , tte_3d, rhum_3d, xlm1_3d, xlte_3d &
                                 , xim1_3d, xite_3d, albedo, aclc, qte_3d &
                                 , xim1, xlm1                       &
                                 , o3h, v3h, pressh
  USE messy_main_grid_def_bi,ONLY: grmass, grmassdry, grvol         & 
                                 , philat_2d, philon_2d, gboxarea_2d&
                                 , deltaZ, altitude_gnd, altitudei_gnd
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,      ONLY: PTR_2D_ARRAY 
  USE messy_main_constants_mem, ONLY: g, rho_air

  IMPLICIT NONE

  INTEGER                             :: idx, jt
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_entr => NULL()  ! tracer entrainment fluxes
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_emis => NULL()  ! tracer emission fluxes

  TYPE T_INIT
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_name     = '' ! name of tracer
     REAL(DP)                     :: BL               ! initial mixing ratio in boundary layer (mol mol-1)
     REAL(DP)                     :: FT               ! initial mixing ratio in free troposphere (mol mol-1)
  END TYPE T_INIT
  
  TYPE T_EMIS
     CHARACTER(LEN=STRLEN_MEDIUM) :: tr_name     = '' ! name of tracer
     CHARACTER(LEN=STRLEN_MEDIUM) :: f_emisfunc  = '' ! name of emission function
     REAL(DP)                     :: maxflux          ! maximum emission flux
  END TYPE T_EMIS

  INTEGER, PUBLIC, PARAMETER :: NMAXTRAC = 200        !total number of tracers
  TYPE(T_INIT), DIMENSION(NMAXTRAC), SAVE :: INIT_TRAC 
  TYPE(T_EMIS), DIMENSION(NMAXTRAC), SAVE :: EMIS
  INTEGER, PUBLIC, SAVE      :: ntracer = 0           !number of tracers  
  INTEGER, PUBLIC, SAVE      :: nemis = 0             !number of simple emis tracers

  ! OZONE profile related
  INTEGER :: GP_3D_LEV_O3
  INTEGER :: DIMID_LEV_O3
  INTEGER, PARAMETER :: O3_climlev = 19 ! number of levels in our O3 climatology above/over the MXL levels 

  PUBLIC :: mxl_initialize
  PUBLIC :: mxl_init_memory
  PUBLIC :: mxl_init_tracer
  PUBLIC :: mxl_global_start
  PUBLIC :: mxl_free_memory

CONTAINS
  ! --------------------------------------------------------------------
  SUBROUTINE mxl_initialize
      
    ! BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tools,         ONLY: find_next_free_unit
    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mxl_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    CALL start_message_bi(modstr,'INITIALIZATION', substr)

    ! INITIALIZE CTRL-NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mxl_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF
    CALL p_bcast(lat, p_io)
    CALL p_bcast(lon, p_io)
    CALL p_bcast(l_verbose, p_io) 
    CALL p_bcast(l_chem_ft, p_io) 
    CALL p_bcast(l_ustconst, p_io) 

   ! INITIALIZE IC_MXL-NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mxl_read_nml_ic_mxl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    ENDIF

    CALL p_bcast(hbl_ic,  p_io)
    CALL p_bcast(psurf,  p_io)
    CALL p_bcast(wthetasmax,  p_io)
    CALL p_bcast(wqsmax,  p_io)
    CALL p_bcast(beta,  p_io)
    CALL p_bcast(omega,  p_io)
    CALL p_bcast(thetam_ic,  p_io)
    CALL p_bcast(dtheta_ic,  p_io)
    CALL p_bcast(gammatheta,  p_io)
    CALL p_bcast(l_gamma,  p_io)
    CALL p_bcast(hcrit,  p_io)
    CALL p_bcast(gammath2,  p_io)
    CALL p_bcast(advtheta,  p_io)
    CALL p_bcast(qm_ic,  p_io)
    CALL p_bcast(dq_ic,  p_io)
    CALL p_bcast(gammaq,  p_io)
    CALL p_bcast(advq,  p_io)
    CALL p_bcast(um_ic, p_io)
    CALL p_bcast(ug, p_io)
    CALL p_bcast(gammau, p_io)
    CALL p_bcast(vm_ic, p_io)
    CALL p_bcast(vg, p_io) 
    CALL p_bcast(gammav, p_io)
    CALL p_bcast(uws_ic, p_io)
    CALL p_bcast(vws_ic, p_io)

    ! add fluxfunctions 

    CALL p_bcast(l_surfacelayer,  p_io)
    CALL p_bcast(z0m,  p_io)
    CALL p_bcast(z0h,  p_io)

    CALL p_bcast(l_radiation,  p_io)
    CALL p_bcast(Cc         ,  p_io)
    CALL p_bcast(salbedo    ,  p_io)

    CALL p_bcast(l_landsurface,  p_io)
    CALL p_bcast(Tsurf    ,  p_io)
    CALL p_bcast(wwilt    ,  p_io)    
    CALL p_bcast(w2       ,  p_io)    
    CALL p_bcast(w1       ,  p_io)    
    CALL p_bcast(wfc      ,  p_io)    
    CALL p_bcast(wsat     ,  p_io)    
    CALL p_bcast(CLa      ,  p_io)    
    CALL p_bcast(CLb      ,  p_io)    
    CALL p_bcast(CLc      ,  p_io)    
    CALL p_bcast(C1sat    ,  p_io)    
    CALL p_bcast(C2ref    ,  p_io)    
    CALL p_bcast(gD       ,  p_io)    
    CALL p_bcast(rsmin    ,  p_io)    
    CALL p_bcast(rssoilmin,  p_io)    
    CALL p_bcast(LAI      ,  p_io)    
    CALL p_bcast(cveg     ,  p_io)    
    CALL p_bcast(Tsoil1   ,  p_io)    
    CALL p_bcast(Tsoil2   ,  p_io)    
    CALL p_bcast(Wl       ,  p_io)    
    CALL p_bcast(Lambda   ,  p_io)    
    CALL p_bcast(CGsat    ,  p_io)

    CALL p_bcast(hc       ,  p_io)
    CALL p_bcast(drag     ,  p_io)
    CALL p_bcast(soilph   ,  p_io)

    call p_bcast(laip     ,  p_io)
    call p_bcast(btr_frac ,  p_io)
    call p_bcast(ntr_frac ,  p_io)
    call p_bcast(shr_frac ,  p_io)
    call p_bcast(hrb_frac ,  p_io)

    call p_bcast(CH4_conc    , p_io)
    call p_bcast(NOemisclass1, p_io)
    call p_bcast(NOemisclass2, p_io)
    call p_bcast(emis_NO_cult, p_io)
    call p_bcast(emis_NO_fert, p_io)
    call p_bcast(OA_bg       , p_io)

    ! INITIALIZE CHEMISTRY
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mxl_read_nml_init_chem(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       ntracer = 0
        DO jt=1, NMAXTRAC
           IF (TRIM(INIT_TRAC(jt)%tr_name) == '') CYCLE
           ntracer = ntracer + 1
           INIT_TRAC(ntracer)%tr_name = TRIM(INIT_TRAC(jt)%tr_name)
           INIT_TRAC(ntracer)%BL      = INIT_TRAC(jt)%BL
           INIT_TRAC(ntracer)%FT      = INIT_TRAC(jt)%FT
        END DO
    END IF

    ! INIT_CHEM namelist
    CALL p_bcast(ntracer, p_io)
    DO jt=1, ntracer
       CALL p_bcast(INIT_TRAC(jt)%tr_name, p_io)
       CALL p_bcast(INIT_TRAC(jt)%BL,     p_io)
       CALL p_bcast(INIT_TRAC(jt)%FT,     p_io)
    END DO

    ! INITIALIZE SIMPLE EMISSIONS
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mxl_read_nml_emis_simple(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       nemis = 0
        DO jt=1, NMAXTRAC
           IF (TRIM(EMIS(jt)%tr_name) == '') CYCLE
           nemis = nemis + 1
           EMIS(nemis)%tr_name    = TRIM(EMIS(jt)%tr_name)
           EMIS(nemis)%f_emisfunc = TRIM(EMIS(jt)%f_emisfunc)
           EMIS(nemis)%maxflux    = EMIS(jt)%maxflux * 1e-9_dp
        END DO
    END IF

    ! EMIS namelist
    CALL p_bcast(nemis, p_io)
    DO jt=1, nemis
       CALL p_bcast(EMIS(jt)%tr_name   , p_io)
       CALL p_bcast(EMIS(jt)%f_emisfunc, p_io)
       CALL p_bcast(EMIS(jt)%maxflux   , p_io)
    END DO

    CALL end_message_bi(modstr,'INITIALIZATION', substr)

  END SUBROUTINE mxl_initialize
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE mxl_init_memory


    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_2D_HORIZONTAL,  &
                                            GP_3D_INT, GP_3D_MID, SCALAR,&
                                            DIMID_LON, DIMID_LAT
    USE messy_main_channel,           ONLY: new_channel, new_channel_object &
                                          , new_attribute, get_channel_object
    USE messy_main_channel_repr,      ONLY: new_representation, AUTO
    USE messy_main_channel_dimensions,   ONLY: new_dimension            &
                                             , write_dimension          &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att
    USE messy_main_tracer_mem_bi,     ONLY: GPTRSTR, ntrac_gp, ti_gp

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mxl_init_memory'
    REAL(DP), ALLOCATABLE, DIMENSION(:) :: array
    INTEGER :: status
    INTEGER :: i,jt,jc

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

       ALLOCATE(flux_entr(ntrac_gp))
    
       CALL new_channel(status,modstr//'_entr', reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)  
   
        DO jt=1,ntrac_gp 
          CALL new_channel_object(status, modstr//'_entr'     &
               , TRIM(ti_gp(jt)%tp%ident%fullname)//'_entrflux'                &
               , p2=flux_entr(jt)%ptr)                            
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_entr'                &
               , TRIM(ti_gp(jt)%tp%ident%fullname)//'_entrflux'                      &
               , 'longname', c='entrainment flux of '//ti_gp(jt)%tp%ident%fullname )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_entr'                &
               , TRIM(ti_gp(jt)%tp%ident%fullname)//'_entrflux'                      &
               , 'units', c='mol mol-1 m s-1')
          CALL channel_halt(substr, status)
        END DO

        IF (l_emis_simple) then
          ALLOCATE(flux_emis(nemis))
        
          CALL new_channel(status,modstr//'_emis', reprid=GP_2D_HORIZONTAL)
          CALL channel_halt(substr, status)  

          DO jt=1,ntrac_gp
            DO jc=1,nemis
              IF (TRIM(ti_gp(jt)%tp%ident%fullname) == TRIM(EMIS(jc)%tr_name) &
                        .and. TRIM(EMIS(jc)%f_emisfunc) /= 'IMPORT') THEN
                CALL new_channel_object(status, modstr//'_emis'     &
                     , TRIM(EMIS(jc)%tr_name)//'_emisflux'                &
                     , p2=flux_emis(jc)%ptr)                            
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_emis'                &
                     , TRIM(EMIS(jc)%tr_name)//'_emisflux'                      &
                     , 'longname', c='emission flux for '//EMIS(jc)%tr_name )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//'_emis'                &
                     , TRIM(EMIS(jc)%tr_name)//'_emisflux'                      &
                     , 'units', c='mol mol-1 m s-1')
                CALL channel_halt(substr, status)
              ENDIF
            END DO
          END DO
          DO jc= 1,nemis
            IF (TRIM(EMIS(jc)%f_emisfunc) == 'IMPORT') then
              CALL new_channel_object(status, modstr//'_emis'     &
                   , TRIM(EMIS(jc)%tr_name)//'_emisflux'                &
                   , p2=flux_emis(jc)%ptr)                            
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//'_emis'                &
                   , TRIM(EMIS(jc)%tr_name)//'_emisflux'                      &
                   , 'longname', c='emission flux for '//EMIS(jc)%tr_name )
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//'_emis'                &
                   , TRIM(EMIS(jc)%tr_name)//'_emisflux'                      &
                   , 'units', c='mol mol-1 m s-1')
              CALL channel_halt(substr, status)
            END IF
          END DO          

        ENDIF 

    ! create new representation
    CALL new_dimension(status, DIMID_LEV_O3, 'O3_climlev', O3_climlev)
    CALL channel_halt(substr, status)
    !
    ! ... dimension variable ...
    ALLOCATE(array(O3_climlev))
    DO i=1, O3_climlev
       array(i) = REAL(i, DP)
    END DO
    CALL add_dimension_variable(status, 'O3_climlev', 'O3_climlev', array)
    CALL channel_halt(substr, status)
    DEALLOCATE(array)

    ! ... with attribute 'long_name' ...
    CALL add_dimension_variable_att(status, 'O3_climlev', 'O3_climlev', &
         'long_name', c='level index')
    CALL channel_halt(substr, status)

    ! ... with attribute 'units' ...
    CALL add_dimension_variable_att(status, 'O3_climlev', 'O3_climlev', &
         'units', c='level')
    CALL channel_halt(substr, status)

    CALL new_representation(status, GP_3D_LEV_O3, 'GP_3D_LEV_O3' &
         , rank = 3, link = 'xxx-', dctype = 0                &
         , dimension_ids = (/ DIMID_LON, DIMID_LAT, DIMID_LEV_O3 /) &
         , ldimlen       = (/ AUTO,  AUTO, AUTO   /) &
         , output_order  = (/ 1,2,3 /)                            &
         , axis = 'XYZ-'                                          &
         )
    CALL channel_halt(substr, status)


    ! create new channel
    CALL new_channel (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!      new channel objects MXL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,'hbl' , &
         p2=hbl_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hbl', &
         'long_name', c='boundary layer height')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hbl', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'hsl' , &
         p2=hsl_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hsl', &
         'long_name', c='surface layer height')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hsl', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wthetasmax' , &
         p2=wthetasmax_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetasmax', &
         'long_name', c='daily maximum surface heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetasmax', 'units', c='K m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wqsmax' , &
         p2=wqsmax_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqsmax', &
         'long_name', c='daily maximum surface moisture flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqsmax', 'units', c='g kg-1 m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'beta' , &
         p0=beta_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'beta', &
         'long_name', c='ratio of entrainment/surface heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'beta', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'omega' , &
         p0=omega_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'omega', &
         'long_name', c='subsidence rate')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'omega', 'units', c='s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'thetam' , &
         p3=thetam_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetam', &
         'long_name', c='potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetam', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'thetate' , &
         p3=thetate_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetate', &
         'long_name', c='potential temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetate', 'units', c='K s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'dtheta' , &
         p2=dtheta_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dtheta', &
         'long_name', c='potential temparature jump')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dtheta', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'gammatheta' , &
         p0=gammatheta_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammatheta', &
         'long_name', c='free troposphere temperature lapse rate')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammatheta', 'units', c='K m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'hcrit' , &
         p2=hcrit_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hcrit', &
         'long_name', c='critical BL height')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hcrit', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'advtheta' , &
         p2=advtheta_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'advtheta', &
         'long_name', c='advection of heat')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'advtheta', 'units', c='K s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'qm' , &
         p3=qm_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm', &
         'long_name', c='specific humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm', 'units', c='g kg-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'qte_mem' , &
         p3=qte_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_mem', &
         'long_name', c='specific humidity tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_mem', 'units', c='g kg-1 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'dq' , &
         p2=dq_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dq', &
         'long_name', c='specific humidity jump')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dq', 'units', c='g kg-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'gammaq' , &
         p0=gammaq_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammaq', &
         'long_name', c='free troposphere specific moisture lapse rate')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammaq', 'units', c='g kg-1 m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'advq_mem' , &
         p2=advq_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'advq_mem', &
         'long_name', c='advection of moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'advq_mem', 'units', c='g kg-1 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'rh_mem' , &
         p3=rh_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_mem', &
         'long_name', c='relative humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_mem', 'units', c='%')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'we' , &
         p2=we_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'we', &
         'long_name', c='entrainment velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'we', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ws' , &
         p2=ws_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws', &
         'long_name', c='subsidence velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ws', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'thetavm' , &
         p3=thetavm_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetavm', &
         'long_name', c='mixed layer virtual potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetavm', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'thetavsurf' , &
         p2=thetavsurf_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetavsurf', &
         'long_name', c='surface virtual potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetavsurf', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'dthetav' , &
         p2=dthetav_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dthetav', &
         'long_name', c='virtual potential temperature jump')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dthetav', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wthetas' , &
         p2=wthetas_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetas', &
         'long_name', c='surface heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetas', 'units', c='K m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wqs' , &
         p2=wqs_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqs', &
         'long_name', c='surface moisture flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqs', 'units', c='g kg-1 m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wthetavs' , &
         p2=wthetavs_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetavs', &
         'long_name', c='surface buoyancy flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetavs', 'units', c='K m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wqe' , &
         p2=wqe_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqe', &
         'long_name', c='entrainment moisture flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wqe', 'units', c='g kg-1 m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'wthetave' , &
         p2=wthetave_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetave', &
         'long_name', c='entrainment buoyancy flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wthetave', 'units', c='K m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'betaq' , &
         p2=betaq_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'betaq', &
         'long_name', c='ratio of entrainment/surface moisture flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'betaq', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'SH' , &
         p2=SH_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SH', &
         'long_name', c='surface sensible heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'SH', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'LE' , &
         p2=LE_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LE', &
         'long_name', c='surface latent heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LE', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'um' , &
         p3=um_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um', &
         'long_name', c='mixed-layer wind velocity in x-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'du' , &
         p2=du_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'du', &
         'long_name', c='wind velocity in x-direction jump')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'du', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'gammau' , &
         p0=gammau_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammau', &
         'long_name', c='free troposphere wind velocity in x-direction lapse rate')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammau', 'units', c='m s-1 m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'vm' , &
         p3=vm_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm', &
         'long_name', c='mixed-layer wind velocity in y-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'dv' , &
         p2=dv_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dv', &
         'long_name', c='wind velocity in y-direction jump')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dv', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'gammav' , &
         p0=gammav_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammav', &
         'long_name', c='free troposphere wind velocity in y-direction lapse rate')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gammav', 'units', c='m s-1 m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'uws' , &
         p2=uws_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uws', &
         'long_name', c='surface momentum in x-direction flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'uws', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'vws' , &
         p2=vws_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vws', &
         'long_name', c='surface momentum in y-direction flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vws', 'units', c='m2 s-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ueff' , &
         p2=ueff_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ueff', &
         'long_name', c='effective wind speed')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ueff', 'units', c='m s-1')
    CALL channel_halt(substr, status)

!-------------------------------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!      new channel objects MXL_radiation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,'Swin' , &
         p2=Swin_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Swin', &
         'long_name', c='incoming short wave radiation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Swin', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Swout' , &
         p2=Swout_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Swout', &
         'long_name', c='outgoing short wave radiation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Swout', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Lwin' , &
         p2=Lwin_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lwin', &
         'long_name', c='incoming long wave radiation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lwin', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Lwout' , &
         p2=Lwout_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lwout', &
         'long_name', c='outgoing long wave radiation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lwout', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Rn' , &
         p2=Rn_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn', &
         'long_name', c='net radiation at the surface')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'cossza' , &
         p0=cossza_mem, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cossza', &
         'long_name', c='cosine of solar zenith angle')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cossza', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'salbedo' , &
         p2=albedo_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'salbedo', &
         'long_name', c='surface albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'salbedo', 'units', c='-')
    CALL channel_halt(substr, status)

!-------------------------------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!      new channel objects MXL_surfacelayer
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,'z0h' , &
         p2=z0h_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'z0h', &
         'long_name', c='roughness length for heat')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'z0h', 'units', c='m')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'z0m' , &
         p2=z0m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'z0m', &
         'long_name', c='roughness length for momentum')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'z0m', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'T2m' , &
         p2=T2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'T2m', &
         'long_name', c='2 meter temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'T2m', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'q2m' , &
         p2=q2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'q2m', &
         'long_name', c='2 meter specific humidyt')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'q2m', 'units', c='g kg-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'u2m' , &
         p2=u2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u2m', &
         'long_name', c='2 meter wind velocity in x-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u2m', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'v2m' , &
         p2=v2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v2m', &
         'long_name', c='2 meter wind velocity in y-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v2m', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'u10m' , &
         p2=u10m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10m', &
         'long_name', c='10 meter wind velocity in x-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10m', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'v10m' , &
         p2=v10m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10m', &
         'long_name', c='10 meter wind velocity in y-direction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10m', 'units', c='m s-1')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'rh2m' , &
         p2=rh2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh2m', &
         'long_name', c='2 meter relative humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh2m', 'units', c='%')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'Rib' , &
         p2=Rib_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rib', &
         'long_name', c='bulk Richardson number')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rib', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ustar' , &
         p2=ustar_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustar', &
         'long_name', c='friction velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ustar', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'thetasurf' , &
         p2=thetasurf_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetasurf', &
         'long_name', c='surface potential temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'thetasurf', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'qsurf' , &
         p2=qsurf_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qsurf', &
         'long_name', c='surface specific humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qsurf', 'units', c='K')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr,'e2m' , &
         p2=e2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'e2m', &
         'long_name', c='vapor pressure at 2m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'e2m', 'units', c='Pa')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr,'esat2m' , &
         p2=esat2m_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'esat2m', &
         'long_name', c='saturation vapor pressure at 2m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'esat2m', 'units', c='Pa')
    CALL channel_halt(substr, status)
        
!-------------------------------------------------------------------------------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!      new channel objects MXL_landsurface
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL new_channel_object(status, modstr,'ra' , &
         p2=ra_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra', &
         'long_name', c='aerodynamic resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ra', 'units', c='s m-1')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'rs' , &
         p2=rs_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rs', &
         'long_name', c='stomatal resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rs', 'units', c='s m-1')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'rssoil' , &
         p2=rssoil_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rssoil', &
         'long_name', c='soil resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rssoil', 'units', c='s m-1')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'Tsurf' , &
         p2=Tsurf_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsurf', &
         'long_name', c='surface temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsurf', 'units', c='K')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr,'wwilt' , &
         p2=wwilt_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wwilt', &
         'long_name', c='volumetric soil moisture at wilting point')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wwilt', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'w1' , &
         p2=w1_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w1', &
         'long_name', c='volumetric soil moisture layer 1')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w1', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr,'w2' , &
         p2=w2_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w2', &
         'long_name', c='volumetric soil moisture layer 2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'w2', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)
        
    CALL new_channel_object(status, modstr,'wfc' , &
         p2=wfc_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wfc', &
         'long_name', c='volumetric soil moisture at field capacity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wfc', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,'wsat' , &
         p2=wsat_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsat', &
         'long_name', c='saturated volumetric water content')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsat', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'CLa' , &
         p2=CLa_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLa', &
         'long_name', c='Clapp-Hornberger retention curve parameter')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLa', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'CLb' , &
         p2=CLb_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLb', &
         'long_name', c='Clapp-Hornberger retention curve parameter')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLb', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'CLc' , &
         p2=CLc_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLc', &
         'long_name', c='Clapp-Hornberger retention curve parameter')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CLc', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'C1sat' , &
         p2=C1sat_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'C1sat', &
         'long_name', c='coefficient force term moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'C1sat', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'C2ref' , &
         p2=C2ref_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'C2ref', &
         'long_name', c='coefficient restore term moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'C2ref', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'gD' , &
         p2=gD_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gD', &
         'long_name', c='correction factor for vapor pressure deficit')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gD', 'units', c='Pa-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'rsmin' , &
         p2=rsmin_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rsmin', &
         'long_name', c='minimum resistance transpiration')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rsmin', 'units', c='s m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'rssoilmin' , &
         p2=rssoilmin_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rssoilmin', &
         'long_name', c='minimum resistance soil evaporation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rssoilmin', 'units', c='s m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'LAI' , &
         p2=LAI_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LAI', &
         'long_name', c='leaf area index')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'LAI', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'cveg' , &
         p2=cveg_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cveg', &
         'long_name', c='vegetation fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cveg', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Tsoil1' , &
         p2=Tsoil1_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsoil1', &
         'long_name', c='soil temperature layer 1 (top)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsoil1', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Tsoil2' , &
         p2=Tsoil2_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsoil2', &
         'long_name', c='soil temperature layer 2')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Tsoil2', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Wl' , &
         p2=Wl_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Wl', &
         'long_name', c='equivalent water layer depth for wet vegetation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Wl', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Lambda' , &
         p2=Lambda_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lambda', &
         'long_name', c='thermal diffusivity of the skin layer')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Lambda', 'units', c='W m-2 K-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'CGsat' , &
         p2=CGsat_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CGsat', &
         'long_name', c='saturated soil conductivity for heat')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CGsat', 'units', c='K m-2 J-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Ch' , &
         p2=Ch_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Ch', &
         'long_name', c='drag coefficient for heat and moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Ch', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Cm' , &
         p2=Cm_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Cm', &
         'long_name', c='drag coefficient for momentum')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Cm', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'GR' , &
         p2=GR_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'GR', &
         'long_name', c='ground heat flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'GR', 'units', c='W m-2')
    CALL channel_halt(substr, status)   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects DDEP
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,  'hc', &
         p2=hc_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hc', &
         'long_name', c='canopy height')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hc', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'drag', &
         p2=drag_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'drag', &
         'long_name', c='drag')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'drag', 'units', c='?')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'soilph', &
         p3=soilpH_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'soilph', &
         'long_name', c='soil pH')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'soilph', 'units', c='-')
    CALL channel_halt(substr, status)
!-------------------------------------------------------------------

    CALL new_channel_object(status, modstr, 'tm1',  &
         p3=tm1, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1' &
         , 'long_name', c = 'temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1' &
         , 'units', c = 'K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cvs', &
         p2=cvs, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs', &
         'long_name', c='snow cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvs', 'units', c='fraction')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cvw', &
         p2=cvw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw', &
         'long_name', c='wet skin fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cvw', 'units', c='fraction')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'alake', &
         p2=alake, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alake', &
         'long_name', c='lake fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'alake', 'units', c='fraction')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'slm', &
         p2=slm, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', &
         'long_name', c='land mask')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slm', 'units', c='fraction')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'slf', &
         p2=slf, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf', &
         'long_name', c='sea-land fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'slf', 'units', c='fraction')
    CALL channel_halt(substr, status)
    
    CALL new_channel_object(status, modstr,  'seaice', &
         p2=seaice, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', &
         'long_name', c='sea ice fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'seaice', 'units', c='fraction')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'vgrat', &
         p2=vgrat, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgrat', &
         'long_name', c='vegetation fraction')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vgrat', 'units', c='fraction')
    CALL channel_halt(substr, status)
 
    CALL new_channel_object(status, modstr,  'cdnl', &
         p2=cdnl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', &
         'long_name', c='neutral drag coefficient., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnl', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfncl', &
         p2=cfncl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', &
         'long_name', c='exchange parameter., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncl', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfml', &
         p2=cfml, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', &
         'long_name', c='momentum drag coefficient., land')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfml', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'um1', &
         p3=um1, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1', &
         'long_name', c='horizontal wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'um1', 'm s-1', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'vm1', &
         p3=vm1, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1', &
         'long_name', c='horizontal wind velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vm1', 'm s-1', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tvir', &
         p2=tvir, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', &
         'long_name', c='surface virtual temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvir', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tvl', &
         p2=tvl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', &
         'long_name', c='surface virtual temperature (land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvl', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'ril', &
         p2=ril, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', &
         'long_name', c='Richardson number (land)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ril', 'units', c='1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cdnw', &
         p2=cdnw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', &
         'long_name', c='neutral drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdnw', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfmw', &
         p2=cfmw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', &
         'long_name', c='momentum drag coeff., water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmw', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfncw', &
         p2=cfncw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', &
         'long_name', c='exchange parameter, water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfncw', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'riw', &
         p2=riw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', &
         'long_name', c='Richardson number, water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'riw', 'units', c='1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tvw', &
         p2=tvw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', &
         'long_name', c='surface virtual temperature (water)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvw', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cdni', &
         p2=cdni, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', &
         'long_name', c='neutral drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdni', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfmi', &
         p2=cfmi, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', &
         'long_name', c='momentum drag coeff., ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfmi', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cfnci', &
         p2=cfnci, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', &
         'long_name', c='exchange parameter, ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cfnci', 'units', c='kg m-2 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'rii', &
         p2=rii, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', &
         'long_name', c='Richardson number, ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rii', 'units', c='1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tvi', &
         p2=tvi, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', &
         'long_name', c='surface virtual temperature (ice)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tvi', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tslm1', &
         p2=tslm1, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1', &
         'long_name', c='ground surface temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tslm1', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'rh_2m', &
         p2=rh_2m, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_2m', &
         'long_name', c='relative humidity at 2m')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rh_2m', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'rco_leaf', &
         p2=rco_leaf, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', &
         'long_name', c='leaf stomatal resistance')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rco_leaf', 'units', c='s m-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'fws', &
         p2=fws, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', &
         'long_name', c='soil moisture stress function')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fws', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'u10', &
         p2=u10, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', &
         'long_name', c='10m u-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'u10', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'v10', &
         p2=v10, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', &
         'long_name', c='10m v-velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v10', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects MEGAN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,  'laip', &
         p2=laip_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'laip', &
         'long_name', c='LAI of previous month')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'laip', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'btr_frac', &
         p2=btr_frac_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'btr_frac', &
         'long_name', c='broadleaf coverage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'btr_frac', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'ntr_frac', &
         p2=ntr_frac_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ntr_frac', &
         'long_name', c='needleleaf coverage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ntr_frac', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'shr_frac', &
         p2=shr_frac_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'shr_frac', &
         'long_name', c='shrub coverage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'shr_frac', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'hrb_frac', &
         p2=hrb_frac_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hrb_frac', &
         'long_name', c='herb/grass/crop coverage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hrb_frac', 'units', c='-')
    CALL channel_halt(substr, status)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects ONEMIS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,  'philat', &
         p2=philat_2d, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat', &
         'long_name', c='geographical latitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philat', 'units', c='degree')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tsw', &
         p2=tsw, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw', &
         'long_name', c='surface temperature of water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsw', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'rho_air_dry_3d', &
         p3=rho_air_dry_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rho_air_dry_3d', &
         'long_name', c='density of dry air')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rho_air_dry_3d', 'units', c='kg m-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'wind10_2d', &
         p2=wind10_2d, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10_2d', &
         'long_name', c='10m wind speed')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wind10_2d', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'srfl', &
         p2=srfl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', &
         'long_name', c='net surface readiative flux')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'srfl', 'units', c='W m-2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'prc', &
         p2=prc, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc', &
         'long_name', c='convective precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prc', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'prl', &
         p2=prl, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl', &
         'long_name', c='large-scale precipitation')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prl', 'units', c='m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'wsoil', &
         p2=wsoil, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsoil', &
         'long_name', c='soil moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsoil', 'units', c='m3 m-3') 
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tsoil', &
         p3=tsoil, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsoil', &
         'long_name', c='deep soil temperatures')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tsoil', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'gboxarea_2d', &
         p2=gboxarea_2d, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gboxarea_2d', &
         'long_name', c='gridbox area')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'gboxarea_2d', 'units', c='m2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'philon', &
         p2=philon_2d, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon', &
         'long_name', c='geographical longitude')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'philon', 'units', c='degree')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'glac', &
         p2=glac, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', &
         'long_name', c='fraction of land covered by glaciers')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'glac', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'wsmx', &
         p2=wsmx, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', &
         'long_name', c='soil moisture at field capacity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'wsmx', 'units', c='m3 m-3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'zust_2d', &
         p2=zust_2d, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust_2d', &
         'long_name', c='surface friction velocity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zust_2d', 'units', c='m s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tm1_3d', &
         p3=tm1_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1_3d', &
         'long_name', c='air temperature')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tm1_3d', 'units', c='K')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'qm1_3d', &
         p3=qm1_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1_3d', &
         'long_name', c='specific humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1_3d', 'units', c='g g-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'CH4_conc', &
         p2=CH4_conc_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CH4_conc', &
         'long_name', c='CH4 climatological concentration') ! check name RJ
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'CH4_conc', 'units', c='mol mol-1 ?') ! check unit RJ
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'NOemisclass1', &
         p3=NOemisclass1_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NOemisclass1', &
         'long_name', c='NO emission class 1') ! check name RJ
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NOemisclass1', 'units', c='?') ! check unit RJ
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'NOemisclass2', &
         p3=NOemisclass2_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NOemisclass2', &
         'long_name', c='NO emission class 2') ! check name RJ
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'NOemisclass2', 'units', c='?') ! check unit RJ
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'emis_NO_cult', &
         p2=emis_NO_cult_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'emis_NO_cult', &
         'long_name', c='') ! check name RJ
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'emis_NO_cult', 'units', c='?') ! check unit RJ
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'emis_NO_fert', &
         p2=emis_NO_fert_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'emis_NO_fert', &
         'long_name', c='') ! check name RJ
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'emis_NO_fert', 'units', c='?') ! check unit RJ
    CALL channel_halt(substr, status)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects JVAL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,  'tte_3d', &
         p3=tte_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte_3d', &
         'long_name', c='air temperature tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tte_3d', 'units', c='K s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'rhum_3d', &
         p3=rhum_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhum_3d', &
         'long_name', c='relative humidity')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhum_3d', 'units', c='%')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'xlm1_3d', &
         p3=xlm1_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1_3d', &
         'long_name', c='cloud water')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlm1_3d', 'units', c='g g-1')
    CALL channel_halt(substr, status)
    xlm1 => xlm1_3d     ! mz_ho_20160412

    CALL new_channel_object(status, modstr,  'xlte_3d', &
         p3=xlte_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte_3d', &
         'long_name', c='cloud water tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xlte_3d', 'units', c='g g-1 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'xim1_3d', &
         p3=xim1_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1_3d', &
         'long_name', c='cloud ice')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xim1_3d', 'units', c='g g-1')
    CALL channel_halt(substr, status)
    xim1 => xim1_3d     ! mz_ho_20160412

    CALL new_channel_object(status, modstr,  'xite_3d', &
         p3=xite_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite_3d', &
         'long_name', c='cloud ice tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xite_3d', 'units', c='g g-1 s-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'albedo', &
         p2=albedo, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', &
         'long_name', c='surface albedo')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'albedo', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'Cc' , &
         p2=Cc_mem, reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Cc', &
         'long_name', c='Cloud coverage')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Cc', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'aclc', &
         p3=aclc, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', &
         'long_name', c='cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', 'units', c='-')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'cdisse', &
         p0=cdisse, reprid=SCALAR, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse', &
         'long_name', c='distance earth-sun')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'cdisse', 'units', c='AU')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'o3h', &
         p3=o3h, reprid=GP_3D_LEV_O3, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'o3h', &
         'long_name', c='column relative ozone')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'o3h', 'units', c='mol/mol')
    CALL channel_halt(substr, status)
 
    CALL new_channel_object(status, modstr,  'v3h', &
         p3=v3h, reprid=GP_3D_LEV_O3, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v3h', &
         'long_name', c='vertical ozone column')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v3h', 'units', c='molec./cm2')
    CALL channel_halt(substr, status)
 
    CALL new_channel_object(status, modstr,  'pressh', &
         p3=pressh, reprid=GP_3D_LEV_O3, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressh', &
         'long_name', c='column pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressh', 'units', c='Pa')
    CALL channel_halt(substr, status)
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects MECCA
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL new_channel_object(status, modstr,  'qm1', &
         p3=qm1, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1', &
         'long_name', c='specific moisture')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qm1', 'units', c='g g-1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'qte_3d', &
         p3=qte_3d, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_3d', &
         'long_name', c='specific moisture tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'qte_3d', 'units', c='g g-1 s-1')
    CALL channel_halt(substr, status)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! new channel objects ORACLE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    CALL new_channel_object(status, modstr,  'OA_bg', &
         p3=OA_bg_mem, reprid=GP_3D_MID, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'OA_bg', &
         'long_name', c='background organic aerosol') 
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'OA_bg', 'units', c='ug m-3') 
    CALL channel_halt(substr, status)

  END SUBROUTINE mxl_init_memory

! --------------------------------------------------------------------

  SUBROUTINE mxl_init_tracer

    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, xtm1
    USE messy_mecca_kpp,          ONLY: NSPEC, SPC_NAMES
    USE messy_main_grid_def_mem_bi, ONLY: nlat, nlon, nlev
    USE messy_main_tracer,        ONLY: TRSET, NSETID, t_trinfo_list    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mxl_init_tracer'
    INTEGER :: jc, jt, zid, ntrac
    TYPE(t_trinfo_list), POINTER :: ti => NULL()

    CALL start_message_bi(modstr,'INITIALIZE TRACER', substr)
    
    set_loop: DO zid=1, NSETID

    ! NUMBER OF TRACERS
    ntrac = trset(zid)%ntrac

      IF (trim(trset(zid)%name)=='gp') THEN
        WRITE(*,*) 'initialize TRACER SET from namelist: gp'

        !LOOP OVER TRACERS IN SET
        ti => trset(zid)%tilist
        tracer_loop: DO 
           IF (.NOT. ASSOCIATED(ti)) EXIT
        
           jt = ti%info%ident%idx
             DO jc = 1, ntrac_gp
               IF (TRIM(ti%info%ident%fullname)==TRIM(INIT_TRAC(jc)%tr_name)) THEN
                 trset(zid)%xtm1(:,:,jt,nlev,:)     = INIT_TRAC(jc)%BL * 1.E-9_dp
                 trset(zid)%xtm1(:,:,jt,nlev-1,:)   = INIT_TRAC(jc)%FT * 1.E-9_dp 
               ENDIF
             ENDDO
           ti => ti%next
        END DO tracer_loop
      END IF
    END DO set_loop
    
    CALL end_message_bi(modstr,'INITIALIZE TRACER', substr)
        
  END SUBROUTINE mxl_init_tracer

  ! --------------------------------------------------------------------
  SUBROUTINE mxl_global_start

    USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, xt, xtte, xtm1, ti_gp
    USE messy_main_timer,         ONLY: delta_time, time_step_len, lstart
    USE messy_main_grid_def_mem_bi, ONLY: nlat, nlon, nlev
    USE messy_main_constants_mem, ONLY: STRLEN_SHORT    

    IMPLICIT NONE

    INTEGER :: jc, ji, jk, jj
    REAL(dp), DIMENSION(nlat,nlon,ntrac_gp,nlev) :: conc
    real(dp) :: hasl
    real(dp) :: htr
    real(dp) :: actflux 
    real(dp), allocatable, dimension(:) :: v3
    real(dp), allocatable, dimension(:) :: relo3
    real(dp), allocatable, dimension(:) :: jpress

    conc(:,:,:,:)        = 0_dp
 
    if (lstart) then
       xt(:,:,:,:)           = xtm1(:,:,:,:)
       ! MXL initial conditions
       hbl_mem(:,:)          = hbl_ic
       hsl_mem               = 0.1_dp * hbl_ic
       pressi_3d(:,:,nlev+1) = psurf
       wthetasmax_mem(:,:)   = wthetasmax           
       wqsmax_mem(:,:)       = wqsmax              
       beta_mem              = beta             
       omega_mem             = omega           
       thetam_mem(:,:,nlev)  = thetam_ic           
       thetam_mem(:,:,nlev-1)= thetam_ic + dtheta_ic           
       dtheta_mem(:,:)       = dtheta_ic           
       gammatheta_mem        = gammatheta
       thetavm_mem           = thetam_mem ! thetav initialized as theta
       thetasurf_mem         = thetam_mem(1:nlon,1:nlat,nlev) ! thetasurf initialized as theta
       hcrit_mem(:,:)        = hcrit           
       advtheta_mem(:,:)     = advtheta           
       qm_mem(:,:,nlev)      = qm_ic           
       qm_mem(:,:,nlev-1)    = qm_ic + dq_ic           
       dq_mem(:,:)           = dq_ic           
       gammaq_mem            = gammaq           
       advq_mem(:,:)         = advq
       um_mem(:,:,nlev)      = um_ic
       du_mem(:,:)           = ug - um_ic
       gammau_mem            = gammau
       vm_mem(:,:,nlev)      = vm_ic
       dv_mem(:,:)           = vg - vm_ic
       gammav_mem            = gammav
       uws_mem(:,:)          = uws_ic
       vws_mem(:,:)          = vws_ic   
       ueff_mem              = sqrt(um_mem(:,:,nlev)**2._dp + vm_mem(:,:,nlev)**2._dp)
       z0h_mem(:,:)          = z0h
       z0m_mem(:,:)          = z0m
       albedo_mem(:,:)       = salbedo
       Cc_mem(:,:)           = Cc
       Tsurf_mem(:,:)        = Tsurf
       wwilt_mem(:,:)        = wwilt
       w2_mem(:,:)           = w2
       w1_mem(:,:)           = w1
       wfc_mem(:,:)          = wfc
       wsat_mem(:,:)         = wsat
       CLa_mem(:,:)          = CLa
       CLb_mem(:,:)          = CLb
       CLc_mem(:,:)          = CLc
       C1sat_mem(:,:)        = C1sat
       C2ref_mem(:,:)        = C2ref
       gD_mem(:,:)           = gD
       rsmin_mem(:,:)        = rsmin
       rssoilmin_mem(:,:)    = rssoilmin
       LAI_mem(:,:)          = LAI
       cveg_mem(:,:)         = cveg
       Tsoil1_mem(:,:)       = Tsoil1
       Tsoil2_mem(:,:)       = Tsoil2
       Wl_mem(:,:)           = Wl
       Lambda_mem(:,:)       = Lambda
       CGsat_mem(:,:)        = CGsat
       Ch_mem(:,:)           = 1e-5_dp
       ! DDEP initial conditions
       hc_mem(:,:)           = hc
       drag_mem(:,:)         = drag
       soilpH_mem(:,:,:)     = soilph
       ! MEGAN initial conditions
       laip_mem(:,:)         = laip
       btr_frac_mem(:,:)     = btr_frac
       ntr_frac_mem(:,:)     = ntr_frac
       shr_frac_mem(:,:)     = shr_frac
       hrb_frac_mem(:,:)     = hrb_frac
       ! ONEMIS initial conditions
       CH4_conc_mem(:,:)     = CH4_conc
       emis_NO_cult_mem(:,:) = emis_NO_cult
       emis_NO_fert_mem(:,:) = emis_NO_fert
       OA_bg_mem(:,:,:)      = OA_bg 

     endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! radiation calculations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (l_radiation) then
    CALL mxl_radiation(hsl_mem, pressi_3d, thetam_mem, &
                       Cc_mem, albedo_mem, Tsurf_mem, &
                       Swin_mem,Swout_mem,Lwin_mem,Lwout_mem,Rn_mem,cossza_mem)  ! out
  !else: rad4all
  endif 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! surface layer calculations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (l_surfacelayer) then
    CALL mxl_surfacelayer(pressi_3d,thetavm_mem,thetam_mem,qm_mem,um_mem,vm_mem,hsl_mem, &
                              wthetas_mem,wqs_mem,ueff_mem, & 
                            rs_mem, & 
                            z0h_mem, z0m_mem, &   
                            thetasurf_mem,thetavsurf_mem,qsurf_mem,uws_mem,vws_mem,T2m_mem, & ! out
                              q2m_mem,u2m_mem,v2m_mem,u10m_mem,v10m_mem,rh2m_mem,Rib_mem,&    ! out
                              e2m_mem,esat2m_mem,ustar_mem,Cm_mem,Ch_mem)                     ! out
  do ji = 1,nlon
    do jj = 1,nlat
      if(ustar_mem(ji,jj) .le. 0) stop "u* has to be greater than 0: increase initial wind or momentum flux"
    enddo
  enddo

  else ! no surfacelayer
    CALL mxl_momentumflux(l_ustconst, um_mem, vm_mem, z0m_mem, &
                          uws_mem, vws_mem, ustar_mem)   ! out                        
  endif ! l_surfacelayer

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! land surface calculations
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (l_landsurface) then
    if ((f_wthetas /= 'INTERACT') .or. (f_wqs /= 'INTERACT')) then
      stop 'land surface activated (l_landsurface = T): set f_wthetas and f_wqs to "INTERACT"'
    endif
  do ji = 1,nlon
    do jj = 1,nlat
      if(ustar_mem(ji,jj) .le. 0) then 
        write(*,*) "u* has to be greater than 0:" 
        write(*,*) "  a) increase initial wind velocity (u,v) and set l_uconst = F, or"
        write(*,*) "  b) increase surface momentum flux (uws, vws) and set l_ustconst = T" 
        stop
      endif
    enddo
  enddo
    CALL mxl_landsurface(l_surfacelayer, l_radiation, &
                     rs_mem,ra_mem,rssoil_mem,SH_mem,LE_mem,GR_mem,Tsurf_mem, & ! out
                     pressi_3d,thetam_mem,qm_mem,ueff_mem, &
                     Swin_mem,Rn_mem, &
                     e2m_mem,esat2m_mem,T2m_mem,ustar_mem, &
                     LAI_mem,wwilt_mem,w2_mem,w1_mem,wfc_mem,wsat_mem, &
                       CLa_mem,CLb_mem,CLc_mem,C1sat_mem,C2ref_mem,gD_mem,rsmin_mem, &
                       rssoilmin_mem,Lambda_mem,cveg_mem,Wl_mem,CGsat_mem, &
                       Tsoil1_mem, Tsoil2_mem)

  endif ! l_landsurface 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! calculation mixed-layer dynamics
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL mxl_algorithm(hbl_mem, hsl_mem, press_3d, pressi_3d, wthetasmax_mem, wqsmax_mem, beta_mem, omega_mem, &
                    thetam_mem, thetate_mem, dtheta_mem, gammatheta_mem, &
                    l_gamma, gammath2, hcrit_mem, advtheta_mem, &
                    qm_mem, qte_mem,   dq_mem,     gammaq_mem,     advq_mem, rh_mem, &
                    um_mem,     du_mem,     gammau_mem, uws_mem, &
                    vm_mem,     dv_mem,     gammav_mem, vws_mem, &
                    we_mem, ws_mem, thetavm_mem, dthetav_mem, wthetas_mem, wqs_mem, wthetavs_mem, &
                    wthetave_mem, wqe_mem, betaq_mem, SH_mem, LE_mem, ueff_mem, &
                    f_wthetas, starttime_wths,stoptime_wths,&
                    f_wqs, starttime_wqs, stoptime_wqs, &
                    starttime_adv, stoptime_adv)
  DO ji = 1,nlon
    DO jj = 1,nlat 
      IF (ueff_mem(ji,jj) .lt. 1.e-2) then 
        stop 'effective wind speed below 1 cm/s: increase initial wind speed'
      ENDIF
    ENDDO
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! end of MXL physics
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! assign values to ddep variables   | description          | used in calculation of.. | get from:
   ! ------------------------------------------------------------------------------------------------------
   vgrat     = cveg_mem                ! vegetation fraction  |                          | mxl.nml
   cdnl      = 0.002                   ! neutral drag coeff ??|         [ra]             |
   cfncl     = 0.001                   ! exch param ??        |         [ra]             |
   cfml      = Cm_mem                  !                      |         [ra],            | mxl_surfacelayer
   az0       = z0m_mem                 !                      |         [ra, vdaer]      | mxl.nml
   um1       = um_mem                  !                      |         [ra, vdaer]      | mxl_algorithm 
   vm1       = vm_mem                  !                      |         [ra, vdaer]      | mxl_algorithm 
   tvir      = thetavsurf_mem          !                      |         [ra]             | mxl_surfacelayer
   tvl       = thetavm_mem(:,:,nlev)   !                      |         [ra]             | mxl_algorithm
   ril       = Rib_mem                 !                      |         [ra]             | mxl_surfacelayer
   tslm1     = Tsurf_mem               !                      |   [vdbl, vdaer]    | mxl_surfacelayer, mxl.nml
   rh_2m     = rh2m_mem                !                      |         [vdbl, vdaer]    | mxl_surfacelayer 
   rco_leaf  = rs_mem                  !                      |         [vdbl]           | mxl_landsurface
   fws       = 0                       ! soil moisture stress function??|[vdbl]          | 
   u10       = u10m_mem                !                      |         [vdaer]          | mxl_surfacelayer
   v10       = v10m_mem                !                      |         [vdaer]          | mxl_surfacelayer
   ! variables related to water and snow surfaces
   cvs       = 0._dp
   cvw       = 0._dp
   alake     = 0._dp
   slm       = 1._dp
   slf       = slm
   seaice    = 0._dp
   cdnw      = 0._dp
   cfmw      = 0._dp
   cfncw     = 0._dp
   riw       = 0._dp
   tvw       = 0._dp
   cdni      = 0._dp
   cfmi      = 0._dp
   cfnci     = 0._dp
   rii       = 0._dp
   tvi       = 0._dp

   ! calculate pressure, geopotential   
   htr                    = 10000 ! tropopause heigth, fixed at 10 km
   hasl                   = 0.    ! height above sea level (m)
   ! ---------------------------------------------------------------------------------------
   ! calculation of geopoti follows a different logic than that of other interface variables
   ! by definition: geopoti at surface == 0 (so difference from sea level not accounted for)
   ! upper boundary not defined, here ad hoc solution chosen: calculate from tropopause height
   ! ---------------------------------------------------------------------------------------
   geopoti_3d(:,:,nlev+1) = 0.  ! below surface 
   geopoti_3d(:,:,nlev  ) = 0.  ! surface level 
   geopoti_3d(:,:,nlev-1) = hbl_mem * g !(htr + hasl) * g
   geopot_3d(:,:,nlev-1)  = (htr-hbl_mem)/2._dp * g ! free troposphere          
   geopot_3d(:,:,nlev)    = hbl_mem/2._dp * g       ! boundary layer

   altitude_gnd(:,:,nlev-1) = (htr-hbl_mem)/2._dp
   altitude_gnd(:,:,nlev)   = hbl_mem/2._dp
   altitudei_gnd(:,:,nlev-1) = htr
   altitudei_gnd(:,:,nlev)   = hbl_mem
   altitudei_gnd(:,:,nlev+1) = 0._dp

   deltaZ(:,:,nlev) = hbl_mem
   deltaZ(:,:,nlev-1) = htr-hbl_mem

   pressi_3d(:,:,nlev  )  = pressi_3d(:,:,nlev+1) - rho_air * hbl_mem * g ! pressure at BL-FT interface
   pressi_3d(:,:,nlev-1)  = 10000 ! pressure at tropopause, fixed at 100 hPa
   do idx = 1,nlev
     press_3d(:,:,idx)               =  (pressi_3d(:,:,idx+1) + pressi_3d(:,:,idx)) / 2.0_DP 
     grmass(:,:,idx)                 = ((pressi_3d(:,:,idx+1) - pressi_3d(:,:,idx)) / g) 
   enddo
   grvol(:,:,nlev)        = hbl_mem * 1_dp * 1_dp
   grvol(:,:,nlev-1)      = (htr-hbl_mem)  * 1_dp * 1_dp

   ! assign values to onemis variables 
   philat_2d      = lat
   philon_2d      = lon 
   tsw            = 0.
   rho_air_dry_3d = grmass/grvol
   wind10_2d      = sqrt(u10m_mem**2 + v10m_mem**2)
   srfl           = Rn_mem
   prc            = 0.
   prl            = 0.
   wsoil          = w1_mem
   tsoil(:,:,1)   = Tsoil1_mem
   tsoil(:,:,2)   = Tsoil2_mem
   gboxarea_2d    = 1._dp
   glac           = 0.
   wsmx           = wfc_mem
   zust_2d        = ustar_mem
   tm1_3d         = thetam_mem 
   qm1_3d         = qm_mem*1e-3     ! use commented out in onemis

   ! assign values to jval variables 
   tte_3d  = thetate_mem !temperature tendency
   rhum_3d = rh_mem
   xlm1_3d = 0.
   xlte_3d = 0. 
   xim1_3d = 0.
   xite_3d = 0.
   albedo  = albedo_mem
   aclc    = Cc
   cdisse  = 1.00_dp ! distance sun-earth (AU)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DEFINE A CLIMATOLOGOGICAL O3 CONCENTRATIONS: copied from messy_jval_box.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CAABA VALUES : TO BE IMPROVED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (O3_climlev == 19) then
  allocate(v3(O3_climlev))
  allocate(relo3(O3_climlev))
  allocate(jpress(O3_climlev))
    ! ozone column (Dobson Unit, [mlc/cm^2]), 
  v3    = 1.0*(/ &
     3.366E+17, 1.437E+18, 4.085E+18, 5.428E+18, 6.157E+18, 6.583E+18, &
     6.860E+18, 7.070E+18, 7.227E+18, 7.343E+18, 7.436E+18, 7.523E+18, &
     7.605E+18, 7.678E+18, 7.740E+18, 7.788E+18, 7.822E+18, 7.844E+18, &
     7.857E+18 /)
   ! relative ozone, i.e. ozone mixing ratio [mol/mol]
  relo3 = (/ &
     7.182E-06, 8.319E-06, 4.172E-06, 2.041E-06, 9.525E-07, 4.334E-07, &
     2.571E-07, 1.514E-07, 9.760E-08, 5.775E-08, 5.064E-08, 4.394E-08, &
     3.980E-08, 3.636E-08, 3.209E-08, 2.807E-08, 2.479E-08, 2.242E-08, &
     2.105E-08 /)
   ! pressure [Pa]
  jpress = (/  &
     1000., 3000., 5040., 7339., 10248., 14053., 18935., 24966., 32107., &
     40212., 49027., 58204., 67317., 75897., 83472., 89631., 94099.,     &
     96838., 98169. /) 
elseif (O3_climlev == 90) then
  allocate(v3(O3_climlev))
  allocate(relo3(O3_climlev))
  allocate(jpress(O3_climlev))
  ! vertical ozone column [mcl/cm2]
  v3   = (/ &
     9.585E+12, 2.593E+13, 2.351E+14, 7.025E+14, 1.322E+15, 2.108E+15, &
     3.092E+15, 4.312E+15, 5.843E+15, 7.789E+15, 1.029E+16, 1.346E+16, &
     1.747E+16, 2.246E+16, 2.861E+16, 3.616E+16, 4.540E+16, 5.665E+16, &
     7.030E+16, 8.678E+16, 1.066E+17, 1.302E+17, 1.581E+17, 1.907E+17, &
     2.288E+17, 2.728E+17, 3.233E+17, 3.810E+17, 4.463E+17, 5.199E+17, &
     6.024E+17, 6.941E+17, 7.955E+17, 9.070E+17, 1.029E+18, 1.163E+18, &
     1.308E+18, 1.466E+18, 1.637E+18, 1.822E+18, 2.020E+18, 2.229E+18, &
     2.448E+18, 2.674E+18, 2.906E+18, 3.143E+18, 3.385E+18, 3.628E+18, &
     3.874E+18, 4.121E+18, 4.366E+18, 4.609E+18, 4.848E+18, 5.081E+18, &
     5.306E+18, 5.524E+18, 5.733E+18, 5.930E+18, 6.112E+18, 6.276E+18, &
     6.422E+18, 6.551E+18, 6.664E+18, 6.765E+18, 6.856E+18, 6.940E+18, &
     7.018E+18, 7.092E+18, 7.161E+18, 7.225E+18, 7.284E+18, 7.334E+18, &
     7.376E+18, 7.412E+18, 7.443E+18, 7.472E+18, 7.501E+18, 7.530E+18, &
     7.560E+18, 7.592E+18, 7.627E+18, 7.664E+18, 7.702E+18, 7.742E+18, &
     7.782E+18, 7.823E+18, 7.862E+18, 7.894E+18, 7.916E+18, 7.929E+18, &
     7.935E+18 /)
  ! relative ozone, i.e. ozone mixing ratio [mol/mol]
  relo3 = (/ &
     4.015E-07, 7.810E-07, 9.022E-07, 9.528E-07, 1.010E-06, 1.083E-06, &
     1.174E-06, 1.290E-06, 1.433E-06, 1.604E-06, 1.804E-06, 2.031E-06, &
     2.281E-06, 2.548E-06, 2.840E-06, 3.152E-06, 3.485E-06, 3.832E-06, &
     4.203E-06, 4.586E-06, 4.970E-06, 5.367E-06, 5.762E-06, 6.135E-06, &
     6.483E-06, 6.810E-06, 7.102E-06, 7.363E-06, 7.588E-06, 7.785E-06, &
     7.948E-06, 8.090E-06, 8.201E-06, 8.280E-06, 8.331E-06, 8.354E-06, &
     8.360E-06, 8.366E-06, 8.316E-06, 8.163E-06, 7.909E-06, 7.582E-06, &
     7.204E-06, 6.810E-06, 6.407E-06, 5.993E-06, 5.585E-06, 5.185E-06, &
     4.781E-06, 4.375E-06, 3.960E-06, 3.548E-06, 3.156E-06, 2.794E-06, &
     2.454E-06, 2.132E-06, 1.816E-06, 1.501E-06, 1.225E-06, 9.909E-07, &
     7.922E-07, 6.395E-07, 5.220E-07, 4.307E-07, 3.607E-07, 3.065E-07, &
     2.615E-07, 2.219E-07, 1.852E-07, 1.493E-07, 1.148E-07, 8.642E-08, &
     6.697E-08, 5.531E-08, 4.807E-08, 4.365E-08, 4.096E-08, 3.948E-08, &
     3.833E-08, 3.733E-08, 3.574E-08, 3.348E-08, 3.141E-08, 2.937E-08, &
     2.732E-08, 2.549E-08, 2.413E-08, 2.293E-08, 2.188E-08, 2.083E-08, &
     2.083E-08 /)
  ! pressure [Pa]
  jpress = (/ &
         0.99,      3.18,      5.81,      8.96,     12.74,     17.17, &
        22.27,     28.14,     34.88,     42.64,     51.43,     61.28, &
        72.21,     84.23,     97.45,    111.99,    127.99,    145.59, &
       164.95,    186.24,    209.56,    234.97,    262.66,    292.85, &
       325.76,    361.63,    400.73,    443.34,    489.79,    540.42, &
       595.44,    655.07,    719.68,    789.69,    865.56,    947.77, &
      1036.86,   1133.41,   1238.02,   1351.39,   1474.24,   1607.37, &
      1751.63,   1907.95,   2077.31,   2260.76,   2459.51,   2674.83, &
      2908.17,   3161.09,   3436.18,   3736.45,   4063.97,   4421.67, &
      4813.36,   5242.81,   5713.73,   6230.96,   6799.36,   7422.58, &
      8104.57,   8852.52,   9675.99,  10584.35,  11587.46,  12697.23, &
     13926.26,  15286.00,  16787.98,  18443.60,  20268.79,  22280.20, &
     24495.96,  26934.59,  29615.81,  32566.47,  35815.20,  39391.63, &
     43328.65,  47656.74,  52414.16,  57639.89,  63371.77,  69660.11, &
     76471.58,  83472.54,  89631.42,  94099.31,  96838.38,  98169.53 /)
else
  stop 'O3_climlev should be either 19 or 90' 
endif      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO ji =1, nlon
       DO jj=1, nlat
        DO jk = 1, O3_climlev
            o3h(ji,jj,jk)      =  relo3(jk)
            v3h(ji,jj,jk)      =  v3(jk)
            pressh(ji,jj,jk)   =  jpress(jk)
          ENDDO
        ENDDO
      ENDDO

   ! assign values to mecca variables
   qte_3d  = qte_mem*1e-3 ! spec. hum. tendency, used in messy_mecca_si
   qm1     = qm_mem*1e-3  ! spec. hum.

   ! assign values to oracle variables 
   tm1     = tm1_3d
 
!-------------------------------------------------------------------------------------------------------

   tracer_loop: DO jc=1,ntrac_gp
     IF (ti_gp(jc)%tp%ident%fullname(1:2) == 'PT')  CYCLE
     !----------------------------------------------------
     !  SIMPLE EMISSION CALCULATION (mol mol-1 m s-1)
     !----------------------------------------------------
     IF (l_emis_simple) THEN
       DO jt=1,nemis
         IF (TRIM(ti_gp(jc)%tp%ident%fullname) == TRIM(EMIS(jt)%tr_name)) THEN
           CALL emis_simple(EMIS(jt)%f_emisfunc,EMIS(jt)%maxflux, &
                         starttime_emis, stoptime_emis, actflux)
            flux_emis(jt)%ptr(:,:) = actflux
            xtte(:,:,jc,nlev) = xtte(:,:,jc,nlev) + flux_emis(jt)%ptr(:,:) / hbl_mem(:,:)
         ENDIF
       ENDDO  
     ENDIF

     ! tracer concentration before entrainment
     conc(:,:,jc,:) = xtm1(:,:,jc,:) + xtte(:,:,jc,:)*delta_time
     !----------------------------------------------------
     ! ENTRAINMENT FLUX CALCULATION (mol mol-1 m s-1)
     !----------------------------------------------------
     flux_entr(jc)%ptr(:,:) = we_mem(:,:) * (conc(:,:,jc,nlev-1) - conc(:,:,jc,nlev))
     !----------------------------------------------------
     ! UPDATE CONCENTRATION tendency
     !----------------------------------------------------
     xtte(:,:,jc,nlev) = xtte(:,:,jc,nlev) + flux_entr(jc)%ptr(:,:) / hbl_mem(:,:)
   END DO tracer_loop

  END SUBROUTINE mxl_global_start

! ===========================================================================
  SUBROUTINE mxl_read_nml_ic_mxl(status, iou)
   
    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /IC_MXL/ hbl_ic, psurf, wthetasmax, wqsmax, beta, omega, &
                      thetam_ic, dtheta_ic, gammatheta, l_gamma, hcrit, gammath2, advtheta, & 
                      qm_ic,     dq_ic,     gammaq,     advq, &
                      f_wthetas, starttime_wths, stoptime_wths, &
                      f_wqs, starttime_wqs, stoptime_wqs, &
                      starttime_adv, stoptime_adv, &
                      um_ic, vm_ic, ug, vg, uws_ic, vws_ic, gammau, gammav, l_ustconst, &  
                      l_surfacelayer, z0m, z0h, &
                      l_radiation, Cc, salbedo, &
                      l_landsurface, Tsurf, wwilt, w2, w1, wfc, wsat, CLa, CLb, CLc, C1sat, C2ref, &
                      gD, rsmin, rssoilmin, LAI, cveg, Tsoil1, Tsoil2, Wl, Lambda, CGsat, &
                      hc, drag, soilph, laip, btr_frac, ntr_frac, shr_frac, hrb_frac, &
                      CH4_conc, NOemisclass1, NOemisclass2, emis_NO_cult, emis_NO_fert, &
                      OA_bg

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mxl_read_nml_ic_mxl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'IC_MXL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=IC_MXL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'IC_MXL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   IF (l_verbose) THEN
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '             INITIAL CONDITIONS MXL                      '
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '  hbl_ic       = ', hbl_ic
     WRITE(*,*) '  psurf        = ', psurf
     WRITE(*,*) '  wthetasmax   = ', wthetasmax
     WRITE(*,*) '  f_wthetas    = ', f_wthetas
     WRITE(*,*) '  wqsmax       = ', wqsmax
     WRITE(*,*) '  f_wqs        = ', f_wqs
     WRITE(*,*) '  beta         = ', beta
     WRITE(*,*) '  omega        = ', omega
     WRITE(*,*) '  thetam_ic    = ', thetam_ic
     WRITE(*,*) '  dtheta_ic    = ', dtheta_ic
     WRITE(*,*) '  gammatheta   = ', gammatheta
     WRITE(*,*) '  l_gamma       = ', l_gamma
     WRITE(*,*) '  gammath2     = ', gammath2
     WRITE(*,*) '  hcrit        = ', hcrit
     WRITE(*,*) '  advtheta     = ', advtheta
     WRITE(*,*) '  qm_ic        = ', qm_ic
     WRITE(*,*) '  dq_ic        = ', dq_ic
     WRITE(*,*) '  gammaq       = ', gammaq
     WRITE(*,*) '  advq         = ', advq
     WRITE(*,*) '  um_ic        = ', um_ic
     WRITE(*,*) '  ug           = ', ug
     WRITE(*,*) '  gammau       = ', gammau
     WRITE(*,*) '  vm_ic        = ', vm_ic
     WRITE(*,*) '  vg           = ', vg
     WRITE(*,*) '  gammav       = ', gammav
     WRITE(*,*) '  uws_ic       = ', uws_ic
     WRITE(*,*) '  vws_ic       = ', vws_ic
     WRITE(*,*) '  l_ustcont    = ', l_ustconst
     WRITE(*,*) '  l_surfacelayer = ', l_surfacelayer
     WRITE(*,*) '  z0h          = ', z0h
     WRITE(*,*) '  z0m          = ', z0m
     WRITE(*,*) '  l_radiation  = ', l_radiation
     WRITE(*,*) '  Cc           = ', Cc
     WRITE(*,*) '  salbedo       = ', salbedo
     WRITE(*,*) '  l_landsurface = ', l_landsurface
     WRITE(*,*) '  Tsurf        = ', Tsurf 
     WRITE(*,*) '  wwilt        = ', wwilt 
     WRITE(*,*) '  w2           = ', w2 
     WRITE(*,*) '  w1           = ', w1 
     WRITE(*,*) '  wfc          = ', wfc 
     WRITE(*,*) '  wsat         = ', wsat 
     WRITE(*,*) '  CLa          = ', CLa 
     WRITE(*,*) '  CLb          = ', CLb 
     WRITE(*,*) '  CLc          = ', CLc 
     WRITE(*,*) '  C1sat        = ', C1sat 
     WRITE(*,*) '  C2ref        = ', C2ref 
     WRITE(*,*) '  gD           = ', gD 
     WRITE(*,*) '  rsmin        = ', rsmin 
     WRITE(*,*) '  rssoilmin    = ', rssoilmin 
     WRITE(*,*) '  LAI          = ', LAI 
     WRITE(*,*) '  cveg         = ', cveg 
     WRITE(*,*) '  Tsoil1       = ', Tsoil1 
     WRITE(*,*) '  Tsoil2       = ', Tsoil2 
     WRITE(*,*) '  Wl           = ', Wl 
     WRITE(*,*) '  Lambda       = ', Lambda 
     WRITE(*,*) '  CGsat        = ', CGsat  
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '             INITIAL CONDITIONS DDEP                     '
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '  hc           = ', hc
     WRITE(*,*) '  drag         = ', drag
     WRITE(*,*) '  soilpH       = ', soilph
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '             INITIAL CONDITIONS MEGAN                    '
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '  laip         = ', laip
     WRITE(*,*) '  btr_frac     = ', btr_frac
     WRITE(*,*) '  ntr_frac     = ', ntr_frac
     WRITE(*,*) '  shr_frac     = ', shr_frac
     WRITE(*,*) '  hrb_frac     = ', hrb_frac
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '             INITIAL CONDITIONS ONEMIS                   '
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '  CH4_conc     = ', CH4_conc
     WRITE(*,*) '  NOemisclass1 = ', NOemisclass1
     WRITE(*,*) '  NOemisclass2 = ', NOemisclass2
     WRITE(*,*) '  emis_NO_cult = ', emis_NO_cult
     WRITE(*,*) '  emis_NO_fert = ', emis_NO_fert
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '             INITIAL CONDITIONS ORACLE                   '
     WRITE(*,*) '.........................................................'
     WRITE(*,*) '  OA_bg        = ', OA_bg
     
   ENDIF

  END SUBROUTINE mxl_read_nml_ic_mxl
! ===========================================================================

  SUBROUTINE mxl_read_nml_init_chem(status, iou)
   
    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
!
    NAMELIST /INIT_CHEM/  INIT_TRAC
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mxl_read_nml_init_chem'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'INIT_CHEM', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=INIT_CHEM, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'INIT_CHEM', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '              TRACER MIXING RATIO INITIALIZATION         '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXTRAC
       IF (TRIM(INIT_TRAC(jt)%tr_name) == '') CYCLE

       WRITE(*,*) '  TRACER NO.               ',jt
       WRITE(*,*) '  TRACER NAME            = ', TRIM(INIT_TRAC(jt)%tr_name)
       WRITE(*,*) '  BL mixing ratio (ppb)  = ', INIT_TRAC(jt)%BL
       WRITE(*,*) '  FT mixing ratio (ppb)  = ', INIT_TRAC(jt)%FT
       WRITE(*,*) '.........................................................'
   END DO
 END SUBROUTINE mxl_read_nml_init_chem 

!! ===========================================================================

  SUBROUTINE mxl_read_nml_emis_simple(status, iou)
   
    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
!
    NAMELIST /EMIS_SIMPLE/ l_emis_simple, starttime_emis, stoptime_emis, EMIS
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mxl_read_nml_init_chem'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'EMIS_SIMPLE', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=EMIS_SIMPLE, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'EMIS_SIMPLE', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '              SIMPLE EMISSIONS                           '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXTRAC
       IF (TRIM(EMIS(jt)%tr_name) == '') CYCLE
       WRITE(*,*) '  TRACER NO.                    = ', jt
       WRITE(*,*) '  TRACER NAME                   = ', TRIM(EMIS(jt)%tr_name)
       WRITE(*,*) '  FLUX FUNCTION                 = ', TRIM(EMIS(jt)%f_emisfunc)
       WRITE(*,*) '  maximum emission (ppb m s-1)  = ', EMIS(jt)%maxflux
       WRITE(*,*) '.........................................................'
   END DO
 END SUBROUTINE mxl_read_nml_emis_simple 

!! ===========================================================================

  SUBROUTINE mxl_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(flux_entr))   DEALLOCATE(flux_entr)
    IF (ASSOCIATED(flux_emis))   DEALLOCATE(flux_emis)

  END SUBROUTINE mxl_free_memory

! **********************************************************************
END MODULE messy_mxl_si
! **********************************************************************
