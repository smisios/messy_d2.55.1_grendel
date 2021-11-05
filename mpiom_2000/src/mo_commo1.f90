!>
!! @ingroup common
!!
      MODULE MO_COMMO1

      USE mo_kind, ONLY: dp, wp
      USE mo_param1, ONLY: ie, je, ke, kep, ie_g, je_g, ill, ilt, ito, jto

      IMPLICIT NONE

      REAL(dp), PARAMETER :: zero   = 0.0_dp
      REAL(dp), PARAMETER :: one    = 1.0_dp
      REAL(dp), PARAMETER :: two    = 2.0_dp
      REAL(dp), PARAMETER :: four   = 4.0_dp

      REAL(dp), PARAMETER :: half   = 0.5_dp
      REAL(dp), PARAMETER :: fourth = 0.25_dp
      REAL(dp), PARAMETER :: eighth = 0.125_dp
      REAL(dp), PARAMETER :: tenm4  = 1.0e-4_dp

      INTEGER IAUFR,IAUFW
      INTEGER I3DREST

      ! set of istart options:
      !> COMPLETELY NEW SETUP,
      !! topography read from anta and written to topo (!!! WARNING !!!)
      !! start from climatology
      INTEGER, PARAMETER :: istart_new_topo_update=0
      !> new run, topography read from topo
      !! start from horizontally uniform ts-profile (taf,saf)
      INTEGER, PARAMETER :: istart_new_run=1
      !> new run, topography read from topo
      !! start from climatology
      INTEGER, PARAMETER :: istart_new_topo=2
      !> continuing run (default)
      INTEGER, PARAMETER :: istart_run_continue=3
      INTEGER :: istart
      !> timestep length in seconds
      REAL(wp) :: dt

      REAL(wp) creltem, crelsal, DTI, rleadclose(3)
      INTEGER KBM,KM,MATR,NMAX,imean &
           ,LYEARS, LYEAR1,LYEAR2,LYEAR,LMONTH &
           ,LMONTS,LMONT1,LMONT2 &
           ,LDAYS,LDAY,LDAY1,LDTDAYC &
           ,IOASISFLUX,nfixYearLen,ndtday

      LOGICAL lforcediag,lhfldiag,lconvdiag,ldiffdiag,lgmdiag, &
           lcalcdifi,ltstranspose,ltswrite,lmpitype,lnonblock, &
           lbounds_exch_tp,ltidal,ltidal_diag,lisopyc, lswr_jerlov, &
           lwith_barotropic_stokes_drift

      INTEGER IOCAD,iocaduv,ICONTRO,IOCONV, ibbl_transport

      LOGICAL :: ladpo, ladfs   ! flags for advection routine (ocadpo or oadfs)
                                ! to be called depending on specified iocad value

      !> ALMost ZERo: used in various places to prevent division by
      !! zero errors
      REAL(wp), PARAMETER :: ALMZER = 1.E-19_wp
      REAL(wp) :: CONO,CONN, &
     &   AULAPUV,AUS, &
     &   CAH00,CSTABEPS,CAVOCON,CDVOCON
      INTEGER :: ibolk
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      REAL(wp), ALLOCATABLE, TARGET :: DLXP(:,:)   !< grid_x_distance_at_pressure_point
      REAL(wp), ALLOCATABLE, TARGET :: DLYP(:,:)   !< grid_y_distance_at_pressure_point
      REAL(wp), ALLOCATABLE, TARGET :: DLXU(:,:)   !< grid_x_distance_at_u_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DLYU(:,:)   !< grid_y_distance_at_u_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DLYV(:,:)   !< grid_y_distance_at_v_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DLXV(:,:)   !< grid_x_distance_at_v_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DLXPSI(:,:) !< grid_x_distance_at_psi_point
      REAL(wp), ALLOCATABLE, TARGET :: DLYPSI(:,:) !< grid_y_distance_at_psi_point
      REAL(wp), ALLOCATABLE, TARGET :: DDUE(:,:,:) !< ocean_level_thickness_at_v_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DDUO(:,:,:) !< ocean_level_thickness_at_u_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DDPO(:,:,:) !< ocean_level_thickness_at_pressure_point
      REAL(wp), ALLOCATABLE, TARGET :: DEPTO(:,:)  !< depth_at_pressure_point
      REAL(wp), ALLOCATABLE, TARGET :: DEUTE(:,:)  !< depth_at_v_vector_point
      REAL(wp), ALLOCATABLE, TARGET :: DEUTO(:,:)  !< depth_at_u_vector_point

      REAL(wp), ALLOCATABLE :: deution(:,:),deutien(:,:) &
                           ,surlon(:,:),surlen(:,:)

      REAL(wp), ALLOCATABLE :: AREA(:,:)
      REAL(wp), ALLOCATABLE :: DTDXUO(:,:)
      REAL(wp), ALLOCATABLE :: DTDXUE(:,:)
      REAL(wp), ALLOCATABLE :: DTDYO(:,:)
      REAL(wp), ALLOCATABLE :: DTDXPO(:,:)
      REAL(wp), ALLOCATABLE :: DTDXPE(:,:)
      REAL(wp), ALLOCATABLE :: DPYE(:,:)
      REAL(wp), ALLOCATABLE :: DPYO(:,:)
      REAL(wp), ALLOCATABLE :: DPIO(:,:,:)
      REAL(wp), ALLOCATABLE :: DDPSIO(:,:,:)
      REAL(wp), ALLOCATABLE :: DWI(:)
      REAL(wp), ALLOCATABLE :: DEUTIE(:,:)
      REAL(wp), ALLOCATABLE :: DEUTIO(:,:)
      REAL(wp), ALLOCATABLE :: AREAIN(:,:)

      REAL(wp), ALLOCATABLE :: dlxpi(:,:)   !< inverse grid_x_distance_at_pressure_point
      REAL(wp), ALLOCATABLE :: dlypi(:,:)   !< inverse grid_y_distance_at_pressure_point
      REAL(wp), ALLOCATABLE :: dlxui(:,:)   !< inverse grid_x_distance_at_u_vector_point
      REAL(wp), ALLOCATABLE :: dlyui(:,:)   !< inverse grid_y_distance_at_u_vector_point
      REAL(wp), ALLOCATABLE :: dlyvi(:,:)   !< inverse grid_y_distance_at_v_vector_point
      REAL(wp), ALLOCATABLE :: dlxvi(:,:)   !< inverse grid_x_distance_at_v_vector_point

      REAL(wp), ALLOCATABLE :: DZ(:),DZW(:),TIESTU(:),TIESTW(:),DI(:)

!:: GEOGRAPHIC POSITION:
      REAL(wp), ALLOCATABLE :: GILA(:,:),GIPH(:,:)
      REAL(wp), ALLOCATABLE, TARGET :: ALAT(:,:),ALON(:,:)
      REAL(wp), ALLOCATABLE, TARGET :: ALATU(:,:),ALONU(:,:)
      REAL(wp), ALLOCATABLE, TARGET :: ALATV(:,:),ALONV(:,:)

      REAL(wp), ALLOCATABLE :: ALAT_G(:,:),ALON_G(:,:)
      REAL(wp), ALLOCATABLE :: ALATPSI_G(:,:),ALONPSI_G(:,:)

!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      REAL(wp), ALLOCATABLE, TARGET :: WETO(:,:,:),AMSUE(:,:,:),AMSUO(:,:,:)
#ifdef WETO_STORE_LOGICAL
      LOGICAL, ALLOCATABLE :: lwetol1_g(:,:)
#elif defined WETO_STORE_BITVECTOR
      INTEGER, ALLOCATABLE :: bvwetol1_g(:)
#else
      REAL(wp), ALLOCATABLE :: WETOL1_G(:,:)
#endif

      REAL(wp), ALLOCATABLE :: DLXP_G(:,:),DLYP_G(:,:)
      REAL(wp), ALLOCATABLE :: DLXU_G(:,:),DLYU_G(:,:)
      REAL(wp), ALLOCATABLE :: DLXV_G(:,:),DLYV_G(:,:)
!
!:: ONE-DIMENSIONAL FIELDS
      REAL(wp), ALLOCATABLE :: TAF(:),SAF(:)  &
      ,PREFF(:),PREFFW(:)                                     &
      ,SKAL(:),B(:),X(:)
!
!:: three-dimensional fields
      REAL(wp), ALLOCATABLE, TARGET ::     rhoo(:,:,:)  !< sea_water_density
      REAL(wp), ALLOCATABLE, TARGET ::     uko(:,:,:)   !< sea_water_x_velocity
      REAL(wp), ALLOCATABLE, TARGET ::     vke(:,:,:)   !< sea_water_y_velocity
      REAL(wp), ALLOCATABLE, TARGET ::     tho(:,:,:)   !< sea_water_potential_temperature
      REAL(wp), ALLOCATABLE, TARGET ::     sao(:,:,:)   !< sea_water_salinity
      REAL(wp), ALLOCATABLE, TARGET ::     wo(:,:,:)    !< upward_sea_water_velocity
      REAL(wp), ALLOCATABLE, TARGET ::     po(:,:,:)    !< sea_water_pressure
      REAL(wp), ALLOCATABLE, TARGET ::     dvo(:,:,:)   !< sea_water_vertical_diffusivity
      REAL(wp), ALLOCATABLE, TARGET ::     avo(:,:,:)   !< sea_water_vertical_viscosity
      REAL(wp), ALLOCATABLE, TARGET ::     uaccel(:,:,:)!< sea_water_x_acceleration
      REAL(wp), ALLOCATABLE, TARGET ::     vaccel(:,:,:)!< sea_water_y_acceleration

      REAL(wp), ALLOCATABLE         ::     vk1e(:,:,:)
      REAL(wp), ALLOCATABLE         ::     uk1o(:,:,:)
      REAL(wp), ALLOCATABLE         ::     voe(:,:,:)
      REAL(wp), ALLOCATABLE         ::     uoo(:,:,:)
      REAL(wp), ALLOCATABLE         ::     t1o(:,:,:)
      REAL(wp), ALLOCATABLE         ::     s1o(:,:,:)
      REAL(wp), ALLOCATABLE         ::     stabio(:,:,:)
!
!:: TWO-DIMENSIONAL FIELDS
      REAL(wp), ALLOCATABLE, TARGET :: ZO(:,:)          !< sea_surface_height_above_sea_level
      REAL(wp), ALLOCATABLE, TARGET :: Z1O(:,:)         !< sea_surface_height_above_sea_level_change

      REAL(wp), ALLOCATABLE :: &
     &   Z1E(:,:), &
     &   U1E(:,:),V1E(:,:),U1O(:,:),V1O(:,:), &
     &   USO(:,:),VSE(:,:),UCOS(:,:),VCOS(:,:), &
     &   TLOW(:,:),SLOW(:,:), &
     &   PXOIN(:,:),PYEIN(:,:),UZO(:,:),VZE(:,:), &
     &   FTWOU(:,:),FTWOV(:,:),CURVAV(:,:),B1E(:,:), &
     &   B1O(:,:)!,SURDIS(:,:)



      INTEGER, ALLOCATABLE :: NUM(:,:),KCONDEP(:,:),KBOT(:,:)
      INTEGER, ALLOCATABLE :: NUM_G(:,:)
!
!:: ICE-MODEL FIELDS
      REAL(wp), ALLOCATABLE,TARGET  :: SICTHO(:,:)   !< sea_ice_thickness
      REAL(wp), ALLOCATABLE,TARGET  :: SICOMO(:,:)   !< sea_ice_area_fraction
      REAL(wp), ALLOCATABLE,TARGET  :: SICUO(:,:)    !< sea_ice_x_velocity
      REAL(wp), ALLOCATABLE,TARGET  :: SICVE(:,:)    !< sea_ice_y_velocity
      REAL(wp), ALLOCATABLE,TARGET  :: TICE(:,:)     !< ice_surface_temperature
      REAL(wp), ALLOCATABLE,TARGET  :: HIBZETO(:,:)  !< fixme
      REAL(wp), ALLOCATABLE,TARGET  :: HIBETO(:,:)   !< fixme
      REAL(wp), ALLOCATABLE,TARGET  :: HIBZETE(:,:)  !< fixme
      REAL(wp), ALLOCATABLE,TARGET  :: HIBETE(:,:)   !< fixme
      REAL(wp), ALLOCATABLE,TARGET  :: SICSNO(:,:)   !< surface_snow_thickness_where_sea_ice

      REAL(wp), ALLOCATABLE         :: SICOMP(:,:)
      REAL(wp), ALLOCATABLE         :: HIBDELO(:,:)
      REAL(wp), ALLOCATABLE         :: HIBDELE(:,:)
      REAL(wp), ALLOCATABLE,TARGET  :: TAUWATU(:,:)
      REAL(wp), ALLOCATABLE,TARGET  :: TAUWATV(:,:)
      REAL(wp), ALLOCATABLE         :: FWO(:,:)


!
!:: FORCING FIELDS
      REAL(wp), ALLOCATABLE,TARGET :: TXO(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: TYE(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: TAFO(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: EMINPO(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: RELSAO(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: RELTHO(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FCLOU(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FSWR(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FU10(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FPREC(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FTDEW(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: RIVRUN(:,:)

      REAL(wp), ALLOCATABLE,TARGET :: GIRIV(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: FSLP(:,:) !< sea level pressure [Pa]


      REAL(wp), ALLOCATABLE,TARGET :: WGO(:,:,:)
      REAL(wp), ALLOCATABLE,TARGET :: BOLX(:,:)
      REAL(wp), ALLOCATABLE,TARGET :: BOLY(:,:)

      CONTAINS

      SUBROUTINE alloc_mem_forcing

        !:: forcing fields
        ALLOCATE( txo(ie,je),tye(ie,je),tafo(ie,je),      &
             fclou(ie,je),fswr(ie,je),fu10(ie,je),        &
             fprec(ie,je),ftdew(ie,je),fslp(ie,je))

      END SUBROUTINE alloc_mem_forcing

      SUBROUTINE alloc_mem_stokes_drift

!     These fields are needed with namelist option lwith_barotropic_stokes_drift=.true.

      USE MO_MPI

        ALLOCATE(deution(ie,je),deutien(ie,je) &
             ,surlon(ie,je), surlen(ie,je))

        deution = 0._wp
        deutien = 0._wp
        surlon = 0._wp
        surlen = 0._wp


      END SUBROUTINE alloc_mem_stokes_drift

      SUBROUTINE alloc_mem_commo1

      USE MO_MPI
#if defined WETO_STORE_BITVECTOR
      INTEGER(i8) :: bvweto_bits, bvweto_ints, bvweto_mod
      INTEGER, PARAMETER :: bits_per_bvweto = bit_size(bvwetol1_g(1))
#endif
      ALLOCATE (DZ(KE+1),DZW(KE),TIESTU(KE+1),TIESTW(KE+1),DI(KE+1))
!
!:: GRID DESCRIPTORS AND RELATED FIELDS
      ALLOCATE( DLXP(IE,JE),DLYP(IE,JE),DLXU(IE,JE),DLYV(IE,JE),        &
     &  DLXV(IE,JE),DLYU(IE,JE),DTDXUO(IE,JE),                          &
     &  DTDXUE(IE,JE),DTDYO(IE,JE),DTDXPO(IE,JE),DTDXPE(IE,JE),         &
     &  DPYE(IE,JE),DPYO(IE,JE),                                        &
     &  DDUE(IE,JE,KE),DDUO(IE,JE,KE),DDPO(IE,JE,KE),DPIO(IE,JE,KE),    &
     &  DDPSIO(IE,JE,KE),DWI(KE),                                       &
     &  DEUTE(IE,JE),DEUTO(IE,JE),DEUTIE(IE,JE),DEUTIO(IE,JE),          &
     &  DEPTO(IE,JE),DLXPSI(IE,JE),DLYPSI(IE,JE), AREA(IE,JE),          &
     &  AREAIN(IE,JE),                                                  &
     &  dlxpi(ie,je), dlypi(ie,je), dlxui(ie,je), dlyui(ie,je),         &
     &  dlxvi(ie,je), dlyvi(ie,je) )

!:: GEOGRAPHIC POSITION:
      ALLOCATE(GILA(ITO,JTO),GIPH(ITO,JTO)      &
              ,ALAT(IE,JE),ALON(IE,JE)          &
              ,ALATU(IE,JE),ALONU(IE,JE)        &
              ,ALATV(IE,JE),ALONV(IE,JE))

!:: LAND SEA MASKS ON SKALAR/VECTOR POINTS
      ALLOCATE(WETO(IE,JE,KE),AMSUE(IE,JE,KE),AMSUO(IE,JE,KE))

#ifdef WETO_STORE_LOGICAL
      ALLOCATE(lwetol1_g(ie_g, je_g))
#elif defined WETO_STORE_BITVECTOR
      bvweto_bits = INT(ie_g, i8) * INT(je_g, i8)
      bvweto_ints = bvweto_bits / INT(bits_per_bvweto, i8)
      bvweto_mod = MODULO(bvweto_bits, INT(bits_per_bvweto, i8))
      IF (bvweto_mod .GT. 0_i8) bvweto_ints = bvweto_ints + 1_i8
      ALLOCATE(bvwetol1_g(0:bvweto_ints - 1))
#else
      ALLOCATE(WETOL1_G(IE_G,JE_G))
#endif


      IF ( p_pe == p_io ) THEN
        ALLOCATE(alat_g(ie_g,je_g),alon_g(ie_g,je_g),                 &
             alatpsi_g(1:ie_g,0:je_g),alonpsi_g(1:ie_g,0:je_g),       &
             dlxp_g(ie_g,je_g),dlyp_g(ie_g,je_g),                     &
             dlxu_g(ie_g,je_g),dlyu_g(ie_g,je_g),                     &
             dlxv_g(ie_g,je_g),dlyv_g(ie_g,je_g))
      ELSE
#ifdef MESSY
         ALLOCATE(alat_g(0,0),alon_g(0,0),                            &
              alatpsi_g(0,0),alonpsi_g(0,0),                          &
              dlxp_g(ie_g,je_g),dlyp_g(ie_g,je_g),                    &
              dlxu_g(ie_g,je_g),dlyu_g(ie_g,je_g),                    &
              dlxv_g(ie_g,je_g),dlyv_g(ie_g,je_g))
#else
         ALLOCATE(alat_g(0,0),alon_g(0,0),                            &
              alatpsi_g(0,0),alonpsi_g(0,0),                          &
              dlxp_g(0,0),dlyp_g(0,0),                                &
              dlxu_g(0,0),dlyu_g(0,0),                                &
              dlxv_g(0,0),dlyv_g(0,0))

#endif
      ENDIF

!
!:: ONE-DIMENSIONAL FIELDS
      ALLOCATE(TAF(KE),SAF(KE),PREFF(KE),       &
     &   PREFFW(KEP),SKAL(ILL),B(ILT),X(ILT))
!
!:: three-dimensional fields
      ALLOCATE( rhoo(ie,je,ke),uko(ie,je,ke),vke(ie,je,ke),             &
     &   tho(ie,je,ke),sao(ie,je,ke),wo(ie,je,kep),                     &
     &   po(ie,je,ke),vk1e(ie,je,ke),                                   &
     &   uk1o(ie,je,ke),voe(ie,je,ke),uoo(ie,je,ke),                    &
     &   t1o(ie,je,ke),s1o(ie,je,ke),dvo(ie,je,kep),avo(ie,je,kep),     &
     &   stabio(ie,je,ke),uaccel(ie,je,ke),vaccel(ie,je,ke))
!
!:: TWO-DIMENSIONAL FIELDS
      ALLOCATE( ZO(IE,JE),Z1E(IE,JE),Z1O(IE,JE),                        &
     & U1E(IE,JE),V1E(IE,JE),U1O(IE,JE),V1O(IE,JE),                     &
     & USO(IE,JE),VSE(IE,JE),UCOS(IE,JE),VCOS(IE,JE),                   &
     & TLOW(IE,JE),SLOW(IE,JE),                                         &
     & PXOIN(IE,JE),PYEIN(IE,JE),UZO(IE,JE),VZE(IE,JE),                 &
     & FTWOU(IE,JE),FTWOV(IE,JE),CURVAV(IE,JE),B1E(IE,JE),              &
     & B1O(IE,JE),NUM(IE,JE)                                &
     &,KCONDEP(IE,JE),KBOT(IE,JE))
!      ALLOCATE(SURDIS(IE,JE))

#ifndef SOR
      ALLOCATE(NUM_G(IE_G,JE_G))
#endif
!
!:: ICE-MODEL FIELDS
       ALLOCATE(SICTHO(IE,JE),SICOMO(IE,JE), SICOMP(IE,JE)              &
               ,SICUO(IE,JE),SICVE(IE,JE),TICE(IE,JE)                   &
               ,HIBDELO(IE,JE),HIBZETO(IE,JE),HIBETO(IE,JE)             &
               ,HIBDELE(IE,JE),HIBZETE(IE,JE),HIBETE(IE,JE)             &
               ,TAUWATU(IE,JE),TAUWATV(IE,JE)                           &
               ,SICSNO(IE,JE),FWO(IE,JE))

!:: RESTORING FIELDS
        ALLOCATE(RELSAO(IE,JE),RELTHO(IE,JE))

!  FIELDS FOR COMMON DIAGNOSTICS
        ALLOCATE(EMINPO(IE,JE),RIVRUN(IE,JE))


!<===END COMMO1
      END SUBROUTINE alloc_mem_commo1

      SUBROUTINE alloc_mem_gmbolus

        INTEGER :: i, j, k

        ALLOCATE(WGO(IE,JE,KEP),BOLX(IE,JE),BOLY(IE,JE))

        DO k=1,kep
          DO j=1,je
            DO i=1,ie
              wgo(i,j,k)=zero
            ENDDO
          ENDDO
        ENDDO
        DO j=1,je
          DO i=1,ie
            bolx(i,j)=zero
            boly(i,j)=zero
          ENDDO
        ENDDO

      END SUBROUTINE alloc_mem_gmbolus

#if defined WETO_STORE_LOGICAL || defined WETO_STORE_BITVECTOR

      ELEMENTAL FUNCTION wetol1_g(i, j)
        REAL(wp) :: wetol1_g
        INTEGER, INTENT(in) :: i, j
        wetol1_g = MERGE(1.0, 0.0, lwetol1_g(i, j))
      END FUNCTION wetol1_g

#endif

#if defined WETO_STORE_LOGICAL

      SUBROUTINE set_wetol1_g(i, j, is_water)
        INTEGER, INTENT(in) :: i, j
        LOGICAL, INTENT(in) :: is_water
        lwetol1_g(i, j) = is_water
      END SUBROUTINE set_wetol1_g

#elif defined WETO_STORE_BITVECTOR

      SUBROUTINE set_wetol1_g(i, j, is_water)
        INTEGER, INTENT(in) :: i, j
        LOGICAL, INTENT(in) :: is_water
        INTEGER(i8) :: int_index, bit_index
        int_index = INT(i - 1, i8) + INT(j - 1, i8) * INT(ie_g, i8)
        bit_index = MOD(int_index, INT(BIT_SIZE(bvwetol1_g(1)), i8))
        int_index = int_index / INT(BIT_SIZE(bvwetol1_g(1)), i8)
        IF (is_water) THEN
          bvwetol1_g(int_index) = IBSET(bvwetol1_g(int_index), bit_index)
        ELSE
          bvwetol1_g(int_index) = IBCLR(bvwetol1_g(int_index), bit_index)
        END IF
      END SUBROUTINE set_wetol1_g

#else

      SUBROUTINE set_wetol1_g(i, j, is_water)
        INTEGER, INTENT(in) :: i, j
        LOGICAL, INTENT(in) :: is_water
        wetol1_g(i, j) = MERGE(1.0_wp, 0.0_wp, is_water)
      END SUBROUTINE set_wetol1_g

#endif


#if defined WETO_STORE_LOGICAL
#elif defined WETO_STORE_BITVECTOR
      ELEMENTAL FUNCTION lwetol1_g(i,j)
        LOGICAL :: lwetol1_g
        INTEGER, INTENT(in) :: i, j
        INTEGER(i8) :: int_index, bit_index
        int_index = INT(i - 1, i8) + INT(j - 1, i8) * INT(ie_g, i8)
        bit_index = MOD(int_index, INT(BIT_SIZE(bvwetol1_g(1)), i8))
        int_index = int_index / INT(BIT_SIZE(bvwetol1_g(1)), i8)
        lwetol1_g = BTEST(bvwetol1_g(int_index), bit_index)
      END FUNCTION lwetol1_g

#else
      ELEMENTAL FUNCTION lwetol1_g(i,j)
        LOGICAL :: lwetol1_g
        INTEGER, INTENT(in) :: i, j
        lwetol1_g = wetol1_g(i, j) .GT. 0.5_wp
      END FUNCTION lwetol1_g
#endif

      ELEMENTAL FUNCTION iwetol1_g(i, j)
        INTEGER :: iwetol1_g
        INTEGER, INTENT(in) :: i, j
        iwetol1_g = MERGE(1, 0, lwetol1_g(i, j))
      END FUNCTION iwetol1_g

      ELEMENTAL FUNCTION lweto(i, j, k)
        LOGICAL :: lweto
        INTEGER, INTENT(in) :: i, j, k
        lweto = weto(i, j, k) .GT. 0.5_wp
      END FUNCTION lweto

      ELEMENTAL FUNCTION iweto(i, j, k)
        INTEGER :: iweto
        INTEGER, INTENT(in) :: i, j, k
        iweto = INT(weto(i, j, k))
      END FUNCTION iweto


      ELEMENTAL FUNCTION istart_in_valid_range(istart) RESULT(valid)
        INTEGER, INTENT(in) :: istart
        LOGICAL :: valid
        valid = istart .GE. istart_new_topo_update &
             .AND. istart .LE. istart_run_continue
      END FUNCTION istart_in_valid_range

    END MODULE MO_COMMO1
