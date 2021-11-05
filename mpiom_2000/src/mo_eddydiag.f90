!> Modul: eddydiag
!> @author H.Haak,JS.vStorch
!> @date 06.10.2008

MODULE mo_eddydiag

  USE mo_kind, ONLY: wp
  USE mo_parallel, ONLY : have_g_js, p_ioff
  USE mo_param1, ONLY : ie,je,ke, ie_g
  USE mo_commo1, ONLY : sictho,sicomo,sicsno,sicuo,sicve,wo,uko,vke  &
                       ,amsuo,amsue,sao,tho,weto,zo,lbounds_exch_tp  &
                       ,po,rhoo
  USE mo_boundsexch

  IMPLICIT NONE

  LOGICAL :: leddydiag

  REAL(wp), ALLOCATABLE, TARGET :: ut0(:,:,:)     !< ut0 : product_of_sea_water_x_velocity_and_sea_water_potential_temperature
  REAL(wp), ALLOCATABLE, TARGET :: us0(:,:,:)     !< us0 : product_of_sea_water_x_velocity_and_sea_water_salinity
  REAL(wp), ALLOCATABLE, TARGET :: uu0(:,:,:)     !< uu0 : square_of_sea_water_x_velocity
  REAL(wp), ALLOCATABLE, TARGET :: uv0(:,:,:)     !< uv0 : product_of_sea_water_x_velocity_and_sea_water_y_velocity
  REAL(wp), ALLOCATABLE, TARGET :: uw0(:,:,:)     !< uw0 : product_of_sea_water_x_velocity_and_upward_sea_water_velocity
  REAL(wp), ALLOCATABLE, TARGET :: vt0(:,:,:)     !< vt0 : product_of_sea_water_y_velocity_and_sea_water_potential_temperature
  REAL(wp), ALLOCATABLE, TARGET :: vs0(:,:,:)     !< vs0 : product_of_sea_water_y_velocity_and_sea_water_salinity
  REAL(wp), ALLOCATABLE, TARGET :: vu0(:,:,:)     !< vu0 : product_of_sea_water_y_velocity_and_sea_water_x_velocity
  REAL(wp), ALLOCATABLE, TARGET :: vv0(:,:,:)     !< vv0 : square_of_sea_water_y_velocity
  REAL(wp), ALLOCATABLE, TARGET :: vw0(:,:,:)     !< vw0 : product_of_sea_water_y_velocity_and_upward_sea_water_velocity
  REAL(wp), ALLOCATABLE, TARGET :: wt0(:,:,:)     !< wt0 : product_of_upward_sea_water_velocity_and_sea_water_potential_temperature
  REAL(wp), ALLOCATABLE, TARGET :: ws0(:,:,:)     !< ws0 : product_of_upward_sea_water_velocity_and_sea_water_salinity
  REAL(wp), ALLOCATABLE, TARGET :: wu0(:,:,:)     !< wu0 : product_of_upward_sea_water_velocity_and_sea_water_x_velocity
  REAL(wp), ALLOCATABLE, TARGET :: wv0(:,:,:)     !< wv0 : product_of_upward_sea_water_velocity_and_sea_water_y_velocity
  REAL(wp), ALLOCATABLE, TARGET :: ww0(:,:,:)     !< ww0 : square_of_upward_sea_water_velocity

  REAL(wp), ALLOCATABLE, TARGET :: tt0(:,:,:)     !< tt0 : square_of_sea_water_potential_temperature
  REAL(wp), ALLOCATABLE, TARGET :: ss0(:,:,:)     !< ss0 : square_of_sea_water_salinity
  REAL(wp), ALLOCATABLE, TARGET :: pp0(:,:,:)     !< pp0 : square_of_sea_water_pressure
  REAL(wp), ALLOCATABLE, TARGET :: rr0(:,:,:)     !< rr0 : square_of_sea_water_density

  REAL(wp), ALLOCATABLE, TARGET :: uic0(:,:)      !< uic0 : product_of_sea_ice_x_velocity_and_sea_ice_concentration
  REAL(wp), ALLOCATABLE, TARGET :: vic0(:,:)      !< vic0 : product_of_sea_ice_y_velocity_and_sea_ice_concentration
  REAL(wp), ALLOCATABLE, TARGET :: uih0(:,:)      !< uih0 : product_of_sea_ice_x_velocity_and_sea_ice_thickness
  REAL(wp), ALLOCATABLE, TARGET :: vih0(:,:)      !< vih0 : product_of_sea_ice_y_velocity_and_sea_ice_thickness
  REAL(wp), ALLOCATABLE, TARGET :: uisn0(:,:)     !< uisn0 : product_of_sea_ice_x_velocity_and_snow_thickness
  REAL(wp), ALLOCATABLE, TARGET :: visn0(:,:)     !< visn0 : product_of_sea_ice_y_velocity_and_snow_thickness



CONTAINS

    !> Calulates eddy statistics
    !> @return eddy statistics

  SUBROUTINE alloc_eddydiag

    ALLOCATE(ut0(ie,je,ke),us0(ie,je,ke),uu0(ie,je,ke))
    ALLOCATE(uv0(ie,je,ke),uw0(ie,je,ke),vt0(ie,je,ke))
    ALLOCATE(vs0(ie,je,ke),vu0(ie,je,ke),vv0(ie,je,ke))
    ALLOCATE(vw0(ie,je,ke),wt0(ie,je,ke+1),ws0(ie,je,ke+1))
    ALLOCATE(wu0(ie,je,ke+1),wv0(ie,je,ke+1),ww0(ie,je,ke+1))

    ALLOCATE(tt0(ie,je,ke),ss0(ie,je,ke),pp0(ie,je,ke),rr0(ie,je,ke))


    ALLOCATE(uic0(ie,je),uih0(ie,je),uisn0(ie,je))
    ALLOCATE(vic0(ie,je),vih0(ie,je),visn0(ie,je))

    ut0(:,:,:) = 0._wp
    us0(:,:,:) = 0._wp
    uu0(:,:,:) = 0._wp
    uv0(:,:,:) = 0._wp
    uw0(:,:,:) = 0._wp
    vt0(:,:,:) = 0._wp
    vs0(:,:,:) = 0._wp
    vu0(:,:,:) = 0._wp
    vv0(:,:,:) = 0._wp
    vw0(:,:,:) = 0._wp
    wt0(:,:,:) = 0._wp
    ws0(:,:,:) = 0._wp
    wu0(:,:,:) = 0._wp
    wv0(:,:,:) = 0._wp
    ww0(:,:,:) = 0._wp

    tt0(:,:,:) = 0._wp
    ss0(:,:,:) = 0._wp
    pp0(:,:,:) = 0._wp
    rr0(:,:,:) = 0._wp

    uic0(:,:) = 0._wp
    vic0(:,:) = 0._wp
    uih0(:,:) = 0._wp
    vih0(:,:) = 0._wp
    uisn0(:,:) = 0._wp
    visn0(:,:) = 0._wp


    leddydiag = .FALSE.

  END SUBROUTINE alloc_eddydiag

  SUBROUTINE calc_eddydiag

    ! horizontal fluxes

    REAL(wp) :: fak(2:ie-1,2:je-1)
    INTEGER :: i,j,k,jb

    jb = MERGE(3,2,lbounds_exch_tp .and. have_g_js )

    uu0(:,:,:)=uko(:,:,:)*uko(:,:,:)*amsuo(:,:,:)
    vv0(:,:,:)=vke(:,:,:)*vke(:,:,:)*amsue(:,:,:)
    ww0(:,:,1:ke)=wo(:,:,1:ke)*wo(:,:,1:ke)*weto(:,:,1:ke)
    tt0(:,:,:)=tho(:,:,:)*tho(:,:,:)*weto(:,:,:)
    ss0(:,:,:)=sao(:,:,:)*sao(:,:,:)*weto(:,:,:)
    pp0(:,:,:)=po(:,:,:)*po(:,:,:)*weto(:,:,:)
    rr0(:,:,:)=rhoo(:,:,:)*rhoo(:,:,:)*weto(:,:,:)

    DO j=2,je-1
       DO i=2,ie-1
          fak(i,j) = MERGE (-1.0_wp, 1.0_wp, (lbounds_exch_tp .AND. have_g_js      &
               .AND. j <= 3 .AND. i+p_ioff >= ie_g/2))
       END DO
    END DO

    DO k=1,ke
       DO j=2,je-1
          DO i=2,ie-1

          ut0(i,j,k)=uko(i,j,k)*0.5_wp*(tho(i,j,k)+tho(i+1,j,k))     &
               *amsuo(i,j,k)
          us0(i,j,k)=uko(i,j,k)*0.5_wp*(sao(i,j,k)+sao(i+1,j,k))     &
               *amsuo(i,j,k)
          uw0(i,j,k)=uko(i,j,k)*0.25_wp*(wo(i,j,k)+wo(i+1,j,k)       &
               +wo(i,j,k+1)+wo(i+1,j,k+1))*amsuo(i,j,k)

          vt0(i,j,k)=vke(i,j,k)*0.5_wp*(tho(i,j,k)+tho(i,j+1,k))     &
               *amsue(i,j,k)
          vs0(i,j,k)=vke(i,j,k)*0.5_wp*(sao(i,j,k)+sao(i,j+1,k))     &
               *amsue(i,j,k)
          vw0(i,j,k)=vke(i,j,k)*0.25_wp*(wo(i,j,k)+wo(i,j+1,k)       &
               +wo(i,j,k+1)+wo(i,j+1,k+1))*amsue(i,j,k)

          vu0(i,j,k)=vke(i,j,k)*0.25_wp*(uko(i-1,j,k)+uko(i,j,k)     &
               +uko(i-1,j-1,k)+uko(i,j-1,k))*amsue(i,j,k)

          uv0(i,j,k)=uko(i,j,k)*0.25_wp*(vke(i-1,j,k)+fak(i,j)*vke(i-1,j-1,k)    &
               +vke(i,j,k)+fak(i,j)*vke(i,j-1,k))*amsuo(i,j,k)

        ENDDO
      ENDDO
    ENDDO

    DO j=2,je-1
       DO i=2,ie-1

        uih0(i,j)=sicuo(i,j)*0.5_wp*(sictho(i,j)+sictho(i+1,j))*amsuo(i,j,1)

        uic0(i,j)=sicuo(i,j)*0.5_wp*(sicomo(i,j)+sicomo(i+1,j))*amsuo(i,j,1)

        uisn0(i,j)=sicuo(i,j)*0.5_wp*(sicsno(i,j)+sicsno(i+1,j))*amsuo(i,j,1)

        vih0(i,j)=sicve(i,j)*0.5_wp*(sictho(i,j)+sictho(i,j+1))*amsue(i,j,1)

        vic0(i,j)=sicve(i,j)*0.5_wp*(sicomo(i,j)+sicomo(i,j+1))*amsue(i,j,1)

        visn0(i,j)=sicve(i,j)*0.5_wp*(sicsno(i,j)+sicsno(i,j+1))*amsue(i,j,1)

      ENDDO
    ENDDO


    DO k=2,ke
       DO j=2,je-1
          DO i=2,ie-1

          wt0(i,j,k)=wo(i,j,k)*0.5_wp*(tho(i,j,k) + tho(i,j,k-1))*weto(i,j,k)
          ws0(i,j,k)=wo(i,J,k)*0.5_wp*(sao(i,j,k) + sao(i,j,k-1))*weto(i,j,k)
          wu0(i,j,k)=wo(i,j,k)*0.25_wp*(uko(i-1,j,k)+uko(i,j,k)              &
               +uko(i-1,j,k-1)+uko(i,j,k-1))*weto(i,j,k)
          wv0(i,j,k)=wo(i,j,k)*0.5_wp*(fak(i,j)*vke(i,j-1,k)+vke(i,j,k)           &
               +fak(i,j)*vke(i,j-1,k-1)+vke(i,j,k-1))*weto(i,j,k)

        ENDDO
      ENDDO
    ENDDO


    CALL bounds_exch(1,'u',uih0,'in calc_eddydiag 1')
    CALL bounds_exch(1,'u',uic0,'in calc_eddydiag 2')
    CALL bounds_exch(1,'u',uisn0,'in calc_eddydiag 3')
    CALL bounds_exch(1,'u',ut0,'in calc_eddydiag 4')
    CALL bounds_exch(1,'u',us0,'in calc_eddydiag 5')
    CALL bounds_exch(1,'u',uw0,'in calc_eddydiag 6')
    CALL bounds_exch(1,'u',uv0,'in calc_eddydiag 7')

    CALL bounds_exch(1,'v',vih0,'in calc_eddydiag 8')
    CALL bounds_exch(1,'v',vic0,'in calc_eddydiag 9')
    CALL bounds_exch(1,'v',visn0,'in calc_eddydiag 10')
    CALL bounds_exch(1,'v',vt0,'in calc_eddydiag 11')
    CALL bounds_exch(1,'v',vs0,'in calc_eddydiag 12')
    CALL bounds_exch(1,'v',vw0,'in calc_eddydiag 13')
    CALL bounds_exch(1,'v',vu0,'in calc_eddydiag 14')

    CALL bounds_exch(1,'p',wt0,'in calc_eddydiag 14')
    CALL bounds_exch(1,'p',ws0,'in calc_eddydiag 15')
    CALL bounds_exch(1,'p',wu0,'in calc_eddydiag 16')
    CALL bounds_exch(1,'p',wv0,'in calc_eddydiag 17')

END SUBROUTINE calc_eddydiag



END MODULE mo_eddydiag
