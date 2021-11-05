MODULE mo_ncar_ocean_fluxes

  USE mo_kind, ONLY: i4, i8, wp
  USE mo_parallel, ONLY: p_pe, p_io, read_slice, global_sum, spool_slice
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_param1, ONLY: ie, je
  USE mo_constants, ONLY: api
  USE mo_units, ONLY: io_stdout
  USE mo_commo1, ONLY: lmont1, lday1, lyear1, zo, sao, eminpo, &
       dlxp, dlyp, sictho, sicsno, ddpo, weto, dt
  USE mo_model_time, ONLY: monlen_sum
  USE mo_commoau1, ONLY: tfreeze, tmelt, clb, con
  USE mo_commoau2, ONLY: prech
  USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa
  USE mo_io_config, ONLY: next_free_unit

  IMPLICIT NONE
  PRIVATE

  REAL(wp),ALLOCATABLE,DIMENSION(:,:) :: coru10,corv10,cort10,       &
                                     corq10,corqdlw,corqdsw,     &
                                     corprec,corrunoff
  INTEGER :: io_in_coru10, io_in_corv10, io_in_cort10, io_in_corq10, &
       io_in_corlw, io_in_corsw, io_in_corpr, io_in_corro

  ! hh tune albedos for ncep with Koch
  REAL(wp), PARAMETER :: albi = 0.76_wp
  REAL(wp), PARAMETER :: albm = 0.71_wp
  !> albedo of sea water
  REAL(wp), PARAMETER :: albw = 0.066_wp
  !> snow < 0
  REAL(wp), PARAMETER :: albsn = 0.86_wp
  !> snow > 0
  REAL(wp), PARAMETER :: albsnm = 0.72_wp

  PUBLIC :: albi, albm, albw, albsn, albsnm

  PUBLIC :: alloc_mem_core, budget_ocean_core, budget_ice_core, &
       open_core, read_core, spool_core, normpem, corrunoff

CONTAINS

  SUBROUTINE alloc_mem_core

    ALLOCATE(coru10(ie,je),corv10(ie,je),cort10(ie,je),corq10(ie,je),    &
         corqdlw(ie,je),corqdsw(ie,je),corprec(ie,je),corrunoff(ie,je))

  END SUBROUTINE alloc_mem_core

  SUBROUTINE open_core

    INTEGER :: i,j

    io_in_coru10 = next_free_unit()
    OPEN(io_in_coru10,file='CORU10',form='unformatted')

    io_in_corv10 = next_free_unit()
    OPEN(io_in_corv10,file='CORV10',form='unformatted')

    io_in_cort10 = next_free_unit()
    OPEN(io_in_cort10,file='CORT10',form='unformatted')

    io_in_corq10 = next_free_unit()
    OPEN(io_in_corq10,file='CORQ10',form='unformatted')

    io_in_corlw = next_free_unit()
    OPEN(io_in_corlw,file='CORQDLW',form='unformatted')

    io_in_corsw = next_free_unit()
    OPEN(io_in_corsw,file='CORQDSW',form='unformatted')

    io_in_corpr = next_free_unit()
    OPEN(io_in_corpr,file='CORPREC',form='unformatted')

    io_in_corro = next_free_unit()
    OPEN(io_in_corro,file='CORRIV',form='unformatted')

    DO j=1,je
       DO i=1,ie
          coru10(i,j) = 0._wp      ! 10m u-wind speed [m s-1]
          corv10(i,j) = 0._wp      ! 10m v-wind speed [m s-1]
          cort10(i,j) = 0._wp      ! 10m air temperature [C]
          corq10(i,j) = 0._wp      ! 10m humidity [??]
          corqdlw(i,j) = 0._wp     ! downward longwave[W m-2]
          corqdsw(i,j) = 0._wp     ! downward shortwave[W m-2]
          corprec(i,j) = 0._wp     ! precipitation [kg m-2]
          corrunoff(i,j) = 0._wp   ! runoff [kg m-2]
       ENDDO
    ENDDO

  END SUBROUTINE open_core

  SUBROUTINE spool_core(nread_per_day)
    INTEGER(i4), INTENT(IN) :: nread_per_day

    INTEGER(i8) :: ii1,ii2,ii3,ii4
    INTEGER     :: nrec_spool, lrec

    !spool the fields to the actual month

    IF ( lmont1 > 1 .OR. lday1 > 1 ) THEN
      nrec_spool=monlen_sum(1,lmont1-1,lyear1) + lday1 - 1
      DO lrec=1,nrec_spool*nread_per_day
             WRITE(IO_STDOUT,*)'in spool'

             IF(p_pe==p_io) THEN
                READ(io_in_coru10)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_coru10)

             IF(p_pe==p_io) THEN
                READ(io_in_corv10)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corv10)

             IF(p_pe==p_io) THEN
                READ(io_in_cort10)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_cort10)

             IF(p_pe==p_io) THEN
                READ(io_in_corq10)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corq10)

             IF(p_pe==p_io) THEN
                READ(io_in_corlw)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corlw)

             IF(p_pe==p_io) THEN
                READ(io_in_corsw)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corsw)

             IF(p_pe==p_io) THEN
                READ(io_in_corpr)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corpr)

             IF(p_pe==p_io) THEN
                READ(io_in_corro)ii1,ii2,ii3,ii4
             ENDIF
             CALL spool_slice(io_in_corro)

             IF(p_pe==p_io) THEN
                WRITE(io_stdout,*)'spool: record=', lrec, 'idate=', ii1
             ENDIF
       ENDDO
       IF(p_pe==p_io) THEN
          WRITE(0,*) 'forcing data is spooled to month ',lmont1, ' and day ',lday1
       ENDIF
    ENDIF

  END SUBROUTINE spool_core


  SUBROUTINE read_core

    INTEGER(i8) :: ii1,ii2,ii3,ii4

             WRITE(io_stdout,*)'in read'

             IF(p_pe==p_io) THEN
                READ(io_in_coru10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_coru10,coru10)

             IF(p_pe==p_io) THEN
                READ(io_in_corv10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corv10,corv10)

             IF(p_pe==p_io) THEN
                READ(io_in_cort10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_cort10,cort10)

             IF(p_pe==p_io) THEN
                READ(io_in_corq10)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corq10,corq10)

             IF(p_pe==p_io) THEN
                READ(io_in_corlw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corlw,corqdlw)

             IF(p_pe==p_io) THEN
                READ(io_in_corsw)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corsw,corqdsw)

             IF(p_pe==p_io) THEN
                READ(io_in_corpr)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corpr,corprec)


             IF(p_pe==p_io) THEN
                READ(io_in_corro)ii1,ii2,ii3,ii4
             ENDIF
             CALL read_slice(io_in_corro,corrunoff)

             WRITE(io_stdout,*)'after read'

  END SUBROUTINE read_core

  SUBROUTINE budget_ocean_core(uo,vo,to,qla,qse,qnsw,qnlw,qpre,qnet,fw,taux,tauy)
    REAL(wp),INTENT(in),DIMENSION(ie,je) :: uo,vo,to

    REAL(wp),INTENT(out),DIMENSION(ie,je) :: qpre,qla,qse,qnsw,qnlw,qnet
    REAL(wp),INTENT(out),DIMENSION(ie,je) :: fw,taux,tauy

    REAL(wp),DIMENSION(ie,je) :: evap,du,dv,uvdel,qo,z,ce,cd,ch,ustar,bstar
    ! FIXME: these should be parameters
    REAL(wp) :: ra,q1,q2,f1,vlamda,cp,flamda,tv,emiss,stebo
    INTEGER :: i,j,n_itts


    n_itts=2                    ! number of iterations in ncar_ocean_fluxes
                                ! ocean: n_itts=2 ; sea ice: n_itts=5
    ra = 1.22_wp                ! near surface air density [kg/m3]
    q1 = 640380._wp             ! coefficiant of q saturation function [kg/m3]
    q2 = -5107.4_wp             ! coefficiant of q saturation function [k]
    f1 = 0.98_wp                ! "saturation effect" applied over sea water
    flamda = 3.337e5_wp         ! latent heat of fusion [J kg-1]
    vlamda = 2.5e6_wp           ! latent heat of vaporisation [J kg-1]
    cp = 1005._wp               ! specific heat of air [J kg-1 K-1]
    stebo = 5.67e-8_wp          ! Stefan-Boltzmann Constant [W m-2 K-4]
    emiss = 1.0_wp              ! emissivity of sea water

    qpre(:,:) = 0._wp
    qla(:,:) = 0._wp
    qse(:,:) = 0._wp
    qnsw(:,:) = 0._wp
    qnlw(:,:) = 0._wp
    qnet(:,:) = 0._wp
    fw(:,:) = 0._wp
    taux(:,:) = 0._wp
    tauy(:,:) = 0._wp
    evap(:,:) = 0._wp
    du(:,:) = 0._wp
    dv(:,:) = 0._wp
    uvdel(:,:) = 0._wp
    qo(:,:) = 0._wp
    z(:,:) = 0._wp
    ce(:,:) = 0._wp
    cd(:,:) = 0._wp
    ch(:,:) = 0._wp
    ustar(:,:) = 0._wp
    bstar(:,:) = 0._wp

    DO j=2,je-1
       DO i=2,ie-1
          du(i,j)=(coru10(i,j)-uo(I,j))
          dv(i,j)=(corv10(i,j)-vo(i,j))
          uvdel(i,j)=SQRT(du(i,j)**2+dv(i,j)**2)
          ! L-Y eqn. 5
          qo(i,j) = (1._wp / ra) * f1 * q1 * EXP(q2 / to(i, j))
          z(i,j) = 10._wp
       ENDDO
    ENDDO

    CALL bounds_exch(1,'u',du,'mo_ncar_ocean_fluxes 1')
    CALL bounds_exch(1,'v',dv,'mo_ncar_ocean_fluxes 2')
    CALL bounds_exch(1,'p',uvdel,'mo_ncar_ocean_fluxes 3')
    CALL bounds_exch(1,'p',qo,'mo_ncar_ocean_fluxes 4')
    CALL bounds_exch(1,'p',z,'mo_ncar_ocean_fluxes 5')


    CALL ncar_ocean_fluxes(n_itts,uvdel,cort10,to,corq10,qo,z,cd,ch,ce,ustar,bstar)

    call bounds_exch(1,'p',ce,'mo_ncar_ocean_fluxes 6')
    call bounds_exch(1,'p',cd,'mo_ncar_ocean_fluxes 7')
    call bounds_exch(1,'p',ch,'mo_ncar_ocean_fluxes 8')

    DO j=2,je-1
       DO i=2,ie-1

!          ce(i,j)=1.75e-3

          ! calculate evaporation
          evap(i,j)=ra*ce(i,j)*(corq10(i,j)-qo(i,j))*uvdel(i,j)    ! L-Y eqn. 4a

          ! calculate precipitaion heat flux
          IF (cort10(i,j) > tmelt) THEN
             qpre(i,j)=flamda*corprec(i,j)                        ! L-Y eqn. 14
          ENDIF

          ! calculate latent heat flux
          qla(i,j)=vlamda*evap(i,j)                           ! L-Y eqn. 4b

          ! calculate sensible heat flux
          tv = cort10(i,j)

          qse(i,j)=ra*cp*ch(i,j)*(tv-to(i,j))*uvdel(i,j)  ! L-Y eqn. 4c

          ! calculate net short wave heat flux

          qnsw(i, j) = corqdsw(i, j) * (1._wp - albw)     ! L-Y eqn. 11


          ! calculate net long wave heat flux


          qnlw(i,j)=corqdlw(i,j)-emiss*stebo*to(i,j)**4             ! L-Y eqn. 12

          ! calculate net heat flux

          qnet(i,j)=-(qnsw(i,j)+qnlw(i,j) + qla(i,j)                      &
               +qse(i,j)  + qpre(i,j)) /clb          ! L-Y eqn. 3a

          ! calculate net fw flux
          fw(i,j)=corprec(i,j)
!                +evap(i,j)+corrunoff(i,j)

          ! calculate wind stress
          taux(i,j)=ra*cd(i,j)*uvdel(i,j)*du(i,j)             ! L-Y eqn. 4d
          tauy(i,j)=ra*cd(i,j)*uvdel(i,j)*dv(i,j)             ! L-Y eqn. 4d
       ENDDO
    ENDDO

    call  bounds_exch(1,'u',taux,'mo_ncar_ocean_fluxes 9')
    call  bounds_exch(1,'v',tauy,'mo_ncar_ocean_fluxes 10')
    call  bounds_exch(1,'p',fw,'mo_ncar_ocean_fluxes 11')
    call  bounds_exch(1,'p',qnet,'mo_ncar_ocean_fluxes 12')
    call  bounds_exch(1,'p',qnlw,'mo_ncar_ocean_fluxes 13')
    call  bounds_exch(1,'p',qnsw,'mo_ncar_ocean_fluxes 14')
    call  bounds_exch(1,'p',qse,'mo_ncar_ocean_fluxes 15')
    call  bounds_exch(1,'p',qla,'mo_ncar_ocean_fluxes 16')
    call  bounds_exch(1,'p',qpre,'mo_ncar_ocean_fluxes 17')
    call  bounds_exch(1,'p',evap,'mo_ncar_ocean_fluxes 18')

  END SUBROUTINE budget_ocean_core


  SUBROUTINE budget_ice_core(ui,vi,ti,sictho,          &
                             qla,qse,qnsw,qnlw,qres,qcon,fw,taux,tauy)

    IMPLICIT NONE

    REAL(wp),INTENT(in),DIMENSION(ie,je) :: ui,vi,ti,sictho

    REAL(wp),INTENT(out),DIMENSION(ie,je) :: qla,qse,qnsw,qnlw,qres,qcon
    REAL(wp),INTENT(out),DIMENSION(ie,je) :: fw,taux,tauy

    REAL(wp),DIMENSION(ie,je) :: evap,qnet,du,dv,uvdel,qo,z,ce,cd,ch!,ustar,bstar
    REAL(wp) :: ra,q1,q2,f1,vlamda,cp,flamda,tv,emiss,stebo
    REAL(wp) :: albvdr, albndr, albvdf, albndf
    INTEGER :: i,j,n_itts

    REAL(wp) :: tb

    n_itts=5                    ! number of iterations in ncar_ocean_fluxes
    ! ocean: n_itts=2 ; sea ice: n_itts=5 see L&Y page 10


    ra = 1.22_wp                ! near surface air density [kg/m3]
    q1 = 640380._wp             ! coefficiant of q saturation function [kg/m3]
    q2 = -5107.4_wp             ! coefficiant of q saturation function [k]
    f1 = 0.98_wp                ! "saturation effect" applied over sea water
    flamda = 3.337e5_wp         ! latent heat of fusion [J kg-1]
    vlamda = 2.5e6_wp           ! latent heat of vaporisation [J kg-1]
    cp = 1005._wp               ! specific heat of air [J kg-1 K-1]

    albvdr = 0.95_wp            ! Fraction of dsw 0.29 ; visible direct albedo of snow/sea ice
    albvdf = 0.85_wp            ! Fraction of dsw 0.31 ; visible difuse albedo of snow/sea ice
    albndr = 0.5_wp             ! Fraction of dsw 0.24 ; near infrared direct albedo of snow/sea ice
    albndf = 0.4_wp             ! Fraction of dsw 0.16 ; near infrared difuse albedo of snow/sea ice

    stebo = 5.67e-8_wp            ! Stefan-Boltzmann Constant [W m-2 K-4]
    emiss = 0.95_wp             ! emissivity of sea ice

    tb = tfreeze + tmelt

!    call flx_blk_albedo(ti,sictho,sicsno,palb,palcn,palbp,palcnp)
!           palb         , & !  albedo of ice under overcast sky
!           palcn        , & !  albedo of ocean under overcast sky
!           palbp        , & !  albedo of ice under clear sky
!           palcnp           !  albedo of ocean under clear sky


    DO j=2,je-1
       DO i=2,ie-1
          du(i,j)=(coru10(i,j)-ui(I,j))
          dv(i,j)=(corv10(i,j)-vi(i,j))
          uvdel(i,j)=SQRT(du(i,j)**2+dv(i,j)**2)
          qo(i,j) = (1._wp / ra) * f1 * q1 * EXP(q2 / ti(i, j))  ! L-Y eqn. 5
          z(i, j) = 10._wp
       ENDDO
    ENDDO

!    CALL ncar_ocean_fluxes(n_itts,uvdel,cort10,ti,corq10,qo,z,cd,ch,ce,ustar,bstar)

    DO j=2,je-1
       DO i=2,ie-1

          ch(i,j) = 1.63e-3_wp                                ! L-Y eqn. 21
          cd(i,j) = 1.63e-3_wp                                ! L-Y eqn. 21
          ce(i,j) = 1.63e-3_wp                                ! L-Y eqn. 21

          ! calculate evaporation
          evap(i,j)=ra*ce(i,j)*(corq10(i,j)-qo(i,j))*uvdel(i,j)   ! L-Y eqn. 20a

          ! calculate latent heat flux
          qla(i,j)=vlamda*evap(i,j)                           ! L-Y eqn. 20b

          ! calculate sensible heat flux
          tv = cort10(i,j)

          qse(i,j)=ra*cp*ch(i,j)*(tv-ti(i,j))*uvdel(i,j) ! L-Y eqn. 20c

          ! calculate net short wave heat flux                ! L-Y eqn. 22
          qnsw(i, j) = corqdsw(i, j) * (  0.29_wp * (1._wp - albvdr) &
                              &         + 0.31_wp * (1._wp - albvdf) &
                              &         + 0.24_wp * (1._wp - albndr) &
                              &         + 0.16_wp * (1._wp - albndf))


          ! calculate net long wave heat flux

          qnlw(i,j)=emiss*corqdlw(i,j)-emiss*stebo*ti(i,j)**4  ! L-Y eqn. 23

          ! calculate net heat flux
          qnet(i,j)=qnsw(i,j)+qnlw(i,j)+qla(i,j)       &
               +qse(i,j)                                          ! L-Y eqn. 19a

          ! calculate net fw flux
          fw(i,j)=corprec(i,j)
!                             +evap(i,j)                         ! L-Y eqn. 19b

          ! calculate wind stress
          taux(i,j)=ra*cd(i,j)*uvdel(i,j)*du(i,j)                ! L-Y eqn. 20d
          tauy(i,j)=ra*cd(i,j)*uvdel(i,j)*dv(i,j)                ! L-Y eqn. 20d

          qcon(i,j)=((tb-ti(i,j))/sictho(i,j)*con)/clb           ! conductive hflx

          qres(i,j)=(-qnet(i,j)/clb)-qcon(i,j)                   ! residual hflx


       ENDDO
    ENDDO

    call  bounds_exch(1,'u',taux,'mo_ncar_ocean_fluxes 17')
    call  bounds_exch(1,'v',tauy,'mo_ncar_ocean_fluxes 18')
    call  bounds_exch(1,'p',fw,'mo_ncar_ocean_fluxes 19')
    call  bounds_exch(1,'p',qnet,'mo_ncar_ocean_fluxes 20')
    call  bounds_exch(1,'p',qnlw,'mo_ncar_ocean_fluxes 21')
    call  bounds_exch(1,'p',qnsw,'mo_ncar_ocean_fluxes 22')
    call  bounds_exch(1,'p',qse,'mo_ncar_ocean_fluxes 23')
    call  bounds_exch(1,'p',qla,'mo_ncar_ocean_fluxes 24')
    call  bounds_exch(1,'p',evap,'mo_ncar_ocean_fluxes 25')
    call  bounds_exch(1,'p',qres,'mo_ncar_ocean_fluxes 26')
    call  bounds_exch(1,'p',qcon,'mo_ncar_ocean_fluxes 27')

  END SUBROUTINE budget_ice_core


  SUBROUTINE ncar_ocean_fluxes(n_itts,udel,t,ts,q,qs,z,cd,ch,ce,ustar,bstar)

    IMPLICIT NONE

    REAL(wp), PARAMETER :: GRAV   = 9.80_wp
    REAL(wp), PARAMETER :: VONKARM = 0.40_wp

    REAL(wp),INTENT(in) :: udel(ie,je), t(ie,je), ts(ie,je),        &
         q(ie,je), qs(ie,je), z(ie,je)
    REAL(wp),INTENT(inout) :: cd(ie,je), ch(ie,je), ce(ie,je),       &
         ustar(ie,je), bstar(ie,je)

    REAL(wp) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
    REAL(wp) :: cd_rt                                ! full drag coefficients @ z
    REAL(wp) :: zeta, x2, x, psi_m, psi_h            ! stability parameters
    REAL(wp) :: u, u10, tv, tstar, qstar, z0, xx, stab

!    INTEGER, PARAMETER :: n_itts = 2
    INTEGER               i, j, jj,n_itts

    DO i=2,ie-1
       DO j=2,je-1
          !          if (avail(i,j)) then
          tv = t(i, j) * (1._wp + 0.608_wp * q(i, j));
          ! 0.5 m/s floor on wind (undocumented NCAR)
          u = MAX(udel(i,j), 0.5_wp);
          u10 = u;                                                ! first guess 10m wind
          ! L-Y eqn. 6a
          cd_n10 = (2.7_wp / u10 + 0.142_wp + 0.0764_wp * u10)/1.e3_wp;
          cd_n10_rt = SQRT(cd_n10);
          ! L-Y eqn. 6b
          ce_n10 = 34.6_wp *cd_n10_rt / 1.e3_wp;
          stab = 0.5_wp + SIGN(0.5_wp, t(i,j) - ts(i,j))
          ch_n10 = (18.0_wp * stab + 32.7_wp * (1._wp - stab)) &
               * cd_n10_rt / 1.e3_wp;       ! L-Y eqn. 6c

          cd(i,j) = cd_n10;                                         ! first guess for exchange coeff's at z
          ch(i,j) = ch_n10;
          ce(i,j) = ce_n10;
          DO jj=1,n_itts                                           ! Monin-Obukhov iteration
             cd_rt = SQRT(cd(i,j));
             ustar(i,j) = cd_rt*u;                                   ! L-Y eqn. 7a
             tstar    = (ch(i,j)/cd_rt)*(t(i,j)-ts(i,j));                ! L-Y eqn. 7b
             qstar    = (ce(i,j)/cd_rt)*(q(i,j)-qs(i,j));                ! L-Y eqn. 7c
             bstar(i,j) = grav*(tstar/tv+qstar/(q(i,j) + 1._wp / 0.608_wp));
             zeta     = vonkarm*bstar(i,j)*z(i,j)/(ustar(i,j)*ustar(i,j)); ! L-Y eqn. 8a
             ! undocumented NCAR
             zeta = SIGN(MIN(ABS(zeta), 10.0_wp), zeta );
             x2 = SQRT(ABS(1._wp - 16._wp * zeta));                            ! L-Y eqn. 8b
             ! undocumented NCAR
             x2 = MAX(x2, 1.0_wp);
             x = SQRT(x2);

             IF (zeta > 0._wp) THEN
               ! L-Y eqn. 8c
               psi_m = -5._wp * zeta;
               ! L-Y eqn. 8c
               psi_h = -5._wp * zeta;
             ELSE
               ! L-Y eqn. 8d
               psi_m = LOG((1._wp + 2._wp * x + x2)*(1._wp + x2)/8._wp) &
                    - 2._wp * (ATAN(x) - api/4._wp);
               ! L-Y eqn. 8e
               psi_h = 2._wp * LOG((1._wp + x2)/2._wp);
             END IF

             ! L-Y eqn. 9
             u10 = u/(1._wp + cd_n10_rt * (LOG(z(i,j)/10._wp) - psi_m)/ vonkarm);
             ! L-Y eqn. 6a again
             cd_n10 = (2.7_wp / u10 + 0.142_wp + 0.0764_wp * u10)/1.e3_wp;
             cd_n10_rt = SQRT(cd_n10);
             ! L-Y eqn. 6b again
             ce_n10 = 34.6_wp * cd_n10_rt / 1.e3_wp;
             stab = 0.5_wp + SIGN(0.5_wp, zeta)
             ch_n10 = (18.0_wp * stab + 32.7_wp * (1._wp - stab)) &
                  * cd_n10_rt / 1.e3_wp;         ! L-Y eqn. 6c again
             z0 = 10._wp * EXP(-vonkarm/cd_n10_rt);             ! diagnostic

             xx = (LOG(z(i,j)/10._wp)-psi_m)/vonkarm;
             cd(i,j) = cd_n10 / (1._wp + cd_n10_rt * xx)**2;       ! L-Y 10a
             xx = (LOG(z(i,j)/10._wp)-psi_h)/vonkarm;
             ch(i,j) = ch_n10 / (1._wp + ch_n10 * xx/cd_n10_rt)**2;                !       b
             ce(i,j) = ce_n10 / (1._wp + ce_n10 * xx/cd_n10_rt)**2;                !       c
          END DO
          !          end if
       END DO
    END DO




  END SUBROUTINE ncar_ocean_fluxes

  SUBROUTINE flx_blk_albedo( tice, sictho, sicsno, palb , palcn , palbp , palcnp )

    !!----------------------------------------------------------------------
    !!               ***  ROUTINE flx_blk_albedo  ***
    !!
    !! ** Purpose :   Computation of the albedo of the snow/ice system
    !!      as well as the ocean one
    !!
    !! ** Method  : - Computation of the albedo of snow or ice (choose the
    !!      right one by a large number of tests
    !!              - Computation of the albedo of the ocean
    !!
    !! References :
    !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
    !!
    !! History :
    !!  8.0   !  01-04  (LIM 1.0)
    !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
    !!----------------------------------------------------------------------
    !! * Modules used

    !! * Arguments
    REAL(wp), DIMENSION(ie,je), INTENT(out) ::  &
           palb         , & !  albedo of ice under overcast sky
           palcn        , & !  albedo of ocean under overcast sky
           palbp        , & !  albedo of ice under clear sky
           palcnp                  !  albedo of ocean under clear sky

    !! * Local variables
    INTEGER ::   i,j        ! dummy loop indices
    ! FIXME: these should probably be parameters,
    ! currently they are in constant storage
    REAL(wp) ::   &
         c1     = 0.05_wp, & ! constants values
         c2     = 0.1_wp, &
         albice = 0.50_wp, & !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
         cgren  = 0.06_wp, & !  correction of the snow or ice albedo to take into account
         !  effects of cloudiness (Grenfell & Perovich, 1984)
         alphd  = 0.80_wp, & ! coefficients for linear interpolation used to compute
         alphdi = 0.72_wp, & !  albedo between two extremes values (Pyane, 1972)
         alphc  = 0.65_wp, &
         zmue   = 0.4_wp, &     !  cosine of local solar altitude
         zzero   = 0.0_wp,  &
         zone    = 1.0_wp,  &
         rt0_snow = 0.0_wp,  &
         rt0_ice = 0.0_wp


    REAL(wp) ::   &
         zmue14         , & !  zmue**1.4
         zalbpsnm       , & !  albedo of ice under clear sky when snow is melting
         zalbpsnf       , & !  albedo of ice under clear sky when snow is freezing
         zalbpsn        , & !  albedo of snow/ice system when ice is coverd by snow
         zalbpic        , & !  albedo of snow/ice system when ice is free of snow
         zithsn         , & !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
         zitmlsn        , & !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
         zihsc1         , & !  = 1 hsn <= c1 ; = 0 hsn > c1
         zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2

    REAL(wp), DIMENSION(ie,je) ::  &
         zalbfz         , & !  ( = alphdi for freezing ice ; = albice for melting ice )
         zficeth                  !  function of ice thickness
    LOGICAL , DIMENSION(ie,je) ::  &
         llmask

    REAL(wp), DIMENSION(ie,je) ::      &
          tice   ,                 &   !: Sea-Ice Surface Temperature (Kelvin )
          sicsno  ,                &   !: Snow thickness
          sictho                       !: Ice thickness


    !-------------------------
    !  Computation of  zficeth
    !--------------------------

    llmask = (sicsno == 0.0_wp) .AND. ( tice >= rt0_ice )
    WHERE ( llmask )   !  ice free of snow and melts
       zalbfz = albice
    ELSEWHERE
       zalbfz = alphdi
    END WHERE

    DO i = 1, ie
       DO j = 1, je
          IF (sictho(i,j) > 1.5_wp) THEN
             zficeth(i,j) = zalbfz(i,j)
          ELSEIF( sictho(i,j) > 1.0_wp .AND. sictho(i,j) <= 1.5_wp) THEN
             zficeth(i, j) = 0.472_wp + 2.0_wp * (zalbfz(i, j) - 0.472_wp) &
                  * (sictho(i, j) - 1.0_wp)
          ELSEIF (sictho(i, j) > 0.05_wp .AND. sictho(i,j) <= 1.0_wp) THEN
             zficeth(i, j) = 0.2467_wp + 0.7049_wp * sictho(i, j) &
                          - 0.8608_wp * sictho(i, j) * sictho(i, j) &
                          + 0.3812_wp * sictho(i, j) * sictho(i, j) * sictho (i, j)
          ELSE
             zficeth(i,j) = 0.1_wp + 3.6_wp * sictho(i, j)
          ENDIF
       END DO
    END DO

    !-----------------------------------------------
    !    Computation of the snow/ice albedo system
    !-------------------------- ---------------------

    !    Albedo of snow-ice for clear sky.
    !-----------------------------------------------
    DO i = 1, ie
       DO j = 1, je
          !  Case of ice covered by snow.
          !  melting snow
         ! FIXME: replace zzero with zero from mo_commo1
          zihsc1 = 1.0_wp - MAX(zzero, SIGN(zone, -(sicsno(i, j) - c1)))
          zalbpsnm = (1.0_wp - zihsc1) &
               * ( zficeth(i,j) + sicsno(i,j) * ( alphd - zficeth(i,j) ) / c1 ) &
               &                 + zihsc1   * alphd
          !  freezing snow
          zihsc2       = MAX ( zzero , SIGN ( zone , sicsno(i,j) - c2 ) )
          zalbpsnf = (1.0_wp - zihsc2 ) * &
               ( albice + sicsno(i,j) * ( alphc - albice ) / c2 ) &
               &                 + zihsc2   * alphc

          zitmlsn      =  MAX ( zzero , SIGN ( zone , tice(i,j) - rt0_snow ) )
          zalbpsn = zitmlsn * zalbpsnf + (1.0_wp - zitmlsn) * zalbpsnm

          !  Case of ice free of snow.
          zalbpic      = zficeth(i,j)

          ! albedo of the system
          zithsn = 1.0_wp - MAX(zzero, SIGN(zone, -sicsno(i, j)))
          palbp(i,j) = zithsn * zalbpsn + (1.0_wp - zithsn) * zalbpic
       END DO
    END DO

    !    Albedo of snow-ice for overcast sky.
    !----------------------------------------------
    palb(:,:)   = palbp(:,:) + cgren

    !--------------------------------------------
    !    Computation of the albedo of the ocean
    !-------------------------- -----------------


    !  Parameterization of Briegled and Ramanathan, 1982
    zmue14      = zmue**1.4_wp
    palcnp(:,:) = 0.05_wp / (1.1_wp * zmue14 + 0.15_wp )

    !  Parameterization of Kondratyev, 1969 and Payne, 1972
    palcn(:,:)  = 0.06_wp

  END SUBROUTINE flx_blk_albedo

  SUBROUTINE normpem

! close fw budget for each timestep as suggested in the core release notes

  REAL(wp) :: rimbl,rarea, re, zold, znew, delz, draft
  INTEGER :: i,j

  rimbl = 0._wp
  rarea = 0._wp
  DO j=2,je-1
     DO i=2,ie-1
        rimbl=rimbl+(prech(i,j)+eminpo(i,j))          &
             *dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
        rarea=rarea+dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
     ENDDO
  ENDDO

  CALL global_sum(rimbl)
  CALL global_sum(rarea)

  re=rimbl/rarea

  WRITE(io_stdout, *) 'core:  before normalisation P-E+R+r [Sv]: ', &
       re * rarea * 1.e-6_wp

    DO j=2,je-1
     DO i=2,ie-1
        IF (weto(i,j,1) .ge. 0.5_wp) THEN
            draft=ddpo(i,j,1)+zo(i,j)-((sictho(i,j)*rhoicwa) &
                 -(sicsno(i,j)*rhosnwa))
            zold=zo(i,j)
            znew=zold-(re*dt)
            zo(i,j)=znew
            eminpo(i,j)=eminpo(i,j)-re
            delz=znew-zold
            sao(i,j,1)=sao(i,j,1)*(draft/(delz+draft))
        ENDIF
     ENDDO
    ENDDO

  rimbl = 0._wp
  rarea = 0._wp
  DO j=2,je-1
     DO i=2,ie-1
        rimbl=rimbl+(prech(i,j)+eminpo(i,j))          &
             *dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
        rarea=rarea+dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
     ENDDO
  ENDDO

  CALL global_sum(rimbl)
  CALL global_sum(rarea)

  re=rimbl/rarea

  WRITE(io_stdout, *) 'core: after normalisation P-E+R+r [Sv]: ', &
       re * rarea * 1.e-6_wp

  CALL bounds_exch(1,'p',zo,'mo_ncar_ocean_fluxes 28')
  CALL bounds_exch(1,'p',sao,'mo_ncar_ocean_fluxes 29')
  CALL bounds_exch(1,'p',eminpo,'mo_ncar_ocean_fluxes 39')

  END SUBROUTINE normpem

END MODULE mo_ncar_ocean_fluxes















