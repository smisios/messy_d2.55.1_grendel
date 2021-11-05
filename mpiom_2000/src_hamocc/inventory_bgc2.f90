MODULE mo_inventory

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1_bgc

      USE MO_COMMO1, ONLY: area,weto,lbounds_exch_tp
      USE mo_grid, ONLY: thkcello
      USE mo_planetary_constants, ONLY : rhoicwa,rhosnwa

      USE mo_parallel

      IMPLICIT NONE

      PRIVATE

      REAL(wp) :: zpowtrato(npowtra)
      REAL(wp) :: zsedlayto(nsedtra),zburial(nsedtra)
      REAL(wp) :: zocetrato(nocetra)
      REAL(wp) :: zprorca,zprcaca,zsilpro

      REAL(wp),TARGET :: totalcarbon     !< total carbon [kmolC] code 301
      REAL(wp),TARGET :: zco2flux        !< global carbon flux [kmolC/day] code 303
      REAL(wp),TARGET :: totalcarboce    !< total carbon of water column (incl. pore water) [kmolC] code 305
      REAL(wp),TARGET :: totalcarbsed    !< total carbon of sediment code [kmolC] 306
      REAL(wp),TARGET :: calcitesediment !< total calcite of sediment [kmolC] code 307
      REAL(wp),TARGET :: organsed        !< organic carbon of sediment [kmolC] code 308
      REAL(wp),TARGET :: totalphos       !< total phosphate  [kmolP] code 309
      REAL(wp),TARGET :: watphos         !< total phosphate of water column (incl. pore water) [kmolP] code 310
      REAL(wp),TARGET :: sedphos         !< total phospate of sediment code 311
      REAL(wp),TARGET :: totalsil        !< total silicate  [kmolSi(OH)4] code 312
      REAL(wp),TARGET :: watsil          !< total silicate of water column (incl. pore water) [kmolSi(OH)4] code 313
      REAL(wp),TARGET :: sedsil          !< total silicate of sediment [kmolSi(OH)4] code 314

      PUBLIC :: totalcarbon, zco2flux, totalcarboce, totalcarbsed, &
           calcitesediment, organsed, totalphos, watphos, sedphos, &
           totalsil, watsil, sedsil, inventory_bgc2

    CONTAINS

      SUBROUTINE inventory_bgc2(kpie,kpje,kpke,kplyear,kplmon,kplday)

      integer :: kpie,kpje,kpke,i,j,k,l,nnn
      integer :: kplyear,kplmon,kplday

      INTEGER :: ib1,ib2,ib3,ib4,jb
      REAL(wp) :: r1


        jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

!----------------------------------------------------------------------

        DO k=1,ks
          ! aqueous sediment tracer
          DO l=1,npowtra
            zpowtrato(l) = 0.0_wp
            zpowtrato(l)=zpowtrato(l)+SUM(WETO(2:kpie-1,jb:kpje-1,1)  &
                 *powtra(2:kpie-1,jb:kpje-1,k,l)                      &
                 *seddw(k)*area(2:kpie-1,jb:kpje-1)*porwat(k))
          ENDDO

          ! non aqueous sediment tracer
          DO l=1,nsedtra
            zsedlayto(l) = 0.0_wp
            zsedlayto(l)=zsedlayto(l)+SUM(WETO(2:kpie-1,jb:kpje-1,1)    &
                 *sedlay(2:kpie-1,jb:kpje-1,k,l)                        &
                 *porsol(k)*seddw(k)*area(2:kpie-1,jb:kpje-1))
          ENDDO

          DO l=1,nsedtra
            zburial(l) = 0.0_wp
            zburial(l)=zburial(l)+SUM(burial(2:kpie-1,jb:kpje-1,l)  &
                 *WETO(2:kpie-1,jb:kpje-1,1)*area(2:kpie-1,jb:kpje-1))
          ENDDO

        ENDDO

        CALL global_sum(zpowtrato)
        CALL global_sum(zsedlayto)
        CALL global_sum(zburial)


        !  oceanic tracers
!----------------------------------------------------------------------

        zocetrato(:) = 0.0_wp

        DO l=1,nocetra

          zocetrato(l)=SUM(WETO(2:kpie-1,jb:kpje-1,1)*                &
               ocetra(2:kpie-1,jb:kpje-1,1,l)*                        &
               area(2:kpie-1,jb:kpje-1)*thkcello(2:kpie-1,jb:kpje-1,1))

          DO k=2,kpke
            zocetrato(l)=zocetrato(l)+SUM(WETO(2:kpie-1,jb:kpje-1,K)* &
                 ocetra(2:kpie-1,jb:kpje-1,k,l)*                    &
                 area(2:kpie-1,jb:kpje-1)*thkcello(2:kpie-1,jb:kpje-1,k))
          ENDDO

        ENDDO

        CALL global_sum(zocetrato)

        zprorca=SUM(prorca(2:kpie-1,jb:kpje-1)                        &
             *weto(2:kpie-1,jb:kpje-1,1)*area(2:kpie-1,jb:kpje-1))

        zprcaca=SUM(prcaca(2:kpie-1,jb:kpje-1)                        &
             *weto(2:kpie-1,jb:kpje-1,1)*area(2:kpie-1,jb:kpje-1))

        zsilpro=SUM(silpro(2:kpie-1,jb:kpje-1)                        &
             *weto(2:kpie-1,jb:kpje-1,1)*area(2:kpie-1,jb:kpje-1))


! ks   : co2flux_cpl is co2-flux at actual time step
!      if the coupling interval is one day, co2flux_cpl can be used
!      given a higher coupling interval the daily integral of
!      co2flux_cpl has to be introduced
!      co2flux_cpl has different units, convert to kmol/day

#ifdef __cpl_co2
        zco2flux=SUM(co2flux_cpl(2:kpie-1,jb:kpje)                    &
             *WETO(2:kpie-1,jb:kpje,1)*area(2:kpie-1,jb:kpje))

        ! FIXME: replace 44.011 with constant
        zco2flux = zco2flux * 86400._wp / 44.011_wp

        CALL global_sum(zco2flux)  ! code 303
#endif
        CALL global_sum(zprorca,zprcaca,zsilpro)

! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2

! Sum of inventory
!----------------------------------------------------------------------
! Units in P have a C:P Ratio of 122:1 (rcar=122.)
! #slo: prorca and prcaca are fluxes into sediment
!       set to zero after sediment is filled
! ks: at this stage of the model, prorca, prcaca and silpro must be zero, could be removed

      totalcarboce=                          &   ! total carbon of water column (incl. pore water) code 305
           (zocetrato(idet)                  &   !
         +  zocetrato(idoc)                  &   !  dissolved organic carbon
         +  zocetrato(iphy)                  &   !  phytoplancton
         +  zocetrato(izoo))*rcar            &   !  zooplancton
         + zocetrato(isco212)                &   !  disolved inorganic carbon
         + zocetrato(icalc)                  &
         + zpowtrato(ipowaic)

      calcitesediment=                       &   !  total calcite of sediment code 307
           zsedlayto(isssc12)                &
         + zburial(isssc12)                  &   !  burial reservoir layer
          + zprcaca - calcinpglint

      organsed=                              &   !  organic carbon of sediment code 308
           (zsedlayto(issso12)               &
         +  zburial(issso12)                 &
         +  zprorca - orginpglint ) * rcar

      totalcarbsed=calcitesediment+organsed   ! code 306

      totalcarbon = totalcarboce +  totalcarbsed ! code 301

      watphos=                               &  !code 310
           zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)              &
           + zocetrato(izoo)+zocetrato(iphosph)                           &
           + zpowtrato(ipowaph)

      sedphos=zsedlayto(issso12)+zburial(issso12) & !code 311
           +zprorca-orginpglint

      totalphos=watphos+sedphos  ! code 309

      watsil=                                & !code 313
           zocetrato(isilica)                &
         + zocetrato(iopal)                  &
         + zsilpro                           &
         + zpowtrato(ipowasi)

      sedsil=                                & !code 314
         + zsedlayto(issssil)                &
         + zburial(issssil) - silinpglint

      totalsil=watsil+sedsil ! code 312

 
    END SUBROUTINE INVENTORY_BGC2

  END MODULE mo_inventory
