MODULE mo_bgc_varlist

#ifndef NO_NEW_IO

  USE mo_kind,     ONLY : dp, sp, wp
  USE mo_varlist, ONLY : varlist, new_var, generate_gridid, generate_zaxisid, print_varlist, varlist_add

  IMPLICIT NONE

#include "pointer_intent_macro.inc"
  INCLUDE 'cdi.inc'

  TYPE (varlist), POINTER :: bgc_varlist

CONTAINS


  SUBROUTINE generate_sed_zaxisid(this)

#ifndef NOCDI
    USE mo_parallel, ONLY : p_io, p_pe
#endif

    USE mo_sedmnt, ONLY: z_sed
    USE mo_param1_bgc, ONLY: ks

    TYPE (varlist), POINTERINTENT(in) :: this
    TYPE (varlist), POINTER :: current

#ifndef NOCDI
    INTEGER :: sed_zaxisid

    REAL(wp) :: sed_zaxis(ks)

!   sediment
    sed_zaxis(:) = 10.0_wp * z_sed(:)
    sed_zaxisid = zaxiscreate(ZAXIS_DEPTH_BELOW_LAND , ks)

    CALL zaxisdeflevels(sed_zaxisid, sed_zaxis)

    IF ( p_pe == p_io ) THEN

      current => this

      DO WHILE ( ASSOCIATED(current) )

        IF (current%vardata%zaxis == 'sed') THEN

          current%vardata%zaxisid  = sed_zaxisid

        ENDIF

        current => current%next ! make current alias of next varnode

      END DO

    ENDIF

#endif

  END SUBROUTINE generate_sed_zaxisid

  SUBROUTINE build_bgc_varlist

    USE mo_parallel, ONLY : p_pe, p_io

    USE mo_param1_bgc,ONLY : isco212,ialkali,iphosph,iano3,isilica,iiron      &
                            ,ioxygen,idet,icalc,iphy,igasnit,iopal,izoo,idoc  &
                            ,ian2o,idms,ifdust,iiron,issso12,isssc12,issssil  &
                            ,issster,ipowaic,ipowaal,ipowaph,ipowaox,ipown2   &
                            ,isco213,isco214,idet13,idet14,icalc13,icalc14    &
                            ,issso13,issso14,isssc13,isssc14,ipowc13,ipowc14  &
                            ,inos,iadust                                      &
                            ,icfc11,icfc12                                    &
                            ,ipowno3,ipowasi,iatmco2,iatmo2,iatmn2
    USE mo_biomod, ONLY: rcar,riron

    USE mo_carbch,ONLY : ocetra,hi,co3,aksp,atm    &
           ,silinpglint,orginpglint,calcinpglint
    USE mo_sedmnt,ONLY : sedlay,burial,sedhpl,powtra
    USE mo_bgc_diagnostic

    USE mo_inventory, ONLY : totalcarbon, zco2flux, totalcarboce,        &
         totalcarbsed, calcitesediment, organsed, totalphos, watphos,    &
         sedphos, totalsil, watsil, sedsil

 IMPLICIT NONE
    REAL(wp), PARAMETER :: kilo = 1.e3_wp

    ! build up the list

    NULLIFY(bgc_varlist) ! initially nullify list (empty)

    CALL varlist_add(bgc_varlist, &
         new_var('sco212', 'dissolved_inorganic_carbon', &
         'kmol C m-3', 7, ocetra, isco212, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('alkali', 'total_alkalinity', &
         'kmol m-3', 10, ocetra, ialkali, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('phosph', 'dissolved_phosphate', &
         'kmol P m-3', 11, ocetra, iphosph, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('ano3', 'dissolved_nitrate', &
         'kmol N m-3', 14, ocetra, iano3, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('silica', 'dissolved_silicate', &
         'kmol Si m-3', 15, ocetra, isilica, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('iron', 'iron_concentration', &
         'kmol Fe m-3', 31, ocetra, iiron, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('oxygen', 'oxygen_concentration', &
         'kmol O2 m-3', 12, ocetra, ioxygen, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('gasnit', 'gaseous_nitrogen', &
         'kmol N2 m-3', 13, ocetra, igasnit, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('doc', 'dissolved_organic_carbon', &
         'kmol P m-3', 16, ocetra, idoc, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('det', 'detrital_organic_carbon', &
         'kmol P m-3', 17, ocetra, idet, 'p', 'c'))

    CALL varlist_add(bgc_varlist, &
         new_var('hi', 'hydrogen ion_concentration', &
         'kmol m-3', 20, hi, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('co3', 'dissolved carbonate', &
         'kmol C m-3', 21, co3, 'p', 'c'))

    CALL varlist_add(bgc_varlist, &
         new_var('phy', 'phytoplankton_concentration', &
         'kmol P m-3', 22, ocetra, iphy, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('zoo', 'zooplankton_concentration', &
         'kmol P m-3', 23, ocetra, izoo, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('calc', 'calcium carbonate', &
         'kmol C m-3', 24, ocetra, icalc, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('opal', 'biogenic silica', &
         'kmol Si m-3', 27, ocetra, iopal, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('an2o', 'laughing gas', &
         'kmol N2 m-3', 28, ocetra, ian2o, 'p', 'c'))

    CALL varlist_add(bgc_varlist, &
         new_var('dms', 'DiMethylSulfide', &
         'kmol S m-3', 29, ocetra, idms, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('fdust', 'Non-aggregated dust', &
         'kg m-3', 30, ocetra, ifdust, 'p', 'c'))

    CALL varlist_add(bgc_varlist, &
         new_var('aksp', 'apparent solubility product for calcite', &
         '(kmol m-3)2', 37, aksp, 'p', 'c'))

    CALL varlist_add(bgc_varlist, &
         new_var('ssso12', 'Sediment accumulated organic carbon', &
         'kmol P m-3', 38, sedlay, issso12, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('sssc12', 'Sediment accumulated calcium carbonate', &
         'kmol C m-3', 41, sedlay, isssc12, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('ssssil', 'Sediment accumulated opal', &
         'kmol Si m-3', 44, sedlay, issssil, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('ssster', 'Sediment accumulated clay', &
         'kmol m-3', 45, sedlay, issster, 'p', 'sed'))

    CALL varlist_add(bgc_varlist, &
         new_var('bsso12', 'Burial layer of organic carbon', &
         'kmol P m-3', 46, burial, issso12, 'p', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('bssc12', 'Burial layer of calcium carbonate', &
         'kmol C m-3', 47, burial, isssc12, 'p', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('bsssil', 'Burial layer of opal', &
         'kmol Si m-3', 48, burial, issssil, 'p', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('bsster', 'Burial layer of clay', &
         'kmol m-3', 49, burial, issster, 'p', 'g'))

    CALL varlist_add(bgc_varlist, &
         new_var('sedhpl', 'Sediment accumulated hydrogen ions', &
         'kmol m-3', 50, sedhpl, 'p', 'sed'))

    CALL varlist_add(bgc_varlist, &
         new_var('powaic', 'Sediment pore water DIC', &
         'kmol C m-3', 51, powtra, ipowaic, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powaal', 'Sediment pore water alkalinity', &
         'kmol m-3', 54, powtra, ipowaal, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powaph', 'Sediment pore water phosphate', &
         'kmol P m-3', 55, powtra, ipowaph, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powaox', 'Sediment pore water oxygen', &
         'kmol O2 m-3', 56, powtra, ipowaox, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('pown2', 'Sediment pore water gaseous nitrogen', &
         'kmol N2 m-3', 57, powtra, ipown2, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powno3', 'Sediment pore water nitrate (NO3)', &
         'kmol N m-3', 58, powtra, ipowno3, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powasi', 'Sediment pore water silicid acid (Si(OH)4)', &
         'kmol Si m-3', 59, powtra, ipowasi, 'p', 'sed'))

    CALL varlist_add(bgc_varlist, &
         new_var('atmco2', 'atmospheric CO2', &
         'ppm', 61, atm, iatmco2, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('atmo2', 'atmospheric O2', &
         'ppm', 62, atm, iatmo2, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('atmn2', 'atmospheric N2', &
         'ppm', 64, atm, iatmn2, 'p', 's'))

!#ifdef __c_isotopes
    CALL varlist_add(bgc_varlist, &
         new_var('sco213', 'dissolved_inorganic_carbon13', &
         'kmol m-3', 8, ocetra, isco213, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('sco214', 'dissolved_inorganic_carbon14', &
         'kmol m-3', 9, ocetra, isco214, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('det13', 'particulate organic carbon13', &
         'kmolP m-3', 18, ocetra, idet13, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('det14', 'particulate organic carbon14', &
         'kmolP m-3', 19, ocetra, idet14, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('calc13', 'calcium carbonate13', &
         'kmolP m-3', 25, ocetra, icalc13, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('calc14', 'calcium carbonate14', &
         'kmolP m-3', 26, ocetra, icalc14, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('ssso13', 'Sediment accumulated organic carbon13', &
         'kmol m-2', 39, sedlay, issso13, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('ssso14', 'Sediment accumulated organic carbon14', &
         'kmol m-2', 40, sedlay, issso14, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('sssc13', 'Sediment accumulated calcium carbonate13', &
         'kmol m-2', 42, sedlay, isssc13, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('sssc14', 'Sediment accumulated calcium carbonate14' &
         , 'kmol m-2', 43, sedlay, isssc14, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powc13', 'Sediment pore water DIC13', &
         'kmol m-3', 52, powtra, ipowc13, 'p', 'sed'))
    CALL varlist_add(bgc_varlist, &
         new_var('powc14', 'Sediment pore water DIC14', &
         'kmol m-3', 53, powtra, ipowc14, 'p', 'sed'))
!#endif
!#ifdef AGG
    CALL varlist_add(bgc_varlist, &
         new_var('nos', 'marine snow aggregates per cm3', &
         'cm-3', 32, ocetra, inos, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('adust', 'Aggregated dust', &
         'kg m-3', 33, ocetra, iadust, 'p', 'c'))
!#endif
!#ifdef PCFC
    CALL varlist_add(bgc_varlist, &
         new_var('cfc11', 'CFC11', &
         'kmol m-3', 35, ocetra, icfc11, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('cfc12', 'CFC12', &
         'kmol m-3', 36, ocetra, icfc12, 'p', 'c'))
!#endif
#ifdef AMMO
    CALL varlist_add(bgc_varlist, &
         new_var('ammo', 'dissolved NH4', &
         'kmol N m-3', 150, ocetra, iammo, 'p', 'c'))
#endif

! variables from bgcmean_bioz, vertical resolvled
    CALL varlist_add(bgc_varlist, &
         new_var('phosy', 'photosynthesis', &
         'mol C m-3 s-1', 100, bgcprod, kphosy, 'p', 'c', &
         factor = rcar * kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('graz', 'zoo_grazing', &
         'mol Fe m-3 s-1', 101, bgcprod, kgraz, 'p', 'c', &
         factor = riron * kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('export', 'detritus_production', &
         'kmol C m-3 s-1', 84, bgcprod, kexport, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('delsil', 'opal_production', &
         'kmol Si m-3 s-1', 86, bgcprod, kdelsil, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('delcar', 'calcium_production', &
         'kmol C m-3 s-1', 85, bgcprod, kdelcar, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('dmsprod', 'DMS_Production', &
         'kmol S m-3 s-1', 69, bgcprod, kdmsprod, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('dms_bac', 'DMS_bacterial_consumption', &
         'kmol S m-3 s-1', 70, bgcprod, kdms_bac, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('dms_uv', 'DMS_uv_light_destruction', &
         'kmol S m-3 s-1', 71, bgcprod, kdms_uv, 'p', 'c'))

! variables from bgcmean_2d
    CALL varlist_add(bgc_varlist, &
         new_var('co2flux', 'co2flux', &
         'kg C m-2 s-1', 92, bgcflux, kco2flux, 'p', 's', factor=-12._wp))
    CALL varlist_add(bgc_varlist, &
         new_var('co214f', 'co214flux', &
         'kmol m-2 s-1', 91, bgcflux, kco214f, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('o2flux', 'oxygen_flux', &
         'mol O2 m-2 s-1', 72, bgcflux, ko2flux, 'p', 's', factor=-kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('n2flux', 'nitrogen_flux', &
         'kmol N2 m-2 s-1', 74, bgcflux, kn2flux, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('n2oflux', 'n2oflux', &
         'kmol N2 m-2 s-1', 93, bgcflux, kn2oflux, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('dmsflux', 'DMS_Flux', &
         'mol S m-2 s-1', 68, bgcflux, kdmsflux, 'p', 's', factor=kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('prorca', 'detritus sediment flux', &
         'mol P m-2 s-1', 94, bgcflux, kprorca, 'p', 'g', factor=kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('prcaca', 'CaCO3 sediment flux', &
         'mol C m-2 s-1', 95, bgcflux, kprcaca, 'p', 'g', factor=kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('silpro', 'opal sediment flux', &
         'kmol Si m-2 s-1', 96, bgcflux, ksilpro, 'p', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('produs', 'dust flux to sediment (free+aggregated)', &
         'kmol m-2 s-1', 97, bgcflux, kprodus, 'p', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('kwco2', 'co2_exchange_coefficient', &
         'kmol m-2 s-1 ppm-1', 73, bgcflux, kkwco2, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('pco2', 'CO2_partial_pressure', &
         'ppm', 67, bgcflux, kpco2, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('co2fxd', 'CO2_Flux_down', &
         'kmol C m-2 s-1', 65, bgcflux, kco2fxd, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('co2fxu', 'CO2_Flux_Up', &
         'kmol C m-2 s-1', 66, bgcflux, kco2fxu, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('opex90', 'opal flux in 90 m', &
         'mol Si m-2 s-1', 75, bgcflux, kopex90, 'p', 's=90', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('opex1000', 'opal flux in 1000 m', &
         'mol Si m-2 s-1', 76, bgcflux, kopex1000, 'p', 's=1000', &
         factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('opex2000', 'opal flux in 2000 m', &
         'mol Si m-2 s-1', 77, bgcflux, kopex2000, 'p', 's=2000', &
         factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('caex90', 'calc flux in 90 m', &
         'mol C m-2 s-1', 78, bgcflux, kcaex90, 'p', 's=90', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('caex1000', 'calc flux in 1000 m', &
         'mol C m-2 s-1', 79, bgcflux, kcaex1000, 'p', 's=1000', &
         factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('caex2000', 'calc flux in 2000 m', &
         'mol C m-2 s-1', 80, bgcflux, kcaex2000, 'p', 's=2000', &
         factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('coex90', 'OM flux in 90 m', &
         'mol C m-2 s-1', 81, bgcflux, kcoex90, 'p', 's=90', &
         factor = rcar * kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('coex1000', 'OM flux in 1000 m', &
         'mol C m-2 s-1', 82, bgcflux, kcoex1000, 'p', 's=1000', &
         factor = rcar * kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('coex2000', 'OM flux in 2000 m', &
         'mol C m-2 s-1', 83, bgcflux, kcoex2000, 'p', 's=2000', &
         factor = rcar * kilo))

! all optional variables from bgcmean files

!#ifdef ANTC14
    CALL varlist_add(bgc_varlist, &
         new_var('ac14fx', 'anthro. C14 air-sea-flux', &
         'kmol m-2 s-1', 106, bgcflux, kac14fx, 'p', 's'))
!#ifdef __c_isotopes
    CALL varlist_add(bgc_varlist, &
         new_var('c13flux', 'C13 air-sea-flux', &
         'kmol m-2 s-1', 103, bgcflux, kc13flux, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('c14flux', 'C14 air-sea-flux', &
         'kmol m-2 s-1', 104, bgcflux, kc14flux, 'p', 's'))
!#ifdef PCFC
    CALL varlist_add(bgc_varlist, &
         new_var('cfc11fx', 'CFC11 sea_air flux', &
         'kmol m-2 s-1', 87, bgcflux, kcfc11fx, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('cfc12fx', 'CFC12 sea_air flux', &
         'kmol m-2 s-1', 88, bgcflux, kcfc12fx, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('pcfc11', 'CFC11 atmos. partial pressure', &
         'ppm', 89, bgcflux, kpcfc11, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('pcfc12', 'CFC12 atmos. partial pressure', &
         'ppm', 90, bgcflux, kpcfc12, 'p', 's'))
!#ifdef AMMO
    CALL varlist_add(bgc_varlist, &
         new_var('nh3flux', 'NH3 sea-air-flux', &
         'kmol N m-2 s-1', 105, bgcflux, knh3flux, 'p', 's'))
!#endif
    CALL varlist_add(bgc_varlist, &
         new_var('n2fix', 'nitrogen_fixation_rate', &
         'mol N m-2 s-1', 157, bgcflux, kn2fix, 'p', 's', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('lysokl', 'depth_of_lysokline', &
         'm', 160, bgcflux, klysokl, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('featm', 'surface_flux_of_iron', &
         'mol Fe m-2 s-1', 161, bgcflux, kfeatm, 'p', 's', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('o2min', 'oxygen_minimum_concentration', &
         'mol O2 m-3', 162, bgcomz, komz, 'p', 's', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('zo2min', 'depth_of_oxygen_minimum_concentration', &
         'm', 163, bgcomz, komz_depth, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('dpco2', 'delta_CO2_partial_pressure', &
         'ppm',164, bgcflux, kdpco2, 'p', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('dpo2', 'delta_O2_partial_pressure', &
         'ppm', 165, bgcflux, kdpo2, 'p', 's'))

    CALL varlist_add(bgc_varlist, &
         new_var('caco3_diss', 'dissolution of caco3', &
         'kmol C m-3 s-1', 158, bgc_o_pro, kdissol, 'p', 'c'))
    CALL varlist_add(bgc_varlist, &
         new_var('denitr', 'denitrification', &
         'mol N m-2 s-1', 159, bgc_o_pro, kdenit, 'p', 'si', factor = kilo))

    ! CMIP5 surface tracer concentrations [units mol/m*3]
    CALL varlist_add(bgc_varlist, &
         new_var('surf_sco212', 'surface_dissolved_inorganic_carbon', &
         'mol C m-3', 107, ocetra, isco212, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_alkali', 'surface_total_alkalinity', &
         'mol m-3', 110, ocetra, ialkali, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_phosph', 'surface_dissolved_phosphate', &
         'mol P m-3', 111, ocetra, iphosph, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_ano3', 'surface_dissolved_nitrate', &
         'mol N m-3', 114, ocetra, iano3, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_silica', 'surface_dissolved_silicate', &
         'mol Si m-3', 115, ocetra, isilica, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_iron', 'surface_iron_concentration', &
         'mol Fe m-3', 131, ocetra, iiron, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_oxygen', 'surface_oxygen_concentration', &
         'mol O2 m-3', 112, ocetra, ioxygen, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_doc', 'surface_dissolved_organic_carbon', &
         'mol C m-3', 116, ocetra, idoc, 'p', 's=0', factor = rcar*kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_det', 'surface_detrital_organic_carbon', &
         'mol C m-3', 117, ocetra, idet, 'p', 's=0', factor = rcar*kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_hi', 'surface_hydrogen_ion_concentration', &
         'kmol m-3', 120, hi, 'p', 's=0'))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_phy', 'surface_phytoplankton_concentration', &
         'mol C m-3', 122, ocetra, iphy, 'p', 's=0', factor = rcar*kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_zoo', 'surface_zooplankton_concentration', &
         'mol C m-3', 123, ocetra, izoo, 'p', 's=0', factor = rcar*kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_calc', 'surface_calcium_carbonate', &
         'mol C m-3', 124, ocetra, icalc, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_opal', 'surface_biogenic_silica', &
         'mol Si m-3', 127, ocetra, iopal, 'p', 's=0', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('surf_dms', 'surface_DiMethylSulfide', &
         'mol S m-3', 129, ocetra, idms, 'p', 's=0', factor = kilo))

    !CMIP5 integrated  variables
    CALL varlist_add(bgc_varlist, &
         new_var('intpp', 'integrated_primary_production', &
         'mol C m-2 s-1', 200, bgcprod, kphosy, 'p', 'si', factor = rcar*kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('intdelsil', 'integrated_opal_production', &
         'mol Si m-2 s-1', 186, bgcprod, kdelsil, 'p', 'si', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('intdelcar', 'integrated_calcium_production', &
         'mol C m-2 s-1', 185, bgcprod, kdelcar, 'p', 'si', factor = kilo))


    ! carbon cycle diagnostic
    CALL varlist_add(bgc_varlist, &
         new_var('calcinpglint', 'total_inorganic_calcium_input_to_sediment', &
         'kmol', 203, calcinpglint, 'g', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('orginpglint', 'total_organic_calcium_input_to_sediment', &
         'kmol', 204, orginpglint, 'g', 'g'))
    CALL varlist_add(bgc_varlist, &
         new_var('silinpglint', 'total_silicat_input_to_sediment', &
         'kmol', 205, silinpglint, 'g', 'g'))

    CALL varlist_add(bgc_varlist, &
         new_var('intsco212', 'integrated dissolved_inorganic_carbon', &
         'kg C m-2', 206, ocetra, isco212, 'p', 'si', factor = 12._wp))

    !CMIP5 integrated variables top 100 m
    CALL varlist_add(bgc_varlist, &
         new_var('fddtdic', 'change in dissolved_inorganic_carbon', &
         'mol C m-2 s-1', 207, bgcrate, ksco212, 'p', 'si<100', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('fddtalk', 'change in total_alkalinity', &
         'mol m-2 s-1', 210, bgcrate, kalkali, 'p', 'si<100', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('fddtdip', 'change in dissolved_phosphate', &
         'mol P m-2 s-1', 211, bgcrate, kphosph, 'p', 'si<100', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('fddtdin', 'change in dissolved_nitrate', &
         'mol N m-3', 214, bgcrate, kano3, 'p', 'si<100', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('fddtdisi', 'change in dissolved_silicate', &
         'mol Si m-2 s-1', 215, bgcrate, ksilica, 'p', 'si<100', factor = kilo))
    CALL varlist_add(bgc_varlist, &
         new_var('fddtdife', 'change in iron_concentration', &
         'mol Fe m-2 s-1', 231, bgcrate, kiron, 'p', 'si<100', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('totalcarbon', 'total_carbon', &
         'mol C', 301, totalcarbon, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('zco2flux', 'global_carbon_flux', &
         'mol C/day', 303, zco2flux, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('totalcarboce', 'total_carbon_of_water_column_incl_pore_water', &
         'mol C', 305, totalcarboce, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('totalcarbsed', 'total_carbon_of_sediment', &
         'mol C', 306, totalcarbsed, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('calcitesediment', 'total_calcite_of_sediment', &
         'mol C', 307, calcitesediment, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('organsed', 'organic_carbon_of_sediment', &
         'mol C', 308, organsed, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('totalphos', 'total_phosphate', &
         'mol P', 309, totalphos, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('watphos', 'total_phosphate_of_water_column_incl_pore_water', &
         'mol P', 310, watphos, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('sedphos', 'total_phosphate_of_sediment', &
         'mol P', 311, sedphos, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('totalsil', 'total_silicate', &
         'mol Si', 312, totalsil, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('watsil', 'total_silicate_of_water_column_incl_pore_water', &
         'mol Si', 313, watsil, 'g', 's', factor = kilo))

    CALL varlist_add(bgc_varlist, &
         new_var('sedsil', 'total_silicate_of_sediment', &
         'mol Si', 314, sedsil, 'g', 's', factor = kilo))


    !Time series (output in netcdf allows higher code numbers)
    CALL varlist_add(bgc_varlist, &
         new_var('global_primary_production', 'global_primary_production', &
         'kmol P s-1', 500, bgcprod, kphosy, 'g', 'si'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_zooplankton_grazing', 'global_zooplankton_grazing', &
         'kmol P s-1', 501, bgcprod, kgraz, 'g', 'si'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_OM_export_at_90m', 'global_OM_export_at_90m', &
         'kmol P s-1', 502, bgcflux, kcoex90, 'g', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_calc_export_at_90m', 'global_calc_export_at_90m', &
         'kmol C s-1', 503, bgcflux, kcaex90, 'g', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_opal_export_at_90m', 'global_opal_export_at_90m', &
         'kmol Si s-1', 504, bgcflux, kopex90, 'g', 's'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_opal_production', 'global_opal_production', &
         'kmol Si s-1', 505, bgcprod, kdelsil, 'g', 'si'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_caco3_production', 'global_caco3_production', &
         'kmol C s-1', 506, bgcprod, kdelcar, 'g', 'si'))
    CALL varlist_add(bgc_varlist, &
         new_var('global_net_co2_flux', 'global_net_co2_flux', &
         'kmol C s-1', 507, bgcflux, kco2flux, 'g', 's'))

    IF ( p_pe == p_io ) CALL print_varlist('hamocc.partab', bgc_varlist)

    CALL generate_gridid(bgc_varlist)
    CALL generate_zaxisid(bgc_varlist)
    CALL generate_sed_zaxisid(bgc_varlist)

  END SUBROUTINE build_bgc_varlist

#endif/*ndef NO_NEW_IO */

END MODULE mo_bgc_varlist
