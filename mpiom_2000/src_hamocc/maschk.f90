      SUBROUTINE maschk(kpie,kpje,kpke,imas)
!*******************************************************************
!
!**** *INVENTORY_BGC* - calculate the BGC inventory.
!
!     P.Wetzel,              *MPI-Met, HH*    29.07.02
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     - calculate the BGC inventory.
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *INVENTORY_BGC*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_parallel, ONLY : have_g_js
      USE mo_commo1, ONLY : lbounds_exch_tp
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc

      USE MO_COMMO1
      USE mo_planetary_constants, ONLY : rhoicwa, rhosnwa

      USE mo_parallel

      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l,imas,jb

      REAL(wp) :: zpowtrato(npowtra)
      REAL(wp) :: zsedlayto(nsedtra),zburial(nsedtra)
      REAL(wp) :: zocetrato(nocetra)
      REAL(wp) :: ztotvol,ztotarea,vol,zsedhplto
      REAL(wp) :: zhito,zco3to,sum,sum2,sum3,zprorca,zprcaca,zsilpro
      REAL(wp) :: zatmco2,zatmo2,zatmn2
      REAL(wp) :: sco2flux,so2flux,sn2flux,sn2oflux
      REAL(wp) :: totalcarbon,totalphos,totalsil,totalnitr,totaloxy
      REAL(wp) :: totalalk,zalkn2,zh2o
      REAL(wp) :: ozkohl

#ifdef __c_isotopes
      REAL(wp) :: sc13flux,sc14flux
      REAL(wp) :: zatm13c,zatm14c
      REAL(wp) :: total13c,total14c,zpror13,zpror14,zprca13,zprca14
#endif
     return

      jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

!      if (imas.le.100) return !----------------------------------------------

! aqueous sediment tracer
!----------------------------------------------------------------------
      ztotvol = 0._wp
      DO k=1,ks
      DO j=jb,kpje-1
      DO i=2,kpie-1
        ztotvol=ztotvol+WETO(I,J,1)*seddw(k)                          &
     &         *dlxp(i,j)*dlyp(i,j)*porwat(k)
      ENDDO
      ENDDO
      ENDDO

      CALL global_sum(ztotvol)

      DO l=1,npowtra
         zpowtrato(l) = 0.0_wp
         DO k=1,ks
         DO j=jb,kpje-1
         DO i=2,kpie-1
            vol    = seddw(k)*dlxp(i,j)*dlyp(i,j)*porwat(k)
            zpowtrato(l)=zpowtrato(l)+WETO(I,J,1)*powtra(i,j,k,l)*vol
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      CALL global_sum(zpowtrato)

!
        zalkn2 = 0._wp
        zh2o = 0._wp
         DO k=1,ks
         DO j=jb,kpje-1
         DO i=2,kpie-1
            vol    = seddw(k)*dlxp(i,j)*dlyp(i,j)*porwat(k)
            zalkn2=zalkn2+WETO(I,J,1)*pown2bud(i,j,k)*vol  ! calculated implicit effect of N2 changes in sediment
            zh2o = zh2o+WETO(I,J,1)*powh2obud(i,j,k)*vol ! budget for O2 when water is "produced"
         ENDDO
         ENDDO
         ENDDO


      WRITE(io_stdo_bgc,*) 'aufruf maschk ',imas
      WRITE(io_stdo_bgc,*)'Global inventory of aqueous sediment tracer'
      WRITE(io_stdo_bgc,*)'-------------------------------------------'
      WRITE(io_stdo_bgc,*) '       total[kmol]    concentration[mol/L]'
      DO l=1,npowtra
      WRITE(io_stdo_bgc,*)'No. ',l,' ', &
     &       zpowtrato(l),'  ',zpowtrato(l)/ztotvol
      ENDDO
      WRITE(io_stdo_bgc,*) ' '

! non aqueous sediment tracer
!----------------------------------------------------------------------
      DO l=1,nsedtra
        zsedlayto(l) = 0.0_wp
        DO k=1,ks
        DO j=jb,kpje-1
        DO i=2,kpie-1
          vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
          zsedlayto(l)=zsedlayto(l)+WETO(I,J,1)*sedlay(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zsedlayto)

      DO l=1,nsedtra
        zburial(l)   = 0.0_wp
        DO j=jb,kpje-1
        DO i=2,kpie-1
           zburial(l)=zburial(l)+burial(i,j,l)                         &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zburial)

         zsedhplto = 0.0_wp
         DO k=1,ks
         DO j=jb,kpje-1
         DO i=2,kpie-1
           vol    = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
           zsedhplto= zsedhplto+WETO(I,J,1)*sedhpl(i,j,k)*vol
         ENDDO
         ENDDO
         ENDDO
         CALL global_sum(zsedhplto)




      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*)                                            &
     &     'Global inventory of solid sediment constituents'
      WRITE(io_stdo_bgc,*)                                            &
     &     '----------------------------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

      DO l=1,nsedtra
        WRITE(io_stdo_bgc,*) 'Sediment No. ',l,' ',                   &
     &                        zsedlayto(l)
        WRITE(io_stdo_bgc,*) 'Burial No. ',l,' ',                     &
     &                        zburial(l)
      ENDDO
      WRITE(io_stdo_bgc,*) 'hpl ',                                    &
     &       zsedhplto
      WRITE(io_stdo_bgc,*) ' '

!
!  oceanic tracers
!----------------------------------------------------------------------

      ztotvol = 0._wp
      k=1
      DO j=jb,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k) &
     &                     +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                     -SICSNO(I,J)*RHOSNWA)
      ENDDO
      ENDDO


      DO k=2,kpke
      DO j=jb,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
      ENDDO
      ENDDO
      ENDDO

      CALL global_sum(ztotvol)

      DO l=1,nocetra
        zocetrato(l) = 0.0_wp
        k=1
        DO j=jb,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zocetrato(l)=zocetrato(l)+WETO(I,J,K)*ocetra(i,j,k,l)*vol
        ENDDO
        ENDDO

        DO k=2,kpke
        DO j=jb,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zocetrato(l)=zocetrato(l)+WETO(I,J,K)*ocetra(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zocetrato)

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global inventory of advected ocean tracers'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
      WRITE(io_stdo_bgc,*) ' '
      DO l=1,nocetra
      WRITE(io_stdo_bgc,*) 'No. ',l,    &
     &         zocetrato(l),zocetrato(l)/ztotvol
      ENDDO

! additional ocean tracer
!----------------------------------------------------------------------

        zhito = 0._wp
        zco3to = 0._wp

        k=1
        DO j=jb,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zhito     = zhito + WETO(I,J,K)*hi(i,j,k) *vol
          zco3to    = zco3to+ WETO(I,J,K)*co3(i,j,k)*vol
          zalkn2=zalkn2+WETO(I,J,K)*n2budget(i,j,k)*vol      ! calculated implicit effect of N2 changes
          zh2o = zh2o+WETO(I,J,K)*h2obudget(i,j,k)*vol ! budget for O2 when water is "produced"
        ENDDO
        ENDDO

        DO k=2,kpke
        DO j=jb,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zhito     = zhito + WETO(I,J,K)*hi(i,j,k) *vol
          zco3to    = zco3to+ WETO(I,J,K)*co3(i,j,k)*vol
          zalkn2=zalkn2+WETO(I,J,K)*n2budget(i,j,k)*vol
          zh2o = zh2o+WETO(I,J,K)*h2obudget(i,j,k)*vol
        ENDDO
        ENDDO
        ENDDO

      CALL global_sum(zhito,zco3to,zalkn2,zh2o)

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Glob. inventory of additional ocean tracer'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) '      total[kmol]  concentration[kmol/m^3]'
      WRITE(io_stdo_bgc,*) ' '

      WRITE(io_stdo_bgc,*) ' hi',            &
     &       zhito,zhito/ztotvol
      WRITE(io_stdo_bgc,*) ' co3',           &
     &       zco3to,zco3to/ztotvol
      WRITE(io_stdo_bgc,*) ' '

! atmosphere flux and atmospheric CO2
!--------------------------------------------------------------------

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global fluxes into atmosphere'
      WRITE(io_stdo_bgc,*) '-----------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

      sco2flux = 0._wp
      so2flux  = 0._wp
      sn2flux  = 0._wp
      sn2oflux = 0._wp
      zatmco2  = 0._wp
      zatmo2   = 0._wp
#ifdef __c_isotopes
      sc13flux = 0._wp
      sc14flux = 0._wp
      zatm13c = 0._wp
      zatm14c = 0._wp
#endif
      zatmn2 = 0._wp
      ztotarea = 0._wp

      do i=2,kpie-1
      do j=jb,kpje-1
        sco2flux =sco2flux +co2flux(i,j)*dlxp(i,j)*dlyp(i,j)
#ifdef __c_isotopes
        sc13flux =sc13flux +c13flux(i,j)*dlxp(i,j)*dlyp(i,j)
        sc14flux =sc14flux +c14flux(i,j)*dlxp(i,j)*dlyp(i,j)
#endif
        so2flux =so2flux +o2flux(i,j) *dlxp(i,j)*dlyp(i,j)
        sn2flux =sn2flux +n2flux(i,j) *dlxp(i,j)*dlyp(i,j)
        sn2oflux =sn2oflux +n2oflux(i,j) *dlxp(i,j)*dlyp(i,j)
        ztotarea = ztotarea + dlxp(i,j)*dlyp(i,j)
        if(diffat)then
         zatmco2=zatmco2 + atm(i,j,iatmco2)*dlxp(i,j)*dlyp(i,j)
         zatmo2= zatmo2  + atm(i,j,iatmo2) *dlxp(i,j)*dlyp(i,j)
         zatmn2= zatmn2  + atm(i,j,iatmn2) *dlxp(i,j)*dlyp(i,j)
        endif

#ifdef __c_isotopes
        zatm13c= zatm13c  + atm(i,j,iatmc13) *dlxp(i,j)*dlyp(i,j)
        zatm14c= zatm14c  + atm(i,j,iatmc14) *dlxp(i,j)*dlyp(i,j)
#endif
      enddo
      enddo

      CALL global_sum(sco2flux,so2flux,sn2flux,sn2oflux,ztotarea)
#ifdef __c_isotopes
      CALL global_sum(sc13flux,sc14flux)
#endif


      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'CO2Flux  :',sco2flux
      WRITE(io_stdo_bgc,*) 'O2 Flux  :',so2flux
      WRITE(io_stdo_bgc,*) 'N2 Flux  :',sn2flux
      WRITE(io_stdo_bgc,*) 'N2O Flux  :',sn2oflux
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'DIC+FL :',zocetrato(isco212)+sco2flux,imas
! only when diffusive atm
      if(diffat) then
      CALL global_sum(zatmco2,zatmo2,zatmn2)
      WRITE(io_stdo_bgc,*) 'global atm. CO2[ppm] / kmol: ',          &
     &                               zatmco2/ztotarea,zatmco2*ppm2con
      WRITE(io_stdo_bgc,*) 'global atm. O2[ppm] / kmol : ',          &
     &                               zatmo2/ztotarea,zatmo2*ppm2con
      WRITE(io_stdo_bgc,*) 'global atm. N2[ppm] / kmol : ',          &
     &                               zatmn2/ztotarea,zatmn2*ppm2con
      endif

! Complete sum of inventory in between bgc.f90

        zprorca = 0._wp
        zprcaca = 0._wp
#ifdef __c_isotopes
        zpror13 = 0._wp
        zprca13 = 0._wp
        zpror14 = 0._wp
        zprca14 = 0._wp
#endif
        zsilpro = 0._wp
        DO j=jb,kpje-1
        DO i=2,kpie-1
           zprorca=zprorca+prorca(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zprcaca=zprcaca+prcaca(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
#ifdef __c_isotopes
           zpror13=zpror13+pror13(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zpror14=zpror14+pror14(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zprca13=zprca13+prca13(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zprca14=zprca14+prca14(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
#endif
           zsilpro=zsilpro+silpro(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO

      CALL global_sum(zprorca,zprcaca,zsilpro)
#ifdef __c_isotopes
      CALL global_sum(zpror13,zpror14,zprca13,zprca14)
#endif

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Should be zero at the end: '
      WRITE(io_stdo_bgc,*) 'prorca, prcaca, silpro  ',                &
     &       zprorca, zprcaca, zsilpro
      WRITE(io_stdo_bgc,*) ' '


! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2

! Sum of inventory
!----------------------------------------------------------------------
! Units in P have a C:P Ratio of 122:1
!e      write(0,*)'geloestes DIC',zocetrato(isco212),imas
      totalcarbon=                                                    &
     & (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)               &
     & +zocetrato(izoo))*rcar+zocetrato(isco212)+zocetrato(icalc)     &
     & +zpowtrato(ipowaic)+zsedlayto(isssc12)+zsedlayto(issso12)*rcar &
     & +zburial(isssc12)+zburial(issso12)*rcar+zprorca*rcar+zprcaca   &
     & +sco2flux

      ozkohl=                                                         &
     & (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)               &
     & +zocetrato(izoo))*rcar+zocetrato(isco212)+zocetrato(icalc)     &
     & +zpowtrato(ipowaic)+zsedlayto(isssc12)+zsedlayto(issso12)*rcar &
     & +zburial(isssc12)+zburial(issso12)*rcar+zprorca*rcar+zprcaca

#ifdef __c_isotopes
      total13c=                                                    &
     & zocetrato(idet13)+(zocetrato(idoc)+zocetrato(iphy)             &
     & +zocetrato(izoo))*rcar*bifr13+zocetrato(isco213)               &
     & +zocetrato(icalc13)     &
     & +zpowtrato(ipowc13)+zsedlayto(isssc13)+zsedlayto(issso13) &
     & +zburial(isssc13)+zburial(issso13)*rcar+zpror13*rcar+zprca13   &
     & +sc13flux

      total14c=                                                    &
     & zocetrato(idet14)+zocetrato(isco214)               &
     & +zocetrato(icalc14)     &
     & +zpowtrato(ipowc14)+zsedlayto(isssc14)+zsedlayto(issso14) &
     & +zburial(isssc14)+zburial(issso14)*rcar+zpror14*rcar+zprca14   &
     & +sc14flux
#endif /*__c_isotopes*/

     totalalk=                                                        &
     & zocetrato(ialkali)+ zpowtrato(ipowaal)                   & ! alk in water/pore water
     & -rnit*(zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)   & ! org. compounds water
     &        +zocetrato(izoo)                                  &
     &        +zprorca+zsedlayto(issso12)+zburial(issso12))     & ! org. compounds sed.
     & - zalkn2                                                 & ! compensation for biogenic induced change in H+
                                                                  ! due to N2 fixation and N2 production
     & + 2._wp * (zocetrato(icalc)                              & ! calcerous compounds water
     &      +zsedlayto(isssc12)+zburial(isssc12)+zprcaca)         ! calcerous compounds sed.




      totalnitr=                                                      &
     &   (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)             &
     &  +zocetrato(izoo))*rnit+zocetrato(iano3)+zocetrato(igasnit)*2._wp &
     &  +zpowtrato(ipowno3)+zpowtrato(ipown2)*2._wp                   &
     &  +zsedlayto(issso12)*rnit+zburial(issso12)*rnit                &
     &  +zocetrato(ian2o)*2._wp + zprorca*rnit                        &
     & +(sn2flux+sn2oflux) * 2._wp

      totalphos=                                                      &
     &   zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)              &
     &  +zocetrato(izoo)+zocetrato(iphosph)                           &
     &  +zpowtrato(ipowaph)+zsedlayto(issso12)+zburial(issso12)       &
     &  +zprorca

      totalsil=                                                       &
     &   zocetrato(isilica)+zocetrato(iopal)                          &
     &  +zpowtrato(ipowasi)+zsedlayto(issssil)+zburial(issssil)       &
     &  +zsilpro

      totaloxy=                                                       &
     &  (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)              &
     &  +zocetrato(izoo))*(-24._wp)+zocetrato(ioxygen)                   &
     &  +zocetrato(iphosph)*2._wp +zocetrato(isco212)+zocetrato(icalc)   &
     &  +zocetrato(iano3)*1.5_wp+zocetrato(ian2o)*0.5_wp                    &
     &  +zsedlayto(issso12)*(-24._wp) + zsedlayto(isssc12)               &
     &  +zburial(issso12)*(-24._wp)   +   zburial(isssc12)              &
     &  +zpowtrato(ipowno3)*1.5_wp+zpowtrato(ipowaic)                    &
     &  +zpowtrato(ipowaox)+zpowtrato(ipowaph)*2._wp                     &
     &  +zprorca*(-24._wp)+zprcaca                                       &
     & +so2flux+sn2oflux*0.5_wp+sco2flux                                 &
     & +zh2o

#ifdef __c_isotopes
      if(imas.eq.100) then
      total14c=                                                    &
     & zocetrato(idet14)+zocetrato(isco214)                     &
     & +zocetrato(icalc14)     &
     & +zpowtrato(ipowc14)+zsedlayto(isssc14)+zsedlayto(issso14) &
     & +zburial(isssc14)+zburial(issso14)*rcar+zpror14*rcar+zprca14   &
     & +atm_c14*atmacmol
      c14inva=total14c
       c14prosta=total14c*c14dec/atmacmol
       c14prod=c14prosta
      write(io_stdo_bgc,*)'prodc14',atmacmol*c14prosta,atmacmol*c14prod
      write(io_stdo_bgc,*)'prodc14',atmacmol*(c14prosta-c14prod)
      endif
#endif
      WRITE(io_stdo_bgc,*) ' '

#ifdef __c_isotopes
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of kohlen   : ',       &
     & totalcarbon,ozkohl,atm_co2*atmacmol,total13c,total14c,imas
#else
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of kohlen   : ',       &
     & totalcarbon,ozkohl,sco2flux,atm_co2*atmacmol,imas
#endif


      atcoa=atm_co2*atmacmol
      ozkoa=ozkohl
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of phosph.  : ',       &
     & totalphos,imas
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of silicate : ',       &
     & totalsil,imas
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of nitrogen.  : ',     &
     & totalnitr,imas
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of oxygen.  : ',     &
     & totaloxy,imas
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of alkalinity  : ',    &
     & totalalk,zalkn2,imas
      WRITE(io_stdo_bgc,*) ' '

! Write sediment fluxes (integrated over time)

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global fluxes into and out of the sediment'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

        zprorca = 0._wp
        zprcaca = 0._wp
        zsilpro = 0._wp
        DO j=jb,kpje-1
        DO i=2,kpie-1
           zprorca=zprorca+bgct2d(i,j,jprorca)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zprcaca=zprcaca+bgct2d(i,j,jprcaca)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zsilpro=zsilpro+bgct2d(i,j,jsilpro)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO
        CALL global_sum(zprorca,zprcaca,zsilpro)
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Detritus, Calcium Carbonate, Silicate  ', &
     &       zprorca, zprcaca, zsilpro
      WRITE(io_stdo_bgc,*) ' '

      DO l=1,npowtra
        sum = 0._wp
      do i=2,kpie-1
      do j=jb,kpje-1
        sum=sum+sedfluxo(i,j,l)*dlxp(i,j)*dlyp(i,j)
      enddo
      enddo
      CALL global_sum(sum)
      WRITE(io_stdo_bgc,*) 'diffusive flux from sediment to ocean'
      WRITE(io_stdo_bgc,*) 'No. ',l,' ',sum
      ENDDO

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total detritus und shell production'
      WRITE(io_stdo_bgc,*) '------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'
      WRITE(io_stdo_bgc,*) 'integrated over all model timesteps '


        sum = 0._wp
        sum2 = 0._wp
        sum3 = 0._wp
      do i=2,kpie-1
      do j=jb,kpje-1
        sum=sum+expoor(i,j)*dlxp(i,j)*dlyp(i,j)
        sum2=sum2+expoca(i,j)*dlxp(i,j)*dlyp(i,j)
        sum3=sum3+exposi(i,j)*dlxp(i,j)*dlyp(i,j)
      enddo
      enddo
      CALL global_sum(sum,sum2,sum3)

      zprorca = sum * 12._wp * 1.e-12_wp
      zprcaca = sum2 * 12._wp * 1.e-12_wp
      WRITE(io_stdo_bgc,*) 'organic carbon   : ',sum,' kmol ',zprorca,' GtC'
      WRITE(io_stdo_bgc,*) 'calcium carbonate: ',sum2,' kmol ',zprcaca,' GtC'
      WRITE(io_stdo_bgc,*) 'silicate         : ',sum3,' kmol '
      WRITE(io_stdo_bgc,*) ' '

      RETURN
      END
