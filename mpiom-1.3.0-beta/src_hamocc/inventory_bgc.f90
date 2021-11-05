      SUBROUTINE INVENTORY_BGC(kpie,kpje,kpke)
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
!     GLOBAL_SUM (../src_oce/mo_parallel.f90)
!
!**********************************************************************
    
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 
      
      USE MO_COMMO1
      USE MO_COMMOAU1

      USE mo_parallel

      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l

      REAL :: zpowtrato(npowtra)
      REAL :: zsedlayto(nsedtra),zburial(nsedtra)
      REAL :: zocetrato(nocetra)
      REAL :: zalkali
      REAL :: ztotvol,ztotarea,vol,zsedhplto
      REAL :: zhito,zco3to,sum,sum2,sum3,zprorca,zprcaca,zsilpro
      REAL :: zatmco2,zatmo2,zatmn2
      REAL :: sco2flux,so2flux,sn2flux,sn2oflux
      REAL :: totalcarbon,totalphos,totalsil,totalnitr,totaloxy
      REAL :: ppm2con, co2atm

      
! aqueous sediment tracer
!----------------------------------------------------------------------
      ztotvol=0.
      DO k=1,ks
      DO j=2,kpje-1
      DO i=2,kpie-1
        ztotvol=ztotvol+WETO(I,J,1)*seddw(k)                          &
     &              *dlxp(i,j)*dlyp(i,j)*porwat(k)
      ENDDO
      ENDDO
      ENDDO

      CALL global_sum(ztotvol)

      DO l=1,npowtra
         zpowtrato(l) = 0.0
         DO k=1,ks
         DO j=2,kpje-1
         DO i=2,kpie-1
            vol    = seddw(k)*dlxp(i,j)*dlyp(i,j)*porwat(k)
            zpowtrato(l)=zpowtrato(l)+WETO(I,J,1)*powtra(i,j,k,l)*vol
         ENDDO
         ENDDO
         ENDDO
      ENDDO
      CALL global_sum(zpowtrato)


      WRITE(io_stdo_bgc,*) ' '
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
        zsedlayto(l) = 0.0
        DO k=1,ks
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol = porsol(k)*seddw(k)*dlxp(i,j)*dlyp(i,j)
          zsedlayto(l)=zsedlayto(l)+WETO(I,J,1)*sedlay(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zsedlayto)

      DO l=1,nsedtra
      zburial(l)   = 0.0
        DO j=2,kpje-1
        DO i=2,kpie-1
         zburial(l)=zburial(l)+burial(i,j,l)                         &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zburial)

         zsedhplto = 0.0
         DO k=1,ks
         DO j=2,kpje-1
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

      ztotvol=0.
      k=1
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k) &
     &                     +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                     -SICSNO(I,J)*RHOSNWA)       
      ENDDO
      ENDDO     
      
      
      DO k=2,kpke
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k) 
      ENDDO
      ENDDO
      ENDDO

      CALL global_sum(ztotvol)

      DO l=1,nocetra
        zocetrato(l) = 0.0
        k=1
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zocetrato(l)=zocetrato(l)+WETO(I,J,K)*ocetra(i,j,k,l)*vol
        ENDDO
        ENDDO      
      
        DO k=2,kpke
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zocetrato(l)=zocetrato(l)+WETO(I,J,K)*ocetra(i,j,k,l)*vol
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      CALL global_sum(zocetrato)
    
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global inventory of ocean tracers'
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

        zhito =0.
        zco3to=0.
        k=1
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zhito     = zhito + WETO(I,J,K)*hi(i,j,k) *vol
          zco3to    = zco3to+ WETO(I,J,K)*co3(i,j,k)*vol
        ENDDO
        ENDDO      
      
        DO k=2,kpke
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*DDPO(i,j,k)
          zhito     = zhito + WETO(I,J,K)*hi(i,j,k) *vol
          zco3to    = zco3to+ WETO(I,J,K)*co3(i,j,k)*vol        
        ENDDO
        ENDDO
        ENDDO

      CALL global_sum(zhito,zco3to)
             
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


! alkalinity of the surface layer
!-------------------------------------------------------------------- 

      ztotvol=0.
      k=1
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotvol=ztotvol+WETO(I,J,K)*dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k) &
     &                     +ZO(I,J)-SICTHO(I,J)*RHOICWA               &
     &                     -SICSNO(I,J)*RHOSNWA)       
      ENDDO
      ENDDO     

        zalkali = 0.0
        k=1
        DO j=2,kpje-1
        DO i=2,kpie-1
          vol    = dlxp(i,j)*dlyp(i,j)*(DDPO(i,j,k)                   &
     &                   +ZO(I,J)-SICTHO(I,J)*RHOICWA                 &
     &                   -SICSNO(I,J)*RHOSNWA)
          zalkali=zalkali+WETO(I,J,K)*ocetra(i,j,k,ialkali)*vol
        ENDDO
        ENDDO      

      CALL global_sum(ztotvol,zalkali)
    
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global inventory of surface layer alkalinity'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) '       total[kmol]  concentration[kmol/m^3]'
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) zalkali,zalkali/ztotvol
      


! atmosphere flux and atmospheric CO2
!-------------------------------------------------------------------- 

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global fluxes into atmosphere'
      WRITE(io_stdo_bgc,*) '-----------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

      sco2flux  =0.
      so2flux  =0.
      sn2flux  =0.
      sn2oflux =0.
      zatmco2  =0.
      zatmo2   =0.
      zatmn2   =0.
      ztotarea =0.
      ppm2con=0.35e-3
      do i=2,kpie-1
      do j=2,kpje-1
        sco2flux =sco2flux +bgct2d(i,j,jco2flux)*dlxp(i,j)*dlyp(i,j)
      so2flux =so2flux +bgct2d(i,j,jo2flux) *dlxp(i,j)*dlyp(i,j)
      sn2flux =sn2flux +bgct2d(i,j,jn2flux) *dlxp(i,j)*dlyp(i,j)
      sn2oflux=sn2oflux+bgct2d(i,j,jn2oflux)*dlxp(i,j)*dlyp(i,j)
      ztotarea = ztotarea + dlxp(i,j)*dlyp(i,j)
#ifdef DIFFAT      
      zatmco2=zatmco2 + atm(i,j,iatmco2)*dlxp(i,j)*dlyp(i,j)
      zatmo2= zatmo2  + atm(i,j,iatmo2) *dlxp(i,j)*dlyp(i,j)
      zatmn2= zatmn2  + atm(i,j,iatmn2) *dlxp(i,j)*dlyp(i,j)      
#endif
      enddo
      enddo

      CALL global_sum(sco2flux,so2flux,sn2flux,sn2oflux,ztotarea)
#ifdef DIFFAT
      CALL global_sum(zatmco2,zatmo2,zatmn2)
#endif

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'CO2Flux  :',sco2flux
      WRITE(io_stdo_bgc,*) 'O2 Flux  :',so2flux
      WRITE(io_stdo_bgc,*) 'N2 Flux  :',sn2flux
      WRITE(io_stdo_bgc,*) 'N2O Flux :',sn2oflux
      WRITE(io_stdo_bgc,*) ' '
#ifdef DIFFAT            
      WRITE(io_stdo_bgc,*) 'global atm. CO2[ppm] / kmol: ',          &
     &                               zatmco2/ztotarea,zatmco2*ppm2con       
      WRITE(io_stdo_bgc,*) 'global atm. O2[ppm] / kmol : ',          &
     &                               zatmo2/ztotarea,zatmo2*ppm2con 
      WRITE(io_stdo_bgc,*) 'global atm. N2[ppm] / kmol : ',          &
     &                               zatmn2/ztotarea,zatmn2*ppm2con 
     
#endif /*DIFFAT*/

! Complete sum of inventory in between bgc.f90 (js: don't understand this comment)

        zprorca=0.
      zprcaca=0.
      zsilpro=0.
        DO j=2,kpje-1
        DO i=2,kpie-1
         zprorca=zprorca+prorca(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
              zprcaca=zprcaca+prcaca(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
           zsilpro=zsilpro+silpro(i,j)                        &
     &               *WETO(I,J,1)*dlxp(i,j)*dlyp(i,j)
        ENDDO
        ENDDO

      CALL global_sum(zprorca,zprcaca,zsilpro)

!js: the print statement 'should be zero in the end' makes no sense to me
!     because it is set to zero in powach.f90? so why sum it up?
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
  
      totalcarbon=                                                    &
     & (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)               &
     & +zocetrato(izoo))*rcar+zocetrato(isco212)+zocetrato(icalc)     &
     & +zpowtrato(ipowaic)+zsedlayto(isssc12)+zsedlayto(issso12)*rcar &
     & +zburial(isssc12)+zburial(issso12)*rcar+zprorca*rcar+zprcaca   &
#ifdef DIFFAT     
     & +zatmco2*ppm2con  
#else
     & +sco2flux
#endif
      totalnitr=                                                      &
     &   (zocetrato(idet)+zocetrato(idoc)+zocetrato(iphy)             &
     &  +zocetrato(izoo))*rnit+zocetrato(iano3)+zocetrato(igasnit)*2  &
     &  +zpowtrato(ipowno3)+zpowtrato(ipown2)*2                       &
     &  +zsedlayto(issso12)*rnit+zburial(issso12)*rnit                &
     &  +zocetrato(ian2o)*2+zprorca*rnit                              &
#ifdef DIFFAT     
     &  +zatmn2*ppm2con*2
#else
     & +sn2flux*2+sn2oflux*2
#endif     
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
     &  +zocetrato(izoo))*(-24.)+zocetrato(ioxygen)                   &
     &  +zocetrato(iphosph)*2 +zocetrato(isco212)+zocetrato(icalc)    &
     &  +zocetrato(iano3)*1.5+zocetrato(ian2o)*0.5                    &
     &  +zsedlayto(issso12)*(-24.) + zsedlayto(isssc12)               &
!     &  +zburial(issso12)*(-24.)   +   zburial(isssc12)               &
     &  +zpowtrato(ipowno3)*1.5+zpowtrato(ipowaic)                    &
     &  +zpowtrato(ipowaox)+zpowtrato(ipowaph)*2                      &
     &  +zprorca*(-24.)+zprcaca                                       &
#ifdef DIFFAT     
     &  +zatmo2*ppm2con+zatmco2*ppm2con
#else
     & +so2flux+sn2oflux*0.5+sco2flux
#endif     

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of carbon   : ',       &
     & totalcarbon
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of phosph.  : ',       &
     & totalphos
      WRITE(io_stdo_bgc,*) ' '          
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of silicate : ',       &
     & totalsil
      WRITE(io_stdo_bgc,*) ' '          
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of nitrogen.  : ',     &
     & totalnitr
      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total[kmol] of oxygen.  : ',     &
     & totaloxy
      WRITE(io_stdo_bgc,*) ' '

! Write sediment fluxes

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global fluxes into and out of the sediment'
      WRITE(io_stdo_bgc,*) '------------------------------------------'
      WRITE(io_stdo_bgc,*) '        [kmol]'

        zprorca=0.
      zprcaca=0.
      zsilpro=0.
        DO j=2,kpje-1
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
        sum=0.
      do i=2,kpie-1
      do j=2,kpje-1
        sum=sum+sedfluxo(i,j,l)*dlxp(i,j)*dlyp(i,j)
      enddo
      enddo      
      CALL global_sum(sum)
      WRITE(io_stdo_bgc,*) 'No. ',l,' ',sum
      ENDDO

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Global total export production (from expoor, not coex90)'
      WRITE(io_stdo_bgc,*) '------------------------dims ',kpie,'*',kpje
      WRITE(io_stdo_bgc,*) '        [kmol]', '     n90depth= ',n90depth


      sum=0.
      sum2=0.
      sum3=0.
      do i=2,kpie-1
      do j=2,kpje-1
!       sum=sum+expoor(i,j)*dlxp(i,j)*dlyp(i,j)  *weto(i,j,n90depth)
!       sum2=sum2+expoca(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,n90depth)
!       sum3=sum3+exposi(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,n90depth)      
        sum=sum+expoor(i,j)*dlxp(i,j)*dlyp(i,j)  *weto(i,j,1)
        sum2=sum2+expoca(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)
        sum3=sum3+exposi(i,j)*dlxp(i,j)*dlyp(i,j)*weto(i,j,1)      
      enddo
      enddo
      CALL global_sum(sum,sum2,sum3)

      WRITE(io_stdo_bgc,*) 'organic carbon   : ',sum
      WRITE(io_stdo_bgc,*) 'calcium carbonate: ',sum2
      WRITE(io_stdo_bgc,*) 'silicate         : ',sum3
      WRITE(io_stdo_bgc,*) ' '

      RETURN
      END
