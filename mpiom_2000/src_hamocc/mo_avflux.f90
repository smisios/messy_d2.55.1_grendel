      MODULE MO_AVFLUX

#if defined __cpl_co2
      USE mo_kind, ONLY: wp
      USE mo_param1,            ONLY: ie,je,ie_g,je_g
      USE mo_parallel

      IMPLICIT NONE

      REAL(wp), ALLOCATABLE :: areas_g(:,:)
      REAL(wp), ALLOCATABLE :: co2flux_act_g(:,:)
      REAL(wp), ALLOCATABLE :: corrfrac(:,:)
      REAL(wp), ALLOCATABLE :: co2flux_cpl_g(:,:)
      REAL(wp), ALLOCATABLE :: pdlxp_g(:,:)
      REAL(wp), ALLOCATABLE :: pdlyp_g(:,:)
      INTEGER, ALLOCATABLE :: avflux_index(:,:)
      REAL(wp) :: co2ctrl_old,co2ctrl_new

      REAL(wp), ALLOCATABLE :: co2flux_act(:,:)  !actual ocean flux pattern
      REAL(wp), ALLOCATABLE :: areas_l(:,:)      ! local surface area of grid box

      REAL(wp), ALLOCATABLE :: corsum(:)
      REAL(wp), ALLOCATABLE :: avfcpl(:)            !averaged ECHAM flux pattern to correct for small scale sturcture
      REAL(wp), ALLOCATABLE :: avfact(:,:),avarea(:,:)  !averaged ocean flux pattern, area of flux correction
      REAL(wp) :: newcpl,diffcpl,quotdir,quottot,totact,totdiv
      INTEGER :: index_num, lmsum
      INTEGER, ALLOCATABLE :: totcor(:)

      CONTAINS

      SUBROUTINE AVFLUX_INI

        use mo_kind
        use mo_commo1, only : dlxp,dlyp,lwetol1_g,alat_g,alon_g

        integer :: i,j,IT,Ic,ila,ilo
        REAL(wp) :: dll,rlats,rlatl,rlons,rlonl
        allocate(areas_l(ie,je))

        do j=2,je-1
           do i=2,ie-1
              areas_l(i,j)=dlxp(i,j)*dlyp(i,j)
           enddo
        enddo

        allocate(co2flux_cpl_g(ie_g,je_g))
        allocate(co2flux_act(ie,je))


        co2flux_act(:,:)=0.0


        if ( p_pe == p_io ) then
           allocate(co2flux_act_g(ie_g,je_g))
           allocate(areas_g(ie_g,je_g))
           allocate(avflux_index(ie_g,je_g))
           allocate(corrfrac(ie_g,je_g))
           allocate(corsum(ie_g*je_g))
           allocate(avfcpl(ie_g*je_g))
           allocate(avfact(ie_g*je_g,2))
           allocate(avarea(ie_g*je_g,2))
           allocate(totcor(ie_g*je_g))
        else
           allocate(co2flux_act_g(0,0))
           allocate(areas_g(0,0))
           allocate(avflux_index(0,0))
           allocate(corrfrac(0,0))
           allocate(corsum(0))
           allocate(avfcpl(0))
           allocate(avfact(0,0))
           allocate(avarea(0,0))
           allocate(totcor(0))
        endif


        !collect grid box area
        call gather(areas_l,areas_g,p_io)

#ifdef  AVFLUX_KS
        IF ( p_pe == p_io ) THEN

           open(178,file='AVFLUX_INDEX.txt',status='unknown',    &
                access='sequential',form='formatted')
           rewind (178)
           read(178,'(10I8)') ((avflux_index(i,j),i=1,ie_g),j=1,je_g)
           close(178)

!           do i=1,ie_g
!          do j=1,je_g
!          Write(788,*)i,j,avflux_index(i,j)
!           enddo
!          enddo

        ENDIF
#else

        IF (p_pe==p_io) THEN

           dll=5.

           it=1
           ic=0

           avflux_index(:,:)=-99999

           do ila=-90,90,int(dll)
              do ilo=0,360,int(dll)

                 rlats=real(ila,wp)
                 rlons=real(ilo,wp)

                 rlatl=min(90.,rlats+dll)
                 rlonl=min(360.,rlons+dll)

                 do j=1,JE_G
                    do i=2,IE_G-1
                       if (alat_g(i,j).ge.rlats.and.alat_g(i,j).lt.rlatl) then
                          if (alon_g(i,j).ge.rlons.and.alon_g(i,j).lt.rlonl) then
                             if(.not. lwetol1_g(i,j)) then
                                avflux_index(i,j)=-99999
                             else
                                avflux_index(i,j)=it
                                ic=ic+1
                             endif
                          endif
                       endif
                    enddo
                 enddo
                 if (ic.gt. 0) then
                   it=it+1
                   ic=0
                endif
              enddo
           enddo

!           do i=1,ie_g
!          do j=1,je_g
!          Write(789,*)i,j,avflux_index(i,j)
!           enddo
!          enddo
#endif

        endif


        index_num=ie_g*je_g               ! amount for averaging areas for flux redistribution


      end subroutine avflux_ini

      subroutine avflux_redis

        use mo_commo1, only : ddpo,zo,sictho,sicsno,almzer
        USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa
        use mo_param1_bgc, only : isco212
        use mo_carbch, only : ocetra, co2flux_cpl
        use mo_bgcmean, only: bgcm2d, jco2fxd
        use mo_control_bgc, only : dtbgc, io_stdo_bgc

        REAL(wp) :: thickness
        integer :: i, j, lm, lmtest,lmsum1

         !collect flux fields : co2flux_act, co2flux_cpl

         call gather(co2flux_act,co2flux_act_g,p_io)
         call gather_arr(co2flux_cpl,co2flux_cpl_g,p_io)



         IF (p_pe==p_io) THEN

!            lmtest=avflux_index(50,76)

            avarea(:,:)=0.                 ! mean areas - calculate only once in beleg?
            avfcpl(:)=0.                 ! mean coupled flux
            avfact(:,:)=0.                 ! mean flux of actual ocean pattern

            !sum up cpl_fluxes within averaging areas according to index-file



            co2ctrl_old=0.
            DO j=2,je_g-1
               DO i=2,ie_g-1
                  lm=avflux_index(i,j)
                  if(lm .gt. 0) then
                     avfcpl(lm)=avfcpl(lm)+co2flux_cpl_g(i,j)*areas_g(i,j)
                     co2ctrl_old=co2ctrl_old+co2flux_cpl_g(i,j)*areas_g(i,j)
                  end if
               enddo
            enddo


            ! sum up act_fluxes within averaging areas according to index-file and with
            ! respect to the flux direction 1: flux up (+), 2: flux down (-)

            DO j=2,je_g-1
               DO i=2,ie_g-1

                  lm=avflux_index(i,j)

                  if(lm .gt. 0) then
                     if(avfcpl(lm) .gt.0.) then ! flux out of ocean (lm,1) relevant

                        if(co2flux_act_g(i,j).gt.0.) then
                           avfact(lm,1)=avfact(lm,1)+co2flux_act_g(i,j)*areas_g(i,j)
                           avarea(lm,1)=avarea(lm,1)+areas_g(i,j)
                           corrfrac(i,j)=1.
                        else
                           avfact(lm,2)=avfact(lm,2)+co2flux_act_g(i,j)*areas_g(i,j)
                           avarea(lm,2)=avarea(lm,2)+areas_g(i,j)
                           corrfrac(i,j)=0.
                        endif

                     else          ! flux into ocean (lm,2) relevant

                        if(co2flux_act_g(i,j).gt.0.) then
                           avfact(lm,1)=avfact(lm,1)+co2flux_act_g(i,j)*areas_g(i,j)
                           avarea(lm,1)=avarea(lm,1)+areas_g(i,j)
                           corrfrac(i,j)=0.
                        else
                           avfact(lm,2)=avfact(lm,2)+co2flux_act_g(i,j)*areas_g(i,j)
                           avarea(lm,2)=avarea(lm,2)+areas_g(i,j)
                           corrfrac(i,j)=1.
                        endif
                     endif

                  endif

               enddo
            enddo



            ! check if cpl and act have same direction
            lmsum=0
            do lm=1,index_num
               totact=avfact(lm,1)+avfact(lm,2)

               !        if(avfact(lm,1).gt.0. .and. avfact(lm,2).lt.0.)  then
               !          write(0,*) 'both not zero',lm,avfact(lm,1),avfact(lm,2) &
               !      &   ,totact,avfcpl(lm)
               !        end if

               if(abs(avfact(lm,1)).lt.almzer.and.abs(avfact(lm,2)).lt.almzer) then
                  totcor(lm)=-1   ! ocean completely covered with ice, no correction
               else
                  totdiv=avfcpl(lm)/totact
                  if(totdiv.lt.0.) then
                     totcor(lm)=-1                      ! ocean fluxes in different
                     ! direction, no correction
                  else
                     totcor(lm)=1                       !correction of area lm
                     lmsum=lmsum+1
                  endif
               endif
            enddo
         !       write(0,*) 'totcor ',lmsum

         ! dont do a redistribution if all the coupled C02 fluxes are equal to zero
           lmsum1=0
           do lm=1,index_num
                if(abs(avfcpl(lm)).gt.almzer) lmsum1=lmsum1+1
           enddo
           if (lmsum1 .eq. 0) totcor(:)=-1

! calculated correction factor as 2d field(kpie,kpje)


            corsum(:)=0.
            do j=2,je_g-1
               do i=2,ie_g-1
                  lm=avflux_index(i,j)
                  diffcpl=0.
                  if(lm.gt.0) then                          ! wet point

                     if(totcor(lm).gt.0)  then                 !correction

                        if(avfcpl(lm).gt.0.) then
                           !          if(avfact(lm,1) .eq.0.)               &
                           !     &     write(0,*) 'conflict 1',lm,i,j,avfact(lm,1)   &
                           !     &    ,avfact(lm,2),avfcpl(lm)

                           quotdir=abs(avfact(lm,2)/avfact(lm,1))    ! distribution of flux in opposite direction
                           quottot=abs(avfcpl(lm)/avfact(lm,1) )     ! correction of difference in flux magnitude
                        else
                           !          if(avfact(lm,2) .eq.0.)               &
                           !     &     write(0,*) 'conflict 2',lm,i,j,avfact(lm,1)   &
                           !     &    ,avfact(lm,2),avfcpl(lm)

                           quotdir=abs(avfact(lm,1)/avfact(lm,2))     ! distribution of flux in opposite direction
                           quottot=abs(avfcpl(lm)/avfact(lm,2))       ! correction of difference in flux magnitude
                        end if
                        !         write(0,*) 'quoti ',lm,i,j,quotdir,quottot,avfact(lm,1) &
                        !     &   ,avfact(lm,2),corrfrac(i,j),co2flux_act_g(i,j)

                        newcpl =corrfrac(i,j)*co2flux_act_g(i,j)*(quotdir+quottot) &
                             &       +(1.-corrfrac(i,j))*co2flux_act_g(i,j)
                        diffcpl = co2flux_cpl_g(i,j)*areas_g(i,j)-newcpl*areas_g(i,j)
                        !         if(lm.eq.lmtest)                                    &
                        !     &      write(0,*)' patho yes',i,j,co2flux_cpl_g(i,j),newcpl,      &
                        !     &               co2flux_act_g(i,j),avfcpl(lm),avfact(lm,1),avfact(lm,2)

                        co2flux_cpl_g(i,j)=newcpl
                        corsum(lm)=corsum(lm)+diffcpl
                     else
                        !         if(lm.eq.lmtest)                                    &
                        !     &      write(0,*)' patho no',i,j,co2flux_cpl_g(i,j),co2flux_act_g(i,j)

                     endif
                  endif
               enddo
            enddo



            co2ctrl_new=0.
            do j=2,je_g-1
               do i=2,ie_g-1
                  lm=avflux_index(i,j)
                  if(lm .gt. 0) then
                     co2ctrl_new=co2ctrl_new+co2flux_cpl_g(i,j)*areas_g(i,j)
                     if(abs(corsum(lm)).gt.1.E-6) then
                        write(io_stdo_bgc,*) 'corsum ',lm,i,j,co2flux_cpl_g(i,j)
                     endif
                  endif
               enddo
            enddo

            diffcpl= co2ctrl_new-co2ctrl_old
            write(io_stdo_bgc,*) ' total flux ',co2ctrl_new,co2ctrl_old,            &
                 & co2ctrl_new-co2ctrl_old,lmsum

            !      if(abs(diffcpl).gt. 1.) then
            !         write(0,*) ' mass defect '
            !          do 1152 j=2,je_g-1
            !          do 1152 i=2,ie_g-1
            !          lm=avflux_index(i,j)
            !          if(lm.gt.0 ) then
            !           if (corsum(lm) .gt.1.E-6) then
            !            write(0,*) 'mass at ',i,j,lm
            !            write(0,*) 'mass ',co2flux_cpl_g(i,j),co2flux_act_g(i,j)
            !           end if
            !          end if
            !1152      continue
            !       end if

         END IF  ! end calculation on p_io

         CALL scatter_2d(co2flux_cpl_g,co2flux_cpl,p_io)

         DO j=2,je-1
            DO i=2,ie-1
               IF (ddpo(i,j,1).GT.0.5) THEN
                  thickness = DDPO(i,j,1) + ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
                  ocetra(i,j,1,isco212) = ocetra(i,j,1,isco212) &
                       - ( (co2flux_cpl(i,j)*dtbgc)/(44.011*thickness) )
                  bgcm2d(i,j,jco2fxd)=bgcm2d(i,j,jco2fxd) + co2flux_cpl(i,j)*dtbgc/44.011
               ENDIF  ! end wet cell
            END DO
         END DO

         ! end redistribution of ECHAM fluxes (not for isotopes!)


       END SUBROUTINE AVFLUX_REDIS
#endif /* __cpl_co2 */

     END MODULE MO_AVFLUX
