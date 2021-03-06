MODULE MESSY_SATSIMS_ISCCP


  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: isccp_cloud_types_v3_4

CONTAINS
!-------------------------------------------------------------------------------

  SUBROUTINE isccp_cloud_types_v3_4(&
          debug,&
          debugcol,&
          npoints,&
          sunlit,&
          nlev,&
          ncol,&
          pfull,&
          phalf,&
          qv,&
          cc,&
          conv,&
          dtau_s,&
          dtau_c,&
          top_height,&
          overlap,&
          tautab,&
          invtau,&
          skt,&
          emsfc_lw,&
          at,&
          dem_s,&
          dem_c,&
          fq_isccp,&
          totalcldarea,&
          meanptop,&
          meantaucld &!,&
!          boxtau,&
!          boxptop&
     )
        
!$Id: isccp_cloud_types.f,v 3.4 2003/05/06 08:45:28 hadmw Exp $

! Copyright Steve Klein and Mark Webb 2002 - all rights reserved.
!
! This code is available without charge with the following conditions:
!
!  1. The code is available for scientific purposes and is not for 
!     commercial use.
!  2. Any improvements you make to the code should be made available 
!     to the to the authors for incorporation into a future release.
!  3. The code should not be used in any way that brings the authors 
!     or their employers into disrepute.

    IMPLICIT NONE

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

    INTEGER ncolprint
      
!     -----
!     Input 
!     -----

    INTEGER npoints                   !  number of model points in the horizontal
    INTEGER nlev                      !  number of model levels in column
    INTEGER ncol                      !  number of subcolumns

    INTEGER sunlit(npoints)           !  1 for day points, 0 for night time


    REAL(dp) ::  pfull(npoints,nlev)  !  pressure of full model levels (Pascals)
                                      !  pfull(npoints,1) is top level of model
                                      !  pfull(npoints,nlev) is bottom 
                                      !  level of model

    REAL(dp) ::  phalf(npoints,nlev+1)  !  pressure of half model levels (Pa)
                                        !  phalf(npoints,1) is top of model
                                        !  phalf(npoints,nlev+1) is the 
                                        ! surface pressure

    !  water vapor specific humidity (kg vapor/ kg air) on full model levels
    REAL(dp) ::  qv(npoints,nlev)

    !  input cloud cover in each model level (fraction) 
    !         grid box covered by clouds
    !  NOTE:  This is the HORIZONTAL area of each
    REAL(dp) ::  cc(npoints,nlev)

    !  input convective cloud cover in each model level (fraction) 
    !  NOTE:  This is the HORIZONTAL area of each
    !         grid box covered by convective clouds
    REAL(dp) ::  conv(npoints,nlev)

    !  mean 0.67 micron optical depth of stratiform
    !  clouds in each model level
    !  NOTE:  this the cloud optical depth of only the
    !         cloudy part of the grid box, it is not weighted
    !         with the 0 cloud optical depth of the clear
    !         part of the grid box
    REAL(dp) ::  dtau_s(npoints,nlev)         

    !  mean 0.67 micron optical depth of convective
    !  clouds in each  model level.  Same note applies as in dtau_s.
    REAL(dp) ::  dtau_c(npoints,nlev)  

    INTEGER overlap                   !  overlap type
                                      !  1=max
                                      !  2=rand
                                      !  3=max/rand

    INTEGER top_height    
    !  1 = adjust top height using both a computed
    !  infrared brightness temperature and the visible
    !  optical depth to adjust cloud top pressure. Note
    !  that this calculation is most appropriate to compare
    !  to ISCCP data during sunlit hours.
    !  2 = do not adjust top height, that is cloud top
    !  pressure is the actual cloud top pressure
    !  in the model
    !  3 = adjust top height using only the computed
    !  infrared brightness temperature. Note that this
    !  calculation is most appropriate to compare to ISCCP
    !  IR only algortihm (i.e. you can compare to nighttime
    !  ISCCP data with this option)

    !  ISCCP table for converting count value to optical thickness
    REAL(dp) ::  tautab(0:255)

    !  ISCCP table for converting optical thickness to count value
    INTEGER invtau(-20:45000)
!
!     The following input variables are used only if top_height = 1 
!     or top_height = 3
!
    REAL(dp) ::  skt(npoints)                 !  skin Temperature (K)

    !  10.5 micron emissivity of surface (fraction)
    REAL(dp) ::  emsfc_lw           

    REAL(dp) ::  at(npoints,nlev)     !  temperature in each model level (K)

    !  10.5 micron longwave emissivity of stratiform
    !  clouds in each
    !  model level.  Same note applies as in dtau_s.
    REAL(dp) ::  dem_s(npoints,nlev)

    !  10.5 micron longwave emissivity of convective
    !  clouds in each
    !  model level.  Same note applies as in dtau_s.
    REAL(dp) ::  dem_c(npoints,nlev)   
!     ------
!     Output
!     ------

    !  the fraction of the model grid box covered by
    !  each of the 49 ISCCP D level cloud types
    REAL(dp) ::  fq_isccp(npoints,7,7)  

    !  the fraction of model grid box columns
    !  with cloud somewhere in them.  This should
    !  equal the sum over all entries of fq_isccp
    REAL(dp) ::  totalcldarea(npoints) 
        
        
    ! The following three means are averages over the cloudy areas only.  If no
    ! clouds are in grid box all three quantities should equal zero.  
                                        
    !  mean cloud top pressure (mb) - linear averaging in cloud top pressure.
    REAL(dp) ::  meanptop(npoints)      
                                        
    !  mean optical thickness linear averaging in albedo performed.
    REAL(dp) ::  meantaucld(npoints)  
      
    !  optical thickness in each column
    REAL(dp) ::  boxtau(npoints,ncol)         
      
    !  cloud top pressure (mb) in each column
    REAL(dp) ::  boxptop(npoints,ncol) 

!     ------
!     Working variables added when program updated to mimic 
!     Mark Webb's PV-Wave code
!     ------

    ! boxes gridbox divided up into
    ! Equivalent of BOX in original version, but
    ! indexed by column then row, rather than
    ! by row then column
    REAL(dp) ::  frac_out(npoints,ncol,nlev) 

    ! total cloud cover in each model level (fraction)
    ! with extra layer of zeroes on top
    ! in this version this just contains the values input
    ! from cc but with an extra level
    REAL(dp) ::  tca(npoints,0:nlev) 
    
    ! convective cloud cover in each model level (fraction) from conv 
    REAL(dp) ::  cca(npoints,nlev)      

    REAL(dp) ::  threshold(npoints,ncol)   ! pointer to position in gridbox
    REAL(dp) ::  maxocc(npoints,ncol)      ! Flag for max overlapped conv cld
    REAL(dp) ::  maxosc(npoints,ncol)      ! Flag for max overlapped strat cld
    
    ! ordered pointer to position in gridbox
    REAL(dp) ::  boxpos(npoints,ncol)      

    ! minimum value to define range in with new threshold is chosen
    REAL(dp) ::  threshold_min(npoints,ncol) 

    !  working variables for 10.5 micron longwave 
    !  emissivity in part of gridbox under consideration
    REAL(dp) ::  dem(npoints,ncol),bb(npoints) 

    REAL(dp) ::  ran(npoints)                 ! vector of random numbers
    REAL(dp) ::  ptrop(npoints)
    REAL(dp) ::  attrop(npoints)
    REAL(dp) ::  attropmin (npoints)
    REAL(dp) ::  atmax(npoints)
    REAL(dp) ::  atmin(npoints)
    REAL(dp) ::  btcmin(npoints)
    REAL(dp) ::  transmax(npoints)

    INTEGER i,j,ilev,ibox,itrop(npoints)
    INTEGER ipres(npoints)
    INTEGER itau(npoints),ilev2
    INTEGER acc(nlev,ncol)
    INTEGER match(npoints,nlev-1)
    INTEGER nmatch(npoints)
    INTEGER levmatch(npoints,ncol)
      
      !variables needed for water vapor continuum absorption
    REAL(dp) :: fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
    REAL(dp) :: taumin(npoints)
    REAL(dp) :: dem_wv(npoints,nlev), wtmair, wtmh20, Navo, grav, pstd, t0
    REAL(dp) :: press(npoints), dpress(npoints), atmden(npoints)
    REAL(dp) :: rvh20(npoints), wk(npoints), rhoave(npoints)
    REAL(dp) :: rh20s(npoints), rfrgn(npoints)
    REAL(dp) :: tmpexp(npoints),tauwv(npoints)
      
    CHARACTER*1 cchar(6),cchar_realtops(6)
    INTEGER icycle
    REAL(dp) ::  tau(npoints,ncol)
    LOGICAL box_cloudy(npoints,ncol)
    REAL(dp) ::  tb(npoints,ncol)
    REAL(dp) ::  ptop(npoints,ncol)
    REAL(dp) ::  emcld(npoints,ncol)
    REAL(dp) ::  fluxtop(npoints,ncol)
    REAL(dp) ::  trans_layers_above(npoints,ncol)
    REAL(dp) :: isccp_taumin,fluxtopinit(npoints),tauir(npoints)
    REAL(dp) :: meanalbedocld(npoints) 
    REAL(dp) ::  albedocld(npoints,ncol)
    REAL(dp) :: boxarea
    INTEGER debug       ! set to non-zero value to print out inputs
                          ! with step debug
    INTEGER debugcol    ! set to non-zero value to print out column
                          ! decomposition with step debugcol

    INTEGER index1(npoints),num1,jj
    REAL rec2p13,tauchk

    CHARACTER*10 ftn09
      
    DATA isccp_taumin / 0.3_dp /
    DATA cchar / ' ','-','1','+','I','+'/
    DATA cchar_realtops / ' ',' ','1','1','I','I'/
    
    tauchk = -1._dp*LOG(0.9999999_dp)
    rec2p13=1._dp/2.13_dp
    
    ncolprint=0

    IF ( debug.NE.0 ) THEN
      j=1
      WRITE(6,'(a10)') 'j='
      WRITE(6,'(8I10)') j
      WRITE(6,'(a10)') 'debug='
      WRITE(6,'(8I10)') debug
      WRITE(6,'(a10)') 'debugcol='
      WRITE(6,'(8I10)') debugcol
      WRITE(6,'(a10)') 'npoints='
      WRITE(6,'(8I10)') npoints
      WRITE(6,'(a10)') 'nlev='
      WRITE(6,'(8I10)') nlev
      WRITE(6,'(a10)') 'ncol='
      WRITE(6,'(8I10)') ncol
      WRITE(6,'(a10)') 'top_height='
      WRITE(6,'(8I10)') top_height
      WRITE(6,'(a10)') 'overlap='
      WRITE(6,'(8I10)') overlap
      WRITE(6,'(a10)') 'emsfc_lw='
      WRITE(6,'(8f10.2)') emsfc_lw
      WRITE(6,'(a10)') 'tautab='
      WRITE(6,'(8f10.2)') tautab
      WRITE(6,'(a10)') 'invtau(1:100)='
      WRITE(6,'(8i10)') (invtau(i),i=1,100)
      DO j=1,npoints,debug
        WRITE(6,'(a10)') 'j='
        WRITE(6,'(8I10)') j
        WRITE(6,'(a10)') 'sunlit='
        WRITE(6,'(8I10)') sunlit(j)
        WRITE(6,'(a10)') 'pfull='
        WRITE(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'phalf='
        WRITE(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
        WRITE(6,'(a10)') 'qv='
        WRITE(6,'(8f10.3)') (qv(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'cc='
        WRITE(6,'(8f10.3)') (cc(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'conv='
        WRITE(6,'(8f10.2)') (conv(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'dtau_s='
        WRITE(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'dtau_c='
        WRITE(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'skt='
        WRITE(6,'(8f10.2)') skt(j)
        WRITE(6,'(a10)') 'at='
        WRITE(6,'(8f10.2)') (at(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'dem_s='
        WRITE(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
        WRITE(6,'(a10)') 'dem_c='
        WRITE(6,'(8f10.2)') (dem_c(j,i),i=1,nlev)
      ENDDO
    ENDIF
    
    !     ---------------------------------------------------!
    
    !     assign 2d tca array using 1d input array cc
    
    DO j=1,npoints
      tca(j,0)=0._dp
    ENDDO
    
    DO ilev=1,nlev
      DO j=1,npoints
        tca(j,ilev)=cc(j,ilev)
      ENDDO
    ENDDO
    
    !     assign 2d cca array using 1d input array conv
    
    DO ilev=1,nlev
      !jq        do ibox=1,ncol
      DO j=1,npoints
        cca(j,ilev)=conv(j,ilev)
      ENDDO
      !jq        enddo
    ENDDO
    
    IF (ncolprint.NE.0) THEN
      DO j=1,npoints,1000
        WRITE(6,'(a10)') 'j='
        WRITE(6,'(8I10)') j
        WRITE (6,'(a)') 'tca_pp_rev:'
        WRITE (6,'(8f5.2)') &
          ((tca(j,ilev)),&
          ilev=1,nlev)
        
        WRITE (6,'(a)') 'cca_pp_rev:'
        WRITE (6,'(8f5.2)') &
          ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
      ENDDO
    ENDIF
    
    IF (top_height .EQ. 1 .OR. top_height .EQ. 3) THEN 
      
      DO j=1,npoints 
        ptrop(j)=5000._dp
        atmin(j) = 400._dp
        attropmin(j) = 400._dp
        atmax(j) = 0._dp
        attrop(j) = 120._dp
        itrop(j) = 1
      ENDDO
      
      DO 12 ilev=1,nlev
        DO j=1,npoints 
          IF (pfull(j,ilev) .LT. 40000._dp .AND.&
            pfull(j,ilev) .GT.  5000._dp .AND.&
            at(j,ilev) .LT. attropmin(j)) THEN
            ptrop(j) = pfull(j,ilev)
            attropmin(j) = at(j,ilev)
            attrop(j) = attropmin(j)
            itrop(j)=ilev
          END IF
          IF (at(j,ilev) .GT. atmax(j)) atmax(j)=at(j,ilev)
          IF (at(j,ilev) .LT. atmin(j)) atmin(j)=at(j,ilev)
        ENDDO
12    ENDDO
        
      END IF
      
      !     -----------------------------------------------------!
      
      !     ---------------------------------------------------!
      
      DO 13 ilev=1,nlev
        num1=0
        DO j=1,npoints
          IF (cc(j,ilev) .LT. 0._dp .OR. cc(j,ilev) .GT. 1._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  ' error = cloud fraction less than zero'
          WRITE(6,*) ' or '
          WRITE(6,*)  ' error = cloud fraction greater than 1'
          WRITE(6,*) 'value at point ',j,' is ',cc(j,ilev)
          WRITE(6,*) 'level ',ilev
          STOP
        ENDDO
        num1=0
        DO j=1,npoints
          IF (conv(j,ilev) .LT. 0._dp .OR. conv(j,ilev) .GT. 1._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  &
            ' error = convective cloud fraction less than zero'
          WRITE(6,*) ' or '
          WRITE(6,*)  &
            ' error = convective cloud fraction greater than 1'
          WRITE(6,*) 'value at point ',j,' is ',conv(j,ilev)
          WRITE(6,*) 'level ',ilev
          STOP
        ENDDO
        
        num1=0
        DO j=1,npoints
          IF (dtau_s(j,ilev) .LT. 0._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  &
            ' error = stratiform cloud opt. depth less than zero'
          WRITE(6,*) 'value at point ',j,' is ',dtau_s(j,ilev)
          WRITE(6,*) 'level ',ilev
          STOP
        ENDDO
        num1=0
        DO j=1,npoints
          IF (dtau_c(j,ilev) .LT. 0._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  &
            ' error = convective cloud opt. depth less than zero'
          WRITE(6,*) 'value at point ',j,' is ',dtau_c(j,ilev)
          WRITE(6,*) 'level ',ilev
          STOP
        ENDDO
        
        num1=0
        DO j=1,npoints
          IF (dem_s(j,ilev) .LT. 0._dp .OR. dem_s(j,ilev) .GT. 1._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  &
            ' error = stratiform cloud emissivity less than zero'
          WRITE(6,*)'or'
          WRITE(6,*)  &
            ' error = stratiform cloud emissivity greater than 1'
          WRITE(6,*) 'value at point ',j,' is ',dem_s(j,ilev)
          WRITE(6,*) 'level ',ilev
          STOP
        ENDDO
        
        num1=0
        DO j=1,npoints
          IF (dem_c(j,ilev) .LT. 0._dp .OR. dem_c(j,ilev) .GT. 1._dp) THEN
            num1=num1+1
            index1(num1)=j
          END IF
        ENDDO
        DO jj=1,num1   
          j=index1(jj)
          WRITE(6,*)  &
            ' error = convective cloud emissivity less than zero'
          WRITE(6,*)'or'
          WRITE(6,*)  &
            ' error = convective cloud emissivity greater than 1'
          WRITE (6,*) &
            'j=',j,'ilev=',ilev,'dem_c(j,ilev) =',dem_c(j,ilev) 
          STOP
        ENDDO
13    ENDDO
        
        
      DO ibox=1,ncol
        DO j=1,npoints 
          boxpos(j,ibox)=(ibox-.5_dp)/ncol
        ENDDO
      ENDDO
      
        !     ---------------------------------------------------!
        !     Initialise working variables
        !     ---------------------------------------------------!
        
        !     Initialised frac_out to zero
        
      DO ibox=1,ncol
        DO ilev=1,nlev
          DO j=1,npoints
            frac_out(j,ibox,ilev)=0.0_dp
          ENDDO
        ENDDO
      ENDDO

      IF (ncolprint.NE.0) THEN
        WRITE (6,'(a)') 'frac_out_pp_rev:'
        DO j=1,npoints,1000
          WRITE(6,'(a10)') 'j='
          WRITE(6,'(8I10)') j
          WRITE (6,'(8f5.2)') & 
            ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)
          
        ENDDO
        WRITE (6,'(a)') 'ncol:'
        WRITE (6,'(I3)') ncol
      ENDIF
      IF (ncolprint.NE.0) THEN
        WRITE (6,'(a)') 'last_frac_pp:'
        DO j=1,npoints,1000
          WRITE(6,'(a10)') 'j='
          WRITE(6,'(8I10)') j
          WRITE (6,'(8f5.2)') (tca(j,0))
        ENDDO
      ENDIF
      
      !     ---------------------------------------------------!
      !     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
      !     frac_out is the array that contains the information 
      !     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
      !     convective cloud
      
      !loop over vertical levels
      DO 200 ilev = 1,nlev
        
        !     Initialise threshold
        
        IF (ilev.EQ.1) THEN
          ! If max overlap 
          IF (overlap.EQ.1) THEN
            ! select pixels spread evenly
            ! across the gridbox
            DO ibox=1,ncol
              DO j=1,npoints
                threshold(j,ibox)=boxpos(j,ibox)
              ENDDO
            ENDDO
          ELSE
            DO ibox=1,ncol
              CALL random_uniform(ran)
              ! select random pixels from the non-convective
              ! part the gridbox ( some will be converted into
              ! convective pixels below )
              DO j=1,npoints
                threshold(j,ibox)= &
                  cca(j,ilev)+(1-cca(j,ilev))*ran(j)
              ENDDO
            ENDDO
          ENDIF
          IF (ncolprint.NE.0) THEN
            WRITE (6,'(a)') 'threshold_nsf2:'
            DO j=1,npoints,1000
              WRITE(6,'(a10)') 'j='
              WRITE(6,'(8I10)') j
              WRITE (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
            ENDDO
          ENDIF
        ENDIF

        IF (ncolprint.NE.0) THEN
          WRITE (6,'(a)') 'ilev:'
          WRITE (6,'(I2)') ilev
        ENDIF

        DO ibox=1,ncol

          ! All versions
          DO j=1,npoints
            IF (boxpos(j,ibox).LE.cca(j,ilev)) THEN
              maxocc(j,ibox) = 1
            ELSE
              maxocc(j,ibox) = 0
            END IF
          ENDDO

          ! Max overlap
          IF (overlap.EQ.1) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=1
            ENDDO
          ENDIF
          
          ! Random overlap
          IF (overlap.EQ.2) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=cca(j,ilev)
              maxosc(j,ibox)=0
            ENDDO
          ENDIF
          
          ! Max/Random overlap
          IF (overlap.EQ.3) THEN 
            DO j=1,npoints
              threshold_min(j,ibox)=MAX(cca(j,ilev), &
                MIN(tca(j,ilev-1),tca(j,ilev)))
              IF (threshold(j,ibox) &
                .LT.MIN(tca(j,ilev-1),tca(j,ilev)) &
                .AND.(threshold(j,ibox).GT.cca(j,ilev))) THEN
                maxosc(j,ibox)= 1
              ELSE
                maxosc(j,ibox)= 0
              END IF
            ENDDO
          ENDIF
          
          ! Reset threshold 
          
          CALL random_uniform(ran)
          
          DO j=1,npoints
            threshold(j,ibox)= &
              !if max overlapped conv cloud
              maxocc(j,ibox) * ( &                                       
              boxpos(j,ibox)  &                                             
              ) +                  &                                    
              !else
              (1-maxocc(j,ibox)) * ( &                                   
              !if max overlapped strat cloud
              (maxosc(j,ibox)) * ( &                                 
              !threshold=boxpos 
              threshold(j,ibox) &                                       
              ) +                    &                              
              !else
              (1-maxosc(j,ibox)) * (  &                             
              !threshold_min=random[thrmin,1]
              threshold_min(j,ibox)+&
              (1-threshold_min(j,ibox))*ran(j)  &
              ) &
              )
          ENDDO
          
        ENDDO ! ibox

        !       Fill frac_out with 1's where tca is greater than the threshold
        
        DO ibox=1,ncol
          DO j=1,npoints 
            IF (tca(j,ilev).GT.threshold(j,ibox)) THEN
              frac_out(j,ibox,ilev)=1._dp
            ELSE
              frac_out(j,ibox,ilev)=0._dp
            END IF
          ENDDO
        ENDDO
        
        !          Code to partition boxes into startiform and convective parts
        !          goes here
        
        DO ibox=1,ncol
          DO j=1,npoints 
            IF (threshold(j,ibox).LE.cca(j,ilev)) THEN
              ! = 2 IF threshold le cca(j)
              frac_out(j,ibox,ilev) = 2 
            ELSE
              ! = the same IF NOT threshold le cca(j) 
              frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
            END IF
          ENDDO
        ENDDO
        
        !         Set last_frac to tca at this level, so as to be tca 
        !         from last level next time round
        
        IF (ncolprint.NE.0) THEN
          
          DO j=1,npoints ,1000
            WRITE(6,'(a10)') 'j='
            WRITE(6,'(8I10)') j
            WRITE (6,'(a)') 'last_frac:'
            WRITE (6,'(8f5.2)') (tca(j,ilev-1))
            
            WRITE (6,'(a)') 'cca:'
            WRITE (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'max_overlap_cc:'
            WRITE (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'max_overlap_sc:'
            WRITE (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'threshold_min_nsf2:'
            WRITE (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'threshold_nsf2:'
            WRITE (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'frac_out_pp_rev:'
            WRITE (6,'(8f5.2)') &
              ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
          ENDDO
        ENDIF
        
200   ENDDO    !loop over nlev
        
!
!     ---------------------------------------------------!

      
!
!     ---------------------------------------------------!
!     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
!     put into vector tau
 
      !initialize tau and albedocld to zero
      DO 15 ibox=1,ncol
        DO j=1,npoints 
          tau(j,ibox)=0._dp
          albedocld(j,ibox)=0._dp
          boxtau(j,ibox)=0._dp
          boxptop(j,ibox)=0._dp
          box_cloudy(j,ibox)=.FALSE.
        ENDDO
15    ENDDO
        
      !compute total cloud optical depth for each column     
      DO ilev=1,nlev
        !increment tau for each of the boxes
        DO ibox=1,ncol
          DO j=1,npoints 
            IF (frac_out(j,ibox,ilev).EQ.1) THEN
              tau(j,ibox)=tau(j,ibox) &
                + dtau_s(j,ilev)
            ENDIF
            IF (frac_out(j,ibox,ilev).EQ.2) THEN
              tau(j,ibox)=tau(j,ibox) &
                + dtau_c(j,ilev)
            END IF
          ENDDO
        ENDDO ! ibox
      ENDDO ! ilev
      IF (ncolprint.NE.0) THEN
        
        DO j=1,npoints ,1000
          WRITE(6,'(a10)') 'j='
          WRITE(6,'(8I10)') j
          WRITE(6,'(i2,1X,8(f7.2,1X))') &
            ilev, &
            (tau(j,ibox),ibox=1,ncolprint)
        ENDDO
      ENDIF
!
!     ---------------------------------------------------!



!     
!     ---------------------------------------------------!
!     COMPUTE INFRARED BRIGHTNESS TEMPERUATRES
!     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
!
!     again this is only done if top_height = 1 or 3
!
!     fluxtop is the 10.5 micron radiance at the top of the
!              atmosphere
!     trans_layers_above is the total transmissivity in the layers
!             above the current layer
!     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
!             sky versions of these quantities.

      IF (top_height .EQ. 1 .OR. top_height .EQ. 3) THEN


        !----------------------------------------------------------------------
        !    
        !             DO CLEAR SKY RADIANCE CALCULATION FIRST
        !
        !compute water vapor continuum emissivity
        !this treatment follows Schwarkzopf and Ramasamy
        !JGR 1999,vol 104, pages 9467-9499.
        !the emissivity is calculated at a wavenumber of 955 cm-1, 
        !or 10.47 microns 
        wtmair = 28.9644_dp
        wtmh20 = 18.01534_dp
        Navo = 6.023E+23_dp
        grav = 9.806650E+02_dp
        pstd = 1.013250E+06_dp
        t0 = 296._dp
        IF (ncolprint .NE. 0) &
          WRITE(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
        DO 125 ilev=1,nlev
          DO j=1,npoints 
            !press and dpress are dyne/cm2 = Pascals *10
            press(j) = pfull(j,ilev)*10.
            dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
            !atmden = g/cm2 = kg/m2 / 10 
            atmden(j) = dpress(j)/grav
            rvh20(j) = qv(j,ilev)*wtmair/wtmh20
            wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
            rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
            rh20s(j) = rvh20(j)*rhoave(j)
            rfrgn(j) = rhoave(j)-rh20s(j)
            tmpexp(j) = EXP(-0.02_dp*(at(j,ilev)-t0))
            tauwv(j) = wk(j)*1.e-20_dp*( &
              (0.0224697_dp*rh20s(j)*tmpexp(j)) + &
              (3.41817e-7_dp*rfrgn(j)) )*0.98_dp
            dem_wv(j,ilev) = 1._dp - EXP( -1._dp * tauwv(j))
          ENDDO
          IF (ncolprint .NE. 0) THEN
            DO j=1,npoints ,1000
              WRITE(6,'(a10)') 'j='
              WRITE(6,'(8I10)') j
              WRITE(6,'(i2,1X,3(f8.3,3X))') ilev,&
                qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))/(grav/100._dp),&
                tauwv(j),dem_wv(j,ilev)
            ENDDO
          ENDIF
125     ENDDO
          
        !initialize variables
        DO j=1,npoints 
          fluxtop_clrsky(j) = 0._dp
          trans_layers_above_clrsky(j)=1._dp
        ENDDO
        
        DO ilev=1,nlev
          DO j=1,npoints 
            
            ! Black body emission at temperature of the layer
            
            bb(j)=1 / ( EXP(1307.27_dp/at(j,ilev)) - 1._dp )
            !bb(j)= 5.67e-8*at(j,ilev)**4
            
            ! increase TOA flux by flux emitted from layer
            ! times total transmittance in layers above
            
            fluxtop_clrsky(j) = fluxtop_clrsky(j) &
              + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j) 
            
            ! update trans_layers_above with transmissivity
            ! from this layer for next time around loop
            
            trans_layers_above_clrsky(j)= &
              trans_layers_above_clrsky(j)*(1._dp-dem_wv(j,ilev))
            
            
          ENDDO
          IF (ncolprint.NE.0) THEN
            DO j=1,npoints ,1000
              WRITE(6,'(a10)') 'j='
              WRITE(6,'(8I10)') j
              WRITE (6,'(a)') 'ilev:'
              WRITE (6,'(I2)') ilev
              
              WRITE (6,'(a)') &
                'emiss_layer,100.*bb(j),100.*f,total_trans:'
              WRITE (6,'(4(f7.2,1X))') dem_wv(j,ilev),100._dp*bb(j),&
                100._dp*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
            ENDDO
          ENDIF
          
        ENDDO   !loop over level
        
        DO j=1,npoints 
          !add in surface emission
          bb(j)=1/( EXP(1307.27_dp/skt(j)) - 1. )
          !bb(j)=5.67e-8*skt(j)**4
          
          fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) & 
            * trans_layers_above_clrsky(j)
        ENDDO
        
        IF (ncolprint.NE.0) THEN
          DO j=1,npoints ,1000
            WRITE(6,'(a10)') 'j='
            WRITE(6,'(8I10)') j
            WRITE (6,'(a)') 'id:'
            WRITE (6,'(a)') 'surface'
            
            WRITE (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
            WRITE (6,'(4(f7.2,1X))') emsfc_lw,100._dp*bb(j), &
              100._dp*fluxtop_clrsky(j), &
              trans_layers_above_clrsky(j)
          ENDDO
        ENDIF
        
        
        !
        !           END OF CLEAR SKY CALCULATION
        !
        !----------------------------------------------------------------



        IF (ncolprint.NE.0) THEN

          DO j=1,npoints ,1000
            WRITE(6,'(a10)') 'j='
            WRITE(6,'(8I10)') j
            WRITE (6,'(a)') 'ts:'
            WRITE (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'ta_rev:'
            WRITE (6,'(8f7.2)') &
              ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
            
          ENDDO
        ENDIF
        !loop over columns 
        DO ibox=1,ncol
          DO j=1,npoints
            fluxtop(j,ibox)=0.
            trans_layers_above(j,ibox)=1.
          ENDDO
        ENDDO
        
        DO ilev=1,nlev
          DO j=1,npoints 
            ! Black body emission at temperature of the layer
            
            bb(j)=1 / ( EXP(1307.27_dp/at(j,ilev)) - 1. )
            !bb(j)= 5.67e-8*at(j,ilev)**4
          ENDDO
          
          DO ibox=1,ncol
            DO j=1,npoints 
              
              ! emissivity for point in this layer
              IF (frac_out(j,ibox,ilev).EQ.1) THEN
                dem(j,ibox)= 1._dp - & 
                  ( (1._dp - dem_wv(j,ilev)) * (1._dp -  dem_s(j,ilev)) )
              ELSE IF (frac_out(j,ibox,ilev).EQ.2) THEN
                dem(j,ibox)= 1._dp - &
                  ( (1._dp - dem_wv(j,ilev)) * (1._dp -  dem_c(j,ilev)) )
              ELSE
                dem(j,ibox)=  dem_wv(j,ilev)
              END IF
              
              
              ! increase TOA flux by flux emitted from layer
              ! times total transmittance in layers above
              
              fluxtop(j,ibox) = fluxtop(j,ibox) &
                + dem(j,ibox) * bb(j) &
                * trans_layers_above(j,ibox) 
              
              ! update trans_layers_above with transmissivity
              ! from this layer for next time around loop
              
              trans_layers_above(j,ibox)= &
                trans_layers_above(j,ibox)*(1._dp-dem(j,ibox))
              
            ENDDO ! j
          ENDDO ! ibox
          
          IF (ncolprint.NE.0) THEN
            DO j=1,npoints,1000
              WRITE (6,'(a)') 'ilev:'
              WRITE (6,'(I2)') ilev
              
              WRITE(6,'(a10)') 'j='
              WRITE(6,'(8I10)') j
              WRITE (6,'(a)') 'emiss_layer:'
              WRITE (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
              
              WRITE (6,'(a)') '100.*bb(j):'
              WRITE (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
              
              WRITE (6,'(a)') '100.*f:'
              WRITE (6,'(8f7.2)') &
                (100.*fluxtop(j,ibox),ibox=1,ncolprint)
              
              WRITE (6,'(a)') 'total_trans:'
              WRITE (6,'(8f7.2)') &
                (trans_layers_above(j,ibox),ibox=1,ncolprint)
            ENDDO
          ENDIF
          
        ENDDO ! ilev
        

        DO j=1,npoints 
          !add in surface emission
          bb(j)=1/( EXP(1307.27_dp/skt(j)) - 1._dp )
          !bb(j)=5.67e-8*skt(j)**4
        END DO
        
        DO ibox=1,ncol
          DO j=1,npoints 
            
            !add in surface emission
            
            fluxtop(j,ibox) = fluxtop(j,ibox)  &
              + emsfc_lw * bb(j) & 
              * trans_layers_above(j,ibox) 
            
          END DO
        END DO
        
        IF (ncolprint.NE.0) THEN
          
          DO j=1,npoints ,1000
            WRITE(6,'(a10)') 'j='
            WRITE(6,'(8I10)') j
            WRITE (6,'(a)') 'id:'
            WRITE (6,'(a)') 'surface'
            
            WRITE (6,'(a)') 'emiss_layer:'
            WRITE (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') '100.*bb(j):'
            WRITE (6,'(8f7.2)') (100._dp*bb(j),ibox=1,ncolprint)
            
            WRITE (6,'(a)') '100.*f:'
            WRITE (6,'(8f7.2)') (100._dp*fluxtop(j,ibox),ibox=1,ncolprint)
          END DO
        ENDIF
        
        !now that you have the top of atmosphere radiance account
        !for ISCCP procedures to determine cloud top temperature
        
        !account for partially transmitting cloud recompute flux 
        !ISCCP would see assuming a single layer cloud
        !note choice here of 2.13, as it is primarily ice
        !clouds which have partial emissivity and need the 
        !adjustment performed in this section
        !
        !If it turns out that the cloud brightness temperature
        !is greater than 260K, then the liquid cloud conversion
        !factor of 2.56 is used.
        !
        !Note that this is discussed on pages 85-87 of 
        !the ISCCP D level documentation (Rossow et al. 1996)
           
        DO j=1,npoints  
          !compute minimum brightness temperature and optical depth
          btcmin(j) = 1._dp /  ( EXP(1307.27_dp/(attrop(j)-5._dp)) - 1._dp ) 
        ENDDO
        DO ibox=1,ncol
          DO j=1,npoints  
            transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
              /(fluxtop_clrsky(j)-btcmin(j))
            !note that the initial setting of tauir(j) is needed so that
            !tauir(j) has a realistic value should the next if block be
            !bypassed
            tauir(j) = tau(j,ibox) * rec2p13
            taumin(j) = -1._dp*LOG(MAX(MIN(transmax(j),0.9999999_dp),0.001_dp))
            
          ENDDO
          
          IF (top_height .EQ. 1) THEN
            DO j=1,npoints  
              IF (transmax(j) .GT. 0.001_dp .AND. &
                transmax(j) .LE. 0.9999999_dp) THEN
                fluxtopinit(j) = fluxtop(j,ibox)
                tauir(j) = tau(j,ibox) *rec2p13
              ENDIF
            ENDDO
            DO icycle=1,2
              DO j=1,npoints  
                IF (tau(j,ibox) .GT. (tauchk            )) THEN 
                  IF (transmax(j) .GT. 0.001_dp .AND. & 
                    transmax(j) .LE. 0.9999999_dp) THEN
                    emcld(j,ibox) = 1._dp - EXP(-1._dp * tauir(j)  )
                    fluxtop(j,ibox) = fluxtopinit(j) -   &
                      ((1._dp-emcld(j,ibox))*fluxtop_clrsky(j))
                    fluxtop(j,ibox)=MAX(1.E-06_dp, &
                      (fluxtop(j,ibox)/emcld(j,ibox)))
                    tb(j,ibox)= 1307.27_dp &
                      / (LOG(1._dp + (1._dp/fluxtop(j,ibox))))
                    IF (tb(j,ibox) .GT. 260._dp) THEN
                      tauir(j) = tau(j,ibox) / 2.56_dp
                    END IF
                  END IF
                END IF
              ENDDO
            ENDDO
            
          ENDIF
          
          DO j=1,npoints
            IF (tau(j,ibox) .GT. (tauchk            )) THEN 
              !cloudy box
              tb(j,ibox)= 1307.27_dp/ (LOG(1._dp + (1._dp/fluxtop(j,ibox))))
              IF (top_height.EQ.1.AND.tauir(j).LT.taumin(j)) THEN
                tb(j,ibox) = attrop(j) - 5. 
                tau(j,ibox) = 2.13_dp*taumin(j)
              END IF
            ELSE
              !clear sky brightness temperature
              tb(j,ibox) = 1307.27_dp/(LOG(1._dp+(1._dp/fluxtop_clrsky(j))))
            END IF
          ENDDO ! j
        ENDDO ! ibox
        
        IF (ncolprint.NE.0) THEN
          
          DO j=1,npoints,1000
            WRITE(6,'(a10)') 'j='
            WRITE(6,'(8I10)') j
            
            WRITE (6,'(a)') 'attrop:'
            WRITE (6,'(8f7.2)') (attrop(j))
            
            WRITE (6,'(a)') 'btcmin:'
            WRITE (6,'(8f7.2)') (btcmin(j))
            
            WRITE (6,'(a)') 'fluxtop_clrsky*100:'
            WRITE (6,'(8f7.2)') & 
              (100.*fluxtop_clrsky(j))
            
            WRITE (6,'(a)') '100.*f_adj:'
            WRITE (6,'(8f7.2)') (100._dp*fluxtop(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'transmax:'
            WRITE (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'tau:'
            WRITE (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'emcld:'
            WRITE (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'total_trans:'
            WRITE (6,'(8f7.2)') &
              (trans_layers_above(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'total_emiss:'
            WRITE (6,'(8f7.2)') &
              (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'total_trans:'
            WRITE (6,'(8f7.2)') &
              (trans_layers_above(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'ppout:'
            WRITE (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
          ENDDO ! j
        ENDIF
        
      END IF

!     ---------------------------------------------------!

!     
!     ---------------------------------------------------!
!     DETERMINE CLOUD TOP PRESSURE
!
!     again the 2 methods differ according to whether
!     or not you use the physical cloud top pressure (top_height = 2)
!     or the radiatively determined cloud top pressure (top_height = 1 or 3)
!

      !compute cloud top pressure
      DO 30 ibox=1,ncol
        !segregate according to optical thickness
        IF (top_height .EQ. 1 .OR. top_height .EQ. 3) THEN  
          !find level whose temperature
          !most closely matches brightness temperature
          DO j=1,npoints 
            nmatch(j)=0
          ENDDO
          DO 29 ilev=1,nlev-1
            !cdir nodep
            DO j=1,npoints 
              IF ((at(j,ilev)   .GE. tb(j,ibox) .AND. &
                at(j,ilev+1) .LT. tb(j,ibox)) .OR. &
                (at(j,ilev) .LE. tb(j,ibox) .AND. &
                at(j,ilev+1) .GT. tb(j,ibox))) THEN 
                
                nmatch(j)=nmatch(j)+1
                IF(ABS(at(j,ilev)-tb(j,ibox)) .LT. &
                  ABS(at(j,ilev+1)-tb(j,ibox))) THEN
                  match(j,nmatch(j))=ilev
                ELSE
                  match(j,nmatch(j))=ilev+1
                END IF
              END IF
            ENDDO
29        END DO
            
          DO j=1,npoints 
            IF (nmatch(j) .GE. 1) THEN
              ptop(j,ibox)=pfull(j,match(j,nmatch(j)))
              levmatch(j,ibox)=match(j,nmatch(j))   
            ELSE
              IF (tb(j,ibox) .LT. atmin(j)) THEN
                ptop(j,ibox)=ptrop(j)
                levmatch(j,ibox)=itrop(j)
              END IF
              IF (tb(j,ibox) .GT. atmax(j)) THEN
                ptop(j,ibox)=pfull(j,nlev)
                levmatch(j,ibox)=nlev
              END IF
            END IF
          ENDDO ! j
          
        ELSE ! if (top_height .eq. 1 .or. top_height .eq. 3) 
          
          DO j=1,npoints     
            ptop(j,ibox)=0._dp
          ENDDO
          DO ilev=1,nlev
            DO j=1,npoints     
              IF ((ptop(j,ibox) .EQ. 0._dp ) &
                .AND.(frac_out(j,ibox,ilev) .NE. 0)) THEN
                ptop(j,ibox)=pfull(j,ilev)
                levmatch(j,ibox)=ilev
              END IF
            END DO
          END DO
        END IF
        
        DO j=1,npoints
          IF (tau(j,ibox) .LE. (tauchk            )) THEN
            ptop(j,ibox)=0._dp
            levmatch(j,ibox)=0      
          ENDIF
        ENDDO
        
30    ENDDO
      
!
!
!     ---------------------------------------------------!


!     
!     ---------------------------------------------------!
!     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
!
!     Now that ptop and tau have been determined, 
!     determine amount of each of the 49 ISCCP cloud
!     types
!
!     Also compute grid box mean cloud top pressure and
!     optical thickness.  The mean cloud top pressure and
!     optical thickness are averages over the cloudy 
!     area only. The mean cloud top pressure is a linear
!     average of the cloud top pressures.  The mean cloud
!     optical thickness is computed by converting optical
!     thickness to an albedo, averaging in albedo units,
!     then converting the average albedo back to a mean
!     optical thickness.  
!

      !compute isccp frequencies

      !reset frequencies
      DO ilev=1,7
        DO ilev2=1,7
          DO j=1,npoints ! 
            fq_isccp(j,ilev,ilev2)=0._dp
          ENDDO
        ENDDO
      ENDDO

      !reset variables need for averaging cloud properties
      DO j=1,npoints 
        totalcldarea(j) = 0._dp
        meanalbedocld(j) = 0._dp
        meanptop(j) = 0._dp
        meantaucld(j) = 0._dp
      ENDDO ! j

      boxarea = 1./REAL(ncol)
     
      DO 39 ibox=1,ncol
        DO j=1,npoints 

          IF (tau(j,ibox) .GT. (tauchk            ) &
            .AND. ptop(j,ibox) .GT. 0._dp) THEN
            box_cloudy(j,ibox)=.TRUE.
          ENDIF
          
          IF (box_cloudy(j,ibox)) THEN
            
              ! totalcldarea always diagnosed day or night
            totalcldarea(j) = totalcldarea(j) + boxarea
            
            IF (sunlit(j).EQ.1) THEN
              
              ! tau diagnostics only with sunlight
              
              boxtau(j,ibox) = tau(j,ibox)

              !convert optical thickness to albedo
              albedocld(j,ibox) &
!!$                =REAL(invtau(MIN(NINT(100._dp*tau(j,ibox)),45000)))
                =REAL(invtau(MAX(-20,MIN(NINT(100._dp*tau(j,ibox)),45000))))
              
              !contribute to averaging
              meanalbedocld(j) = meanalbedocld(j)  &
                +albedocld(j,ibox)*boxarea
              
            ENDIF
            
          ENDIF
          
          IF (sunlit(j).EQ.1 .OR. top_height .EQ. 3) THEN 
            
            !convert ptop to millibars
            ptop(j,ibox)=ptop(j,ibox) / 100._dp
            
            !save for output cloud top pressure and optical thickness
            boxptop(j,ibox) = ptop(j,ibox)
            
            IF (box_cloudy(j,ibox)) THEN
              
              meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea
              
              !reset itau(j), ipres(j)
              itau(j) = 0
              ipres(j) = 0
              
              !determine optical depth category
              IF (tau(j,ibox) .LT. isccp_taumin) THEN
                itau(j)=1
              ELSE IF (tau(j,ibox) .GE. isccp_taumin &
                .AND. tau(j,ibox) .LT. 1.3_dp) THEN
                itau(j)=2
              ELSE IF (tau(j,ibox) .GE. 1.3_dp & 
                .AND. tau(j,ibox) .LT. 3.6_dp) THEN
                itau(j)=3
              ELSE IF (tau(j,ibox) .GE. 3.6_dp &
                .AND. tau(j,ibox) .LT. 9.4_dp) THEN
                itau(j)=4
              ELSE IF (tau(j,ibox) .GE. 9.4_dp &
                .AND. tau(j,ibox) .LT. 23._dp) THEN
                itau(j)=5
              ELSE IF (tau(j,ibox) .GE. 23._dp &
                .AND. tau(j,ibox) .LT. 60._dp) THEN
                itau(j)=6
              ELSE IF (tau(j,ibox) .GE. 60._dp) THEN
                itau(j)=7
              END IF
              
              !determine cloud top pressure category
              IF (    ptop(j,ibox) .GT. 0._dp  &
                .AND.ptop(j,ibox) .LT. 180._dp) THEN
                ipres(j)=1
              ELSE IF(ptop(j,ibox) .GE. 180._dp &
                .AND.ptop(j,ibox) .LT. 310._dp) THEN
                ipres(j)=2
              ELSE IF(ptop(j,ibox) .GE. 310._dp &
                .AND.ptop(j,ibox) .LT. 440._dp) THEN
                ipres(j)=3
              ELSE IF(ptop(j,ibox) .GE. 440._dp &
                .AND.ptop(j,ibox) .LT. 560._dp) THEN
                ipres(j)=4
              ELSE IF(ptop(j,ibox) .GE. 560._dp &
                .AND.ptop(j,ibox) .LT. 680._dp) THEN
                ipres(j)=5
              ELSE IF(ptop(j,ibox) .GE. 680._dp &
                .AND.ptop(j,ibox) .LT. 800._dp) THEN
                ipres(j)=6
              ELSE IF(ptop(j,ibox) .GE. 800._dp) THEN
                ipres(j)=7
              END IF
              
              !update frequencies
              IF(ipres(j) .GT. 0.AND.itau(j) .GT. 0) THEN
                fq_isccp(j,itau(j),ipres(j))= &
                  fq_isccp(j,itau(j),ipres(j))+ boxarea
              END IF
              
            END IF
            
          END IF
          
        ENDDO ! j
39    ENDDO
      
      !compute mean cloud properties
      DO j=1,npoints 
        IF (totalcldarea(j) .GT. 0._dp) THEN
          meanptop(j) = meanptop(j) / totalcldarea(j)
          IF (sunlit(j).EQ.1) THEN
            meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
            meantaucld(j) = tautab(MIN(255,MAX(1,NINT(meanalbedocld(j))))) 
          END IF
        END IF
      ENDDO ! j
      !
!     ---------------------------------------------------!

!     ---------------------------------------------------!
!     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
!
      IF (debugcol.NE.0) THEN
!     
        DO j=1,npoints,debugcol

          !produce character output
          DO ilev=1,nlev
            DO ibox=1,ncol
              acc(ilev,ibox)=0
            ENDDO
          ENDDO
          
          DO ilev=1,nlev
            DO ibox=1,ncol
              acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
              IF (levmatch(j,ibox) .EQ. ilev) &
                acc(ilev,ibox)=acc(ilev,ibox)+1
            ENDDO
          ENDDO
          
          !print test
          
          WRITE(ftn09,11) j
11        FORMAT('ftn09.',i4.4)
          OPEN(9, FILE=ftn09, FORM='FORMATTED')
          
          WRITE(9,'(a1)') ' '
          WRITE(9,'(10i5)') &
            (ilev,ilev=5,nlev,5)
          WRITE(9,'(a1)') ' '
          
          DO ibox=1,ncol
            WRITE(9,'(40(a1),1x,40(a1))') &
              (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev) &
              ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev) 
          END DO
          CLOSE(9)
          
          IF (ncolprint.NE.0) THEN
            WRITE(6,'(a1)') ' '
            WRITE(6,'(a2,1X,5(a7,1X),a50)') & 
              'ilev',&
              'pfull','at',&
              'cc*100','dem_s','dtau_s',&
              'cchar'
            
            WRITE (6,'(a)') 'skt(j):'
            WRITE (6,'(8f7.2)') skt(j)
            
            WRITE (6,'(8I7)') (ibox,ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'tau:'
            WRITE (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'tb:'
            WRITE (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
            
            WRITE (6,'(a)') 'ptop:'
            WRITE (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
          ENDIF
          
        ENDDO
        
      END IF
      
      RETURN
    END SUBROUTINE isccp_cloud_types_v3_4

!==============================================================================

    SUBROUTINE random_uniform (random_numbers)
  ! The  KISS (Keep It Simple Stupid) random number generator
  !
  ! by George Marsaglia, 1999 (latest update 2007-06-23) 
  !
  ! Combines:
  !
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2 ^ 32
  ! (2) A 3-shift shift-register generator, period 2 ^ 32-1
  ! (3) Two 16-bit multiply-with-carry generators, 
  !     period 597273182964842497 > 2^59
  !
  !  Overall period > 2^123. 
  !
  ! L. Kornblueh, MPI, 2008-09-29
  !
      IMPLICIT NONE
  

      INTEGER :: x = 123456789
      INTEGER :: y = 362436069
      INTEGER :: z = 21288629
      INTEGER :: w = 14921776
      INTEGER :: c = 0



      REAL(dp) :: random_numbers(:)
      
      INTEGER :: i, t, kiss
    
      DO i = 1, SIZE(random_numbers)
        x = x + 545925293
        y = m(m(m(y, 13), -17), 5)
        t = z + w + c 
        z = w 
        c = ISHFT(t, -31) 
        w = IAND(t, 2147483647)
        kiss = x + y + w 
        ! move random number to interval 0 to 1:
        ! r = 1/2 * (1 + kiss / 2147483647)
        random_numbers(i) = 0.5_dp+2.3283064376228985E-10_dp*kiss
      END DO
      
    CONTAINS
      
      INTEGER FUNCTION m(y, k) 
        INTEGER, INTENT(in) :: y, k
        
        m = IEOR (y, ISHFT (y, k) )
      END FUNCTION m
      
    END SUBROUTINE random_uniform
    
!===============================================================================

END MODULE MESSY_SATSIMS_ISCCP
