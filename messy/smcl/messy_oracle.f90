 MODULE messy_oracle
!_______________________________________________________________________________
! DESCRIPTION
! -----------
 USE messy_main_blather,         ONLY: start_message, end_message

 USE messy_main_constants_mem,   ONLY: dp

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAX_SPECIES  = 100

   ! SUBMODEL
      CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr  = 'oracle'     ! name of module
      CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver  = '1.0'    ! module version

! SUBROUTINES
      PUBLIC :: oracle_read_nml_ctrl
      PUBLIC :: oracle_soap
      PUBLIC :: oracle_spfcn
      PUBLIC :: oracle_mode

!      INTEGER, PUBLIC ::  NOPOA = 10
      INTEGER, PUBLIC ::  NfPOA  = 4
      INTEGER, PUBLIC ::  NbbPOA  = 4
      INTEGER, PUBLIC ::  NfSOAsv = 1
      INTEGER, PUBLIC ::  NbbSOAsv = 1
      INTEGER, PUBLIC ::  NfSOAiv = 3
      INTEGER, PUBLIC ::  NbbSOAiv = 3
      INTEGER, PUBLIC ::  NSOAv  = 3
      INTEGER, PUBLIC ::  NSOAP   = 20 !NfPOA+NbbPOA+NfSOAsv+NbbSOAsv+NfSOAiv+NbbSOAiv+NaSOAv+NbSOAv+NaOSOAv+NbOSOAv
!
      CHARACTER(LEN=4), PUBLIC ::  aermod = 'gmxe'   !aerosol model for inorganic aerosols
! STRINGS FOR MOM SPECIES 
      CHARACTER(LEN=1000), PUBLIC :: SOGv01,SOGv02,SOGv03,SOGv04,SOGv05,SOGv06,SOGv07,SOGv08,SOGv09

      INTEGER, PUBLIC ::  nmode = 3   !number of modes used for organic aerosols 
      INTEGER, PUBLIC ::  tmode(3)    !type of modes
!
!-----------------------------------------------------------------------
!     Variables that are initialized in soapdat.f
!
!     mwsoap  -- molecular weights of CG/SOA species (g/mol)
!     csat    -- saturation concentrations of CG/SOA species (ug/m3)
!     cstemp  -- temperatures corresponding to saturation concentrations
!                of CG/SOA species (K)
!     deltah  -- enthalpy of vaporization of CG/SOA species (kJ/mol)
!     flagsoap-- set to 1 if CG/SOA species forms solutions; 0 if not
!-----------------------------------------------------------------------
!
! op_pj_20150312+ did not compile with g95 due to data size inconsistency
!!$      DOUBLE PRECISION, PUBLIC    :: mwsoap(50)
!!$      DOUBLE PRECISION, PUBLIC    :: csat(50)
!!$      DOUBLE PRECISION, PUBLIC    :: cstemp(50)
!!$      DOUBLE PRECISION, PUBLIC    :: deltah(50)
!!$      INTEGER, PUBLIC :: flagsoap(50)
!mz_ap_20170622 increase size to accomodate all tracers
      REAL(dp), PUBLIC    :: SOGv_mw(10,MAX_SPECIES) = 0.0_dp
      REAL(dp), PUBLIC    :: mwsoap(30)
      REAL(dp), PUBLIC    :: csat(30)
      REAL(dp), PUBLIC    :: cstemp(30)
      REAL(dp), PUBLIC    :: deltah(30)
      INTEGER, PUBLIC :: flagsoap(30)
! op_pj_20150312-

      data tmode    /2,3,4/   
!
!      data SOGv_mw  /150.d0, 150.d0, 150.d0, 150.d0,&
!                     150.d0, 150.d0, 150.d0, 150.d0,&
!                     150.d0,                        &
!                     150.d0,                        &
!                     150.d0, 150.d0, 150.d0,        &
!                     150.d0, 150.d0, 150.d0,        &
!                     150.d0, 150.d0, 150.d0, 150.d0,&
!                     150.d0, 150.d0, 150.d0,        &
!                     150.d0, 150.d0, 150.d0, 150.d0,&
!                     150.d0, 150.d0, 150.d0         /
      data mwsoap   /250.d0, 250.d0, 250.d0, 250.d0,&
                     250.d0, 250.d0, 250.d0, 250.d0,&
                     250.d0,                        &
                     250.d0,                        &
                     250.d0, 250.d0, 250.d0,        &
                     250.d0, 250.d0, 250.d0,        &
                     180.d0, 180.d0, 180.d0, 180.d0,&
                     180.d0, 180.d0, 180.d0,        &
                     150.d0, 150.d0, 150.d0, 150.d0,&
                     150.d0, 150.d0, 150.d0         /
      data csat     /1.d-1, 1.d1,  1.d3,  1.d5,     &
                     1.d-1, 1.d1,  1.d3,  1.d5,     &
                     1.d-1,                         &
                     1.d-1,                         &
                     1.d-1, 1.d1,  1.d3,            &
                     1.d-1, 1.d1,  1.d3,            &
                     1.d0,  1.d1,  1.d2,  1.d3,     &
                     1.d0,  1.d1,  1.d2,            &
                     1.d0,  1.d1,  1.d2,  1.d3,     &
                     1.d0,  1.d1,  1.d2             /
      data cstemp   /300.d0, 300.d0, 300.d0, 300.d0,&
                     300.d0, 300.d0, 300.d0, 300.d0,&
                     300.d0,                        &    
                     300.d0,                        &    
                     300.d0, 300.d0, 300.d0,        &
                     300.d0, 300.d0, 300.d0,        &
                     300.d0, 300.d0, 300.d0, 300.d0,&
                     300.d0, 300.d0, 300.d0,        &
                     300.d0, 300.d0, 300.d0, 300.d0,&
                     300.d0, 300.d0, 300.d0         /
      data deltah   /106.d3,94.d3, 82.d3, 70.d3,    &
                     106.d3,94.d3, 82.d3, 70.d3,    &
                     106.d3,                        & 
                     106.d3,                        & 
                     106.d3,94.d3, 82.d3,           &
                     106.d3,94.d3, 82.d3,           &
                     30000.d0, 30000.d0, 30000.d0, 30000.d0,&
                     30000.d0, 30000.d0, 30000.d0,          &
                     30000.d0, 30000.d0, 30000.d0, 30000.d0,&
                     30000.d0, 30000.d0, 30000.d0            /
      data flagsoap /30*1/

!-----------------------------------------------------------------------
 CONTAINS

    SUBROUTINE oracle_read_nml_ctrl(status, iou)
   
! oracle MODULE ROUTINE (CORE)
!
! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
!
! Author: Alexandra Tsimpidi, MPIC, 2013
 
     USE messy_main_tools,           ONLY: read_nml_open, read_nml_check, read_nml_close
     USE messy_main_blather,         ONLY: start_message, end_message
     IMPLICIT NONE
 
! I/O
     INTEGER, INTENT(OUT) :: status ! error status
     INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
 
     NAMELIST /CTRL/ NfPOA  &
          , NbbPOA          &
          , NfSOAsv         &
          , NbbSOAsv        &
          , NfSOAiv         &
          , NbbSOAiv        &
          , NSOAv           &
          , NSOAP           &
          , aermod          &
          , nmode           &
          , tmode           &
          , SOGv_mw         &
          , mwsoap          &
          , csat            &
          , cstemp          &
          , deltah          &
          , flagsoap        &
          , SOGv01          &
          , SOGv02          &
          , SOGv03          &
          , SOGv04          &
          , SOGv05          &
          , SOGv06          &
          , SOGv07          &
          , SOGv08          &
          , SOGv09          
 
! LOCAL
     CHARACTER(LEN=*), PARAMETER :: substr = 'oracle_read_nml_ctrl'
     LOGICAL                     :: lex          ! file exists ?
     INTEGER                     :: fstat        ! file status
 

    CALL start_message(TRIM(modstr),'INITIALISATION', substr)

! INITIALIZE


     SOGv01 = ''
     SOGv02 = ''
     SOGv03 = ''
     SOGv04 = ''
     SOGv05 = ''
     SOGv06 = ''
     SOGv07 = ''
     SOGv08 = ''
     SOGv09 = ''

!
     status = 1 ! ERROR
 
     ! INITIALIZE GLOBAL CONTROL VARIABLES
     ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE
 
     CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
     write (6,*) 'read_nml_open', lex, substr, iou, modstr
     IF (.not.lex) RETURN    ! <modstr>.nml does not exist
 
     READ(iou, NML=CTRL, IOSTAT=fstat)
     CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
     write (6,*) 'read_nml_check', fstat, substr, iou, modstr
     IF (fstat /= 0) RETURN  ! error while reading namelist
 
! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
! CHECK TRACER INITIALIZATION

     WRITE(*,*) '    Initialization of organic aerosol module oracle'
     WRITE(*,*) ' '
     WRITE(*,*) 'DIAGNOSE NAMELISTS'
     WRITE(*,*) ' '
     WRITE(*,*) ' CTRL - oracle general options:'
     WRITE(*,*) ' '
     WRITE(*,*) '              NfPOA      = ', NfPOA
     WRITE(*,*) '              NbbPOA     = ', NbbPOA
     WRITE(*,*) '              NfSOAsv    = ', NfSOAsv
     WRITE(*,*) '              NbbSOAsv   = ', NbbSOAsv
     WRITE(*,*) '              NfSOAiv    = ', NfSOAiv
     WRITE(*,*) '              NbbSOAiv   = ', NbbSOAiv
     WRITE(*,*) '              NSOAv      = ', NSOAv
     WRITE(*,*) '              NSOAP      = ', NSOAP
     WRITE(*,*) '              aermod     = ', aermod
     WRITE(*,*) '              nmode      = ', nmode
     WRITE(*,*) '              tmode      = ', tmode
     WRITE(*,*) '              SOGv_mw    = ', SOGv_mw
     WRITE(*,*) '              mwsoap     = ', mwsoap
     WRITE(*,*) '              csat       = ', csat
     WRITE(*,*) '              cstemp     = ', cstemp
     WRITE(*,*) '              deltah     = ', deltah
     WRITE(*,*) '              flagsoap   = ', flagsoap

    CALL read_nml_close(substr, iou, modstr)
 
    status = 0  ! no ERROR
    CALL end_message(TRIM(modstr),'INITIALISATION', substr)

   END SUBROUTINE oracle_read_nml_ctrl
 
   subroutine oracle_soap(caer,cgas,tempk,p_pe,kproma,klev,csatT)
        implicit none

! DEFINITION OF VARIABLES:
!
!  INPUTS
!
!     ntot    - total number of CG/SOA species pairs
!     caer    - aerosol-phase concentrations of SOA species (ug/m3)
!     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
!     tempk   - temperature (K)
!     convfac - conversion factor: umol/m3 = ppm * convfac
!     iout    - standard output file unit
!     igrdchm - index for grid containing the grid cell
!     ichm    - i index of grid cell
!     jchm    - j index of grid cell
!     kchm    - k index of grid cell
!     mwpre   - molecular weight of pre-existing organic aerosol (g/mol)
!
!   OUTPUTS
!
!     caer    - aerosol-phase concentrations of SOA species (ug/m3)
!     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
!     csatT   - saturation concentrations of CG/SOA species at current T
!               (ug/m3)
!
!   VARIABLES USED WITHIN oracle_SOAP
!
!     i       - counter
!     icont   - counter
!     sumsc     - counter
!     nsol    - total number of solution-forming SOA species
!     cstemp  - temperatures corresponding to saturation concentrations
!               of CG/SOA species (K)
!     csat    - saturation concentrations of CG/SOA species (ug/m3)
!     deltah  - enthalpy of vaporization of CG/SOA species (J/mol)
!     flagsoap- set to 1 if CG/SOA species forms solutions; 0 if not
!     SOGv_mw - molecular weights of CG/SOA species from VOCs (g/mol)
!     mwsoap  - molecular weights of CG/SOA species (g/mol)
!     scaer   - aerosol-phase concentrations of solution-forming
!               SOA species (ug/m3)
!     scgas   - gas-phase concentrations of solution-forming
!               SOA species (ug/m3)
!     scsat   - saturation concentrations of solution-forming
!               SOA species (ug/m3)
!     sctot   - total concentrations of solution-forming SOA species
!               (ug/m3)
!     smw     - molecular weights of solution-forming SOA species
!               (g/mol)
!     znum    - counter for number of iterations
!     conmin  - use simple solution for species below this level (ug/m3)
!     cpremin - no pre-existing organics if cpre < cpremin (ug/m3)
!     xtol    - error tolerance for bi-section method
!
!***********************************************************************
!
! VARIABLE DECLARATION
!
      real (dp)        conmin, cpremin, xtol
      parameter ( conmin  = 1.d-6 )
      parameter ( cpremin = 1.01d-9 )
      parameter ( xtol    = 5.0d-5 )
!
      integer     ntot
      real (dp)   caer(NSOAP), cgas(NSOAP), ctot(NSOAP), csatT(NSOAP)
      real (dp)   smw(NSOAP), scsat(NSOAP)
      real (dp)   sctot(NSOAP), scaer(NSOAP), scgas(NSOAP)
      integer     idx(NSOAP)
!
      real (dp)   mwpre, cpre, tempk, sumsc,convfac
      integer     iout, p_pe,kproma, klev
      integer     i, icont, nsol, znum

      real (dp)   cpx, bb, cc, xend, fend, xmid, fmid, dx
!
!***********************************************************************
!
! Entry point
!
      iout=6 
      ntot=NSOAP
      do i=1,ntot
        ctot(i) = caer(i) + cgas(i)
      enddo
!
!
      mwpre = 220.d0
      cpre = 0.d0
      cpx = cpre/mwpre
!
! CHANGE SATURATION CONCENTRATIONS ACCORDING TO CURRENT TEMPERATURE
!
      do i=1,ntot
         csatT(i)=csat(i)*(cstemp(i)/tempk)*dexp((deltah(i)/8.314d0)&
     &                *(1.d0/cstemp(i)-1.d0/tempk))
      enddo
!
! CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE
! CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS
! COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN conmin ARE IGN0RED
! MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS
!
      icont=0
      do i=1,ntot
         if (flagsoap(i).eq.0) then
            cgas(i) = dmin1(ctot(i), csatT(i))
            caer(i) = ctot(i) - cgas(i)
         elseif (ctot(i).lt.conmin) then
            cgas(i) = ctot(i)
            caer(i) = 0.d0
         else
            icont=icont+1
            idx(icont) = i
            smw(icont)=mwsoap(i)
            scsat(icont)=csatT(i)
            sctot(icont)=ctot(i)
            scaer(icont)=caer(i)
         endif
      enddo
      nsol=icont
!
! Check for a trivial solution
!
      if (nsol.eq.0) goto 1000
      if (nsol.eq.1) then
         if (cpre.lt.cpremin) then
           scgas(1) = dmin1(sctot(1), scsat(1))
           scaer(1) = sctot(1) - scgas(1)
         else ! This case has an analytical solution
           bb = scsat(1)-sctot(1)+cpx*smw(1)
           cc = -sctot(1)*cpx*smw(1)
           scaer(1) = dmin1( sctot(1), .5*(-bb+DSQRT(bb*bb-4.d0*cc)) )
           scgas(1) = sctot(1) - scaer(1)
         endif
         goto 900
      endif
      sumsc=0.d0
      do i=1,nsol
         sumsc = sumsc + sctot(i)/scsat(i)
      enddo
      if (cpre.lt.cpremin .and. sumsc.le.1.d0) then
         do i=1,nsol
            scgas(i)=sctot(i)
            scaer(i)=0.d0
         enddo
         goto 900
      endif
!
! Find the solution using a bi-section method (approach from max)
!
      xend = 0.d0
      do i = 1, nsol
        xend = xend + sctot(i)/smw(i)
      enddo
      xend = xend + cpx
      call oracle_spfcn (nsol,sctot,scsat,scaer,smw,cpx,xend,fend)
      if (dabs(fend).le.xtol*xend) goto 99
      if (fend.gt.0.d0) then
        write (iout,'(//,a)') ' ERROR in oracle_SOAP:'
        write (iout,'(/,a)') ' ERROR: positive end point'
        goto 50
      endif
      dx = xend - cpx
      do znum = 1, 200
        dx = 0.5d0 * dx
        xmid = xend - dx
        call oracle_spfcn (nsol,sctot,scsat,scaer,smw,cpx,xmid,fmid)
        if (dabs(fmid).le.xtol*xmid .or. dx.le.xtol*xmid) goto 100
        if (fmid.lt.0.d0) xend = xmid
      enddo
       write (iout,'(//,a)') ' ERROR in oracle_SOAP:'
       write (iout,'(/,a)') ' ERROR: max number of iterations reached'
 50    write (iout,'(a,i3,i4)') &
     &                 ' cell(kproma,klev) = ', kproma,klev
       write (iout,'(a5,2a15)') ' spec','total [ug/m3]','c* [ug/m3]'
       write (iout,'(i5,1p2e15.6)') (idx(i),sctot(i),scsat(i),i=1,nsol)
       write (iout,'(i5,1p1e15.6)') (idx(i),mwsoap(i),i=1,nsol)
       write (iout,'(a5,e15.6)') ' cpre',cpre
      STOP
!
! Converged
!
  99  xmid = xend
 100  continue
      do i=1,nsol
         scaer(i) = dmin1( sctot(i), scaer(i) )
         scgas(i) = sctot(i) - scaer(i)
      enddo

!
! REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
!
 900  continue
      do i=1,nsol
         caer(idx(i))=scaer(i)
         cgas(idx(i))=scgas(i)
      enddo
!
! Convert to ppm if inputs in ppm
!
 1000 continue

!

      return
      end subroutine oracle_soap
    subroutine oracle_spfcn (n,ct,cs,ca,mw,cpx,tom,fval)
      implicit none
!
! oracle_SPFCN calculates the objective function for the bi-section solver in oracle_SOAP
!     Total Organics in Mole (TOM) = sumsc_i(C_i(aer)/MW_i) + C_pre/MW_pre
!     C_i(aer) = C_i(tot) - x_i * Cstar_i
!              = C_i(tot) - (C_i(aer)/MW_i/TOM) * Cstar_i
!  => C_i(aer) = C_i(tot) * TOM / (TOM + Cstar_i/MW_i)
!  => sumsc_i(C_i(tot) * TOM / (TOM*MW_i + Cstar_i)) + C_pre/MW_pre - TOM = 0
!
! Called by oracle_SOAP
!
      integer     n,i
      real (dp)       ct(n),cs(n),ca(n),mw(n),cpx,tom,fval
!
      fval = 0.d0
      do i = 1, n
        ca(i) = ct(i) * tom / ( tom + cs(i) / mw(i) )
        fval  = fval + ca(i) / mw(i)
      enddo
      fval = fval + cpx - tom
!
      return
      end subroutine oracle_spfcn
     
   subroutine oracle_mode(caer,cgas,qaer,qgas,csatT,rsec,p_pe,kproma,klev) !gmxelink 
   USE messy_main_constants_mem,   ONLY: dp
        implicit none

    integer flag,iout, p_pe,kproma, klev
    integer i,j,iter,itermax
    real(dp) caer(NSOAP), cgas(NSOAP), csatT(NSOAP),qaer(NSOAP,nmode),qgas(NSOAP)
    real(dp) accom           ! accomodation coefficient of condensible gases
    real(dp) rlambda         ! Mean free path of condensible gases
    real(dp) dsec(nmode)     ! mode diameters - gmxelink
    REAL(dp), DIMENSION(nmode) :: rsec     ! mode radius - gmxelink
    real(dp) qt(nmode)      ! Total aerosol concentration at each mode
    real(dp) qn(3)       ! 0th moment (number) of total mass for mode (qt/dsec**3)
    real(dp) tinys           ! minimum non-zero concentration (ug/m3)
    real(dp) dq, frqtot,frt,mnk,mxk,frq(nmode),xsum(nmode)

! gmxelink
    do i=1,nmode
    dsec(i)= MAX(2.d0*rsec(i),1.d-9)
    end do
!    dsec(1)= 2.d0*0.258d-7
!    dsec(2)= 2.d0*0.258d-6
!    dsec(3)= 2.d0*0.258d-5
!gmxelink
    itermax = 1000
    tinys=1.0d-9
    accom=1.0d-1
    rlambda=0.065d0

       do j = 1,nmode                                                      
         xsum(j)=0.0d0
         qt(j) = 0.0d0
         do i = 1,NSOAP                                                        
           xsum(j)=xsum(j)+qaer(i,j)/mwsoap(i)                  
           qt(j)=qt(j)+qaer(i,j)                  
         enddo                                                               
!           qn(j)=qt(j)/dsec(tmode(j)-1)**3 ! calculate 0th moment(number of particles)
            qn(j)=qt(j)/dsec(j)**3 ! calculate 0th moment(number of particles) - gmxelink
         ! give any non-zero value to xsum if it is zero in the mode
         if(xsum(j).lt.tinys) xsum(j) = tinys
       enddo
      
       do i = 1,NSOAP
         ! calculate DQ 
         dq = qgas(i) - cgas(i) ! ug/m3
    
         ! calculate frq
         frqtot = 0.d0
         mnk = 0.d0 
         mxk = 0.d0 
         do j = 1,nmode
           frq(j) = qn(j)*( qgas(i)/mwsoap(i)- qaer(i,j)/mwsoap(i)/xsum(j)*csatT(i))* &                                      
!      &        dsec(tmode(j)-1)/(1.d0+rlambda/(accom*dsec(tmode(j)-1)))
       &        dsec(j)/(1.d0+rlambda/(accom*dsec(j)))  !gmxelink
           mnk = dmin1(mnk,frq(j)) 
           mxk = dmax1(mxk,frq(j)) 
         enddo

         if(dq.gt.0.d0 .and. mnk.lt.0.d0 .and. mxk.gt.0.d0) then
           do j = 1,nmode                                  
             frq(j)=dmax1(frq(j)-mnk,0.d0)                 
           enddo                                               
         elseif(dq.lt.0.d0 .and. mxk.gt.0.d0 .and. mnk.lt.0.d0) then 
           do j = 1,nmode                                  
             frq(j)=dmin1(frq(j)-mxk,0.d0)                 
           enddo                                               
         endif                                                 
         do j = 1,nmode                                    
           frqtot = frqtot + frq(j)                         
         enddo                                                 
         ! normalize frq
         do j = 1,nmode
           frq(j) = frq(j) / frqtot
         enddo

         ! condense all condensing species
         if(dq.gt.0.d0) then
           do j = 1,nmode
             qaer(i,j) = qaer(i,j) + dq * frq(j)
           enddo
           ! set the gas species concentration to the equilibrium value
           qgas(i) = cgas(i)


         ! evaporate all evaporating species
         elseif(dq.lt.0.d0) then
           iter = 0
  100      frt = 1.d0
           do j = 1,nmode
             if(frq(j).gt.0.d0) frt = DMAX1( DMIN1(qaer(i,j)/(-dq*frq(j)),frt), 0.d0)
           enddo
           frqtot = 0.d0
           do j = 1,nmode
             qaer(i,j) = DMAX1(qaer(i,j)+frt*dq*frq(j),0.d0)
             if(qaer(i,j).lt.tinys) frq(j) = 0.d0
             frqtot = frqtot + frq(j)
           enddo

           ! check if we should evaporate more
           dq = (1.d0 - frt) * dq
           if(dq.lt.-1.d-7) then
             if(frqtot.gt.tinys) then ! we have sections which are not empty
               if(iter.le.itermax) then ! check infinite loop
                 iter = iter + 1
                 do j = 1,nmode
                   frq(j) = frq(j) / frqtot
                 enddo
                 goto 100
               endif
             endif
             ! we need to evaporate more to achieve equilibrium
             ! but we completely evaporate the species in all sections
             ! or exceed itermax
             if(dq.lt.-1.d-3) then                                       
               write(6,'(A9,E15.5,4I4)')'oracle_mode:',dq,i,p_pe,kproma,klev
             endif                                                      
           endif
           ! now set the gas species concentration conservatively
           qgas(i) = cgas(i)+caer(i)
           do j = 1,nmode
             qgas(i) = qgas(i) - qaer(i,j)
           enddo
         endif
       enddo
       return
   end subroutine oracle_mode

 END MODULE messy_oracle
