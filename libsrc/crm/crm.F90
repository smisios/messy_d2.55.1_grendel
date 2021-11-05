!---------------------------------------------------------------
!  Super-parameterization's main driver 
!  Marat Khairoutdinov, 2001-2009
!---------------------------------------------------------------

subroutine crm        (lchnk, icol, &
                       tl, ql, qccl, qiil, ul, vl, &
                       ps, pmid, pdel, phis, &
                       zmid, zint, dt_gl, plev, &
                       ultend, vltend, qltend, qcltend, qiltend, sltend, &
                       crm_buffer, qrad_crm, &
                       qv_crm, qc_crm, qi_crm, qpc_crm, qpi_crm, prec_crm, &  !!! added qv_crm um_hr_20140813
                       cld_crm, cdnc_crm, radlp_crm, radip_crm,  &            !!! added for cloudopt um_hr_20150123 
                       t_rad, qv_rad, qc_rad, qi_rad, &
                       precc, precl, precsc, precsl, &
                       cltot, clhgh, clmed, cllow, cld, cdnc, &
                       cldtop, cldbot, convtoppdf, convbotpdf, &   !!! added cldbot, convtoppdf, convbotpdf um_hr_20140918
                       gicewp, gliqwp, &
                       mc, mcup, mcdn, mcuup, mcudn, &
                       crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
                       tkez, tkesgsz, flux_u, flux_v, flux_qt, fluxsgs_qt,flux_qp, &
                       pflx, qt_ls, qt_trans, qp_trans, qp_fall, &
                       qp_evp, qp_src, t_ls, prectend, precstend, &
                       ocnfrac, wnd, tau00, bflxls, &
                       taux_crm, tauy_crm, z0m, timing_factor, conv_counter, &
                       massfu, massfd, u_entr, u_detr, d_entr, d_detr, &
                       pxtp_crm ) !!! added new up-/downdraft massfluxes and entr-/detrainment and tracer variable

!            dolong, doshort, nrad0, &
!            latitude00, longitude00, day00, pres00, tabs_s0, case0, &
!            radlwup0, radlwdn0, radswup0, radswdn0, radqrlw0, radqrsw0, &
!            lwnsxy,swnsxy,lwntxy,swntxy,solinxy,lwnscxy,swnscxy,lwntcxy,swntcxy,lwdsxy,swdsxy)


!---------------------------------------------------------------

!        use shr_kind_mod, only: r8 => shr_kind_r8
!        use buffer, only: crmvars
        use vars
        use params
        use microphysics
        use crmtracers

        ! um_hr_20190315+
        use micro_sam1mom
        use micro_m2005
        ! um_hr_20190315-

        implicit none

! um_hr_20131119+
        ! added to get rid of external module uses
!        integer, parameter :: crmvars = 7
        integer, parameter :: R8 = selected_real_kind(12) 
! um_hr_20131911-

!        integer, parameter :: r8 = 8

!  Input:

         integer, intent(in) :: lchnk    ! chunk identifier                           !
         integer, intent(in) :: icol     ! column identifier                          !
         integer, intent(in) :: plev     ! number of levels                           !
         real(r8), intent(in) :: ps ! Global grid surface pressure (Pa)               !
         real(r8), intent(in) :: pmid(plev) ! Global grid pressure (Pa)               !
         real(r8), intent(in) :: pdel(plev) ! Layer's pressure thickness (Pa)         !
         real(r8), intent(in) :: phis ! Global grid surface geopotential (m2/s2)      !  is zero?!
         real(r8), intent(in) :: zmid(plev) ! Global grid height (m)                  !
         real(r8), intent(in) :: zint(plev+1)! Global grid interface height (m)       !
         real(r8), intent(in) :: qrad_crm(crm_nx, crm_ny, crm_nz) ! CRM rad. heating
         real(r8), intent(in) :: dt_gl ! global model's time step
         real(r8), intent(in) :: ocnfrac ! area fraction of the ocean
         real(r8), intent(in) :: tau00  ! large-scale surface stress (N/m2)
         real(r8), intent(in) :: wnd  ! large-scale surface wind (m/s)
         real(r8), intent(in) :: bflxls  ! large-scale surface buoyancy flux (K m/s)
!         logical, intent(in)  :: doshort ! compute shortwave radiation
!         logical, intent(in)  :: dolong ! compute longwave radiation
!         real(r8), intent(in) :: day00 ! initial day
!         real(r8), intent(in) :: latitude00
!         real(r8), intent(in) :: longitude00
!         real(r8), intent(in) :: pres00
!         real(r8), intent(in) :: tabs_s0
!         integer , intent(in) :: nrad0
!         character *40 case0  ! 8-symbol id-string to identify a case-name

!  Input/Output:
         
         real(r8), intent(inout) :: tl(plev) ! Global grid temperature (K)            !
         real(r8), intent(inout) :: ql(plev) ! Global grid water vapor (g/g)          !
         real(r8), intent(inout) :: qccl(plev)! Global grid cloud liquid water (g/g)  !
         real(r8), intent(inout) :: qiil(plev)! Global grid cloud ice (g/g)           !
         real(r8), intent(inout) :: ul(plev) ! Global grid u (m/s)                    !
         real(r8), intent(inout) :: vl(plev) ! Global grid v (m/s)                    !
         real(r8), intent(inout), target :: crm_buffer(crm_nx, crm_ny, crm_nz,1:crmvars)
         real(r8), intent(inout) :: cltot ! shaded cloud fraction
         real(r8), intent(inout) :: clhgh ! shaded cloud fraction
         real(r8), intent(inout) :: clmed ! shaded cloud fraction
         real(r8), intent(inout) :: cllow ! shaded cloud fraction
         real(r8), intent(inout) :: pxtp_crm(crm_nx, crm_ny, crm_nz, ntracers) ! crm tracer fields  ! um_hr20141210
         
!  Output
         
         real(r8), intent(out) :: ultend(plev) ! tendency of ul
         real(r8), intent(out) :: vltend(plev) ! tendency of vl
         real(r8), intent(out) :: sltend(plev) ! tendency of static energy
         real(r8), intent(out) :: qltend(plev) ! tendency of water vapor
         real(r8), intent(out) :: qcltend(plev)! tendency of cloud liquid water
         real(r8), intent(out) :: qiltend(plev)! tendency of cloud ice
         real(r8), intent(out) :: t_rad (crm_nx, crm_ny, crm_nz) ! rad temperature
         real(r8), intent(out) :: qv_rad(crm_nx, crm_ny, crm_nz) ! rad vapor
         real(r8), intent(out) :: qc_rad(crm_nx, crm_ny, crm_nz) ! rad cloud water
         real(r8), intent(out) :: qi_rad(crm_nx, crm_ny, crm_nz) ! rad cloud ice
         real(r8), intent(out) :: precc ! convective precip rate (m/s)
         real(r8), intent(out) :: precl ! stratiform precip rate (m/s)
         real(r8), intent(out) :: cld(plev)  ! cloud fraction
         real(r8), intent(out) :: cdnc(plev) ! mean cloud droplet number concentration
         real(r8), intent(out) :: cldtop(plev)  ! cloud top pdf
         real(r8), intent(out) :: cldbot(plev)  ! cloud bot pdf           ! um_hr_20140909 added variable cldbot
         real(r8), intent(out) :: convtoppdf(plev) ! conv. cloud top pdf  ! um_hr_20140918 added variable convtoppdf
         real(r8), intent(out) :: convbotpdf(plev) ! conv. cloud bot pdf  ! um_hr_20140918 added variable convbotpdf

         real(r8), intent(out) :: gicewp(plev)  ! ice water path
         real(r8), intent(out) :: gliqwp(plev)  ! liquid water path
         real(r8), intent(out) :: mc(plev)   ! cloud mass flux
         real(r8), intent(out) :: mcup(plev) ! updraft cloud mass flux
         real(r8), intent(out) :: mcdn(plev) ! downdraft cloud mass flux
         real(r8), intent(out) :: mcuup(plev) ! unsat updraft cloud mass flux
         real(r8), intent(out) :: mcudn(plev) ! unsat downdraft cloud mass flux

         ! um_hr_20141024 added variables+
         real(r8), intent(out) :: massfu(plev) ! convective updraft cloud mass flux ! 
         real(r8), intent(out) :: massfd(plev) ! convective downdraft cloud mass flux !
         real(r8), intent(out) :: u_entr(plev) ! convective updraft entrainment !
         real(r8), intent(out) :: u_detr(plev) ! convective updraft detrainment !
         real(r8), intent(out) :: d_entr(plev) ! convective downdraft entrainment !
         real(r8), intent(out) :: d_detr(plev) ! convective downdraft detrainment !
         ! um_hr_20141024 added variables-

         real(r8), intent(out) :: crm_qc(plev)  ! mean cloud water
         real(r8), intent(out) :: crm_qi(plev)  ! mean cloud ice
         real(r8), intent(out) :: crm_qs(plev)  ! mean snow
         real(r8), intent(out) :: crm_qg(plev)  ! mean graupel
         real(r8), intent(out) :: crm_qr(plev)  ! mean rain
         real(r8), intent(out) :: flux_qt(plev) ! nonprecipitating water flux
         real(r8), intent(out) :: fluxsgs_qt(plev) ! sgs nonprecipitating water flux
         real(r8), intent(out) :: tkez(plev) ! tke profile
         real(r8), intent(out) :: tkesgsz(plev) ! sgs tke profile
         real(r8), intent(out) :: flux_u(plev) ! x-momentum flux
         real(r8), intent(out) :: flux_v(plev) ! y-momentum flux
         real(r8), intent(out) :: flux_qp(plev) ! precipitating water flux
         real(r8), intent(out) :: pflx(plev)    ! precipitation flux
         real(r8), intent(out) :: qt_ls(plev) ! tendency of nonprec water due to large-scale
         real(r8), intent(out) :: qt_trans(plev)! tendency of nonprec water due to transport
         real(r8), intent(out) :: qp_trans(plev) ! tendency of prec water due to transport
         real(r8), intent(out) :: qp_fall(plev) ! tendency of prec water due to fall-out
         real(r8), intent(out) :: qp_src(plev) ! tendency of prec water due to conversion
         real(r8), intent(out) :: qp_evp(plev) ! tendency of prec water due to evp
         real(r8), intent(out) :: t_ls(plev) ! tendency of lwse  due to large-scale
         real(r8), intent(out) :: prectend ! column integrated tendency in precipitating water+ice (kg/m2/s)
         real(r8), intent(out) :: precstend ! column integrated tendency in precipitating ice (kg/m2/s)
         real(r8), intent(out) :: precsc ! convective snow rate (m/s)
         real(r8), intent(out) :: precsl ! stratiform snow rate (m/s)
         real(r8), intent(out) :: taux_crm  ! zonal CRM surface stress perturbation (N/m2)
         real(r8), intent(out) :: tauy_crm  ! merid CRM surface stress perturbation (N/m2)
         real(r8), intent(out) :: z0m ! surface stress (N/m2)
         real(r8), intent(out) :: timing_factor ! crm cpu efficiency
         real(r8), intent(out) :: conv_counter  ! no. of convective active CRM cells ! add conv_counter um_hr_20140925

         real(r8), intent(inout)   :: qv_crm (crm_nx, crm_ny, crm_nz)! CRM water vapor ! add qv_crm um_hr_20140813
         real(r8), intent(inout)   :: qc_crm (crm_nx, crm_ny, crm_nz)! CRM cloud water
         real(r8), intent(inout)   :: qi_crm (crm_nx, crm_ny, crm_nz)! CRM cloud ice
         real(r8), intent(inout)   :: qpc_crm(crm_nx, crm_ny, crm_nz)! CRM precip water
         real(r8), intent(inout)   :: qpi_crm(crm_nx, crm_ny, crm_nz)! CRM precip ice
         real(r8), intent(out)     :: prec_crm(crm_nx, crm_ny)       ! CRM precipiation rate

         ! um_hr_20150123+
         real(r8), intent(out)   :: cld_crm(crm_nx, crm_ny, crm_nz)   ! CRM cloud cover (3D)
         real(r8), intent(out)   :: cdnc_crm(crm_nx, crm_ny, crm_nz)  ! CRM cloud droplet number concentration
         real(r8), intent(out)   :: radlp_crm(crm_nx, crm_ny, crm_nz) ! CRM effective radii for liquid droplets
         real(r8), intent(out)   :: radip_crm(crm_nx, crm_ny, crm_nz) ! CRM effective radii for ice droplets
         ! um_hr_20150123-

!         real(r8), intent(out) :: radlwup0(crm_nz)
!         real(r8), intent(out) :: radlwdn0(crm_nz)
!         real(r8), intent(out) :: radswup0(crm_nz)
!         real(r8), intent(out) :: radswdn0(crm_nz)
!         real(r8), intent(out) :: radqrlw0(crm_nz)
!         real(r8), intent(out) :: radqrsw0(crm_nz)
!         double precision, intent(out) :: lwnsxy,swnsxy,lwntxy,swntxy,solinxy
!         double precision, intent(out) :: lwnscxy,swnscxy,lwntcxy,swntcxy,lwdsxy,swdsxy

!  Local space:

        real(r8), pointer :: u_crm  (:,:,:) ! CRM v-wind component
        real(r8), pointer :: v_crm  (:,:,:) ! CRM v-wind component
        real(r8), pointer :: w_crm  (:,:,:) ! CRM w-wind component
        real(r8), pointer :: t_crm  (:,:,:) ! CRM temperature
        real dummy(nz), t00(nz)
        real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
        real tln(plev), qln(plev), qccln(plev), qiiln(plev), uln(plev), vln(plev)
        real cwp(nx,ny), cwph(nx,ny), cwpm(nx,ny), cwpl(nx,ny)
        real(r8) factor_xy, idt_gl
        real tmp1, tmp2
        real u2z,v2z,w2z
        integer i,j,k,l,ptop,nn,icyc, nstatsteps
        ! um_ht_20140414+
!        real(r8), parameter :: umax = 0.5*crm_dx/crm_dt ! maximum amplitude of the l.s. wind
        real(r8) umax
        ! um_ht_20140414-
        real(r8), parameter :: wmin = 2.   ! minimum up/downdraft velocity for stat
        real, parameter :: cwp_threshold = 0.001 ! threshold for cloud condensate for shaded fraction calculation
        logical flag_top(nx,ny), flag_bot(nx,ny), flag_dd(nx,ny) ! um_hr_20140918 added flag_bot and flag_dd
        integer convbot(nx,ny), convtop(nx,ny) ! um_hr_20141024 added local convbot, convtop variables
        real ustar, bflx, z0_est, qsat, omg, delta ! um_hr_20141024 added delta
        real umflx(nx,ny,crm_nz), dmflx(nx,ny,crm_nz) !um_hr_20141024 added local massflux variables
        real uentr(nx,ny,crm_nz), udetr(nx,ny,crm_nz), dentr(nx,ny,crm_nz), ddetr(nx,ny,crm_nz) ! um_hr_20141024 added local entrainment/detrainment variables
        real colprec,colprecs
        real(r8) zs ! surface elevation

!-----------------------------------------------

        umax = 0.5*crm_dx/crm_dt !um_ht_20140414+
        idt_gl = 1._r8/dt_gl
        ptop = plev-nzm+1
        factor_xy = 1._r8/dble(nx*ny)
        dummy = 0.
        t_rad = 0.
        qv_rad = 0.
        qc_rad = 0.
        qi_rad = 0.
        zs=phis/ggr
        bflx = bflxls

!-----------------------------------------
!        i = get_nstep() 
!        if(i.ge.40.and.lchnk.eq.708.and.icol.eq.3) then
!         write(i) lchnk, icol, &
!                       tl, ql, qccl= 0.0, qiil=0.0, ul, vl, &
!                       ps, pmid, pdel,
! phis = 0.0000000 -> zs=0.000 surface elevation -> not used
!                       zmid, zint, dt_gl, plev, &
! ultend, vltend, qltend, qcltend, qiltend, sltend, & 
!                       crm_buffer, qrad_crm, &
!                       qc_crm, qi_crm, qpc_crm, qpi_crm, prec_crm, &
!                       t_rad, qv_rad, qc_rad, qi_rad, &
!                       precc, precl, precsc, precsl, &
!                       cltot, clhgh, clmed, cllow, cld, cldtop, &
!                       gicewp, gliqwp, &
!                       mc, mcup, mcdn, mcuup, mcudn, &
!                       crm_qc, crm_qi, crm_qs, crm_qg, crm_qr, &
!                       tkez, tkesgsz, flux_u, flux_v, flux_qt, fluxsgs_qt,flux_qp, &
!                       pflx, qt_ls, qt_trans, qp_trans, qp_fall, &
!                       qp_evp, qp_src, t_ls, prectend, precstend, &
!                       ocnfrac, wnd, tau00, bflxls, &
!                       taux_crm, tauy_crm, z0m, timing_factor
!        close(i)
!        endif
!-----------------------------------------
      
        call setparm_crm()

!        doshortwave = doshort
!        dolongwave = dolong
!        day0 = day00-dt_gl/86400.
!        latitude = latitude00
!        longitude = longitude00
!        pres0 = pres00
!        tabs_s = tabs_s0
!        case = case0

        if(ocnfrac.gt.0.5) then
           OCEAN = .true.
        else
           LAND = .true.
        end if

        ! create CRM vertical grid and initialize some vertical reference arrays:
        do k = 1, nzm

           z(k) = zmid(plev-k+1) - zint(plev+1)
           zi(k) = zint(plev-k+2)- zint(plev+1)
           pres(k) = pmid(plev-k+1)/100.
           bet(k) = ggr/tl(plev-k+1)
           gamaz(k)=ggr/cp*z(k)

        end do
        zi(nz) =  zint(plev-nz+2)

        dz = 0.5*(z(1)+z(2))
        do k=2,nzm
           adzw(k) = (z(k)-z(k-1))/dz
        end do
        adzw(1) = 1.
        adzw(nz) = adzw(nzm)
        adz(1) = 1.
        do k=2,nzm-1
          adz(k) = 0.5*(z(k+1)-z(k-1))/dz
        end do
        adz(nzm) = adzw(nzm)

        do k=1,nzm
           grdf_x(k) = min(16.,dx**2/(adz(k)*dz)**2)
           grdf_y(k) = min(16.,dy**2/(adz(k)*dz)**2)
           grdf_z(k) = 1.
        end do
        
        do k = 1,nzm
          rho(k) = pdel(plev-k+1)/ggr/(adz(k)*dz)
        end do
        do k=2,nzm
          rhow(k) = 0.5*(rho(k)+rho(k-1))
        end do
        rhow(1) = 2*rhow(2) - rhow(3)
        rhow(nz)= 2*rhow(nzm) - rhow(nzm-1)
        colprec=0
        colprecs=0

        !  Initialize:
        u_crm => crm_buffer(1:nx,1:ny,1:nzm,1)
        v_crm => crm_buffer(1:nx,1:ny,1:nzm,2)
        w_crm => crm_buffer(1:nx,1:ny,1:nzm,3)
        t_crm => crm_buffer(1:nx,1:ny,1:nzm,4)

        ! limit the velocity at the very first step:
        if(u_crm(1,1,1).eq.u_crm(2,1,1).and.u_crm(3,1,2).eq.u_crm(4,1,2)) then
         do k=1,nzm
          do j=1,ny
           do i=1,nx
             u_crm(i,j,k) = min( umax, max(-umax,u_crm(i,j,k)) )
             v_crm(i,j,k) = min( umax, max(-umax,v_crm(i,j,k)) )*YES3D
           end do
          end do
         end do
        
        end if

        u(1:nx,1:ny,1:nzm)    = u_crm(1:nx,1:ny,1:nzm)
        v(1:nx,1:ny,1:nzm)    = v_crm(1:nx,1:ny,1:nzm)*YES3D
        w(1:nx,1:ny,1:nzm)    = w_crm(1:nx,1:ny,1:nzm)
        tabs(1:nx,1:ny,1:nzm) = t_crm(1:nx,1:ny,1:nzm)

        micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields) = crm_buffer(1:nx,1:ny,1:nzm,5:4+nmicro_fields)

        ! um_hr_20190315+
        IF (domicro_sam1mom) THEN
           qn(1:nx,1:ny,1:nzm) = crm_buffer(1:nx,1:ny,1:nzm,7)
        END IF
        ! um_hr_20190315-

        w(:,:,nz)                 = 0.
        dudt(1:nx,1:ny,1:nzm,1:3) = 0.
        dvdt(1:nx,1:ny,1:nzm,1:3) = 0.
        dwdt(1:nx,1:ny,1:nz,1:3)  = 0.
        tke(1:nx,1:ny,1:nzm)      = 0.
        tk(1:nx,1:ny,1:nzm)       = 0.
        tkh(1:nx,1:ny,1:nzm)      = 0.
        p(1:nx,1:ny,1:nzm)        = 0.

        ! um_hr_20190315+
        IF (domicro_sam1mom) call micro_init_sam1mom
        IF (domicro_m2005)   THEN
           ! um_hr_20190321+ 
           ! set micro_fields and temperature 
           ! in order to account for satadj in micro init...
           qv(1:nx,1:ny,1:nzm)  = micro_field(1:nx,1:ny,1:nzm,iqv)
           qcl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqcl)
           qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)
           do k=1,nzm
              tabs0(k) = 0.
              q0(k)    = 0.
              do j=1,ny
                 do i=1,nx
                    tabs0(k) = tabs0(k)+tabs(i,j,k)
                    q0(k)    = q0(k)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
                 end do
              end do
              tabs0(k) = tabs0(k) * factor_xy
              q0(k)    = q0(k)    * factor_xy
           end do
           ! um_hr_20190321-
           call micro_init_m2005
        END IF
        ! um_hr_20190315-

        do k=1,nzm
          
          u0(k)    = 0.
          v0(k)    = 0.
          t0(k)    = 0.
          t00(k)   = 0.
          tabs0(k) = 0.
          q0(k)    = 0.
          qv0(k)   = 0.
          qn0(k)   = 0.
          qp0(k)   = 0.
          tke0(k)  = 0.
          
          do j=1,ny
           do i=1,nx
            
            t(i,j,k) = tabs(i,j,k)+gamaz(k) &
                        -fac_cond*qcl(i,j,k)-fac_sub*qci(i,j,k) &
                        -fac_cond*qpl(i,j,k)-fac_sub*qpi(i,j,k)

            colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(plev-k+1)
            colprecs=colprecs+qpi(i,j,k)*pdel(plev-k+1)
            u0(k)=u0(k)+u(i,j,k)
            v0(k)=v0(k)+v(i,j,k)
            t0(k)=t0(k)+t(i,j,k)
            t00(k)=t00(k)+t(i,j,k)+fac_cond*qpl(i,j,k)+fac_sub*qpi(i,j,k)
            tabs0(k)=tabs0(k)+tabs(i,j,k)
            q0(k)=q0(k)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
            qv0(k) = qv0(k) + qv(i,j,k)
            qn0(k) = qn0(k) + qcl(i,j,k) + qci(i,j,k)
            qp0(k) = qp0(k) + qpl(i,j,k) + qpi(i,j,k)
            tke0(k)=tke0(k)+tke(i,j,k)

           end do
          end do

          u0(k) = u0(k) * factor_xy
          v0(k) = v0(k) * factor_xy
          t0(k) = t0(k) * factor_xy
          t00(k) = t00(k) * factor_xy
          tabs0(k) = tabs0(k) * factor_xy
          q0(k) = q0(k) * factor_xy
          qv0(k) = qv0(k) * factor_xy
          qn0(k) = qn0(k) * factor_xy
          qp0(k) = qp0(k) * factor_xy
          tke0(k) = tke0(k) * factor_xy

          l = plev-k+1
          uln(l) = min( umax, max(-umax,ul(l)) )
          vln(l) = min( umax, max(-umax,vl(l)) )*YES3D
          ttend(k) = (tl(l)+gamaz(k)- &
               fac_cond*(qccl(l)+qiil(l))-fac_fus*qiil(l)-t00(k))*idt_gl  !!!! l.-s. temp. tend. are equally distributed for every layer
          qtend(k) = (ql(l)+qccl(l)+qiil(l)-q0(k))*idt_gl
          utend(k) = (uln(l)-u0(k))*idt_gl
          vtend(k) = (vln(l)-v0(k))*idt_gl
          ug0(k) = uln(l)
          vg0(k) = vln(l)
          tg0(k) = tl(l)+gamaz(k)-fac_cond*qccl(l)-fac_sub*qiil(l)
          qg0(k) = ql(l)+qccl(l)+qiil(l)

        end do ! k

        uhl = u0(1)
        vhl = v0(1)

        ! estimate roughness length assuming logarithmic profile of velocity near the surface:
        ustar = sqrt(tau00/rho(1))
        z0    = z0_est(z(1),bflx,wnd,ustar)
        z0    = max(0.00001,min(1.,z0))

        timing_factor = 0.

        prectend  = colprec
        precstend = colprecs

        fluxbu   = 0.
        fluxbv   = 0.
        fluxbt   = 0.
        fluxbq   = 0.
        fluxtu   = 0.
        fluxtv   = 0.
        fluxtt   = 0.
        fluxtq   = 0.
        fzero    = 0.
        precsfc  = 0.
        precssfc = 0.

!---------------------------------------------------
        cld    = 0.
        cdnc   = 0.
        cldtop = 0.
        cldbot = 0.       ! um_hr_20140909
        convtoppdf = 0.   ! um_hr_20140918
        convbotpdf = 0.   ! um_hr_20140918
        convbot = 0
        convtop = 0
        gicewp  = 0
        gliqwp  = 0
        mc      = 0.
        mcup    = 0.
        mcdn    = 0.
        mcuup   = 0.
        mcudn   = 0.
        
        ! um_hr_20141024+
        cld_crm   = 0.
        cdnc_crm  = 0.
        radlp_crm = 0.
        radip_crm = 0.

        umflx = 0.
        dmflx = 0.
        uentr = 0.
        udetr = 0.
        dentr = 0.
        ddetr = 0.
        delta = 0.
        ! um_hr_20141024-

        crm_qc     = 0.
        crm_qi     = 0.
        crm_qs     = 0.
        crm_qg     = 0.
        crm_qr     = 0.
        flux_qt    = 0.
        flux_u     = 0.
        flux_v     = 0.
        fluxsgs_qt = 0.
        tkez       = 0.
        tkesgsz    = 0.
        flux_qp    = 0.
        pflx       = 0.
        qt_trans   = 0.
        qp_trans   = 0.
        qp_fall    = 0.
        qp_evp     = 0.
        qp_src     = 0.
        qt_ls      = 0.
        t_ls       = 0.

        uwle     = 0.
        uwsb     = 0.
        vwle     = 0.
        vwsb     = 0.
        qpsrc    = 0.
        qpfall   = 0.
        qpevp    = 0.
        precflux = 0.

!        radlwup0 = 0.
!        radlwdn0 = 0.
!        radswup0 = 0.
!        radswdn0 = 0.
!        radqrlw0 = 0.
!        radqrsw0 = 0.
!        lwnsxy = 0.
!        swnsxy = 0.
!        lwntxy = 0.
!        swntxy = 0.
!        solinxy = 0.
!        lwnscxy = 0.
!        swnscxy = 0.
!        lwntcxy = 0.
!        swntcxy = 0.
!        lwdsxy = 0.
!        swdsxy = 0.

!--------------------------------------------------

     ! um_hr_20190315+
     IF (domicro_sam1mom) THEN
        if(doprecip) call precip_init()
     END IF
     ! um_hr_20190315-
     
     tracer(1:nx,1:ny,1:nzm,1:ntracers) = pxtp_crm(1:nx,1:ny,1:nzm,1:ntracers)
     
     if(u(1,1,1).eq.u(2,1,1).and.u(3,1,2).eq.u(4,1,2)) THEN
        call setperturb()
     endif

        nstop = dt_gl/dt
        dt = dt_gl/nstop
        nsave3D = nint(60/dt)
!       if(nint(nsave3D*dt).ne.60)then
!          print *,'CRM: time step=',dt,' is not divisible by 60 seconds'
!          print *,'this is needed for output every 60 seconds'
!          stop
!       endif
        nstep = 0
        nprint = 1
        ncycle = 0
!        nrad = nstop/nrad0
        day=day0

!------------------------------------------------------------------
!   Main time loop    
!------------------------------------------------------------------

do while(nstep.lt.nstop) 
  nstep = nstep + 1
  time  = time + dt
  day   = day0 + time/86400.
!------------------------------------------------------------------
!  Check if the dynamical time step should be decreased 
!  to handle the cases when the flow being locally linearly unstable
!------------------------------------------------------------------

  ncycle = 1

  call kurant()       !!! ncycle is probably increased, for local unstable flow

  timing_factor = timing_factor+ncycle

  do icyc=1,ncycle

     icycle   = icyc
     dtn      = dt/ncycle
     dt3(na)  = dtn
     dtfactor = dtn/dt

!---------------------------------------------
!       the Adams-Bashforth scheme in time

     call abcoefs()    !!! calculation of coefficients for the Adam-Bashforth scheme
 
!---------------------------------------------
!       initialize stuff: 
        
     call zero()       !!! wind tendencies are set to zero (dudt, dvdt, dwdt)

!-----------------------------------------------------------
!       Buoyancy term:
             
     call buoyancy()   !!! calculation of dwdt

!-----------------------------------------------------------
!       Large-scale and surface forcing:

     call forcing()    !!! update crm prognostic variables with large-scale tendencies
     do k=1,nzm
      do j=1,ny
        do i=1,nx
          t(i,j,k) = t(i,j,k) + qrad_crm(i,j,k)*dtn
        end do
      end do
     end do

!----------------------------------------------------------
!       suppress turbulence near the upper boundary (spange):

     if(dodamping) call damping()   ! damping at domain top regions

!----------------------------------------------------------
!      Update the subdomain's boundaries for velocity

     call boundaries(0)             ! periodic boundary conditions for wind velocity

!---------------------------------------------------------
!       SGS TKE equation:       
           
     if(dosgs) call tke_full()      ! dosgs=.true.   (switch in setparm)

!---------------------------------------------------------
!   Ice fall-out

     if(docloud) then
        call ice_fall()             ! docloud=.true. (switch in setparm)
     end if

!---------------------------------------------------------
!        Update boundaries for scalars, sst,  SGS exchange coefficients 

     call boundaries(2)

!-----------------------------------------------
!       advection of momentum:

     call advect_mom()

!-----------------------------------------------
!       surface fluxes:

     if(dosurface) call crmsurface(bflx)  ! dosurface=.true. (switch in setparm)

!----------------------------------------------------------
!       SGS diffusion of momentum:

     if(dosgs) call diffuse_mom()   ! dosgs=.true.     (switch in setparm)

!-----------------------------------------------------------
!       Coriolis force:
             
     if(docoriolis) call coriolis() ! docoriolis=.false. (switch in setparm)
         
!---------------------------------------------------------
!       compute rhs of the Poisson equation and solve it for pressure. 

     call pressure()

!---------------------------------------------------------
!       find velocity field at n+1/2 timestep needed for advection of scalars:
         
     call adams()                   !

!----------------------------------------------------------
!     Update boundaries for velocity fields to use for advection of scalars:

     call boundaries(1)

!---------------------------------------------------------
!      advection of scalars :

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
     
     if(dosgs.and..not.dosmagor) then
      call advect_scalar(tke,dummy,tkewle,dummy,dummy,dummy,.false.)
     else if(doscalar) then
      call advect_scalar(tke,dummy,tkewle,s2leadv,s2legrad,swleadv,.true.)
     end if

!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no matter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
        end if
     end do

!   Precipitation fallout:
!

    if(doprecip) then
       ! um_hr_20190315+
       IF (domicro_sam1mom) call micro_precip_fall_sam1mom
       IF (domicro_M2005)   call micro_precip_fall_m2005
!      call micro_precip_fall()
       ! um_hr_20190315+
    end if

!---------------------------------------------------------
!      diffusion of scalars :

!        Update boundaries for scalars:

    if(dosgs) call boundaries(3)

      call diffuse_scalar(t,fluxbt,fluxtt,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
     
      if(.not.dosmagor) then
          call diffuse_scalar(tke,fzero,fzero,dummy,tkewsb, &
                                    dummy,dummy,dummy,.false.)
      else if(doscalar) then
          call diffuse_scalar(tke,fluxbq,fluxtq,dummy,tkewsb, &
                           s2lediff,s2lediss,swlediff,.true.)
      end if

!
!    diffusion of microphysics prognostics:
!
      ! um_hr_20190315+
      IF (domicro_sam1mom) call micro_flux_sam1mom()
      IF (domicro_m2005)   call micro_flux_m2005()
!      call micro_flux()
      ! um_hr_20190315-

      do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no matter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k)
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(micro_field(:,:,:,k),fluxbmk(:,:,k),fluxtmk(:,:,k), &
!!$                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
       end if
      end do

 ! diffusion of tracers:

      if(dotracers) then

        call tracers_flux()

        do k = 1,ntracers
          fluxbtmp = fluxbtr(:,:,k)
          fluxttmp = fluxttr(:,:,k)
          call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, &
               trdiff(:,k),trwsb(:,k), &
               dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k),trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)
 
        end do

      end if

!-----------------------------------------------------------
!    Cloud condensation/evaporation and precipitation processes:

      ! um_hr_20190315+
      if(docloud.or.dosmoke) THEN
      ! call micro_proc()
         IF (domicro_sam1mom) call micro_proc_sam1mom()
         IF (domicro_m2005)   call micro_proc_m2005()
      END IF
      ! um_hr_20190315-

!-----------------------------------------------------------
!    Compute field diagnostics and update the velocity field:

      call uvw()

      call diagnose()

!----------------------------------------------------------
!    Rotate the dynamic tendency arrays for Adams-bashforth scheme:

      nn=na
      na=nc
      nc=nb
      nb=nn

   end do ! icycle      
          
!----------------------------------------------------------
!----------------------------------------------------------

        cwp  = 0.
        cwph = 0.
        cwpm = 0.
        cwpl = 0.

        flag_top(:,:)  =  .true.
        !!!! um_hr_20140909+++
        flag_bot(:,:)  =  .true.

        do k=1,nzm
           l = plev-k+1
           do j=1,ny
              do i=1,nx
                 !!!! added cldbot to crm variables !!!!
                 tmp1 = rho(k)*adz(k)*dz*(qcl(i,j,k)+qci(i,j,k))
                 cwp(i,j) = cwp(i,j)+tmp1
                 if(cwp(i,j).gt.cwp_threshold.and.flag_bot(i,j)) then
                    cldbot(l) = cldbot(l) + 1.     
                    flag_bot(i,j) = .false.
                 end if
              end do
           end do
        end do
        
        cwp = 0.
        !!!! um_hr_20140909---

        do k=1,nzm
           l = plev-k+1
           do j=1,ny
              do i=1,nx

                 crm_qc(l) = crm_qc(l) + qcl(i,j,k)
                 crm_qi(l) = crm_qi(l) + qci(i,j,k)
                 crm_qr(l) = crm_qr(l) + qpl(i,j,k)
                 if (domicro_sam1mom) then
                    omg = max(0.,min(1.,(tabs(i,j,k)-tgrmin)*a_gr))
                    crm_qg(l) = crm_qg(l) + qpi(i,j,k)*omg
                    crm_qs(l) = crm_qs(l) + qpi(i,j,k)*(1.-omg)
                 elseif (domicro_m2005) then
                    crm_qg(l) = crm_qg(l) + qpi(i,j,k)
                    crm_qs(l) = crm_qs(l) + 0.     ! temporary solution
                    if (dopredictNc) then
                       cdnc(l) = cdnc(l) + micro_field(i,j,k,incl)
                    end if
                 end if
                                  
                 tmp1 = rho(nz-k)*adz(nz-k)*dz*(qcl(i,j,nz-k)+qci(i,j,nz-k))
                 cwp(i,j) = cwp(i,j)+tmp1
                 if(cwp(i,j).gt.cwp_threshold.and.flag_top(i,j)) then
                    cldtop(plev-(nz-k)+1) = cldtop(plev-(nz-k)+1) + 1.  !!! um_hr_20140908 changed to plev-(nz-k)+1 -> top to bottom GCM levels (cldtop(k) makes no sense)
                    flag_top(i,j) = .false.
                 end if
                 if(pres(nz-k).ge.700.) then
                    cwpl(i,j) = cwpl(i,j)+tmp1
                 else if(pres(nz-k).lt.400.) then
                    cwph(i,j) = cwph(i,j)+tmp1
                 else
                    cwpm(i,j) = cwpm(i,j)+tmp1
                 end if
                 ! qsat = qsatw_crm(tabs(i,j,k),pres(k))
                 ! if(qcl(i,j,k)+qci(i,j,k).gt.min(1.e-5,0.01*qsat)) then
                 
                 ! include calculation of overall mass flux
                 tmp1 = rho(k)*adz(k)*dz
                 if(tmp1*(qcl(i,j,k)+qci(i,j,k)).gt.cwp_threshold) then
                    cld_crm(i,j,k) = cld_crm(i,j,k) + 1.   ! um_hr_20150123 (include subgrid cloud cover for output)
                    cld(l) = cld(l) + 1.
                    if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                       mcup(l) = mcup(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end if
                    if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                       mcdn(l) = mcdn(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end if
                 else 
                    if(w(i,j,k+1)+w(i,j,k).gt.2*wmin) then
                       mcuup(l) = mcuup(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end if
                    if(w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then
                       mcudn(l) = mcudn(l) + rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end if
                 end if
                 
                 t_rad (i,j,k) = t_rad (i,j,k)+tabs(i,j,k)
                 qv_rad(i,j,k) = qv_rad(i,j,k)+max(0.,qv(i,j,k))
                 qc_rad(i,j,k) = qc_rad(i,j,k)+qcl(i,j,k)
                 qi_rad(i,j,k) = qi_rad(i,j,k)+qci(i,j,k)
                 gliqwp(l)=gliqwp(l)+qcl(i,j,k)
                 gicewp(l)=gicewp(l)+qci(i,j,k)
                 
              end do
           end do
        end do
        

! um_hr_20142310+
! CALCULATION OF ENTRAINMENT / DETRAINMENT AND UPDRAFT / DOWNDRAFT MASSFLUXES 
! CRITERION FOR CLOUDY/CONVECTIVE CRM CELLS: w > wmin, cwp > cwpmin
! CRITERION FOR DOWNDRAFTS: w < -wmin

!UPDRAFTS:
           flag_bot(:,:) = .false.
           flag_top(:,:) = .false.
           cwp = 0. ! reset cloud water path

           do j=1,ny
              do i=1,nx
                 do k=1,nzm
                    tmp1 = rho(k)*adz(k)*dz*(qcl(i,j,k)+qci(i,j,k))               ! density is not distinguished between different CRM cells!!
                    cwp(i,j) = cwp(i,j)+tmp1
                    if (.not. flag_bot(i,j) .and. w(i,j,k+1)+w(i,j,k).gt.2*wmin .and. cwp(i,j).gt.cwp_threshold) then        ! convective cloud bottom criterion
                       convbot(i,j)  = k
                       flag_bot(i,j) = .true.
                       convbotpdf(l) = convbotpdf(l) + 1.
                       EXIT
                    end if
                 end do
              end do
           end do
           
           cwp = 0.  ! reset cloud water path

           do j=1,ny
              do i=1,nx
                 do k=nzm,1,-1
                    tmp1 = rho(k)*adz(k)*dz*(qcl(i,j,k)+qci(i,j,k))
                    cwp(i,j) = cwp(i,j)+tmp1
                    if(.not. flag_top(i,j) .and. cwp(i,j).gt.cwp_threshold) then         
                       convtop(i,j)  = k
                       flag_top(i,j) = .true.
                       convtoppdf(l) = convtoppdf(l) + 1.
                       EXIT
                    end if
                 end do
              end do
           end do

           do j=1,ny
              do i=1,nx
                 if (flag_bot(i,j) .and. convtop(i,j) .gt. convbot(i,j)) then   ! flag_bot identifies CRM column as convective!
                    conv_counter = conv_counter + 1
                    do k=1,convtop(i,j)
                       ! calculation updraft massflux
                       umflx(i,j,k) = rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end do
                    do k=2,convtop(i,j)
                       ! calculation of updraft entrainment and detrainment
                       delta = umflx(i,j,k) - umflx(i,j,k-1)
                       if (delta .gt. 0.) then
                          uentr(i,j,k) = delta
                       else
                          udetr(i,j,k) = delta
                       end if
                    end do
                    uentr(i,j,1) = umflx(i,j,1)
                    udetr(i,j,convtop(i,j)) = umflx(i,j,convtop(i,j))
                 end if
              end do
           end do

!DOWNDRAFTS
           flag_dd(:,:) = .false.

           do j=1,ny
              do i=1,nx
                 do k=1,nzm
                    if (.not. flag_dd(i,j) .and. w(i,j,k+1)+w(i,j,k).lt.-2*wmin) then        ! downdraft criterion
                       flag_dd(i,j) = .true.
                       EXIT
                    end if
                 end do
              end do
           end do

           do j=1,ny
              do i=1,nx
                 if (flag_dd(i,j) ) then
                    do k=1,nzm
                       ! calculation downdraft massflux
                       dmflx(i,j,k) = rho(k)*0.5*(w(i,j,k+1)+w(i,j,k))
                    end do
                    do k=2,nzm
                       ! calculation of downdraft entrainment and detrainment
                       delta = dmflx(i,j,k) - dmflx(i,j,k-1)
                       if (delta .gt. 0.) then
                          dentr(i,j,k) = delta
                       else
                          ddetr(i,j,k) = delta
                       end if
                    end do
                    dentr(i,j,1)   = dmflx(i,j,1)
                    ddetr(i,j,nzm) = dmflx(i,j,nzm)
                 end if
              end do
           end do
           
           ! Modifiy massflux calculation...??? -> loop from convbot to convtop
           ! set upper and lower boundaries to the following: 
           ! uentr(i,j,convbot-1)=umflx(i,j,convbot-1) if convbot==1 then uentr(i,j,1)=umflx(i,j,1)
           ! uentr(i,j,convtop+1)=0 if convtop=29 then uentr(i,j,convtop)=0
           ! udetr(i,j,convbot-1)=0 if convbot=1 then udetr(i,j,convbot=0)
           ! udetr(i,j,convtop+1)=umflx(i,j,convtop+1) if convtop=29 then udetr(i,j,convtop)=umflx(i,j,convtop)
           !------------------------------
           ! dentr(i,j,convbot-1)=0 if convbot==1 then uentr(i,j,1)=0
           ! dentr(i,j,convtop+1)=dmflx(i,j,convtop+1) if convtop=29 then uentr(i,j,convtop)=dmflx(i,j,convtop)
           ! ddetr(i,j,convbot-1)=dmflx(i,j,convbot-1) if convbot=1 then udetr(i,j,convbot)=dmflx(i,j,conbot)
           ! ddetr(i,j,convtop+1)=0 if convtop=29 then udetr(i,j,convtop)=0

           ! save variables after crm timestep on GCM grid
           do k=1,nzm
              l = plev-k+1
              do j=1,ny
                 do i=1,nx
                    massfu(l) = massfu(l) + umflx(i,j,k)
                    massfd(l) = massfd(l) + dmflx(i,j,k)
                    u_entr(l) = u_entr(l) + uentr(i,j,k)
                    u_detr(l) = u_detr(l) + udetr(i,j,k)
                    d_entr(l) = d_entr(l) + dentr(i,j,k)
                    d_detr(l) = d_detr(l) + ddetr(i,j,k)
                 end do
              end do
           end do

! um_hr_20142310-

!        do k=1,nzm
!         radlwup0(k)=radlwup0(k)+radlwup(k)
!         radlwdn0(k)=radlwdn0(k)+radlwdn(k)
!         radqrlw0(k)=radqrlw0(k)+radqrlw(k)
!         radswup0(k)=radswup0(k)+radswup(k)
!         radswdn0(k)=radswdn0(k)+radswdn(k)
!         radqrsw0(k)=radqrsw0(k)+radqrsw(k)
!        end do
        
        do j=1,ny
         do i=1,nx
           if(cwp(i,j).gt.cwp_threshold)  cltot = cltot + 1.
           if(cwph(i,j).gt.cwp_threshold) clhgh = clhgh + 1.
           if(cwpm(i,j).gt.cwp_threshold) clmed = clmed + 1.
           if(cwpl(i,j).gt.cwp_threshold) cllow = cllow + 1.
         end do
        end do

!        call stepout()
!----------------------------------------------------------
        end do ! main loop
!----------------------------------------------------------

        tmp1 = 1._r8/ dble(nstop)
        t_rad = t_rad * tmp1
        qv_rad = qv_rad * tmp1
        qc_rad = qc_rad * tmp1
        qi_rad = qi_rad * tmp1

! no CRM tendencies above its top
        
        tln(1:ptop-1) = tl(1:ptop-1)
        qln(1:ptop-1) = ql(1:ptop-1)
        qccln(1:ptop-1)= qccl(1:ptop-1)
        qiiln(1:ptop-1)= qiil(1:ptop-1)
        uln(1:ptop-1) = ul(1:ptop-1)
        vln(1:ptop-1) = vl(1:ptop-1)

!  Compute tendencies due to CRM:
        
        tln(ptop:plev) = 0.
        qln(ptop:plev) = 0.
        qccln(ptop:plev)= 0.
        qiiln(ptop:plev)= 0.
        uln(ptop:plev) = 0.
        vln(ptop:plev) = 0.
        
        colprec=0
        colprecs=0
        do k = 1,nzm
         l = plev-k+1
         do i=1,nx
          do j=1,ny
             colprec=colprec+(qpl(i,j,k)+qpi(i,j,k))*pdel(plev-k+1)
             colprecs=colprecs+qpi(i,j,k)*pdel(plev-k+1)
             tln(l) = tln(l)+tabs(i,j,k)
             qln(l) = qln(l)+qv(i,j,k)
             qccln(l)= qccln(l)+qcl(i,j,k)
             qiiln(l)= qiiln(l)+qci(i,j,k)
             uln(l) = uln(l)+u(i,j,k)
             vln(l) = vln(l)+v(i,j,k)
          end do ! k
         end do
        end do ! i

        tln(ptop:plev) = tln(ptop:plev) * factor_xy
        qln(ptop:plev) = qln(ptop:plev) * factor_xy
        qccln(ptop:plev) = qccln(ptop:plev) * factor_xy
        qiiln(ptop:plev) = qiiln(ptop:plev) * factor_xy
        uln(ptop:plev) = uln(ptop:plev) * factor_xy
        vln(ptop:plev) = vln(ptop:plev) * factor_xy

        sltend = cp * (tln - tl) * idt_gl
        qltend = (qln - ql) * idt_gl
        qcltend = (qccln - qccl) * idt_gl
        qiltend = (qiiln - qiil) * idt_gl
        ultend = (uln - ul ) * idt_gl    ! mz_hr_20140625
        vltend = (vln - vl ) * idt_gl    ! mz_hr_20140625

        prectend=(colprec-prectend)/ggr*factor_xy * idt_gl
        precstend=(colprecs-precstend)/ggr*factor_xy * idt_gl

! don't use CRM tendencies from two crm top levels
!        sltend(ptop:ptop+1) = 0.
!        qltend(ptop:ptop+1) = 0.
!        qcltend(ptop:ptop+1) = 0.
!        qiltend(ptop:ptop+1) = 0.
!-------------------------------------------------------------
! 
! Save the last step to the permanent core:
        
        u_crm  (1:nx,1:ny,1:nzm) = u   (1:nx,1:ny,1:nzm)
        v_crm  (1:nx,1:ny,1:nzm) = v   (1:nx,1:ny,1:nzm)
        w_crm  (1:nx,1:ny,1:nzm) = w   (1:nx,1:ny,1:nzm)
        t_crm  (1:nx,1:ny,1:nzm) = tabs(1:nx,1:ny,1:nzm)

        crm_buffer(1:nx,1:ny,1:nzm,5:4+nmicro_fields) = micro_field(1:nx,1:ny,1:nzm,1:nmicro_fields)

        pxtp_crm(1:nx,1:ny,1:nzm,1:ntracers) = tracer(1:nx,1:ny,1:nzm,1:ntracers)   !!! um_hr_20141217 
        
        qv_crm(1:nx,1:ny,1:nzm)  = qv(1:nx,1:ny,1:nzm)  !!! added qv_crm um_hr_20140813
        qc_crm(1:nx,1:ny,1:nzm)  = qcl(1:nx,1:ny,1:nzm)
        qi_crm(1:nx,1:ny,1:nzm)  = qci(1:nx,1:ny,1:nzm)
        qpc_crm(1:nx,1:ny,1:nzm) = qpl(1:nx,1:ny,1:nzm)
        qpi_crm(1:nx,1:ny,1:nzm) = qpi(1:nx,1:ny,1:nzm)
        if (domicro_sam1mom) then
           crm_buffer(1:nx,1:ny,1:nzm,7) = qn(1:nx,1:ny,1:nzm)
        ! um_hr_20190321+
        elseif (domicro_m2005) then
           if (dopredictNc) then
              do k=1,nzm
                 cdnc_crm(1:nx,1:ny,k) = micro_field(1:nx,1:ny,k,incl)*rho(k) ! convert to number concentration
              end do
            end if
         end if
        ! um_hr_20190321-

        z0m = z0 
        taux_crm = taux0 / dble(nstop)
        tauy_crm = tauy0 / dble(nstop)

!---------------------------------------------------------------
!
!  Diagnostics:
        
        cld_crm = min(1._r8,cld_crm/float(nstop)) ! um_hr_20150123
        cld     = min(1._r8,cld/float(nstop)*factor_xy)
        cldtop  = min(1._r8,cldtop/float(nstop)*factor_xy)
        cldbot  = min(1._r8,cldbot/float(nstop)*factor_xy)
        
        !!! um_hr_20140918+
        conv_counter = min(1._r8,conv_counter/float(nstop)*factor_xy)
        convtoppdf   = min(1._r8,convtoppdf/float(nstop)*factor_xy)
        convbotpdf   = min(1._r8,convbotpdf/float(nstop)*factor_xy)

        massfu = massfu / float(nstop) * factor_xy
        massfd = massfd / float(nstop) * factor_xy
        u_entr = u_entr / float(nstop) * factor_xy
        u_detr = u_detr / float(nstop) * factor_xy
        d_entr = d_entr / float(nstop) * factor_xy
        d_detr = d_detr / float(nstop) * factor_xy
        !!! um_hr_20140918-

        gicewp(:) = gicewp*pdel(:)*1000./ggr/float(nstop)*factor_xy
        gliqwp(:) = gliqwp*pdel(:)*1000./ggr/float(nstop)*factor_xy
        mcup      = mcup / float(nstop) * factor_xy
        mcdn      = mcdn / float(nstop) * factor_xy
        mcuup     = mcuup / float(nstop) * factor_xy
        mcudn     = mcudn / float(nstop) * factor_xy
        mc        = mcup + mcdn + mcuup + mcudn
             
        crm_qc = crm_qc / float(nstop) * factor_xy
        crm_qi = crm_qi / float(nstop) * factor_xy
        crm_qs = crm_qs / float(nstop) * factor_xy
        crm_qg = crm_qg / float(nstop) * factor_xy
        crm_qr = crm_qr / float(nstop) * factor_xy
        cdnc   = cdnc   / float(nstop) * factor_xy

        precc  = 0.
        precl  = 0.
        precsc = 0.
        precsl = 0.
        do j=1,ny
         do i=1,nx
          precsfc(i,j) = precsfc(i,j)*dz/dt/dble(nstop)
          precssfc(i,j) = precssfc(i,j)*dz/dt/dble(nstop)
          if(precsfc(i,j).gt.10./86400.) then
             precc = precc + precsfc(i,j)
             precsc = precsc + precssfc(i,j)
          else
             precl = precl + precsfc(i,j)
             precsl = precsl + precssfc(i,j)
          end if
         end do
        end do

        prec_crm = precsfc/1000.
        precc = precc*factor_xy/1000.
        precl = precl*factor_xy/1000.
        precsc = precsc*factor_xy/1000.
        precsl = precsl*factor_xy/1000.
        
        cltot = cltot *factor_xy/nstop
        clhgh = clhgh *factor_xy/nstop
        clmed = clmed *factor_xy/nstop
        cllow = cllow *factor_xy/nstop

!-------------------------------------------------------------
!       Fluxes and other stat:
!-------------------------------------------------------------
        do k=1,nzm
          u2z = 0.
          v2z = 0.
          w2z = 0.
          do j=1,ny
           do i=1,nx
             u2z = u2z+(u(i,j,k)-u0(k))**2
             v2z = v2z+(v(i,j,k)-v0(k))**2
             w2z = w2z+0.5*(w(i,j,k+1)**2+w(i,j,k)**2)
           end do
          end do
          tmp1 = dz/rhow(k)
          tmp2 = tmp1/dtn
          mkwsb(k,:) = mkwsb(k,:) * tmp1*rhow(k) * factor_xy/nstop
          mkwle(k,:) = mkwle(k,:) * tmp2*rhow(k) * factor_xy/nstop
          mkadv(k,:) = mkadv(k,:) * factor_xy*idt_gl
          mkdiff(k,:) = mkdiff(k,:) * factor_xy*idt_gl
          qpsrc(k) = qpsrc(k) * factor_xy*idt_gl
          qpevp(k) = qpevp(k) * factor_xy*idt_gl
          qpfall(k) = qpfall(k) * factor_xy*idt_gl
          precflux(k) = precflux(k) * factor_xy*dz/dt/nstop
          l = plev-k+1
          flux_u(l) = (uwle(k) + uwsb(k))*tmp1*factor_xy/nstop
          flux_v(l) = (vwle(k) + vwsb(k))*tmp1*factor_xy/nstop

          if (domicro_sam1mom) then
             flux_qt(l) = mkwle(k,1) + mkwsb(k,1)
             fluxsgs_qt(l) =  mkwsb(k,1)
             flux_qp(l) = mkwle(k,2) + mkwsb(k,2)
             qt_trans(l) = mkadv(k,1) + mkdiff(k,1)
             qp_trans(l) = mkadv(k,2) + mkdiff(k,2)
          elseif (domicro_m2005) then
             flux_qt(l) = mkwle(k,1) + mkwsb(k,1) +  &
                  mkwle(k,iqcl) + mkwsb(k,iqcl) + mkwle(k,iqci) + mkwsb(k,iqci)
             fluxsgs_qt(l) =  mkwsb(k,1) + mkwsb(k,iqcl) + mkwsb(k,iqci)
             flux_qp(l) = mkwle(k,iqr) + mkwsb(k,iqr) +  &
                  mkwle(k,iqs) + mkwsb(k,iqs) + mkwle(k,iqg) + mkwsb(k,iqg)
             qt_trans(l) = mkadv(k,1) + mkadv(k,iqcl) + mkadv(k,iqci) + &
                  mkdiff(k,1) + mkdiff(k,iqcl) + mkdiff(k,iqci) 
             qp_trans(l) = mkadv(k,iqr) + mkadv(k,iqs) + mkadv(k,iqg) + &
                  mkdiff(k,iqr) + mkdiff(k,iqs) + mkdiff(k,iqg) 
          end if

          tkesgsz(l) = rho(k)*sum(tke(1:nx,1:ny,k))*factor_xy
          tkez(l)    = rho(k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + tkesgsz(l)
          pflx(l)    = precflux(k)/1000.
          qp_fall(l) = qpfall(k)
          qp_evp(l)  = qpevp(k)
          qp_src(l)  = qpsrc(k)
          qt_ls(l)   = qtend(k)
          t_ls(l)    = ttend(k)

!          radlwup0(k)=radlwup0(k)* factor_xy/nstop
!          radlwdn0(k)=radlwdn0(k)* factor_xy/nstop
!          radqrlw0(k)=radqrlw0(k)* factor_xy/nstop
!          radswup0(k)=radswup0(k)* factor_xy/nstop
!          radswdn0(k)=radswdn0(k)* factor_xy/nstop
!          radqrsw0(k)=radqrsw0(k)* factor_xy/nstop
!          lwnsxy = lwnsxy * factor_xy/nstop
!          swnsxy = swnsxy * factor_xy/nstop
!          lwntxy = lwntxy * factor_xy/nstop
!          swntxy = swntxy * factor_xy/nstop
!          lwnscxy = lwnscxy * factor_xy/nstop
!          swnscxy = swnscxy * factor_xy/nstop
!          lwntcxy = lwntcxy * factor_xy/nstop
!          swntcxy = swntcxy * factor_xy/nstop
!          solinxy = solinxy * factor_xy/nstop
!          lwdsxy = lwdsxy * factor_xy/nstop
!          swdsxy = swdsxy * factor_xy/nstop
        
        end do
        
        timing_factor = timing_factor / nstop
        
end subroutine crm

! ====================================================================

  SUBROUTINE t_startf(event)
    
    character(len=*), intent(in) :: event  ! performance timer event name

    RETURN

  END SUBROUTINE T_STARTF

! ====================================================================

  SUBROUTINE t_stopf(event)

    character(len=*), intent(in) :: event  ! performance timer event name
    
    RETURN
    
  END SUBROUTINE T_STOPF

! ====================================================================
