MODULE MESSY_CONVECT_ZHANG

! This module contains the subroutines for the Zhang/Hack Convection
! Scheme. They are called from the zhang_cumastr from 
! MESSY_convection or other subroutines within this module

! Author:  H.Tost,   MPICH, May 2004

! changes are mostly due to adaption to f90 Code and transfer of grid
! information as well as elimination of include statements and common 
! blocks

! original code from CCM group from MATCH delopers

USE messy_convect_zhang_param

IMPLICIT NONE

PRIVATE
PUBLIC :: conv_ccm
PUBLIC :: cmfmca
PUBLIC :: convtran
PUBLIC :: conv_ccm_pjr

INTRINSIC :: ABS, MAX, MIN, LOG, EXP, NINT, REAL

!==========================================================================
CONTAINS
!==========================================================================


! $Id: conv_zhang.F,v 1.2.2.4 1999/04/03 18:54:30 eaton Exp $
  subroutine conv_ccm(plev    ,plevp   ,plon    ,plond   ,          &  ! nlev, nlevp1, kproma, kbdim
                      t       ,qh      ,pcpc    ,jctop   ,jcbot   , &
                      pblh    ,zm      ,geos    ,zi      ,qtg     , &
                      ttg     ,pap     ,paph    ,dpp     ,ts      , &
                      delt    ,mcon    ,cme                       , &
                      tpert   ,dlf     ,pflx    ,zdu              , &
!+bee
                      cmfdqr  ,                                     &
!-bee
!+scyc
                      mu2     ,md2     ,du2     ,eu2     ,ed2     , &
                      wdp     ,dsubcld ,jt      ,maxg    ,ideep   , &
                      lengath ,ql      ,cape,    lwc,     iwc,      &
                      rform   ,sform )
!-scyc

!- MAIN DRIVER FOR ZHANG-MCFARLANE CONVECTION SCHEME ------

!----------------------------------------------------------------------
!This is contributed code not fully standardized by the CCM core group.
!All variables have been typed, where most are identified in comments
!The current procedure will be reimplemented in a subsequent version 
!of the CCM where it will include a more straightforward formulation 
!and will make use of the standard CCM nomenclature
!----------------------------------------------------------------------

!$Id: conv_zhang.F,v 1.2.2.4 1999/04/03 18:54:30 eaton Exp $

!******************************************************************

!******************************************************************

!same as conv.up except saturation vapor pressure is calculated
!in a different way.

!jul 17/92 - guang jun zhang, m.lazare. calls new buoyan, q1q2
!            and moment (several new work fields added for later).

!nov 21/91 - m.lazare. like previous conv except calls new
!                      clpdprp.
!feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.
!            previous version conv.
!performs deep convective adjustment based on mass-flux closure
!algorithm.

!************************ index of variables **********************
! i      => input arrays.
! i/o    => input/output arrays.
! w      => work arrays.
! wg     => work arrays operating only on gathered points.
! ic     => input data constants.
! c      => data constants pertaining to subroutine itself.

! wg * alpha    array of vertical differencing used (=1. for upstream).
! wg * betad    downward mass flux at cloud base.
! wg * betau    upward   mass flux at cloud base.
! w  * cape     convective available potential energy.
! wg * capeg    gathered convective available potential energy.
! c  * capelmt  threshold value for cape for deep convection.
! ic  * cpres    specific heat at constant pressure in j/kg-degk.
! i  * dpp      local sigma half-level thickness (i.e. dshj).
! ic  * delt     length of model time-step in seconds.
! wg * dp       layer thickness in mbs (between upper/lower interface).
! wg * dqdt     mixing ratio tendency at gathered points.
! wg * dsdt     dry static energy ("temp") tendency at gathered points.
! wg * dudt     u-wind tendency at gathered points.
! wg * dvdt     v-wind tendency at gathered points.
! wg * dsubcld  layer thickness in mbs between lcl and maxi.
! ic  * grav     acceleration due to gravity in m/sec2.
! wg * du       detrainment in updraft. specified in mid-layer
! wg * ed       entrainment in downdraft.
! wg * eu       entrainment in updraft.
! wg * hmn      moist static energy.
! wg * hsat     saturated moist static energy.
! w  * ideep    holds position of gathered points vs longitude index.
! ic  * plev     number of model levels.
! ic  * ilg      lon+2 = size of grid slice.
! wg * j0       detrainment initiation level index.
! wg * jd       downdraft   initiation level index.
! ic  * jlat     gaussian latitude index.
! ic  * jlatpr   gaussian latitude index for printing grids (if needed).
! wg * jt       top  level index of deep cumulus convection.
! ic  * kount    current model timestep number.
! w  * lcl      base level index of deep cumulus convection.
! wg * lclg     gathered values of lcl.
! w  * lel      index of highest theoretical convective plume.
! wg * lelg     gathered values of lel.
! w  * lon      index of onset level for deep convection.
! wg * long     gathered values of lon.
! ic  * lev      plev+1.
! w  * maxi     index of level with largest moist static energy.
! wg * maxg     gathered values of maxi.
! wg * mb       cloud base mass flux.
! wg * mc       net upward (scaled by mb) cloud mass flux.
! wg * md       downward cloud mass flux (positive up).
! wg * mu       upward   cloud mass flux (positive up). specified 
!               at interface
! ic  * msg      number of missing moisture levels at the top of model.
! c  * kups     number of points undergoing deep convection.
! w  * p        grid slice of ambient mid-layer pressure in mbs.
! i  * pblt     row of pbl top indices.
! i/o * pcp      row of precipitable water in metres.
! w  * pcpdh    scaled surface pressure.
! w  * pf       grid slice of ambient interface pressure in mbs.
! wg * pg       grid slice of gathered values of p.
! i  * pressg   row of surface pressure in pa.
! w  * q        grid slice of mixing ratio.
! wg * qd       grid slice of mixing ratio in downdraft.
! wg * qdb      row of qd at cloud base.
! wg * qg       grid slice of gathered values of q.
! i/o * qh       grid slice of specific humidity.
! w  * qh0      grid slice of initial specific humidity.
! wg * qhat     grid slice of upper interface mixing ratio.
! wg * ql       grid slice of cloud liquid water.
! wg * qs       grid slice of saturation mixing ratio.
! w  * qstp     grid slice of parcel temp. saturation mixing ratio.
! wg * qstpg    grid slice of gathered values of qstp.
! wg * qu       grid slice of mixing ratio in updraft.
! ic  * rgas     dry air gas constant.
! wg * rl       latent heat of vaporization.
! w  * s        grid slice of scaled dry static energy (t+gz/cp).
! wg * sd       grid slice of dry static energy in downdraft.
! wg * sdb      row of sd at cloud base.
! wg * sg       grid slice of gathered values of s.
! wg * shat     grid slice of upper interface dry static energy.
! i  * shbj     grid slice of local bottom interface sigma values.
! i  * shj      grid slice of local half-level sigma values.
! i  * shtj     row of local top interfaces of first level.
! wg * su       grid slice of dry static energy in updraft.
! wg * sumde    row of vertically-integrated moist static energy 
!               change.
! wg * sumdq    row of vertically-integrated scaled mixing ratio 
!               change.
! wg * sumdt    row of vertically-integrated dry static energy change.
! wg * sumq     row of vertically-integrated mixing ratio change.
! i/o * t        grid slice of temperature at mid-layer.
! o  * jctop    row of top-of-deep-convection indices passed out.
! o  * jcbot    row of base of cloud indices passed out.
! w  * tf       grid slice of temperature at interface.
! wg * tg       grid slice of gathered values of t.
! w  * tl       row of parcel temperature at lcl.
! wg * tlg      grid slice of gathered values of tl.
! w  * tp       grid slice of parcel temperatures.
! wg * tpg      grid slice of gathered values of tp.
! i/o * u        grid slice of u-wind (real).
! wg * ug       grid slice of gathered values of u.
! i/o * utg      grid slice of u-wind tendency (real).
! i/o * v        grid slice of v-wind (real).
! w  * va       work array re-used by called subroutines.
! wg * vg       grid slice of gathered values of v.
! i/o * vtg      grid slice of v-wind tendency (real).
! i  * w        grid slice of diagnosed large-scale vertical velocity.
! w  * z        grid slice of ambient mid-layer height in metres.
! w  * zf       grid slice of ambient interface height in metres.
! wg * zfg      grid slice of gathered values of zf.
! wg * zg       grid slice of gathered values of z.

!----------------------------------------------------------------------
        
    IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plevp,              &    ! number of levels + 1
                 plon,               &    ! "longitude index", kproma
                 plond

!multi-level i/o fields:

!input/output arguments:

      REAL(dp) :: t(plond,plev)

      real(dp) :: qh(plond,plev,1) 
      real(dp) :: u(plond,plev) 
      real(dp) :: v(plond,plev) 
!     real(dp) :: utg(plond,plev) 
!     real(dp) :: vtg(plond,plev) 
      real(dp) :: qtg(plond,plev) 
      real(dp) :: ttg(plond,plev)

!      input arguments
      real(dp) :: pap(plond,plev) 
      real(dp) :: paph(plond,plev+1) 
      real(dp) :: dpp(plond,plev) 
      real(dp) :: zm(plond,plev) 
      real(dp) :: geos(plond) 
      real(dp) :: zi(plond,plev+1)
      real(dp) :: pblh(plond) 
      real(dp) :: zs(plond) 
      real(dp) :: tpert(plond) 

!output arguments

      real(dp) :: pcpck(plond,plev)
      real(dp) :: mcon(plond,plev) 
      real(dp) :: dlg(plond,plev)     ! gathrd version of the detraining cld h2o tend
      real(dp) :: dlf(plond,plev)     ! scattrd version of the detraining cld h2o tend
      real(dp) :: pflx(plond,plevp)   ! scattered precip flux at each level
      real(dp) :: pflxg(plond,plevp)  ! gather precip flux at each level
      real(dp) :: cug(plond,plev)     ! gathered condensation rate 
      real(dp) :: evpg(plond,plev)    ! gathered evap rate of rain in downdraft
      real(dp) :: mumax(plond) 
      real(dp) :: cme(plond,plev)
      real(dp) :: zdu(plond,plev)
!+bee
      real(dp) :: cmfdqr(plond,plev)
!-bee
!+scyc, move these vars from local storage to output so that convective
!       transports can be done in outside of conv_ccm.
      real(dp) :: mu2(plond,plev) 
      real(dp) :: eu2(plond,plev) 
      real(dp) :: du2(plond,plev) 
      real(dp) :: md2(plond,plev) 
      real(dp) :: ed2(plond,plev)
      real(dp) :: wdp(plond,plev) 
      real(dp) :: dsubcld(plond) 
      INTEGER  :: jt(plond) 
      INTEGER  :: maxg(plond)
      INTEGER  :: ideep(plond) 
      INTEGER  :: lengath
!    diagnostic field used by chem/wetdep codes
      real(dp) :: ql(plond,plev) 
!-scyc

!single-level i/o fields:

!      input arguments
      real(dp) :: ts(plond) 
      real(dp) :: pblt(plond)

!      input/output arguments:

      real(dp) :: paprc(plond) 
      real(dp) :: paprs(plond) 

!      output arguments:

      real(dp) :: jctop(plond) 
      real(dp) :: jcbot(plond)
      real(dp) :: pcpr(plond) 
      real(dp) :: pcps(plond) 
      real(dp) :: pcpc(plond)
    ! mz_ht_20070918+
      REAL(dp) :: lwc(plond,plev)
      REAL(dp) :: iwc(plond,plev)
      REAL(dp) :: rform(plond,plev)
      REAL(dp) :: sform(plond,plev)
      REAL(dp) :: lwcg(plond,plev)
      REAL(dp) :: iwcg(plond,plev)
      REAL(dp) :: rformg(plond,plev)
      REAL(dp) :: sformg(plond,plev)
      ! mz_ht_20070918-
!----------------------------------------------------------------------

!general work fields (local variables):

      real(dp) :: q(plond,plev) 
      real(dp) :: p(plond,plev) 
      real(dp) :: z(plond,plev) 
      real(dp) :: s(plond,plev) 
      real(dp) :: qh0(plond,plev) 
      real(dp) :: tp(plond,plev) 
      real(dp) :: zf(plond,plev+1) 
      real(dp) :: pf(plond,plev+1) 
      real(dp) :: qstp(plond,plev) 

      real(dp) :: cape(plond) 
      real(dp) :: tl(plond) 
      real(dp) :: sumq(plond) 
      real(dp) :: pcpdh(plond)
!     real(dp) :: sumdt(plond) 
!     real(dp) :: sumdq(plond) 
!     real(dp) :: sumde(plond)

      INTEGER  :: lcl(plond) 
      INTEGER  :: lel(plond) 
      INTEGER  :: lon(plond) 
      INTEGER  :: maxi(plond) 
      INTEGER  :: indecs(plond)
      real(dp) :: precip

!gathered work fields:

      real(dp) :: qg(plond,plev) 
      real(dp) :: tg(plond,plev) 
      real(dp) :: pg(plond,plev) 
      real(dp) :: zg(plond,plev) 
      real(dp) :: sg(plond,plev) 
      real(dp) :: tpg(plond,plev) 
      real(dp) :: zfg(plond,plev+1) 
      real(dp) :: qstpg(plond,plev) 
      real(dp) :: ug(plond,plev) 
      real(dp) :: vg(plond,plev) 
      real(dp) :: cmeg(plond,plev)

!+bee
      real(dp) :: cmfdqrg(plond,plev)
!-bee
      real(dp) :: capeg(plond) 
      real(dp) :: tlg(plond)

      INTEGER  :: lclg(plond) 
      INTEGER  :: lelg(plond) 

!work fields arising from gathered calculations.

      real(dp) :: mu(plond,plev) 
      real(dp) :: eu(plond,plev) 
      real(dp) :: dqdt(plond,plev) 
      real(dp) :: dsdt(plond,plev) 
      real(dp) :: du(plond,plev) 
      real(dp) :: md(plond,plev) 
      real(dp) :: ed(plond,plev) 
      real(dp) :: sd(plond,plev) 
      real(dp) :: qd(plond,plev) 
      real(dp) :: mc(plond,plev) 
      real(dp) :: qhat(plond,plev) 
      real(dp) :: qu(plond,plev) 
      real(dp) :: su(plond,plev) 
      real(dp) :: qs(plond,plev) 
      real(dp) :: shat(plond,plev) 
      real(dp) :: hmn(plond,plev) 
      real(dp) :: hsat(plond,plev) 
!+scyc
      real(dp) :: qlg(plond,plev)
!-scyc
      real(dp) :: dudt(plond,plev) 
      real(dp) :: dvdt(plond,plev) 
      real(dp) :: ud(plond,plev) 
      real(dp) :: vd(plond,plev)

!     real(dp) :: deltat(plond,plev) 
!     real(dp) :: deltaq(plond,plev)

      real(dp) :: betau(plond) 
      real(dp) :: betad(plond) 
      real(dp) :: mb(plond) 
      real(dp) :: totpcp(plond) 
      real(dp) :: totevp(plond)

      INTEGER  :: jlcl(plond) 
      INTEGER  :: j0(plond) 
      INTEGER  :: jd(plond)

      real(dp) :: capelmt
      real(dp) :: cpres
      real(dp) :: delt

      INTEGER  :: i
      INTEGER  :: ii
      INTEGER  :: k
      INTEGER  :: msg
      real(dp) :: psdiss
      real(dp) :: psevap
      real(dp) :: psheat
      real(dp) :: psrain
      real(dp) :: qdifr
      real(dp) :: qeff
      real(dp) :: qmin
      real(dp) :: rl
      real(dp) :: sdifr

      logical  :: momentm
!-------------------------Data statements------------------------------

      capelmt = 70._dp
      momentm =.FALSE.

!Set internal variable "msg" (convection limit) to "limcnv-1"

      msg = limcnv - 1

!initialize necessary arrays.
!zero out variables not used in ccm

      do i = 1,plon
         paprc(i) = 0.
         paprs(i) = 0.
      end do
      psdiss = 0.
      psheat = 0.
      psevap = 0.
      psrain = 0.
!     jlatpr = 32
      cpres = 1004.64_dp
      a = 21.656_dp
      b = 5418._dp
      c1 = 6.112_dp
      c2 = 17.67_dp
      c3 = 243.5_dp
      eps1 = 0.622_dp
      qmin = 1.E-20_dp
      tfreez = 273.16_dp
      rl = 2.5104E6_dp

!initialize convective tendencies

      do k = 1,plev
         do i = 1,plon
            dqdt(i,k) = 0.
            dsdt(i,k) = 0.
            dudt(i,k) = 0.
            dvdt(i,k) = 0.
!           deltaq(i,k) = qh(i,k,1)
!           deltat(i,k) = t(i,k)
            pcpck(i,k) = 0.
            pflx(i,k) = 0.
            pflxg(i,k) = 0.
            cme(i,k) = 0.
!+bee
            cmfdqr(i,k) = 0.
!-bee
            zdu(i,k) = 0.
!+scyc
            ql(i,k) = 0.
            qlg(i,k) = 0.
!-scyc
            lwcg(i,k) = 0._dp
            iwcg(i,k) = 0._dp
            rformg(i,k) = 0._dp
            sformg(i,k) = 0._dp
         end do
      end do
      do i = 1,plon
         pflx(i,plevp) = 0.
         pflxg(i,plevp) = 0.
      end do
      if (.not.momentm) then
!         do k = msg + 1,plev
        do k=1,plev
            do i = 1,plon
               u(i,k) = 0.
               v(i,k) = 0.
            end do
         end do
      end if

      do i = 1,plon
         pblt(i) = real(plev,dp)
         pcpr(i) = 0.
         pcps(i) = 0.
         dsubcld(i) = 0.
         sumq(i) = 0.
!        sumdt(i) = 0.
!        sumdq(i) = 0.
         pcpdh(i) = rgrav
         jctop(i) = real(plev,dp)
         jcbot(i) = 1._dp
      end do

!calculate local pressure (mbs) and height (m) for both interface
!and mid-layer locations.

      do i = 1,plon
         zs(i) = geos(i)*rgrav
         pf(i,plev+1) = paph(i,plev+1)*0.01
         zf(i,plev+1) = zi(i,plev+1) + zs(i)
      end do
      do k = 1,plev
         do i = 1,plon
            p(i,k) = pap(i,k)*0.01
            pf(i,k) = paph(i,k)*0.01
            z(i,k) = zm(i,k) + zs(i)
            zf(i,k) = zi(i,k) + zs(i)
         end do
      end do

      do k = plev - 1,msg + 1,-1
         do i = 1,plon
            if (abs(z(i,k)-zs(i)-pblh(i)).lt.         &
              (zf(i,k)-zf(i,k+1))*0.5) pblt(i) = real(k,dp)
         end do
      end do

!store incoming specific humidity field for subsequent calculation
!of precipitation (through change in storage).
!convert from specific humidity (bounded by qmin) to mixing ratio.
!define dry static energy (normalized by cp).

      do k = 1,plev
         do i = 1,plon
            qh0(i,k) = qh(i,k,1)
            qeff = max(qh(i,k,1),qmin)
            q(i,k) = qeff
            s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
            tp(i,k)=0.0
            shat(i,k) = s(i,k)
            qhat(i,k) = q(i,k)
            wdp(i,k) = dpp(i,k)*0.01
            qg(i,k) = q(i,k)
            tg(i,k) = t(i,k)
            pg(i,k) = p(i,k)
            zg(i,k) = z(i,k)
            sg(i,k) = s(i,k)
            tpg(i,k) = tp(i,k)
            zfg(i,k) = zf(i,k)
            qstpg(i,k) = q(i,k)
            ug(i,k) = u(i,k)
            vg(i,k) = v(i,k)
            dlg(i,k) = 0.
            dlf(i,k) = 0.
         end do
      end do
      do i = 1,plon
         zfg(i,plev+1) = zf(i,plev+1)
         capeg(i) = 0.
         lclg(i) = 1
         lelg(i) = plev
         maxg(i) = 1
         tlg(i) = 400._dp
         dsubcld(i) = 0.
         betau(i) = 0.
         betad(i) = 0.
      end do

!evaluate covective available potential energy (cape).

      call buoyan(plev    ,plon    ,plond             ,   &  !nlev, kproma, kbdim
                  q       ,t       ,p       ,z       ,pf       ,   &
                  tp      ,qstp    ,tl      ,rl      ,cape     ,   &
                  pblt    ,lcl     ,lel     ,lon     ,maxi     ,   &
                  rgas    ,grav    ,cpres   ,msg               ,   &
                  tpert   )

!determine whether grid points will undergo some deep convection
!(ideep=1) or not (ideep=0), based on values of cape,lcl,lel
!(require cape.gt. 0 and lel<lcl as minimum conditions).

         call whenfgt(plon,cape,1,capelmt,indecs,lengath)
         if (lengath.eq.0) return

         do ii=1,lengath
            i=indecs(ii)
            ideep(ii)=i
         end do
!       jyes = 0
!       jno = plon - 1 + 2
!       do il = 1,plon
!          if (cape(il).gt.capelmt) then
!             jyes = jyes + 1
!             ideep(jyes) = il
!          else
!             jno = jno - 1
!             ideep(jno) = il
!          end if
!       end do
!       lengath = jyes
!       if (lengath.eq.0) return

!obtain gathered arrays necessary for ensuing calculations.

      do k = 1,plev
         do i = 1,lengath
            wdp(i,k) = 0.01*dpp(ideep(i),k)
            qg(i,k) = q(ideep(i),k)
            tg(i,k) = t(ideep(i),k)
            pg(i,k) = p(ideep(i),k)
            zg(i,k) = z(ideep(i),k)
            sg(i,k) = s(ideep(i),k)
            tpg(i,k) = tp(ideep(i),k)
            zfg(i,k) = zf(ideep(i),k)
            qstpg(i,k) = qstp(ideep(i),k)
            ug(i,k) = u(ideep(i),k)
            vg(i,k) = v(ideep(i),k)
         end do
      end do

      do i = 1,lengath
         zfg(i,plev+1) = zf(ideep(i),plev+1)
      end do
      do i = 1,lengath
         capeg(i) = cape(ideep(i))
         lclg(i) = lcl(ideep(i))
         lelg(i) = lel(ideep(i))
         maxg(i) = maxi(ideep(i))
         tlg(i) = tl(ideep(i))
      end do

!calculate sub-cloud layer pressure "thickness" for use in
!closure and tendency routines.

      do k = msg + 1,plev
         do i = 1,lengath
            if (k.ge.maxg(i)) then
               dsubcld(i) = dsubcld(i) + wdp(i,k)
            end if
         end do
      end do

!define array of factors (alpha) which defines interfacial
!values, as well as interfacial values for (q,s) used in
!subsequent routines.

      do k = msg + 2,plev
         do i = 1,lengath
            sdifr = 0.
            qdifr = 0.
            if (sg(i,k).gt.0. .or. sg(i,k-1).gt. 0.)     &
                sdifr = abs((sg(i,k)-sg(i,k-1))/         &
                max(sg(i,k-1),sg(i,k)))
            if (qg(i,k).gt.0. .or. qg(i,k-1).gt.0.)      &
                qdifr = abs((qg(i,k)-qg(i,k-1))/         &
                max(qg(i,k-1),qg(i,k)))
            if (sdifr.gt.1.E-6_dp) then
              shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/      &
                (sg(i,k-1)-sg(i,k))
            else
              shat(i,k) = 0.5_dp* (sg(i,k)+sg(i,k-1))
            end if
            if (qdifr.gt.1.E-6_dp) then
              qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/      &
                (qg(i,k-1)-qg(i,k))
            else
              qhat(i,k) = 0.5_dp* (qg(i,k)+qg(i,k-1))
            end if
          end do
        end do

!obtain cloud properties.

      call cldprp(plev    ,plevp   ,plond                     , &  !nlev, nlevp1, kbdim
                  qg      ,tg      ,ug      ,vg      ,pg      , &
                  zg      ,sg      ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,ud      , &
                  vd      ,mc      ,qu      ,su      ,zfg     , &
                  qs      ,hmn     ,hsat    ,shat             , &
!+scyc
                  qlg     ,totpcp  ,totevp  ,cmeg    ,maxg    , &
!-scyc
                  lelg    ,jt      ,jlcl    ,maxg    ,j0      , &
                  jd      ,rl      ,lengath ,rgas             , &
                  grav    ,cpres   ,msg                       , &
                  pflxg   ,evpg    ,cug     ,mu2     ,eu2     , &
!+bee
! jr added limcnv
                  du2     ,md2     ,ed2     ,cmfdqrg,           &
                  lwcg    ,iwcg    ,rformg  ,sformg  ) !,limcnv  )
!-bee

!determine cloud base mass flux.


      do i = 1,lengath
         betad(i) = md(i,maxg(i))
         betau(i) = mu(i,maxg(i))
      end do

!convert detrainment from units of "1/m" to "1/mb".

      do k = msg + 1,plev
         do i = 1,lengath
            du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
!+bee, bug fix from /home/pjr/ccm2/omega0.10.1/fspj01/.
            cmeg(i,k) = cmeg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            cmfdqrg(i,k) = cmfdqrg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
!-bee
            evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
         end do
      end do

      call closure(plev    ,plond                              ,&  !nlev, kbdim
                   qg      ,tg      ,pg      ,sg               ,&
                   tpg     ,qu      ,su      ,mc               ,&
                   du      ,mu      ,md      ,qd      ,sd      ,&
                   qhat    ,shat    ,wdp     ,qstpg            ,&
!+scyc
                   zfg     ,qlg     ,dsubcld ,mb      ,capeg   ,&
!-scyc
                   tlg     ,lclg    ,lelg    ,jt      ,maxg    ,&
                   1       ,lengath ,rgas    ,grav    ,cpres   ,&
                   rl      ,msg     ,capelmt )

!limit cloud base mass flux to theoretical upper bound.

      do i=1,lengath
        mumax(i) = 0.
      end do
      do k=msg + 2,plev
        do i=1,lengath
          mumax(i) = max(mumax(i), mu(i,k)/wdp(i,k))
        end do
      end do
      do i=1,lengath
        if (mumax(i).gt.0.) then
          mb(i) = min(mb(i),0.5/(delt*mumax(i)))
        else
          mb(i) = 0.
        endif
      end do
      do k=msg+1,plev
        do i=1,lengath
          mu(i,k) = mu(i,k)*mb(i)
          md(i,k) = md(i,k)*mb(i)
          mc(i,k) = mc(i,k)*mb(i)
          du(i,k) = du(i,k)*mb(i)
          eu(i,k) = eu(i,k)*mb(i)
          ed(i,k) = ed(i,k)*mb(i)
          cmeg(i,k) = cmeg(i,k)*mb(i)
!+bee
!-bee
          cug(i,k) = cug(i,k)*mb(i)
          evpg(i,k) = evpg(i,k)*mb(i)
          pflxg(i,k+1) = pflxg(i,k+1)*mb(i)*100./grav
          cmfdqrg(i,k) = cmfdqrg(i,k)*mb(i)
          mu2(i,k) = mu2(i,k)*mb(i)
          md2(i,k) = md2(i,k)*mb(i)
          du2(i,k) = du2(i,k)*mb(i)
          eu2(i,k) = eu2(i,k)*mb(i)
          ed2(i,k) = ed2(i,k)*mb(i)
        end do
      end do
      do i = 1,lengath
         betau(i) = betau(i)*mb(i)
         betad(i) = betad(i)*mb(i)

!totpcp from cldprp has the dimension of kg/kg, here it is 
!converted to kg/(m^2*s), the precipitation rate

         totpcp(i) = totpcp(i)*mb(i)*100./grav
         totevp(i) = totevp(i)*mb(i)*100./grav
      end do

!compute temperature and moisture changes due to convection.

      call q1q2_pjr(plev    ,plond                              ,&  !nlev, kbdim
                    dqdt    ,dsdt                               ,&
                    qu      ,su      ,du                        ,&
                    qhat    ,shat    ,wdp     ,mu      ,md      ,&
!+scyc
                    sd      ,qd      ,qlg     ,dsubcld          ,&
!-scyc
                    jt      ,maxg    ,1       ,lengath          ,&
                    cpres   ,rl      ,msg                       ,&
                    dlg     ,evpg    ,cug     )

!compute momentum changes due to convection, if desired (i.e
!if logical switch set).

!     if(momentm)                                                   
!      then
!       call moment(dudt,dvdt,du,alpha,wdp,ed,eu,mc,md,mu,
!    1             pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,
!    2             dsubcld,maxg,jd,jt,rl,
!    3             msg,2.*delt,grav,cpres,rgas,plev,1,lengath,plond,lat)
!     endif

!+scyc, move the convective transport outside of conv_ccm.

!gather back temperature and mixing ratio.

      do k = msg + 1,plev
         do i = 1,lengath
            psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+                   &
                     dvdt(i,k)*v(ideep(i),k))*dpp(ideep(i),k)/grav

!q is updated to compute net precip, and then reset to old value.
!the last line is overwritten. so the input basi!variables, i.e.
!q, t, u and v are updated to include convective increments. 
!(5/2/95)

            q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
            t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)
            u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
            v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
            cme(ideep(i),k) = cmeg(i,k)
!+bee
            cmfdqr(ideep(i),k) = cmfdqrg(i,k)
!-bee
            zdu(ideep(i),k) = du2(i,k)
            mcon(ideep(i),k) = mc(i,k)
            qtg(ideep(i),k) = dqdt(i,k)
            ttg(ideep(i),k) = dsdt(i,k)
!           utg(ideep(i),k) = dudt(i,k)
!           vtg(ideep(i),k) = dvdt(i,k)
            dlf(ideep(i),k) = dlg(i,k)
            pflx(ideep(i),k) = pflxg(i,k)
!+scyc
            ql(ideep(i),k) = qlg(i,k)
!-scyc       
            ! mz_ht_20070918+            
            lwc(ideep(i),k) = lwcg(i,k)
            iwc(ideep(i),k) = iwcg(i,k)
            rform(ideep(i),k) = rformg(i,k)
            sform(ideep(i),k) = sformg(i,k)
            ! mz_ht_20070918-
         end do
      end do

      do i = 1,lengath
         jctop(ideep(i)) = real(jt(i),dp)
!+bee
         jcbot(ideep(i)) = real(maxg(i),dp)
!-bee
         pflx(ideep(i),plevp) = pflxg(i,plevp)
         psevap = psevap + totevp(i)
         psrain = psrain + totpcp(i)
      end do

!convert back to specific humidity from mixing ratio.
!take into account any moisture added to ensure positiveness
!of specific humidity at start of routine.

      do k = msg + 1,plev
         do i = 1,plon
            qh(i,k,1) = q(i,k)
            qh(i,k,1) = qh(i,k,1) - max((qmin-qh0(i,k)),0._dp)
         end do
      end do
      do k = plev,msg + 1,-1
         do i = 1,plon
            sumq(i) = sumq(i) - dpp(i,k)* (qh(i,k,1)-qh0(i,k))

!account for the detraining cloud water in the precip 

            sumq(i) = sumq(i) - dpp(i,k)*dlf(i,k)*2._dp*delt
            pcpck(i,k) = max(0._dp,sumq(i))
         end do
      end do

!obtain final precipitation rate.

      do i = 1,plon
!        llo1 = ts(i) .ge. tfreez

!here pcpr and pcps are in units of kg/m^2, ie. precip per
!time step

!        pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),0.),0.,llo1)
!        pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),0.),llo1)
         precip = pcpdh(i)*max(sumq(i),0._dp)
         if (ts(i) .ge. tfreez) then
           pcpr(i) = precip
           pcps(i) = 0.
         else
           pcpr(i) = 0.
           pcps(i) = precip
         end if
      end do


!accumulate precipitation, the 1000. is the density of water, so
!paprc and paprs are now in units of meters.

      do i = 1,plon
         paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
         paprs(i) = paprs(i) + pcps(i)/1000.
      end do

!convert precipitation to m/s, ie, precip rate.

      do i = 1,plon
         pcpr(i) = pcpr(i)/ (2.*delt)/1000.
         pcps(i) = pcps(i)/ (2.*delt)/1000.
         pcpc(i) = pcpr(i) + pcps(i)
         psheat = psheat + (pcps(i)+pcpr(i))*rl
      end do
      do k = msg + 1,plev
         do i = 1,plon
            pcpck(i,k) = pcpdh(i)*pcpck(i,k)/ (2.*delt)
         end do
      end do

!calculate conservation of quantities.

!      if(lat.eq.jlatpr)then
!       do l=msg+1,plev
!       do i=1,lengath
!         sumdq(i) = sumdq(i) + 2.*delt*(rl/cpres)*dpp(ideep(i),l)*
!    1                          dqdt(i,l)
!         sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
!       end do
!       end do
!       print *,'sumdq,sumdt,sumde in convection subroutine########'
!       do i=1,lengath
!         sumde(i) = sumdt(i) + sumdq(i)
!         write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
!       end do
!       print *,'sumdq,sumdt,sumde ... all points'
!     do i=1,plon
!        sumdq(i) = 0.0
!        sumdt(i) = 0.0
!     end do
!     do l=msg+1,plev
!     do i=1,plon
!       deltaq(i,l) = qh(i,l,1) - deltaq(i,l)
!       deltat(i,l) = t (i,l) - deltat(i,l)
!       sumdq(i) = sumdq(i) + (rl/cpres)*dpp(i,l)*deltaq(i,l)
!       sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
!     end do
!     end do
!     do i=1,plon
!       sumde(i) = sumdt(i) + sumdq(i)
!     end do
!       write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,plon)
! 901   format(1x,3e20.12, i10, i10)
! 902   format(1x,i10, 3e20.12)
!     endif

      return
    end subroutine conv_ccm

!========================================================================
!------------------------subroutines for conv_ccm start -----------------

    subroutine buoyan(plev    ,plon    ,plond                     ,&  !nlev, nlevp1, kproma, kbdim
                      q       ,t       ,p       ,z       ,pf      ,&
                      tp      ,qstp    ,tl      ,rl      ,cape    ,&
                      pblt    ,lcl     ,lel     ,lon     ,mx      ,&
                      rd      ,grav    ,cp      ,msg     ,         &
                      tpert     )
!----------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.

! the documentation has been enhanced to the degree that we are able

! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
!----------------------------------------------------------------------
!  $Id: conv_zhang.F,v 1.2.2.4 1999/04/03 18:54:30 eaton Exp $
!----------------------------------------------------------------------

!----------------------------------------------------------------------

! jul 14/92 - guang jun zhang, m.lazare, n.mcfarlane.  as in
!             previous version buoyan except remove pathalogical
!             cases of "zig-zags" in profiles where lel defined
!             far too high (by use of lelten array, which assumes
!             a maximum of five such crossing points).
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
!             version buoyan.
  IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plon,               &    ! "longitude index", kproma
                 plond                    ! number of longitudes, kbdim
                 
! input arguments

      REAL(dp) :: q(plond,plev)        ! spec. humidity
      REAL(dp) :: t(plond,plev)        ! temperature
      REAL(dp) :: p(plond,plev)        ! pressure
      REAL(dp) :: z(plond,plev)        ! height
      REAL(dp) :: pf(plond,plev+1)     ! pressure at interfaces
      REAL(dp) :: pblt(plond)          ! index of pbl depth
      REAL(dp) :: tpert(plond)         ! perturbation temperature by pbl processes

! output arguments

      REAL(dp) :: tp(plond,plev)       ! parcel temperature
      REAL(dp) :: qstp(plond,plev)     ! saturation mixing ratio of parcel
      REAL(dp) :: tl(plond)            ! parcel temperature at lcl
      REAL(dp) :: cape(plond)          ! convective aval. pot. energy.
      INTEGER  :: lcl(plond)        ! 
      INTEGER  :: lel(plond)        ! 
      INTEGER  :: lon(plond)        ! level of onset of deep convection
      INTEGER  :: mx(plond)         ! level of max moist static energy

!-------------------------Local Variables------------------------------

      REAL(dp) :: capeten(plond,5)     ! provisional value of cape
      REAL(dp) :: tv(plond,plev)       ! 
      REAL(dp) :: tpv(plond,plev)      ! 
      REAL(dp) :: buoy(plond,plev)

      REAL(dp) :: a1(plond) 
      REAL(dp) :: a2(plond) 
      REAL(dp) :: estp(plond) 
      REAL(dp) :: pl(plond) 
      REAL(dp) :: plexp(plond) 
      REAL(dp) :: hmax(plond) 
      REAL(dp) :: hmn(plond) 
      REAL(dp) :: y(plond)

      logical  :: plge600(plond) 
      INTEGER  :: knt(plond) 
      INTEGER  :: lelten(plond,5)

      REAL(dp) :: cp
      REAL(dp) :: e
      REAL(dp) :: grav

      INTEGER  :: i
      INTEGER  :: k
      INTEGER  :: msg
      INTEGER  :: n

      REAL(dp) :: rd
      REAL(dp) :: rl
#ifdef PERGRO
      REAL(dp) :: rhd
#endif
!----------------------------------------------------------------------

      qstp = 0._dp

      do n = 1,5
        do i = 1,plon
          lelten(i,n) = plev
          capeten(i,n) = 0.
        end do
      end do

      do i = 1,plon
        lon(i) = plev
        knt(i) = 0
        lel(i) = plev
        mx(i) = lon(i)
        cape(i) = 0.
        hmax(i) = 0.
      end do

! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.

#ifdef PERGRO
      do k = plev,msg + 1,-1
        do i = 1,plon
          hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)

! Reset max moist static energy level when relative difference exceeds 1.e-4

          rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
          if (k.ge.nint(pblt(i)) .and. k.le.lon(i) .and.          &
             rhd.gt.-1.e-4) then
            hmax(i) = hmn(i)
            mx(i) = k
          end if
        end do
      end do
#else
      do k = plev,msg + 1,-1
        do i = 1,plon
          hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
          if (k.ge.nint(pblt(i)) .and. k.le.lon(i) .and.          & 
             hmn(i).gt.hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
          end if
        end do
      end do
#endif

      do i = 1,plon
        lcl(i) = mx(i)
        e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
        tl(i) = 2840._dp/ (3.5_dp*log(t(i,mx(i)))-log(e)-4.805_dp) + 55._dp
        if (tl(i).lt.t(i,mx(i))) then
          plexp(i) = (1._dp/ (0.2854_dp* (1._dp-0.28_dp*q(i,mx(i)))))
          pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)

        else
          tl(i) = t(i,mx(i))
          pl(i) = p(i,mx(i))
        end if
      end do

! calculate lifting condensation level (lcl).

      do k = plev,msg + 2,-1
        do i = 1,plon
          if (k.le.mx(i) .and. (p(i,k).gt.pl(i).and.           &
              p(i,k-1).le.pl(i))) then
            lcl(i) = k - 1
          end if
        end do
      end do

! if lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).

      do i = 1,plon
        plge600(i) = pl(i).ge.600._dp
      end do

! initialize parcel properties in sub-cloud layer below lcl.

      do k = plev,msg + 1,-1
        do i=1,plon
          if (k.gt.lcl(i) .and. k.le.mx(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._dp+1.608_dp*q(i,k))/ (1._dp+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**         &
                      (0.2854_dp* (1._dp-0.28_dp*q(i,mx(i))))

! buoyancy is increased by 0.5 k as in tiedtke

!jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
!jjh     1                     (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))*                     &
                       (1._dp+1.608_dp*q(i,mx(i)))/ (1._dp+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_dp
          end if
        end do
      end do

! define parcel properties at lcl (i.e. level immediately above pl).

      do k = plev,msg + 1,-1
        do i=1,plon
          if (k.eq.lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._dp+1.608_dp*q(i,k))/ (1._dp+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = tl(i)* (p(i,k)/pl(i))**                &
                      (0.2854_dp* (1._dp-0.28_dp*qstp(i,k)))
!              estp(i)  =exp(a-b/tp(i,k))
! use of different formulas for est has about 1 g/kg difference
! in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
! above giving larger qs.
!
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/            &
                      ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1._dp+qstp(i,k)/eps1)*rl*   &
                    eps1/ (rd*tp(i,k)**2)
            a2(i) = 0.5_dp* (qstp(i,k)* (1._dp+2._dp/eps1*qstp(i,k))*      &
                    (1._dp+qstp(i,k)/eps1)*eps1**2*rl*rl/           &
                    (rd**2*tp(i,k)**4)-qstp(i,k)*                &
                    (1._dp+qstp(i,k)/eps1)*2.*eps1*rl/              &
                    (rd*tp(i,k)**3))
            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = q(i,mx(i)) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/             &  
                    ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))

! buoyancy is increased by 0.5 k in cape calculation.
! dec. 9, 1994
!jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))

            tpv(i,k) = (tp(i,k)+tpert(i))* (1._dp+1.608_dp*qstp(i,k))/ &
                       (1._dp+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_dp
          end if
        end do
      end do

! main buoyancy calculation.

      do k = plev - 1,msg + 1,-1
        do i=1,plon
          if (k.lt.lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1._dp+1.608_dp*q(i,k))/ (1._dp+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**             &
                      (0.2854_dp* (1._dp-0.28_dp*qstp(i,k)))
!          estp(i) = exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/             &
                      ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1._dp+qstp(i,k)/eps1)*rl*   &
                    eps1/ (rd*tp(i,k)**2)
            a2(i) = 0.5_dp* (qstp(i,k)* (1._dp+2._dp/eps1*qstp(i,k))*      &
                    (1._dp+qstp(i,k)/eps1)*eps1**2*rl*rl/           &
                    (rd**2*tp(i,k)**4)-qstp(i,k)*                &
                    (1._dp+qstp(i,k)/eps1)*2.*eps1*rl/              &
                    (rd*tp(i,k)**3))
            a1(i) = 1._dp/a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/             &
                      ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
!jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
!jt            (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))* (1._dp+1.608_dp*qstp(i,k))/ &
                       (1._dp+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5_dp
          end if
        end do
      end do

      do k = msg + 2,plev
        do i = 1,plon
          if (k.lt.lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1).gt.0. .and. buoy(i,k).le.0.) then
              knt(i) = min(5,knt(i) + 1)
              lelten(i,knt(i)) = k
            end if
          end if
        end do
      end do

! calculate convective available potential energy (cape).

      do n = 1,5
        do k = msg + 1,plev
          do i = 1,plon
            if (plge600(i) .and.k.le.mx(i) .and.k.gt.lelten(i,n)) then
              capeten(i,n) = capeten(i,n) +                     &    
                             rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
          end do
        end do
      end do

! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995

      do n = 1,5
        do i = 1,plon
          if (capeten(i,n).gt.cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
          end if
        end do
      end do

! put lower bound on cape for diagnostic purposes.

      do i = 1,plon
        cape(i) = max(cape(i), 0._dp)
      end do

      return
    end subroutine buoyan
!===========================================================================================

    subroutine convtran(q     ,ncnst  ,mu     ,md      ,du     ,eu,   &
                        ed    ,wdp    ,jt     ,mx              ,      &
                        ideep ,il1g   ,il2g   ,                       &
                        delt  ,fracis ,plond   ,plev)
!-----------------------------------------------------------------------

! Convective transport of trace species

! Note that we are still assuming that the tracers are in a moist mixing ratio
! this will change soon

!-------------------------Code History----------------------------------

! Original version:  P. Rasch, Jan 1996 
! Standardized:      L. Buja,  Feb 1996
! Reviewed:          P. Rasch, Feb 1996      
 
!-----------------------------------------------------------------------
      implicit none

!-----------------------------Arguments---------------------------------
 
! Input

! +rvk moved here
      INTEGER, INTENT(IN) :: plond, plev   ! number of columns, levels
      
      integer  :: ncnst             ! number of tracers to transport

      REAL(dp) :: mu(plond,plev)       ! Mass flux up
      REAL(dp) :: md(plond,plev)       ! Mass flux down
      REAL(dp) :: du(plond,plev)       ! Mass detraining from updraft
      REAL(dp) :: eu(plond,plev)       ! Mass entraining from updraft
      REAL(dp) :: ed(plond,plev)       ! Mass entraining from downdraft
      REAL(dp) :: wdp(plond,plev)      ! Delta pressure between interfaces
      REAL(dp) :: fracis(plond,plev,ncnst) ! fraction of tracer that is insoluble

      integer  :: jt(plond)         ! Index of cloud top for each column
      integer  :: mx(plond)         ! Index of cloud top for each column
      integer  :: ideep(plond)      ! Gathering array
      integer  :: il1g              ! Gathered min lon indices over which to operate
      integer  :: il2g              ! Gathered max lon indices over which to operate



      REAL(dp) :: delt                 ! Time step

! input/output

      REAL(dp) :: q(plond,plev,ncnst)  ! Tracer array including moisture

!--------------------------Local Variables------------------------------

      integer  :: i                 ! Work index
      integer  :: k                 ! Work index
      integer  :: kbm               ! Highest altitude index of cloud base
      integer  :: kk                ! Work index
      integer  :: kkp1              ! Work index
      integer  :: km1               ! Work index
      integer  :: kp1               ! Work index
      integer  :: ktm               ! Highest altitude index of cloud top
      integer  :: m                 ! Work index

      REAL(dp) :: cabv                 ! Mix ratio of constituent above
      REAL(dp) :: cbel                 ! Mix ratio of constituent below
      REAL(dp) :: cdifr                ! Normalized diff between cabv and cbel
      REAL(dp) :: chat(plond,plev)     ! Mix ratio in env at interfaces
      REAL(dp) :: cond(plond,plev)     ! Mix ratio in downdraft at interfaces
      REAL(dp) :: const(plond,plev)    ! Gathered tracer array 
      REAL(dp) :: fisg(plond,plev)     ! gathered insoluble fraction of tracer
      REAL(dp) :: conu(plond,plev)     ! Mix ratio in updraft at interfaces
      REAL(dp) :: dcondt(plond,plev)   ! Gathered tend array 
      REAL(dp) :: small                ! A small number
      REAL(dp) :: mbsth                ! Threshold for mass fluxes
      REAL(dp) :: mupdudp              ! A work variable
      REAL(dp) :: minc                 ! A work variable
      REAL(dp) :: maxc                 ! A work variable
      REAL(dp) :: qn                   ! A work variable
      REAL(dp) :: fluxin               ! A work variable
      REAL(dp) :: fluxout              ! A work variable
      REAL(dp) :: netflux              ! A work variable

!-----------------------------------------------------------------------

      small = 1.e-36_dp
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
      mbsth = 1.e-15_dp

! Find the highest level top and bottom levels of convection
      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do

! Loop ever each constituent
      do m = 1, ncnst

! Gather up the constituent and set tend to zero
         do k = 1,plev
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,plev
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc.lt.1.0e-20_dp) then
                  cdifr = 0.
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr.gt.1.E-6_dp) then
                  cabv = max(const(i,km1),maxc*1.e-12)
                  cbel = max(const(i,k),maxc*1.e-12)
                  chat(i,k) = log(cabv/cbel)            &
                                /(cabv-cbel)            &
                                *cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0.

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = plev 
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*wdp(i,kk)
            if (mupdudp.gt.mbsth) then
               conu(i,kk) = (                                          &
                             +eu(i,kk)*fisg(i,kk)*const(i,kk)*wdp(i,kk) &
                             )/mupdudp
            endif
            if (md(i,k).lt.-mbsth) then
               cond(i,k) =  (                                             &
                            -ed(i,km1)*fisg(i,km1)*const(i,km1)*wdp(i,km1) &
                             )/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = plev-1,1,-1
            kkp1 = min(plev,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*wdp(i,kk)
               if (mupdudp.gt.mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)                    &
                               +eu(i,kk)*fisg(i,kk)*const(i,kk)*wdp(i,kk)     &
                               )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,plev
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k).lt.-mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)                       &
                                 -ed(i,km1)*fisg(i,km1)*const(i,km1)          &
                                 *wdp(i,km1))/md(i,k)
               endif
            end do
         end do


         do k = ktm,plev
            km1 = max(1,k-1)
            kp1 = min(plev,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =                                         &
!                        +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))      &
!                          -mu(i,k)*   (conu(i,k)-chat(i,k))          &
!                          +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))      &
!                          -md(i,k)*   (cond(i,k)-chat(i,k))          &
!                          )/wdp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)         &
!                       -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)     &
!                       -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)                           &
                       + mu(i,k)*min(chat(i,k),const(i,km1))             &
                       -(md(i,k)  *cond(i,k)                             &
                       + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k)                               &
                        +mu(i,kp1)*min(chat(i,kp1),const(i,k))           &
                       -(md(i,kp1)*cond(i,kp1)                           &
                       + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                  netflux = 0.
               endif
               dcondt(i,k) = netflux/wdp(i,k)
            end do
         end do

         do k = kbm,plev             
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k.eq.mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*           &
!                    (-mu(i,k)*(conu(i,k)-chat(i,k))        &
!                     -md(i,k)*(cond(i,k)-chat(i,k))        &
!                    )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)

! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1))       &
                         - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k)                         &
                         - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                     netflux = 0.
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/wdp(i,k)
               else if (k.gt.mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0.
               end if
            end do
         end do

! Update and scatter data back to full arrays

         do k = 1,plev
            kp1 = min(plev,k+1)
            do i = il1g,il2g
               qn = const(i,k)+dcondt(i,k)*2.*delt
               q(ideep(i),k,m) = qn

            end do
         end do
         
      end do

      return
    end subroutine convtran

!===============================================================================

    subroutine cldprp(plev    ,plevp   ,plond            ,          &  !nlev, nlevp1, kbdim
                      q       ,t       ,u       ,v       ,p       , &
                      z       ,s       ,mu      ,eu      ,du      , &
                      md      ,ed      ,sd      ,qd      ,ud      , &
                      vd      ,mc      ,qu      ,su      ,zf      , &
                      qst     ,hmn     ,hsat    ,shat             , &
                      ql      ,totpcp  ,totevp  ,cmeg    ,jb      , &
                      lel     ,jt      ,jlcl    ,mx      ,j0      , &
                      jd      ,rl      ,il2g    ,rd               , &
                      grav    ,cp      ,msg                       , &
                      pflx    ,evp     ,cu      ,mu2     ,eu2     , &
!+bee
!jr added limcnv
                      du2     ,md2     ,ed2     ,cmfdqr           , &
                      lwc     ,iwc     ,rform   ,sform   )!  ,limcnv  )
!-bee
! -----------------------------------------------------------------------------
!This is contributed code not fully standardized by the CCM core group.

!this code is very much rougher than virtually anything else in the CCM
!there are debug statements left strewn about and code segments disabled
!these are to facilitate future development. We expect to release a
!cleaner code in a future release

!the documentation has been enhanced to the degree that we are able

!Original version:  G. Zhang and collaborators
!Standardized:      Core group staff, 1994 and 195
!Reviewed:          P. Rasch, April 1996

!**** PLEASE NOTE ****

!we are aware of a specific problem in this code 
!(identified by the string ---> PROBLEM ONE)
!during the calculation of the updraft cloud properties,
!rather than adding a perturbation to the updraft temperature of 
!half a degree, (there was an inadvertant addition of cp*0.5) degrees
!or about 500 degrees. (This problem was in the code prior to its 
!contribution to the NCAR effort)

!Fortunately, the erroneous values
!are overwritten later in the code. The problem is quite subtle.
!The erroneous values would persist between cloud base and the lifting 
!condensation level. The addition of the very high perturbation to the updraft
!temperature causes the saturation mixing ratio to be set to zero, 
!and later the lcl to be set to one level above cloud base.
!There are therefore no levels between cloud base and the lcl. Therefore
!all erroneous values are overwritten.

!The only manifestation we are aware of with respect to this problem
!is that the lifting condensation level is constrained to be one level above
!cloud base.

!We discovered the problem after too much had been invested in
!very long integrations (in terms of computer time)
!to allow for a modification and model retuning. It is our expectation that
!this problem will be fixed in the next release of the model.

!*********** 
! ----------------------------------------------------------------------

!$Id: conv_zhang.F,v 1.2.2.4 1999/04/03 18:54:30 eaton Exp $

! ----------------------------------------------------------------------
!nov 20/92 - guang jun zhang,m.lazare. now has deeper (more
!            realistic) downdrafts.
!jul 14/92 - guang jun zhang,m.lazare. add shallow mixing
!            formulation.
!nov 21/91 - m.lazare. like previous cldprop except minimum "f"
!                      now 0.0004 instead of 0.001 (more
!                      realisti!with more deep).
!may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!            original version cldprop.
! -----------------------------------------------------------------------------
IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plevp,              &    ! number of levels + 1
                 plond                    ! number of longitudes, kbdim
!Input arguments

      REAL(dp) :: q(plond,plev)        ! spec. humidity of env
      REAL(dp) :: t(plond,plev)        ! temp of env
      REAL(dp) :: p(plond,plev)        ! pressure of env
      REAL(dp) :: z(plond,plev)        ! height of env
      REAL(dp) :: s(plond,plev)        ! normalized dry static energy of env
      REAL(dp) :: zf(plond,plevp)      ! height of interfaces
      REAL(dp) :: u(plond,plev)        ! zonal velocity of env
      REAL(dp) :: v(plond,plev)        ! merid. velocity of env

      integer  :: jb(plond)         ! updraft base level
      integer  :: lel(plond)        ! updraft launch level
      integer  :: jt(plond)         ! updraft plume top
      integer  :: jlcl(plond)       ! updraft lifting cond level
      integer  :: mx(plond)         ! updraft base level (same is jb)
      integer  :: j0(plond)         ! level where updraft begins detraining
      integer  :: jd(plond)         ! level of downdraft
!      integer  :: limcnv            ! convection limiting level

!output

      REAL(dp) :: cmfdqr(plond,plev)   ! rate of production of precip at that layer
      REAL(dp) :: du(plond,plev)       ! detrainement rate of updraft
      REAL(dp) :: ed(plond,plev)       ! entrainment rate of downdraft
      REAL(dp) :: eu(plond,plev)       ! entrainment rate of updraft
      REAL(dp) :: hmn(plond,plev)      ! moist stat energy of env
      REAL(dp) :: hsat(plond,plev)     ! sat moist stat energy of env
      REAL(dp) :: mc(plond,plev)       ! net mass flux
      REAL(dp) :: md(plond,plev)       ! downdraft mass flux
      REAL(dp) :: mu(plond,plev)       ! updraft mass flux
      REAL(dp) :: pflx(plond,plevp)    ! precipitation flux thru layer
      REAL(dp) :: qd(plond,plev)       ! spec humidity of downdraft
      REAL(dp) :: ql(plond,plev)       ! liq water of updraft
      REAL(dp) :: qst(plond,plev)      ! saturation spec humidity of env.
      REAL(dp) :: qu(plond,plev)       ! spec hum of updraft
      REAL(dp) :: sd(plond,plev)       ! normalized dry stat energy of downdraft
      REAL(dp) :: shat(plond,plev)     ! interface values of dry stat energy
      REAL(dp) :: su(plond,plev)       ! normalized dry stat energy of updraft
      REAL(dp) :: ud(plond,plev)       ! downdraft u
      REAL(dp) :: vd(plond,plev)       ! downdraft v

!    these version of the mass fluxes conserve mass (used in tracer transport)

      REAL(dp) :: mu2(plond,plev)      ! updraft mass flux
      REAL(dp) :: eu2(plond,plev)      ! updraft entrainment
      REAL(dp) :: du2(plond,plev)      ! updraft detrainment
      REAL(dp) :: md2(plond,plev)      ! downdraft mass flux
      REAL(dp) :: ed2(plond,plev)      ! downdraft entrainment
      REAL(dp) :: rl                   ! latent heat of vap

      ! mz_ht_20070918+
      REAL(dp) :: lwc(plond,plev)
      REAL(dp) :: iwc(plond,plev)
      REAL(dp) :: rform(plond,plev)
      REAL(dp) :: sform(plond,plev)
      ! mz_ht_20070918-

      integer  :: il2g              !CORE GROUP REMOVE

      REAL(dp) :: rd                   ! gas constant for dry air
      REAL(dp) :: grav                 ! gravity
      REAL(dp) :: cp                   ! heat capacity of dry air

      integer  :: msg               ! missing moisture vals (always 0)

!Local workspace

      REAL(dp) :: gamma(plond,plev)  
      REAL(dp) :: dz(plond,plev)  
      REAL(dp) :: iprm(plond,plev)  
      REAL(dp) :: hu(plond,plev)  
      REAL(dp) :: hd(plond,plev)  
      REAL(dp) :: eps(plond,plev)  
      REAL(dp) :: f(plond,plev)  
      REAL(dp) :: k1(plond,plev)  
      REAL(dp) :: i2(plond,plev)  
      REAL(dp) :: ihat(plond,plev)  
      REAL(dp) :: i3(plond,plev)  
      REAL(dp) :: idag(plond,plev)  
      REAL(dp) :: i4(plond,plev)  
      REAL(dp) :: qsthat(plond,plev)  
      REAL(dp) :: hsthat(plond,plev)  
      REAL(dp) :: gamhat(plond,plev)  
      REAL(dp) :: cu(plond,plev)  
      REAL(dp) :: evp(plond,plev)  
      REAL(dp) :: cmeg(plond,plev)  
      REAL(dp) :: qds(plond,plev) 
      REAL(dp) :: hmin(plond)  
      REAL(dp) :: expdif(plond)  
      REAL(dp) :: expnum(plond)  
      REAL(dp) :: ftemp(plond)  
      REAL(dp) :: eps0(plond)  
      REAL(dp) :: rmue(plond)  
      REAL(dp) :: zuef(plond)  
      REAL(dp) :: zdef(plond)  
      REAL(dp) :: epsm(plond)  
      REAL(dp) :: ratmjb(plond)  
      REAL(dp) :: est(plond)  
      REAL(dp) :: totpcp(plond)  
      REAL(dp) :: totevp(plond)  
      REAL(dp) :: alfa(plond) 
      REAL(dp) :: beta
!      REAL(dp) :: c0
      REAL(dp) :: ql1
      REAL(dp) :: weight
      REAL(dp) :: tu
      REAL(dp) :: estu
      REAL(dp) :: qstu

      REAL(dp) :: small
      REAL(dp) :: mdt  

      integer  :: khighest
      integer  :: klowest  
      integer  :: kount 
      integer  :: i,k

      logical  :: doit(plond)
      logical  :: done(plond)

! -----------------------------------------------------------------------------

      do i = 1,il2g
        ftemp(i) = 0.
        expnum(i) = 0.
        expdif(i) = 0.
      end do

!jr Change from msg+1 to 1 to prevent blowup

      do k = 1,plev
        do i = 1,il2g
          dz(i,k) = zf(i,k) - zf(i,k+1)
        end do
      end do


!initialize many output and work variables to zero

      do k = msg + 1,plev
        do i = 1,il2g
          k1(i,k) = 0.
          i2(i,k) = 0.
          i3(i,k) = 0.
          i4(i,k) = 0.
          mu(i,k) = 0.
          f(i,k) = 0.
          eps(i,k) = 0.
          eu(i,k) = 0.
          du(i,k) = 0.
          ql(i,k) = 0.
          cu(i,k) = 0.
          evp(i,k) = 0.
          cmeg(i,k) = 0.
          qds(i,k) = q(i,k)
          md(i,k) = 0.
          ed(i,k) = 0.
          sd(i,k) = s(i,k)
          qd(i,k) = q(i,k)
          ud(i,k) = u(i,k)
          vd(i,k) = v(i,k)
          mc(i,k) = 0.
          qu(i,k) = q(i,k)
          su(i,k) = s(i,k)
!       est(i)=exp(a-b/t(i,k))
          est(i) = c1*exp((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
! +bee
          if ( p(i,k)-est(i) .gt. 0. ) then
             qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
          else
             qst(i,k) = 1.0_dp
          end if
! -bee
          gamma(i,k) = qst(i,k)*(1._dp + qst(i,k)/eps1)*eps1*rl/        &
                       (rd*t(i,k)**2)*rl/cp
          hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
          hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
          hu(i,k) = hmn(i,k)
          hd(i,k) = hmn(i,k)
          mu2(i,k) = 0.
          eu2(i,k) = 0.
          du2(i,k) = 0.
          md2(i,k) = 0.
          ed2(i,k) = 0.
          pflx(i,k) = 0.
          cmfdqr(i,k) = 0.
        end do
      end do

!jr Set to zero things which make this routine blow up

      do k=1,msg
        do i=1,il2g
          cmfdqr(i,k) = 0.
          mu2(i,k) = 0.
          eu2(i,k) = 0.
          du2(i,k) = 0.
          md2(i,k) = 0.
          ed2(i,k) = 0.
        end do
      end do

!interpolate the layer values of qst, hsat and gamma to
!layer interfaces

      do i = 1,il2g
        hsthat(i,msg+1) = hsat(i,msg+1)
        qsthat(i,msg+1) = qst(i,msg+1)
        gamhat(i,msg+1) = gamma(i,msg+1)
        totpcp(i) = 0.
        totevp(i) = 0.
      end do
      do k = msg + 2,plev
        do i = 1,il2g
          if (abs(qst(i,k-1)-qst(i,k)).gt.1.E-6_dp) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*       &
                          qst(i,k)/ (qst(i,k-1)-qst(i,k))
          else
            qsthat(i,k) = qst(i,k)
          end if
          hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
          if (abs(gamma(i,k-1)-gamma(i,k)).gt.1.E-6_dp) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*              &
                          gamma(i,k-1)*gamma(i,k)/                   &
                          (gamma(i,k-1)-gamma(i,k))
          else
            gamhat(i,k) = gamma(i,k)
          end if
        end do
      end do

!initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed plev)

      do i = 1,il2g
        jt(i) = max(lel(i),limcnv+1)
        jt(i) = min(jt(i),plev)
        jd(i) = plev
        jlcl(i) = lel(i)
        hmin(i) = 1.E6_dp
      end do

!find the level of minimum hsat, where detrainment starts

      do k = msg + 1,plev
        do i = 1,il2g
          if (hsat(i,k).le.hmin(i) .and. k.ge.jt(i).and.k.le.jb(i)) then
            hmin(i) = hsat(i,k)
            j0(i) = k
          end if
        end do
      end do
      do i = 1,il2g
        j0(i) = min(j0(i),jb(i)-2)
        j0(i) = max(j0(i),jt(i)+2)

!Fix from Guang Zhang to address out of bounds array reference

        j0(i) = min(j0(i),plev)
      end do

!Initialize certain arrays inside cloud

      do k = msg + 1,plev
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.le.jb(i)) then
            hu(i,k) = hmn(i,mx(i)) + cp*0.5_dp
            su(i,k) = s(i,mx(i)) + 0.5_dp
          end if
        end do
      end do

!*********************************************************
!compute taylor series for approximate eps(z) below
!*********************************************************

      do k = plev - 1,msg + 1,-1
        do i = 1,il2g
          if (k.lt.jb(i) .and. k.ge.jt(i)) then
            k1(i,k) = k1(i,k+1) + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
          end if
        end do
      end do

!re-initialize hmin array for ensuing calculation.

      do i = 1,il2g
        hmin(i) = 1.E6_dp
      end do
      do k = msg + 1,plev
        do i = 1,il2g
          if (k.ge.j0(i).and.k.le.jb(i) .and. hmn(i,k).le.hmin(i)) then
            hmin(i) = hmn(i,k)
            expdif(i) = hmn(i,mx(i)) - hmin(i)
          end if
        end do
      end do

!*********************************************************
!compute approximate eps(z) using above taylor series
!*********************************************************

      do k = msg + 2,plev
        do i = 1,il2g
          expnum(i) = 0.
          ftemp(i) = 0.
          if (k.lt.jt(i) .or. k.ge.jb(i)) then
            k1(i,k) = 0.
            expnum(i) = 0.
          else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) +            &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
          end if
          if ((expdif(i).gt.100._dp.and.expnum(i).gt.0.) .and.                       &
              k1(i,k).gt.expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 +                     &
                     (2.*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2*                  &
                     ftemp(i)**3 + (-5.*k1(i,k)*i2(i,k)*i3(i,k)+                  &
                     5.*i2(i,k)**3+k1(i,k)**2*i4(i,k))/                           &
                     k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._dp)
            f(i,k) = min(f(i,k),0.0002_dp)
          end if
        end do
      end do
      do i = 1,il2g
        if (j0(i).lt.jb(i)) then
          if (f(i,j0(i)).lt.1.E-6_dp .and. f(i,j0(i)+1).gt.f(i,j0(i)))               &
            j0(i) = j0(i) + 1
        end if
      end do
      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.le.j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
          end if
        end do
      end do
      do i = 1,il2g
        eps0(i) = f(i,j0(i))
        eps(i,jb(i)) = eps0(i)
      end do
! +pjr
!right now I have set this to do it the nflux41way as I want to match it
!but it is probably better to disable it soon
#define PJRWAY
#ifdef PJRWAY
      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.ge.j0(i) .and. k.le.jb(i)) then
               eps(i,k) = f(i,j0(i))
            end if
         end do
      end do
      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.lt.j0(i) .and. k.ge.jt(i)) eps(i,k) = f(i,k)
         end do
      end do
#else
      do k=plev,msg+1,-1
        do i=1,il2g
          if (k.ge.j0(i)) then
            if (k.le.jb(i)) eps(i,k) = f(i,j0(i))
          else
            if (k.ge.jt(i)) eps(i,k) = f(i,k)
          end if
        end do
      end do
#endif

!specify the updraft mass flux mu, entrainment eu, detrainment du
!and moist static energy hu.
!here and below mu, eu,du, md and ed are all normalized by mb

      do i = 1,il2g
        if (eps0(i).gt.0.) then
!         mu(i,jb(i)) = 1.
!         eu(i,jb(i)) = eps0(i)/2.
! *pjr NOTE TO CORE GROUP mu and mu2 (and eu and eu2 should have the same vals now)
          mu2(i,jb(i)) = 1._dp
          eu2(i,jb(i)) = mu2(i,jb(i))/dz(i,jb(i))
          mu(i,jb(i)) = mu2(i,jb(i))
          eu(i,jb(i)) = eu2(i,jb(i))
        end if
      end do
      do k = plev,msg + 1,-1
        do i = 1,il2g
           if (eps0(i).gt.0. .and. (k.ge.jt(i).and.k.lt.jb(i))) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1._dp/eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._dp)/zuef(i)
            mu(i,k) = (1._dp/eps0(i))* (exp(eps(i,k)*zuef(i))-1._dp)/zuef(i)
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
            mu2(i,k) = mu(i,k)
            eu2(i,k) = eu(i,k)
            du2(i,k) = du(i,k)
          end if
        end do
      end do

      khighest = plevp
      klowest = 1
      do i=1,il2g
        khighest = min(khighest,lel(i))
        klowest = max(klowest,jb(i))
      end do
      do k = klowest-1,khighest,-1
!dir$ ivdep
        do i = 1,il2g
          if (k.le.jb(i)-1 .and. k.ge.lel(i) .and. eps0(i).gt.0.) then
            if (mu(i,k).lt.0.01_dp) then
              hu(i,k) = hu(i,jb(i))
              mu(i,k) = 0.
              mu2(i,k) = mu(i,k)
              eu2(i,k) = 0.
              du2(i,k) = mu2(i,k+1)/dz(i,k)
              eu(i,k) = eu2(i,k)
              du(i,k) = du2(i,k)
            else
              hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) +                    &
                        dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)-              &
                        du(i,k)*hsat(i,k))
            end if
          end if
        end do
      end do

!reset cloud top index beginning from two layers above the
!cloud base (i.e. if cloud is only one layer thick, top is not reset

! *pjr there are diffs here, and I hope to god they dont matter
      do i=1,il2g
        doit(i) = .true.
      end do
      do k=klowest-2,khighest-1,-1
        do i=1,il2g
          if (doit(i) .and. k.le.jb(i)-2 .and. k.ge.lel(i)-1) then
            if (hu(i,k  ).le.hsthat(i,k) .and.                                 &
!              hu(i,k+1).gt.hsthat(i,k+1) .and. mu(i,k).ge.0.02) then
! +mgl set this to 0.01 to be consistent with computation of hu in the cloud column done above
               hu(i,k+1).gt.hsthat(i,k+1) .and. mu(i,k).ge.0.01_dp) then
              if (hu(i,k)-hsthat(i,k).lt.-2000._dp) then
                jt(i) = k + 1
                doit(i) = .false.
              else
                jt(i) = k
                if (eps0(i).le.0.) doit(i) = .false.
              end if
            else if (hu(i,k).gt.hu(i,jb(i)) .or. mu(i,k).lt.0.01_dp               &
! +mgl this condition causes a lot of overturning of only the surface layer!
!                 .or. (hu(i,k)-hsthat(i,k).lt.-2000.)) then
                      ) then
              jt(i) = k + 1
              doit(i) = .false.
            end if
          end if
        end do
      end do
      do k = plev,msg + 1,-1
!dir$ ivdep
        do i = 1,il2g
          if (k.ge.lel(i) .and. k.le.jt(i) .and. eps0(i).gt.0.) then
            mu(i,k) = 0.
            eu(i,k) = 0.
            du(i,k) = 0.
            mu2(i,k) = 0.
            eu2(i,k) = 0.
            du2(i,k) = 0.
            hu(i,k) = hu(i,jb(i))
          end if
          if (k.eq.jt(i) .and. eps0(i).gt.0.) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            du2(i,k) = mu2(i,k+1)/dz(i,k)
            eu2(i,k) = 0.
            mu2(i,k) = 0.
            eu(i,k) = 0.
            mu(i,k) = 0.
          end if
        end do
      end do

!specify downdraft properties (no downdrafts if jd.ge.jb).
!scale down downward mass flux profile so that net flux
!(up-down) at cloud base in not negative.

      do i = 1,il2g

!in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1

        alfa(i) = 0.1_dp
        jt(i) = min(jt(i),jb(i)-1)
        jd(i) = max(j0(i),jt(i)+1)
        jd(i) = min(jd(i),jb(i))
        hd(i,jd(i)) = hmn(i,jd(i)-1)
        ud(i,jd(i)) = u(i,jd(i)-1)
        vd(i,jd(i)) = v(i,jd(i)-1)
        if (jd(i).lt.jb(i) .and. eps0(i).gt.0.) then
          epsm(i) = eps0(i)
!         alfa(i)=2.*epsm(i)*( zf(i,jd(i))-zf(i,jb(i)) )/
!    1         (  exp(2.*epsm(i)*( zf(i,jd(i))-
!              zf(i,jb(i)) ))-1.  )
          md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
          md2(i,jd(i)) = md(i,jd(i))
        end if
      end do
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2.*eps0(i))*                        &
                      (exp(2.*epsm(i)*zdef(i))-1._dp)/zdef(i)
            md2(i,k) = md(i,k)
          end if
        end do
      end do
      do k = msg + 1,plev
!dir$ ivdep
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.  &
              jd(i).lt.jb(i)) then
            ratmjb(i) = min(abs(mu2(i,jb(i))/md2(i,jb(i))),1._dp)
            md2(i,k) = md2(i,k)*ratmjb(i)
!           ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1.)
!           md(i,k) = md(i,k)*ratmjb(i)
               md(i,k) = md2(i,k)
          end if
        end do
      end do

      small = 1.e-20_dp
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.plev) .and. eps0(i).gt.0.) then
            ed2(i,k-1) = (md2(i,k-1)-md2(i,k))/dz(i,k-1)
               ed(i,k-1) = ed2(i,k-1)
            mdt = min(md2(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) -                          &
                       dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
          end if
        end do
      end do

!calculate updraft and downdraft properties.

      do k = msg + 2,plev
        do i = 1,il2g
          if ((k.ge.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.   &
             jd(i).lt.jb(i)) then
!        sd(i,k) = shat(i,k)                                          &
!                +              (hd(i,k)-hsthat(i,k))/                &
!                  (cp    *(1.+gamhat(i,k)))
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/         &
                       (rl*(1._dp + gamhat(i,k)))
          end if
        end do
      end do

      do i = 1,il2g
         done(i) = .false.
      end do
      kount = 0
      do k = plev,msg + 2,-1
        do i = 1,il2g
          if ((.not.done(i) .and. k.gt.jt(i) .and. k.lt.jb(i)) .and.  &
               eps0(i).gt.0.) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) +                   &
                      dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)    
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) +                   &
                      dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)-               &
                      du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cp*zf(i,k)
            estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
            qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.-estu)
            if (qu(i,k).ge.qstu) then
              jlcl(i) = k
              kount = kount + 1
              done(i) = .true.
            end if
          end if
        end do
        if (kount.ge.il2g) goto 690
      end do
 690  continue
      do k = msg + 2,plev
        do i = 1,il2g
          if (k.eq.jb(i) .and. eps0(i).gt.0.) then
            qu(i,k) = q(i,mx(i))
            su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
          end if
          if ((k.gt.jt(i).and.k.le.jlcl(i)) .and. eps0(i).gt.0.) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/             &  
                      (cp* (1._dp + gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*                     &
                     (hu(i,k)-hsthat(i,k))/                          &
                     (rl* (1._dp + gamhat(i,k)))
          end if
        end do
      end do

      do k = plev,msg + 2,-1
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0.) then
            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/        &
                      dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/            &
                      (rl/cp)
            if (k.eq.jt(i)) cu(i,k) = 0.
!              cu(i,k) = max(0.,cu(i,k))
!              cu2     = max(0.,
!              cu2     = max(-1.e99,                                      &
!                        +(eu(i,k)*q(i,k) - du(i,k)*qst(i,k))             &
!                        -(mu(i,k)*qu(i,k)-mu(i,k+1)*qu(i,k+1))/dz(i,k)   & 
!                        )
!                 
!              if (abs(cu(i,k)-cu2)/(abs(cu(i,k))+abs(cu2)+1.e-50)        &
!                   .gt.0.0000001) then
!                 write (6,*) ' inconsistent condensation rates ',        & 
!                      i, k, lat,                                         &
!                      cu(i,k), cu2, jt(i), jb(i), jlcl(i), lel(i)        &
!                      ,mu(i,k)
!              endif
               cu(i,k) = max(0._dp,cu(i,k))
          end if
        end do
      end do

      beta = 0._dp
!      c0 = 2.E-3_dp
      do k = plev,msg + 2,-1
        do i = 1,il2g
          cmfdqr(i,k) = 0.
!this modification is for test3 run, modified on 6/20/1995
!       if(t(i,jt(i) ).gt.tfreez)    c0=0.
!       if(t(i,jt(i) ).le.tfreez   )    c0=2.e-3
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0. .and.          &
               mu(i,k).ge.0.0) then
            if (mu(i,k).gt.0.) then
              ql1 = 1./mu(i,k)* (mu(i,k+1)*ql(i,k+1)-                        &
                    dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
              ql(i,k) = ql1/ (1._dp + dz(i,k)*c0)
            else
              ql(i,k) = 0.
            end if
            totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*                & 
                       (beta*ql(i,k) + (1._dp - beta)*ql(i,k+1)))
            cmfdqr(i,k) = c0*mu(i,k)*ql(i,k)

            ! mz_ht_20070918+
            IF (T(i,k) > Tmelt) THEN
              lwc(i,k)   = ql(i,k)
              rform(i,k) = c0 * ql(i,k)
            ELSE
              iwc(i,k)   = ql(i,k)
              sform(i,k) = c0 * ql(i,k)
            ENDIF
            ! mz_ht_20070918-
          end if
        end do
      end do

      do i = 1,il2g
        qd(i,jd(i)) = qds(i,jd(i))
        sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
      end do

      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jd(i).and.k.lt.jb(i) .and. eps0(i).gt.0.) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) +                                      &
                       (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._dp)
            mdt = min(md(i,k+1),-small)
            sd(i,k+1) = ((rl/cp*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) +            &
                          md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      end do
      do i = 1,il2g
! *guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -           & 
        totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) -                     & 
                    md(i,jb(i))*qd(i,jb(i))
      end do
      if (.true.) then
        do i = 1,il2g
          k = jb(i)
          if (eps0(i).gt.0.) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._dp)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      endif

      do i = 1,il2g
        totpcp(i) = max(totpcp(i),0._dp)
        totevp(i) = max(totevp(i),0._dp)
      end do

      weight = 1.0_dp
      do k = msg + 2,plev
        do i = 1,il2g
          if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
            md2(i,k) = md2(i,k)*min(1._dp,weight*totpcp(i)/                       &
                      (totevp(i)+weight*totpcp(i)))
            ed2(i,k) = ed2(i,k)*min(1._dp,weight*totpcp(i)/                       &
                      (totevp(i)+weight*totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._dp,                                        &
                           weight*totpcp(i)/ (totevp(i)+                       &
                           weight*totpcp(i)))
          else
            md2(i,k) = 0.
            ed2(i,k) = 0.
            evp(i,k) = 0.
          end if
            md(i,k) = md2(i,k)
            ed(i,k) = ed2(i,k)
            cmeg(i,k) = cu(i,k) - evp(i,k)
!          cmeg is the cloud water condensed - rain water evaporated
!          cmfdqr  is the cloud water converted to rain - (rain evaporated)
            cmfdqr(i,k) = cmfdqr(i,k)-evp(i,k)
        end do
      end do
      do k = 2,plevp
        do i = 1,il2g
          pflx(i,k) = pflx(i,k-1) + cmfdqr(i,k-1)*dz(i,k-1)
        end do
      end do
      do i = 1,il2g
        if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
          totevp(i) = totevp(i)*min(1._dp,                                   &
                      weight*totpcp(i)/(totevp(i) + weight*totpcp(i)))
        else
          totevp(i) = 0.
        end if
      end do

      do k = msg + 1,plev
        do i = 1,il2g
          mc(i,k) = mu(i,k) + md(i,k)
        end do
      end do

      return
    end subroutine cldprp
 
!===============================================================================

      subroutine closure(plev    ,plond                              ,&  !nlev, kbdim
                         q       ,t       ,p       ,s                ,&
                         tp      ,qu      ,su      ,mc               ,&
                         du      ,mu      ,md      ,qd      ,sd      ,&
                         qhat    ,shat    ,wdp     ,qstp             ,&
                         zf      ,ql      ,dsubcld ,mb      ,cape    ,&
                         tl      ,lcl     ,lel     ,jt      ,mx      ,&
                         il1g    ,il2g    ,rd      ,grav    ,cp      ,&
                         rl      ,msg     ,capelmt )
! ----------------------------------------------------------------------

!This is contributed code not fully standardized by the CCM core group.

!this code is very much rougher than virtually anything else in the CCM
!We expect to release cleaner code in a future release

!the documentation has been enhanced to the degree that we are able

!Original version:  G. Zhang and collaborators
!Standardized:      Core group staff, 1994 and 195
!Reviewed:          P. Rasch, April 1996
! ----------------------------------------------------------------------

!$Id: conv_zhang.F,v 1.2.2.4 1999/04/03 18:54:30 eaton Exp $

!may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.

! ----------------------------Arguments---------------------------------
IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plond                    ! number of longitudes, kbdim
      REAL(dp) :: q(plond,plev)        ! spec humidity
      REAL(dp) :: t(plond,plev)        ! temperature
      REAL(dp) :: p(plond,plev)        ! pressure (mb)
      REAL(dp) :: s(plond,plev)        ! normalized dry static energy 
      REAL(dp) :: tp(plond,plev)       ! parcel temp
      REAL(dp) :: qu(plond,plev)       ! updraft spec. humidity
      REAL(dp) :: su(plond,plev)       ! normalized dry stat energy of updraft
      REAL(dp) :: mc(plond,plev)       ! net convective mass flux 
      REAL(dp) :: du(plond,plev)       ! detrainment from updraft
      REAL(dp) :: mu(plond,plev)       ! mass flux of updraft
      REAL(dp) :: md(plond,plev)       ! mass flux of downdraft
      REAL(dp) :: qd(plond,plev)       ! spec. humidity of downdraft
      REAL(dp) :: sd(plond,plev)       ! dry static energy of downdraft
      REAL(dp) :: qhat(plond,plev)     ! environment spec humidity at interfaces
      REAL(dp) :: shat(plond,plev)     ! env. normalized dry static energy at intrfcs
      REAL(dp) :: wdp(plond,plev)      ! pressure thickness of layers
      REAL(dp) :: qstp(plond,plev)     ! spec humidity of parcel
      REAL(dp) :: zf(plond,plev+1)     ! height of interface levels
      REAL(dp) :: ql(plond,plev)       ! liquid water mixing ratio

      REAL(dp) :: mb(plond)            ! cloud base mass flux
      REAL(dp) :: cape(plond)          ! available pot. energy of column
      REAL(dp) :: tl(plond)
      REAL(dp) :: dsubcld(plond)       ! thickness of subcloud layer

      integer  :: lcl(plond)        ! index of lcl
      integer  :: lel(plond)        ! index of launch leve
      integer  :: jt(plond)         ! top of updraft
      integer  :: mx(plond)         ! base of updraft

! -------------------------Local variables------------------------------

      REAL(dp) :: dtpdt(plond,plev)
      REAL(dp) :: dqsdtp(plond,plev)
      REAL(dp) :: dtmdt(plond,plev)
      REAL(dp) :: dqmdt(plond,plev)
      REAL(dp) :: dboydt(plond,plev)
      REAL(dp) :: thetavp(plond,plev)
      REAL(dp) :: thetavm(plond,plev)

      REAL(dp) :: dtbdt(plond),dqbdt(plond),dtldt(plond)
      REAL(dp) :: beta
      REAL(dp) :: capelmt
      REAL(dp) :: cp
      REAL(dp) :: dadt
      REAL(dp) :: debdt
      REAL(dp) :: dltaa
      REAL(dp) :: eb
      REAL(dp) :: grav

      integer  :: i
      integer  :: il1g
      integer  :: il2g
      integer  :: k
      integer  :: msg

      REAL(dp) :: rd
      REAL(dp) :: rl
      REAL(dp) :: tau
      save tau

      INTRINSIC :: TINY
!tau=4800. were used in canadian climate center. however, when it
!is used here in echam3 t42, convection is too weak, thus 
!adjusted to 2400. i.e the e-folding time is 1 hour now.

      data tau/7200._dp/
! ----------------------------------------------------------------------
!change of subcloud layer properties due to convection is
!related to cumulus updrafts and downdrafts.
!mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
!to define betau, betad and f(z).
!note that this implies all time derivatives are in effect
!time derivatives per unit cloud-base mass flux, i.e. they
!have units of 1/mb instead of 1/sec.

      do i = il1g,il2g
         mb(i) = 0.
         eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
         dtbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                   &
                     (shat(i,mx(i))-su(i,mx(i)))+                    &
                      md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
         dqbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                   &
                     (qhat(i,mx(i))-qu(i,mx(i)))+                    &
                     md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
         debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
         dtldt(i) = -2840.* (3.5/t(i,mx(i))*dtbdt(i)-debdt/eb)/      & 
                     (3.5*log(t(i,mx(i)))-log(eb)-4.805_dp)**2
      end do

!  dtmdt and dqmdt are cumulus heating and drying.

      do k = msg + 1,plev
         do i = il1g,il2g
            dtmdt(i,k) = 0.
            dqmdt(i,k) = 0.
         end do
      end do

      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.eq.jt(i)) then
               dtmdt(i,k) = (1./wdp(i,k))*                              &
                             (mu(i,k+1)* (su(i,k+1)-shat(i,k+1)-        &
                             rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-    &
                             shat(i,k+1)))
               dqmdt(i,k) = (1./wdp(i,k))*(mu(i,k+1)* (qu(i,k+1)-       & 
                             qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*          &
                             (qd(i,k+1)-qhat(i,k+1)))
            end if
         end do
      end do

      beta = 0.
      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.gt.jt(i) .and. k.lt.mx(i)) then
               dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+               &
                             mc(i,k+1)* (s(i,k)-shat(i,k+1)))/          &
                             wdp(i,k) - rl/cp*du(i,k)*                  &
                             (beta*ql(i,k)+ (1._dp-beta)*ql(i,k+1))
!         dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))                        &
!                     +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/wdp(i,k)         &
!                     +du(i,k)*(qs(i,k)-q(i,k))                         &
!                     +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

               dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+         &
                             cp/rl* (su(i,k+1)-s(i,k)))-                &
                             mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*         &
                             (su(i,k)-s(i,k)))+md(i,k+1)*               &
                             (qd(i,k+1)-qhat(i,k+1)+cp/rl*              & 
                             (sd(i,k+1)-s(i,k)))-md(i,k)*               &
                             (qd(i,k)-qhat(i,k)+cp/rl*                  &
                             (sd(i,k)-s(i,k))))/wdp(i,k) +              &
                             du(i,k)* (beta*ql(i,k)+                    & 
                             (1._dp-beta)*ql(i,k+1))
            end if
         end do
      end do

      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.ge.lel(i) .and. k.le.lcl(i)) then
               thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*        & 
                               (1._dp + 1.608*qstp(i,k)-q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*         &
                               (1._dp + 0.608*q(i,k))
               dqsdtp(i,k) = qstp(i,k)* (1._dp + qstp(i,k)/eps1)*eps1*rl/    &
                               (rd*tp(i,k)**2)

!dtpdt is the parcel temperature change due to change of
!subcloud layer properties during convection.

               dtpdt(i,k) = tp(i,k)/ (1._dp +                           &
                             rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))*   &
                              (dtbdt(i)/t(i,mx(i))+                     &
                             rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/         &
                             tl(i)**2*dtldt(i)))

!dboydt is the integrand of cape change.

               dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1./                   &
                              (1._dp + 1.608*qstp(i,k)-q(i,mx(i)))*     &
                              (1.608 * dqsdtp(i,k) * dtpdt(i,k) -       &
                              dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608/    &
                              (1._dp + 0.608*q(i,k))*dqmdt(i,k)))*grav* &
                              thetavp(i,k)/thetavm(i,k)
            end if
         end do
      end do

      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.gt.lcl(i) .and. k.lt.mx(i)) then
               thetavp(i,k) = tp(i,k)* (1000./p(i,k))** (rd/cp)*        &
                              (1._dp + 0.608*q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000./p(i,k))** (rd/cp)*         &
                              (1._dp + 0.608*q(i,k))

!dboydt is the integrand of cape change.

               dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+                      & 
                              0.608/ (1._dp + 0.608*q(i,mx(i)))*dqbdt(i)-    &
                              dtmdt(i,k)/t(i,k)-                        &
                              0.608/ (1._dp + 0.608*q(i,k))*dqmdt(i,k))*     &
                              grav*thetavp(i,k)/thetavm(i,k) 
            end if
         end do
      end do


!buoyant energy change is set to 2/3*excess cape per 3 hours

      do i = il1g,il2g
         dadt = 0.
         do k = lel(i),mx(i) - 1
            dadt = dadt + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         end do

         dltaa = -1.* (cape(i)-capelmt)
!!         if (dadt.ne.0.) mb(i) = max(dltaa/tau/dadt,0._dp)
         if (ABS(dadt).gt.TINY(0.0_dp)) mb(i) = max(dltaa/tau/dadt,0._dp)
      end do

      return
    end subroutine closure
 
!===========================================================================

    subroutine q1q2_pjr                                            & 
                     (plev    ,plond                              ,&  !nlev, kbdim
                      dqdt    ,dsdt                               ,&
                      qu      ,su      ,du                        ,&
                      qhat    ,shat    ,wdp     ,mu      ,md      ,&
                      sd      ,qd      ,ql      ,dsubcld          ,&
                      jt      ,mx      ,il1g    ,il2g             ,&
                      cp      ,rl      ,msg                       ,&
                      dl      ,evp     ,cu )

      

IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plond                    ! number of longitudes, kbdim
!rewritten by phil rasch dec 19 1995

      REAL(dp) :: dqdt(plond,plev),           &
           dsdt(plond,plev)

      REAL(dp) :: cp
      REAL(dp) :: fact
      integer  :: il1g
      integer  :: il2g
      integer  :: i
      integer  :: k
      integer  :: msg
      REAL(dp) :: &
           qu(plond,plev),             &
           su(plond,plev),             & 
           du(plond,plev),             &
           qhat(plond,plev),           & 
           shat(plond,plev),           &
           wdp(plond,plev),            &
           mu(plond,plev),             &
           md(plond,plev),             &
           sd(plond,plev),             &
           qd(plond,plev),             &
           ql(plond,plev)
      REAL(dp) ::                      &
           dl(plond,plev),             &
           evp(plond,plev),            &
           cu(plond,plev)
      integer  :: kbm
      integer  :: ktm
      REAL(dp) :: emc
      REAL(dp) :: dsubcld(plond)


      integer  ::  jt(plond), mx(plond)

!work fields:

      REAL(dp) :: rl
! ------------------------------------------------------------------
      do k = msg + 1,plev
         do i = il1g,il2g
            dsdt(i,k) = 0.
            dqdt(i,k) = 0.
            dl(i,k) = 0.
         end do
      end do

!find the highest level top and bottom levels of convection

      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do

      fact = 0.
!     fact = 1.

      do k = ktm,plev-1
         do i = il1g,il2g
!           fact = 1.
!           if (q(i,k).gt.0.8*qs(i,k) .and. k.lt.plev-3) fact = 0.

            emc = +fact*du(i,k)*ql(i,k+1) &! evaporating cloud detraining to env
                  -cu(i,k)                &! condensation in updraft
                  +evp(i,k)                ! evaporating rain in downdraft
!           emc = 0

            dsdt(i,k) = -rl/cp*emc                              &
                       + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1))   &
                          -mu(i,k)*   (su(i,k)-shat(i,k))       &
                          +md(i,k+1)* (sd(i,k+1)-shat(i,k+1))   &
                          -md(i,k)*   (sd(i,k)-shat(i,k))       &
                         )/wdp(i,k)
          

            dqdt(i,k) = emc                                     &
                        +(+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1))   &
                          -mu(i,k)*   (qu(i,k)-qhat(i,k))       &
                          +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1))   &
                          -md(i,k)*   (qd(i,k)-qhat(i,k))       &
                          )/wdp(i,k)
                       

            dl(i,k) = (1._dp - fact)*du(i,k)*ql(i,k+1)

         end do
      end do


      do k = kbm,plev             
         do i = il1g,il2g
            if (k.eq.mx(i)) then
               dsdt(i,k) = (1./dsubcld(i))*                &
                    (-mu(i,k)* (su(i,k)-shat(i,k))         &
                     -md(i,k)* (sd(i,k)-shat(i,k))         &
                    )
               dqdt(i,k) = (1./dsubcld(i))*                &
                    (-mu(i,k)*(qu(i,k)-qhat(i,k))          &
                     -md(i,k)*(qd(i,k)-qhat(i,k))          &
                    )
            else if (k.gt.mx(i)) then
               dsdt(i,k) = dsdt(i,k-1)
               dqdt(i,k) = dqdt(i,k-1)
            end if
         end do
      end do

      return
    end subroutine q1q2_pjr
!------------------------------ subroutines for conv_ccm end -------------------- 
!=================================================================================

    subroutine arconvtran(plev  ,plon    ,plond      ,&  !nlev, kproma, kbdim
                          wdp, mu, md, eu,            &
                          mug, mdg, dug, eug,                    &
                          edg, dpg, dsubcld, jtg, jbg,           &
                          ideep, lengath )


!This is a setup routine for convective transport using archived mass
!fluxes from the Zhang scheme.  The setup involves:
!1. Gather mass flux arrays.
!2. Calc the mass fluxes that are determined by mass balance.
!3. Determine top and bottom of convection.

  
IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plon,               &    ! "longitude index", kproma
                 plond                    ! number of longitudes, kbdim
!Input arguments:

      REAL(dp) ::                &
       wdp(plond,plev)    &  ! delta pressure between interfaces (Pa)
     , mu(plond,plev)     &  ! mass flux up (kg/m2/s)
     , md(plond,plev)     &  ! mass flux down (kg/m2/s)
     , eu(plond,plev)        ! mass entraining from updraft (1/s)

!Output arguments:

      REAL(dp) ::                &
       mug(plond,plev)    &  ! gathered mu
     , mdg(plond,plev)    &  ! gathered md
     , dug(plond,plev)    &  ! mass detraining from updraft (gathered)
     , eug(plond,plev)    &  ! gathered eu
     , edg(plond,plev)    &  ! mass entraining from downdraft (gathered)
     , dpg(plond,plev)    &  ! gathered .01*wdp
     , dsubcld(plond)        ! delta pressure from cloud base to sfc (gathered)
      integer  ::         &
       jtg(plond)         &
     , jbg(plond)         &
     , ideep(plond)       &
     , lengath

!Local variables:

      integer  ::       &
       i, k             &
     , lindex(plond)    &
     , lenpos

      REAL(dp) ::               &
       lsum(plond)        &
     , rdpg(plond,plev)      ! gathered 1./(.01*dp)

      REAL(dp) :: gravit
      parameter( gravit = 9.80616_dp )
! ----------------------------------------------------------------------

!    Gathered array contains all columns with a updraft.
      do i = 1, plon
         lsum(i) = 0.
      end do
      do k = 1, plev
         do i = 1, plon
            lsum(i) = lsum(i) + mu(i,k)
         end do
      end do
      call whenfgt( plon, lsum, 1, 0.0_dp, ideep, lengath )

      if ( lengath .eq. 0 ) return

!    Gather input mass fluxes.
      do k = 1, plev
         do i = 1, lengath
            dpg(i,k) = .01*wdp(ideep(i),k)         ! convert Pa -> mb
            rdpg(i,k) = 1./dpg(i,k)
            mug(i,k) = mu(ideep(i),k)*gravit*.01  ! convert kg/m2/s -> mb/s
            mdg(i,k) = md(ideep(i),k)*gravit*.01
            eug(i,k) = eu(ideep(i),k)
         end do
      end do

!    Calc du and ed.
      do k = 1, plev-1
         do i = 1, lengath
            dug(i,k) = eug(i,k) - (mug(i,k)-mug(i,k+1))*rdpg(i,k)
            edg(i,k) = (mdg(i,k)-mdg(i,k+1))*rdpg(i,k)
         end do
      end do
      do i = 1, lengath
         dug(i,plev) = eug(i,plev) - mug(i,plev)*rdpg(i,plev)
         edg(i,plev) = 0.0
      end do
      do k = 1, plev
         do i = 1, lengath
            if ( dug(i,k) .lt. 1.e-7*eug(i,k) ) dug(i,k) = 0.0
         end do
      end do

!    Find top and bottom layers with updrafts.
      do i = 1, lengath
         jtg(i) = plev
         jbg(i) = 1
      end do
      do k = 2, plev
         call whenfgt( lengath, mug(1,k), 1, 0.0_dp, lindex, lenpos )
         do i = 1, lenpos
            jtg(lindex(i)) = min( k-1, jtg(lindex(i)) )
            jbg(lindex(i)) = max( k, jbg(lindex(i)) )
         end do
      end do

!    Calc delta p between srfc and cloud base.
      do i = 1, lengath
         dsubcld(i) = dpg(i,plev)
      end do
      do k = plev-1, 2, -1
         do i = 1, lengath
            if ( jbg(i) .le. k ) then
               dsubcld(i) = dsubcld(i) + dpg(i,k)
            end if
         end do
      end do

      return
    end subroutine arconvtran

!===========================================================================================

!==========================================================================

      subroutine cmfmca(plev    ,plevp   ,plon    ,plond   ,pcnst   , &  !nlev, nlevp1, kproma, kbdim, ntrac
                        tdt     ,pmid    ,pdel                      , &
                        rpdel   ,zm      ,tpert   ,qpert   ,phis    , &
                        pblht   ,t       ,q1      ,cmfdt   ,cmfdq   , &
                        cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   , &
                        qc      ,cnt     ,cnb     ,                   & 
! +scyc
                        icwmr   ,                                     &
! -scyc
! +match
                        q, ncnst, hketa, hkbeta   ,                   &
                        rflx2, sflx2, lwc, iwc,    rform, sform   )
! -match
! ----------------------------------------------------------------------

! Moist convective mass flux procedure:
! If stratification is unstable to nonentraining parcel ascent,
! complete an adjustment making successive use of a simple cloud model
! consisting of three layers (sometimes referred to as a triplet)
! 
! Code generalized to allow specification of parcel ("updraft")
! properties, as well as convective transport of an arbitrary
! number of passive constituents (see q array).  The code
! is written so the water vapor field is passed independently
! in the calling list from the block of other transported
! constituents, even though as currently designed, it is the
! first component in the constituents field. 

! ---------------------------Code History-------------------------------

! Original version:  J. J. Hack, March 22, 1990
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, G. Taylor, August 1992

!#######################################################################
!#                                                                     #
!# Debugging blocks are marked this way for easy identification        #
!#                                                                     #
!#######################################################################

! ----------------------------------------------------------------------

! $Id: conv_hack.F,v 1.1.1.1.2.1 1998/11/06 16:16:43 eaton Exp $

    IMPLICIT NONE
    
    INTRINSIC :: TINY
      

      REAL(dp) :: ssfac               ! supersaturation bound (detrained air)
      parameter (ssfac = 1.001_dp)

! -----------------------------Arguments--------------------------------

     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plevp,              &    ! number of levels + 1
                 plon,               &    ! "longitude index", kproma
                 plond,              &    ! number of longitudes, kbdim
                 pcnst                    ! number of tracers
! Input arguments

! +match
      INTEGER  :: ncnst            ! number of passive tracers
! -match

      REAL(dp) :: tdt                 ! 2 delta-t (seconds)
      REAL(dp) :: pmid(plond,plev)    ! pressure
      REAL(dp) :: pdel(plond,plev)    ! delta-p
      REAL(dp) :: rpdel(plond,plev)   ! 1./pdel
      REAL(dp) :: zm(plond,plev)      ! height abv sfc at midpoints
      REAL(dp) :: tpert(plond)        ! PBL perturbation theta
! +match
      REAL(dp) :: qpert(plond)        ! PBL perturbation specific humidity 
! -match
      REAL(dp) :: phis(plond)         ! surface geopotential
      REAL(dp) :: pblht(plond)        ! PBL height (provided by PBL routine)

! Input/output arguments

      REAL(dp) :: t(plond,plev)       ! temperature (t bar)
! +match
      REAL(dp) :: q1(plond,plev)      ! specific humidity (sh bar)
      REAL(dp) :: q(plond,plev,pcnst) ! passive tracers
! -match

! Output arguments

      REAL(dp) :: cmfdt(plond,plev)   ! dt/dt due to moist convection
      REAL(dp) :: cmfdq(plond,plev)   ! dq/dt due to moist convection
      REAL(dp) :: cmfmc(plond,plev )  ! moist convection cloud mass flux
      REAL(dp) :: cmfdqr(plond,plev)  ! dq/dt due to convective rainout 
      REAL(dp) :: cmfsl(plond,plev )  ! convective lw static energy flux
      REAL(dp) :: cmflq(plond,plev )  ! convective total water flux
      REAL(dp) :: precc(plond)        ! convective precipitation rate
      REAL(dp) :: qc(plond,plev)      ! dq/dt due to rainout terms
      REAL(dp) :: cnt(plond)          ! top level of convective activity   
      REAL(dp) :: cnb(plond)          ! bottom level of convective activity
! +scyc
      REAL(dp) :: icwmr(plond,plev)
! -scyc
! +match
      REAL(dp) ::                     &
                   hketa(plond,plev), &  ! save eta
                   hkbeta(plond,plev)    ! save beta 
! -match

      REAL(dp) :: rflx2(plond, plevp) ! rainflux per level
      REAL(dp) :: sflx2(plond, plevp) ! rainflux per level
      REAL(dp) :: lwc(plond, plev)    ! cloud water per level
      REAL(dp) :: iwc(plond, plev)    ! cloud ice per level
      REAL(dp) :: rform(plond, plev)  ! rain production per level (mixing ratio)
      REAL(dp) :: sform(plond, plev)  ! rain production per level (mixing ratio)

! --------------------------Local workspace-----------------------------

      REAL(dp) :: gam(plond,plev)     ! 1/cp (d(qsat)/dT)
      REAL(dp) :: sb(plond,plev)      ! dry static energy (s bar)
      REAL(dp) :: hb(plond,plev)      ! moist static energy (h bar)
      REAL(dp) :: shbs(plond,plev)    ! sat. specific humidity (sh bar star)
      REAL(dp) :: hbs(plond,plev)     ! sat. moist static energy (h bar star)
      REAL(dp) :: shbh(plond,plevp)   ! specific humidity on interfaces
      REAL(dp) :: sbh(plond,plevp)    ! s bar on interfaces
      REAL(dp) :: hbh(plond,plevp)    ! h bar on interfaces
      REAL(dp) :: cmrh(plond,plevp)   ! interface constituent mixing ratio 
      REAL(dp) :: prec(plond)         ! instantaneous total precipitation
      REAL(dp) :: dzcld(plond)        ! depth of convective layer (m)
      REAL(dp) :: beta(plond)         ! overshoot parameter (fraction)
      REAL(dp) :: betamx(plond)       ! local maximum on overshoot
      REAL(dp) :: eta(plond)          ! convective mass flux (kg/m^2 s)
      REAL(dp) :: etagdt(plond)       ! eta*grav*dt
      REAL(dp) :: cldwtr(plond)       ! cloud water (mass)
      REAL(dp) :: rnwtr(plond)        ! rain water  (mass)
      REAL(dp) :: sc  (plond)         ! dry static energy   ("in-cloud")
      REAL(dp) :: shc (plond)         ! specific humidity   ("in-cloud")
      REAL(dp) :: hc  (plond)         ! moist static energy ("in-cloud")
      REAL(dp) :: cmrc(plond)         ! constituent mix rat ("in-cloud")
      REAL(dp) :: dq1(plond)          ! shb  convective change (lower lvl)
      REAL(dp) :: dq2(plond)          ! shb  convective change (mid level)
      REAL(dp) :: dq3(plond)          ! shb  convective change (upper lvl)
      REAL(dp) :: ds1(plond)          ! sb   convective change (lower lvl)
      REAL(dp) :: ds2(plond)          ! sb   convective change (mid level)
      REAL(dp) :: ds3(plond)          ! sb   convective change (upper lvl)
      REAL(dp) :: dcmr1(plond)        ! q convective change (lower lvl)
      REAL(dp) :: dcmr2(plond)        ! q convective change (mid level)
      REAL(dp) :: dcmr3(plond)        ! q convective change (upper lvl)
      REAL(dp) :: estemp(plond,plev)  ! saturation vapor pressure (scratch)
      REAL(dp) :: vtemp1(2*plond)     ! intermediate scratch vector
      REAL(dp) :: vtemp2(2*plond)     ! intermediate scratch vector
      REAL(dp) :: vtemp3(2*plond)     ! intermediate scratch vector
      REAL(dp) :: vtemp4(2*plond)     ! intermediate scratch vector
      integer  :: indx1(plond)        ! longitude indices for condition true
      logical  :: etagt0              ! true if eta > 0.0
      REAL(dp) :: cats                ! modified characteristic adj. time
      REAL(dp) :: rtdt                ! 1./tdt
      REAL(dp) :: qprime              ! modified specific humidity pert.
      REAL(dp) :: tprime              ! modified thermal perturbation
      REAL(dp) :: pblhgt              ! bounded pbl height (max[pblh,1m])
      REAL(dp) :: fac1                ! intermediate scratch variable
      REAL(dp) :: shprme              ! intermediate specific humidity pert.
      REAL(dp) :: qsattp              ! sat mix rat for thermally pert PBL parcels 
      REAL(dp) :: dz                  ! local layer depth
      REAL(dp) :: temp1               ! intermediate scratch variable
      REAL(dp) :: b1                  ! bouyancy measure in detrainment lvl
      REAL(dp) :: b2                  ! bouyancy measure in condensation lvl
      REAL(dp) :: temp2               ! intermediate scratch variable
      REAL(dp) :: temp3               ! intermediate scratch variable
      REAL(dp) :: g                   ! bounded vertical gradient of hb
      REAL(dp) :: tmass               ! total mass available for convective exch
      REAL(dp) :: denom               ! intermediate scratch variable
      REAL(dp) :: qtest1              ! used in negative q test (middle lvl) 
      REAL(dp) :: qtest2              ! used in negative q test (lower lvl) 
      REAL(dp) :: fslkp               ! flux lw static energy (bot interface)
      REAL(dp) :: fslkm               ! flux lw static energy (top interface)
      REAL(dp) :: fqlkp               ! flux total water (bottom interface)
      REAL(dp) :: fqlkm               ! flux total water (top interface)
      REAL(dp) :: botflx              ! bottom constituent mixing ratio flux
      REAL(dp) :: topflx              ! top constituent mixing ratio flux
      REAL(dp) :: efac1               ! ratio q to convectively induced chg (btm lvl)
      REAL(dp) :: efac2               ! ratio q to convectively induced chg (mid lvl)
      REAL(dp) :: efac3               ! ratio q to convectively induced chg (top lvl)
      REAL(dp) :: tb(plond,plev)      ! working storage for temp (t bar)
      REAL(dp) :: shb(plond,plev)     ! working storage for spec hum (sh bar)
      REAL(dp) :: adjfac              ! adjustment factor (relaxation related)
#if ( defined DIAGNS )

!  Following 7 REAL(dp) :: variables are used in diagnostics calculations   

      REAL(dp) :: rh                  ! relative humidity 
      REAL(dp) :: es                  ! sat vapor pressure 
      REAL(dp) :: hsum1               ! moist static energy integral
      REAL(dp) :: qsum1               ! total water integral 
      REAL(dp) :: hsum2               ! final moist static energy integral
      REAL(dp) :: qsum2               ! final total water integral
      REAL(dp) :: fac                 ! intermediate scratch variable
#endif
      integer  :: i,k              ! longitude, level indices
      integer  :: ii               ! index on "gathered" vectors
      integer  :: len1             ! vector length of "gathered" vectors
      integer  :: m                ! constituent index
      integer  :: ktp              ! tmp indx used to track top of convective layer
#if ( defined DIAGNS )
      integer  :: n                ! vertical index     (diagnostics)
      integer  :: kp               ! vertical index     (diagnostics)
      integer  :: kpp              ! index offset, kp+1 (diagnostics)
      integer  :: kpm1             ! index offset, kp-1 (diagnostics)
#endif

! --------------------------Statement functions-------------------------

!      REAL(dp) :: qhalf
!      qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),                        &
!                                 (shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))


! ----------------------------------------------------------------------

! Ensure that characteristic adjustment time scale (cmftau) assumed
! in estimate of eta isn't smaller than model time scale (tdt)
! The time over which the convection is assumed to act (the adjustment
! time scale) can be applied with each application of the three-level
! cloud model, or applied to the column tendencies after a "hard"
! adjustment (i.e., on a 2-delta t time scale) is evaluated
      rform(:,:) = 0._dp
      sform(:,:) = 0._dp
      lwc(:,:)   = 0._dp
      iwc(:,:)   = 0._dp

      tb(:,:) = 0._dp
      if (rlxclm) then
        cats   = tdt               ! relaxation applied to column
        adjfac = tdt/(max(tdt,cmftau))
      else
        cats   = max(tdt,cmftau) ! relaxation applied to triplet
        adjfac = 1.0_dp
      endif
      rtdt = 1.0/tdt

! Move temperature and moisture into working storage

      do k=limcnv,plev
        do i=1,plon
          tb (i,k) = t(i,k)
! +match
          shb(i,k) = q1(i,k)
! -match
        end do
      end do
! +scyc
      do k=1,plev
         do i=1,plon
            icwmr(i,k) = 0.
! +match
            hketa(i,k) = 0.
            hkbeta(i,k) = 0.
! -match
         end do
      end do
! -scyc

! Compute sb,hb,shbs,hbs

      call aqsatd(tb      ,pmid    ,estemp ,shbs    ,gam     ,&
                  plond   ,plon    ,plev   ,limcnv  ,plev    )

      do k=limcnv,plev
        do i=1,plon
          sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
          hb (i,k) = sb(i,k) + hlat*shb(i,k)
          hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
        end do
      end do

! Compute sbh, shbh

      do k=limcnv+1,plev
        do i=1,plon
          sbh (i,k) = 0.5*(sb(i,k-1) + sb(i,k))
          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
        end do
      end do

! Specify properties at top of model (not used, but filling anyway)

      do i=1,plon
        sbh (i,limcnv) = sb(i,limcnv)
        shbh(i,limcnv) = shb(i,limcnv)
        hbh (i,limcnv) = hb(i,limcnv)
      end do

! Zero vertically independent control, tendency & diagnostic arrays

      do i=1,plon
        prec(i)  = 0.0
        dzcld(i) = 0.0
        cnb(i)   = 0.0
        cnt(i)   = REAL(plev+1,dp)
      end do
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    output initial thermodynamic profile if debug diagnostics        #
!#                                                                     #
      if (lat.eq.jloc .and. nstep.ge.nsloc) then
        i = iloc
!#                                                                     #
!#       approximate vertical integral of moist static energy          #
!#       and total preciptable water                                   #
!#                                                                     #
        hsum1 = 0.0
        qsum1 = 0.0
        do k=limcnv,plev
          hsum1 = hsum1 + pdel(i,k)*rgrav*hb(i,k)
          qsum1 = qsum1 + pdel(i,k)*rgrav*shb(i,k)
        end do
!#                                                                     #
        write (6,8010)
        fac = grav*864.
        do k=limcnv,plev
          rh = shb(i,k)/shbs(i,k)
          write(6,8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k),&
              cmfsl(i,k), cmflq(i,k)
          write(6,8040) tb(i,k),shb(i,k),rh,sb(i,k),hb(i,k),hbs(i,k), &
              tdt*cmfdt(i,k),tdt*cmfdq(i,k),tdt*cmfdqr(i,k)
        end do
        write(6, 8000) prec(i)
      end if
#endif
!#                                                                     #
!#                                                                     #
!#######################################################################

! Begin moist convective mass flux adjustment procedure.
! Formalism ensures that negative cloud liquid water can never occur

      do 70 k=plev-1,limcnv+1,-1
        do 10 i=1,plon
          etagdt(i) = 0.0
          eta   (i) = 0.0
          beta  (i) = 0.0
          ds1   (i) = 0.0
          ds2   (i) = 0.0
          ds3   (i) = 0.0
          dq1   (i) = 0.0
          dq2   (i) = 0.0
          dq3   (i) = 0.0

! Specification of "cloud base" conditions

          qprime    = 0.0
          tprime    = 0.0

! Assign tprime within the PBL to be proportional to the quantity
! tpert (which will be bounded by tpmax), passed to this routine by 
! the PBL routine.  Don't allow perturbation to produce a dry 
! adiabatically unstable parcel.  Assign qprime within the PBL to be 
! an appropriately modified value of the quantity qpert (which will be 
! bounded by shpmax) passed to this routine by the PBL routine.  The 
! quantity qprime should be less than the local saturation value 
! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
! linearly reduced toward zero as the PBL top is approached.

          pblhgt = max(pblht(i),1.0_dp)
          if( (zm(i,k+1) .le. pblhgt) .and.                     &
           ABS(dzcld(i)).lt.TINY(0.0_dp) ) then
!!                                 dzcld(i).eq.0.0 ) then
            fac1   = max(0.0_dp,1.0_dp-zm(i,k+1)/pblhgt)
            tprime = min(tpert(i),tpmax)*fac1
            qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
! +match
            shprme = min(min(qpert(i),shpmax)*fac1,             &
! -match
                          max(qsattp-shb(i,k+1),0.0_dp))
            qprime = max(qprime,shprme)
          else
            tprime = 0.0
            qprime = 0.0
          end if

! Specify "updraft" (in-cloud) thermodynamic properties

          sc (i)    = sb (i,k+1) + cp*tprime
          shc(i)    = shb(i,k+1) + qprime
          hc (i)    = sc (i    ) + hlat*shc(i)
          vtemp4(i) = hc(i) - hbs(i,k)
          dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
          if (vtemp4(i).gt.0.0) then
            dzcld(i) = dzcld(i) + dz
          else
            dzcld(i) = 0.0
          end if
10      enddo
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    output thermodynamic perturbation information                    #
!#                                                                     #
        if (lat.eq.jloc .and. nstep.ge.nsloc) then
          write (6,8090) k+1,sc(iloc),shc(iloc),hc(iloc)
        end if
!#                                                                     #
!#######################################################################
#endif

! Check on moist convective instability
! Build index vector of points where instability exists

        call whenfgt(plon,vtemp4,1,0.0_dp,indx1,len1)
        if (len1.le.0) go to 70

! Current level just below top level => no overshoot

        if (k.le.limcnv+1) then
          do ii=1,len1
            i = indx1(ii)
            temp1     = vtemp4(i)/(1.0_dp + gam(i,k))
            cldwtr(i) = max(0.0_dp,(sb(i,k) - sc(i) + temp1))
            beta(i)   = 0.0_dp
            vtemp3(i) = (1.0_dp + gam(i,k))*(sc(i) - sbh(i,k))
          end do
        else

! First guess at overshoot parameter using crude buoyancy closure
! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
! If pre-existing supersaturation in detrainment layer, beta=0
! cldwtr is temporarily equal to hlat*l (l=> liquid water)

!DIR$ IVDEP
          do ii=1,len1
            i = indx1(ii)
            temp1     = vtemp4(i)/(1.0_dp + gam(i,k))
            cldwtr(i) = max(0.0_dp,(sb(i,k)-sc(i)+temp1))
            betamx(i) = 1.0_dp - c0*max(0.0_dp,(dzcld(i)-dzmin))
            b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
            b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
            beta(i)   = max(betamn,min(betamx(i), 1.0_dp + b1/b2))
            if (hbs(i,k-1).le.hb(i,k-1)) beta(i) = 0.0_dp

! Bound maximum beta to ensure physically realistic solutions

! First check constrains beta so that eta remains positive
! (assuming that eta is already positive for beta equal zero)

            vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+     &
                       (1.0_dp + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
            vtemp2(i) = (1.0_dp + gam(i,k))*(sc(i) - sbh(i,k))
            vtemp3(i) = vtemp2(i)
            if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0._dp) then
              betamx(i) = 0.99*(vtemp1(i)/vtemp2(i))
              beta(i)   = max(0.0_dp,min(betamx(i),beta(i)))
            end if
          end do

! Second check involves supersaturation of "detrainment layer"
! small amount of supersaturation acceptable (by ssfac factor)

!DIR$ IVDEP
          do ii=1,len1
            i = indx1(ii)
            if (hb(i,k-1).lt.hbs(i,k-1)) then
              vtemp1(i) = vtemp1(i)*rpdel(i,k)
              temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -        &
                         hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
              temp3 = vtemp3(i)*rpdel(i,k)
              vtemp2(i) = (tdt/cats)*(hc(i) - hbs(i,k))*temp2/           &
                         (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
              if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0.) then
                betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
                beta(i)   = max(0.0_dp,min(betamx(i),beta(i)))
              end if
            else 
              beta(i) = 0.0
            end if
          end do

! Third check to avoid introducing 2 delta x thermodynamic
! noise in the vertical ... constrain adjusted h (or theta e)
! so that the adjustment doesn't contribute to "kinks" in h

!DIR$ IVDEP
          do ii=1,len1
            i = indx1(ii)
            g = min(0.0_dp,hb(i,k) - hb(i,k-1))
            temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/tdt)/                &
                       (hc(i) - hbs(i,k))
            vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*          &
                           rpdel(i,k)
            vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) +                     &
                           (hc(i) - hbh(i,k) - cldwtr(i))*               &
                           (rpdel(i,k) + rpdel(i,k+1))
            if ((beta(i)*vtemp2(i) - vtemp1(i)).gt.0._dp) then
!!              if (vtemp2(i).ne.0.0) then
              if (ABS(vtemp2(i)).gt.TINY(0.0_dp)) then
                betamx(i) = vtemp1(i)/vtemp2(i)
              else
                betamx(i) = 0.0_dp
              end if
              beta(i) = max(0.0_dp,min(betamx(i),beta(i)))
            end if
          end do
        end if

! Calculate mass flux required for stabilization.

! Ensure that the convective mass flux, eta, is positive by
! setting negative values of eta to zero..
! Ensure that estimated mass flux cannot move more than the
! minimum of total mass contained in either layer k or layer k+1.
! Also test for other pathological cases that result in non-
! physical states and adjust eta accordingly.

!DIR$ IVDEP
        do ii=1,len1
          i = indx1(ii)
          beta(i) = max(0.0_dp,beta(i))
          temp1 = hc(i) - hbs(i,k)
          temp2 = ((1.0_dp + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) -  &
                     beta(i)*vtemp3(i))*rpdel(i,k) -                    &
                    (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
          eta(i) = temp1/(temp2*grav*cats)
          tmass = min(pdel(i,k),pdel(i,k+1))*rgrav
          if (eta(i).gt.tmass*rtdt .or. eta(i).le.0.0_dp) eta(i) = 0.0_dp

! Check on negative q in top layer (bound beta)

!!          if (shc(i)-shbh(i,k).lt.0.0 .and. beta(i)*eta(i).ne.0.0) then
          if (shc(i)-shbh(i,k).lt.0.0 .and. &
           ABS(beta(i)*eta(i)).gt.TINY(0.0_dp)) then
            denom = eta(i)*grav*tdt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
            beta(i) = max(0.0_dp,min(-0.999*shb(i,k-1)/denom, beta(i)))
          end if

! Check on negative q in middle layer (zero eta)

          qtest1 = shb(i,k) + eta(i)*grav*tdt*((shc(i) - shbh(i,k+1)) -  &
                   (1.0_dp - beta(i))*cldwtr(i)*rhlat -                  &
                   beta(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
          if (qtest1.le.0.0) eta(i) = 0.0

! Check on negative q in lower layer (bound eta)

          fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
          qtest2 = shb(i,k+1) - eta(i)*grav*tdt*fac1
          if (qtest2 .lt. 0.0) then
            eta(i) = 0.99*shb(i,k+1)/(grav*tdt*fac1)
          end if
          etagdt(i) = eta(i)*grav*tdt
! +match
          hketa(i,k)  = eta(i) * adjfac
          hkbeta(i,k) = beta(i)
! -match
        end do

#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
        if (lat.eq.jloc .and. nstep.ge.nsloc) then
          write(6,8080) beta(iloc), eta(iloc)
        end if
!#                                                                     #
!#######################################################################
#endif

! Calculate cloud water, rain water, and thermodynamic changes

!DIR$ IVDEP
        do 30 ii=1,len1
          i = indx1(ii)
! +scyc
          icwmr(i,k) = cldwtr(i)*rhlat
! -scyc
          cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav
          rnwtr(i) = (1.0_dp - beta(i))*cldwtr(i)
          ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
          dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
          ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +                      &
                   hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*               &
                   (sc(i) - sbh(i,k)))*rpdel(i,k)
          dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) -                    &
                   grav*rnwtr(i) - beta(i)*etagdt(i)*                     &
                   (shc(i) - shbh(i,k)))*rpdel(i,k)
          ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) -                &
                   hlat*grav*cldwtr(i))*rpdel(i,k-1)
          dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)

! Isolate convective fluxes for later diagnostics

          fslkp = eta(i)*(sc(i) - sbh(i,k+1))
          fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) -                    &
                           hlat*cldwtr(i)*rtdt)
          fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
          fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))

! Update thermodynamic profile (update sb, hb, & hbs later)

          tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
          tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
          tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
          shb(i,k+1) = shb(i,k+1) + dq1(i)
          shb(i,k  ) = shb(i,k  ) + dq2(i)
          shb(i,k-1) = shb(i,k-1) + dq3(i)

! ** Update diagnostic information for final budget **
! Tracking precipitation, temperature & specific humidity tendencies,
! rainout term, convective mass flux, convective liquid
! water static energy flux, and convective total water flux
! The variable afac makes the necessary adjustment to the
! diagnostic fluxes to account for adjustment time scale based on
! how relaxation time scale is to be applied (column vs. triplet)

          prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac

          ! mz_ht_20051129+        
          ! conversion to mm / timestep
          IF (T(i,k) > Tmelt) THEN
            lwc(i,k)   = cldwtr(i)
            rform(i,k) = (rnwtr(i)/rhoh2o)*adjfac
            rflx2(i,k) = prec(i) !+ (rnwtr(i)/rhoh2o)*adjfac*1.e3_dp
          ELSE
            iwc(i,k)   = cldwtr(i)
            sform(i,k) = (rnwtr(i)/rhoh2o)*adjfac
            sflx2(i,k) = prec(i) !+ (rnwtr(i)/rhoh2o)*adjfac*1.e3_dp
          ENDIF
          ! mz_ht_20051129-


! The following variables have units of "units"/second

          cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rcp*rtdt*adjfac
          cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rcp*rtdt*adjfac
          cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rcp*rtdt*adjfac
          cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
          cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
          cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac
          qc    (i,k  ) = (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
          cmfdqr(i,k  ) = cmfdqr(i,k  ) + qc(i,k)
          cmfmc (i,k+1) = cmfmc (i,k+1) + eta(i)*adjfac
          cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)*adjfac

! The following variables have units of w/m**2

          cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
          cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
          cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
          cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
30      enddo

! Next, convectively modify passive constituents
! For now, when applying relaxation time scale to thermal fields after 
! entire column has undergone convective overturning, constituents will 
! be mixed using a "relaxed" value of the mass flux determined above
! Although this will be inconsistant with the treatment of the thermal
! fields, it's computationally much cheaper, no more-or-less justifiable,
! and consistent with how the history tape mass fluxes would be used in
! an off-line mode (i.e., using an off-line transport model)

! +match
!        do 50 m=2,pcnst         ! note: indexing assumes water is first field
        do 50 m=1,ncnst
! -match
!DIR$ IVDEP
          do 40 ii=1,len1
            i = indx1(ii)

! If any of the reported values of the constituent is negative in
! the three adjacent levels, nothing will be done to the profile

            if ((q(i,k+1,m).lt.0.0) .or. (q(i,k,m).lt.0.0) .or.           &
                 (q(i,k-1,m).lt.0.0)) go to 40

! Specify constituent interface values (linear interpolation)

            cmrh(i,k  ) = 0.5*(q(i,k-1,m) + q(i,k  ,m))
            cmrh(i,k+1) = 0.5*(q(i,k  ,m) + q(i,k+1,m))

! Specify perturbation properties of constituents in PBL

! +match
! comment out "if" test since if qpert is zero for constituents then only
! the assignment that occurs in the "else" is necessary.
!            pblhgt = max(pblht(i),1.0)
!            if( (zm(i,k+1) .le. pblhgt)                    &
!                                           .and.           & 
!                                     dzcld(i).eq.0.0 ) then
!              fac1 = max(0.0,1.0-zm(i,k+1)/pblhgt)
!              cmrc(i) = q(i,k+1,m) + qpert(i,m)*fac1
!            else
              cmrc(i) = q(i,k+1,m)
!            end if
! -match

! Determine fluxes, flux divergence => changes due to convection
! Logic must be included to avoid producing negative values. A bit
! messy since there are no a priori assumptions about profiles.
! Tendency is modified (reduced) when pending disaster detected.

            botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
            topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
            dcmr1(i) = -botflx*rpdel(i,k+1)
            efac1    = 1.0_dp
            efac2    = 1.0_dp
            efac3    = 1.0_dp

            if (q(i,k+1,m)+dcmr1(i) .lt. 0.0) then
              efac1 = max(small,abs(q(i,k+1,m)/dcmr1(i)) - eps)
            end if

!!            if (efac1.eq.small .or. efac1.gt.1.0_dp) efac1 = 0.0
            if (ABS(efac1-small).lt.tiny(0.0_dp) .or. efac1.gt.1.0_dp) efac1 = 0.0
            dcmr1(i) = -efac1*botflx*rpdel(i,k+1)
            dcmr2(i) = (efac1*botflx - topflx)*rpdel(i,k)
  
            if (q(i,k,m)+dcmr2(i) .lt. 0.0) then
              efac2 = max(small,abs(q(i,k  ,m)/dcmr2(i)) - eps)
            end if

!!            if (efac2.eq.small .or. efac2.gt.1.0_dp) efac2 = 0.0
            if (ABS(efac2-small).lt.tiny(0.0_dp) .or. efac2.gt.1.0_dp) efac2 = 0.0
            dcmr2(i) = (efac1*botflx - efac2*topflx)*rpdel(i,k)
            dcmr3(i) = efac2*topflx*rpdel(i,k-1)

            if (q(i,k-1,m)+dcmr3(i) .lt. 0.0) then
              efac3 = max(small,abs(q(i,k-1,m)/dcmr3(i)) - eps)
            end if

!!            if (efac3.eq.small .or. efac3.gt.1.0_dp) efac3 = 0.0
            if (ABS(efac3-small).lt.tiny(0.0_dp) .or. efac3.gt.1.0_dp) efac3 = 0.0
            efac3    = min(efac2,efac3)
            dcmr2(i) = (efac1*botflx - efac3*topflx)*rpdel(i,k)
            dcmr3(i) = efac3*topflx*rpdel(i,k-1)

            q(i,k+1,m) = q(i,k+1,m) + dcmr1(i)
            q(i,k  ,m) = q(i,k  ,m) + dcmr2(i)
            q(i,k-1,m) = q(i,k-1,m) + dcmr3(i)
40        enddo
50      enddo                ! end of m=2,pcnst loop

! Constituent modifications complete

        if (k.eq.limcnv+1) go to 60

! Complete update of thermodynamic structure at integer levels
! gather/scatter points that need new values of shbs and gamma

        do ii=1,len1
          i = indx1(ii)
          vtemp1(ii     ) = tb(i,k)
          vtemp1(ii+len1) = tb(i,k-1)
          vtemp2(ii     ) = pmid(i,k)
          vtemp2(ii+len1) = pmid(i,k-1)
        end do
        call vqsatd (vtemp1  ,vtemp2  ,estemp  ,vtemp3  , vtemp4  ,          &
                     2*len1   )    ! using estemp as extra long vector
!DIR$ IVDEP
        do ii=1,len1
          i = indx1(ii)
          shbs(i,k  ) = vtemp3(ii     )
          shbs(i,k-1) = vtemp3(ii+len1)
          gam(i,k  ) = vtemp4(ii     )
          gam(i,k-1) = vtemp4(ii+len1)
          sb (i,k  ) = sb(i,k  ) + ds2(i)
          sb (i,k-1) = sb(i,k-1) + ds3(i)
          hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
          hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
          hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
          hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
        end do

! Update thermodynamic information at half (i.e., interface) levels

!DIR$ IVDEP
        do ii=1,len1
          i = indx1(ii)
          sbh (i,k) = 0.5*(sb(i,k) + sb(i,k-1))
          shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
          hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
          sbh (i,k-1) = 0.5*(sb(i,k-1) + sb(i,k-2))
          shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),                       &
                              shbs(i,k-2),shbs(i,k-1))
          hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
        end do

#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    this update necessary, only if debugging diagnostics requested   #
!#                                                                     #
        if (lat.eq.jloc .and. nstep.ge.nsloc) then
          do i=1,plon
            call vqsatd(tb(i,k+1),pmid(i,k+1),es,shbs(i,k+1),gam(i,k+1),  &
                                                                     1)
            sb (i,k+1) = sb(i,k+1) + ds1(i)
            hb (i,k+1) = sb(i,k+1) + hlat*shb(i,k+1)
            hbs(i,k+1) = sb(i,k+1) + hlat*shbs(i,k+1)
            kpp = k + 2
            if(k+1.eq.plev) kpp = k + 1
            do kp=k+1,kpp
              kpm1 = kp-1
              sbh(i,kp)  = 0.5*(sb(i,kpm1) + sb(i,kp))
              shbh(i,kp) = qhalf(shb(i,kpm1),shb(i,kp),shbs(i,kpm1),   &
                                 shbs(i,kp))
              hbh(i,kp)  = sbh(i,kp) + hlat*shbh(i,kp)
            end do
          end do
!#                                                                     #
!#          diagnostic output                                          #
!#                                                                     #
          i = iloc
          write(6, 8060) k
          fac = grav*864.
          do n=limcnv,plev
            rh  = shb(i,n)/shbs(i,n)
            write(6,8020)shbh(i,n),sbh(i,n),hbh(i,n),fac*cmfmc(i,n),  &
                   cmfsl(i,n), cmflq(i,n)
! -------------write(6, 8050)
! -------------write(6, 8030) fac*cmfmc(i,n),cmfsl(i,n), cmflq(i,n)
            write(6, 8040) tb(i,n),shb(i,n),rh,sb(i,n),hb(i,n),       &
                hbs(i,n), tdt*cmfdt(i,n),tdt*cmfdq(i,n),tdt*cmfdqr(i,n)
          end do
          write(6, 8000) prec(i)
        end if
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif

! Ensure that dzcld is reset if convective mass flux zero
! specify the current vertical extent of the convective activity
! top of convective layer determined by size of overshoot param.

   60   do i=1,plon
          etagt0 = eta(i).gt.0.0
          if (.not.etagt0) dzcld(i) = 0.0
          if (etagt0 .and. beta(i).gt.betamn) then
            ktp = k-1
          else
            ktp = k
          end if
          if (etagt0) then
            cnt(i) = min(cnt(i),REAL(ktp,dp))
            cnb(i) = max(cnb(i),REAL(k,dp))
          end if
        end do
70    enddo                 ! end of k loop

! ** apply final thermodynamic tendencies **

      do k=limcnv,plev
        do i=1,plon
          t (i,k) = t (i,k) + cmfdt(i,k)*tdt
!++match
          q1(i,k) = q1(i,k) + cmfdq(i,k)*tdt
! -match
        end do
      end do

! Kludge to prevent cnb-cnt from being zero (in the event
! someone decides that they want to divide by this quantity)

      do i=1,plon
!!        if (cnb(i).ne.0.0 .and. cnb(i).eq.cnt(i)) then
        if (ABS(cnb(i)).gt.tiny(0.0_dp) .and. &
         ABS(cnb(i)-cnt(i)).lt.tiny(0.0_dp)) then
          cnt(i) = cnt(i) - 1.0_dp
        end if
      end do

      do i=1,plon
        precc(i) = prec(i)*rtdt
      end do

#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    we're done ... show final result if debug diagnostics requested  #
!#                                                                     #
      if (lat.eq.jloc .and. nstep.ge.nsloc) then
        i=iloc
        fac = grav*864.
        write(6, 8010)
        do k=limcnv,plev
          rh = shb(i,k)/shbs(i,k)
          write(6, 8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k),   &
              cmfsl(i,k), cmflq(i,k)
          write(6, 8040) tb(i,k),shb(i,k),rh   ,sb(i,k),hb(i,k),       &
              hbs(i,k), tdt*cmfdt(i,k),tdt*cmfdq(i,k),tdt*cmfdqr(i,k)
        end do
        write(6, 8000) prec(i)
!#                                                                     #
!#       approximate vertical integral of moist static energy and      #
!#       total preciptable water after adjustment and output changes   #
!#                                                                     #
        hsum2 = 0.0
        qsum2 = 0.0
        do k=limcnv,plev
          hsum2 = hsum2 + pdel(i,k)*rgrav*hb(i,k)
          qsum2 = qsum2 + pdel(i,k)*rgrav*shb(i,k)
        end do
!#                                                                     #
        write (6,8070) hsum1, hsum2, abs(hsum2-hsum1)/hsum2,          &
             qsum1, qsum2, abs(qsum2-qsum1)/qsum2
      end if
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif
      return                 ! we're all done ... return to calling procedure
#if ( defined DIAGNS )

! Formats

 8000 format(///,10x,'PREC = ',3pf12.6,/)
 8010 format('1**        TB      SHB      RH       SB',                  &
           '       HB      HBS      CAH      CAM       PRECC ',          &
           '     ETA      FSL       FLQ     **', /)
 8020 format(' ----- ',     9x,3p,f7.3,2x,2p,     9x,-3p,f7.3,2x,        &
           f7.3, 37x, 0p,2x,f8.2,  0p,2x,f8.2,2x,f8.2, ' ----- ')
 8030 format(' ----- ',  0p,82x,f8.2,  0p,2x,f8.2,2x,f8.2,               &
             ' ----- ')                      
 8040 format(' - - - ',f7.3,2x,3p,f7.3,2x,2p,f7.3,2x,-3p,f7.3,2x,        &
           f7.3, 2x,f8.3,2x,0p,f7.3,3p,2x,f7.3,2x,f7.3,30x,              &
             ' - - - ')
 8050 format(' ----- ',110x,' ----- ')
 8060 format('1 K =>',  i4,/,                                            &
             '           TB      SHB      RH       SB',                  & 
             '       HB      HBS      CAH      CAM       PREC ',         &
             '     ETA      FSL       FLQ', /)
 8070 format(' VERTICALLY INTEGRATED MOIST STATIC ENERGY BEFORE, AFTER', &
             ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/,    &
             ' VERTICALLY INTEGRATED MOISTURE            BEFORE, AFTER', &
             ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/)
 8080       format(' BETA, ETA => ', 1p,2e12.3)
 8090 format (' k+1, sc, shc, hc => ', 1x, i2, 1p, 3e12.4)
#endif

   end subroutine cmfmca
 
!=========================================================================
!!$ 
!!$      subroutine hkcmfadj(plev    ,plevp   ,plon    ,plond   ,pcnst   , &  !nlev, nlevp1, kproma, kbdim, ntrac
!!$                          tdt, rpdel, eta, beta, q )
!!$
!!$! Compute the convective mass flux adjustment to all tracers using the
!!$! convective mass fluxes and overshoot parameters for the Hack scheme.
!!$! This code was extracted from cmfmca.
!!$
!!$
!!$
!!$! -----------------------------Arguments--------------------------------
!!$IMPLICIT NONE
!!$! Input arguments
!!$  integer, INTENT(IN)  ::         &
!!$                 plev,               &    ! number of levels
!!$                 plevp,              &    ! number of levels + 1
!!$                 plon,               &    ! "longitude index", kproma
!!$                 plond,              &    ! number of longitudes, kbdim
!!$                 pcnst                    ! number of tracers
!!$      real(dp) ::           &
!!$       tdt                  &   ! 2 delta-t (seconds)
!!$     , rpdel(plond,plev)    &   ! 1./pdel
!!$     , eta(plond,plev)      &   ! convective mass flux (kg/m^2 s)
!!$     , beta(plond,plev)         ! overshoot parameter (fraction)
!!$
!!$! Input/output arguments
!!$
!!$      real(dp) ::           &
!!$       q(plond,plev,pcnst)      ! passive tracers
!!$
!!$! --------------------------Local workspace-----------------------------
!!$
!!$      integer  ::         &
!!$       indx1(plond)        ! longitude indices for condition true
!!$      real(dp) ::         &
!!$        adjfac            & ! adjustment factor (relaxation related)
!!$      , etagdt(plond)     & ! eta*grav*dt
!!$      , cmrh(plond,plevp) & ! interface constituent mixing ratio 
!!$      , cmrc(plond)       & ! constituent mix rat ("in-cloud")
!!$      , dcmr1(plond)      & ! q convective change (lower lvl)
!!$      , dcmr2(plond)      & ! q convective change (mid level)
!!$      , dcmr3(plond)      & ! q convective change (upper lvl)
!!$      , botflx            & ! bottom constituent mixing ratio flux
!!$      , topflx            & ! top constituent mixing ratio flux
!!$      , efac1             & ! ratio q to convectively induced chg (btm lvl)
!!$      , efac2             & ! ratio q to convectively induced chg (mid lvl)
!!$      , efac3               ! ratio q to convectively induced chg (top lvl)
!!$
!!$      integer  ::        &
!!$          i,k,           &         ! longitude, level indices
!!$          ii,            &         ! index on "gathered" vectors
!!$          len1,          &         ! vector length of "gathered" vectors
!!$          m                        ! constituent index
!!$! ----------------------------------------------------------------------
!!$
!!$! Ensure that characteristic adjustment time scale (cmftau) assumed
!!$! in estimate of eta isn't smaller than model time scale (tdt)
!!$! The time over which the convection is assumed to act (the adjustment
!!$! time scale) can be applied with each application of the three-level
!!$! cloud model, or applied to the column tendencies after a "hard"
!!$! adjustment (i.e., on a 2-delta t time scale) is evaluated
!!$
!!$      if (rlxclm) then
!!$         adjfac = tdt/(max(tdt,cmftau))
!!$      else
!!$         adjfac = 1.0
!!$      endif
!!$
!!$! Begin moist convective mass flux adjustment procedure.
!!$! Formalism ensures that negative cloud liquid water can never occur
!!$
!!$      do 70 k=plev-1,limcnv+1,-1
!!$        len1 = 0
!!$        do 10 i=1,plon
!!$           if ( eta(i,k) .ne. 0.0 ) then
!!$              etagdt(i) = eta(i,k)*grav*tdt
!!$              len1 = len1 + 1
!!$              indx1(len1) = i
!!$           else
!!$              etagdt(i) = 0.0
!!$           end if
!!$   10   enddo
!!$        if (len1.le.0) go to 70
!!$
!!$! Next, convectively modify passive constituents
!!$! For now, when applying relaxation time scale to thermal fields after 
!!$! entire column has undergone convective overturning, constituents will 
!!$! be mixed using a "relaxed" value of the mass flux determined above
!!$! Although this will be inconsistant with the treatment of the thermal
!!$! fields, it's computationally much cheaper, no more-or-less justifiable,
!!$! and consistent with how the history tape mass fluxes would be used in
!!$! an off-line mode (i.e., using an off-line transport model)
!!$
!!$
!!$! +rvk: this routine used significant portion of the total model run time.
!!$!       When more tracer are treated than the average number of elements in
!!$!       the gathered array (ave(len1) ~ .5 * #levs in troposphere ~ 10-15 )
!!$!       then other order of loops is more efficient
!!$!       for pcnst > 8 vectorization over pcnst is used by CPP otherwise it is
!!$!       checked within with "if" during run
!!$#if ( PCNST > 8 )
!!$        do 41 ii=1,len1
!!$!DIR$ IVDEP
!!$           do 51 m=1,pcnst
!!$#else
!!$        if (len1 .ge. pcnst) then
!!$
!!$        do 50 m=1,pcnst
!!$!DIR$ IVDEP
!!$           do 40 ii=1,len1
!!$
!!$            i = indx1(ii)
!!$
!!$! If any of the reported values of the constituent is negative in
!!$! the three adjacent levels, nothing will be done to the profile
!!$
!!$            if ((q(i,k+1,m).lt.0.0) .or. (q(i,k,m).lt.0.0) .or.         &
!!$                (q(i,k-1,m).lt.0.0)) go to 40
!!$
!!$! Specify constituent interface values (linear interpolation)
!!$
!!$            cmrh(i,k  ) = 0.5*(q(i,k-1,m) + q(i,k  ,m))
!!$            cmrh(i,k+1) = 0.5*(q(i,k  ,m) + q(i,k+1,m))
!!$
!!$            cmrc(i) = q(i,k+1,m)
!!$
!!$! Determine fluxes, flux divergence => changes due to convection
!!$! Logi! must be included to avoid producing negative values. A bit
!!$! messy since there are no a priori assumptions about profiles.
!!$! Tendency is modified (reduced) when pending disaster detected.
!!$
!!$            botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
!!$            topflx   = beta(i,k)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
!!$            dcmr1(i) = -botflx*rpdel(i,k+1)
!!$            efac1    = 1.0
!!$            efac2    = 1.0
!!$            efac3    = 1.0
!!$
!!$            if (q(i,k+1,m)+dcmr1(i) .lt. 0.0) then
!!$              efac1 = max(small,abs(q(i,k+1,m)/dcmr1(i)) - eps)
!!$            end if
!!$
!!$            if (efac1.eq.small .or. efac1.gt.1.0) efac1 = 0.0
!!$            dcmr1(i) = -efac1*botflx*rpdel(i,k+1)
!!$            dcmr2(i) = (efac1*botflx - topflx)*rpdel(i,k)
!!$!  
!!$            if (q(i,k,m)+dcmr2(i) .lt. 0.0) then
!!$              efac2 = max(small,abs(q(i,k  ,m)/dcmr2(i)) - eps)
!!$            end if
!!$
!!$            if (efac2.eq.small .or. efac2.gt.1.0) efac2 = 0.0
!!$            dcmr2(i) = (efac1*botflx - efac2*topflx)*rpdel(i,k)
!!$            dcmr3(i) = efac2*topflx*rpdel(i,k-1)
!!$
!!$            if (q(i,k-1,m)+dcmr3(i) .lt. 0.0) then
!!$              efac3 = max(small,abs(q(i,k-1,m)/dcmr3(i)) - eps)
!!$            end if
!!$
!!$            if (efac3.eq.small .or. efac3.gt.1.0) efac3 = 0.0
!!$            efac3    = min(efac2,efac3)
!!$            dcmr2(i) = (efac1*botflx - efac3*topflx)*rpdel(i,k)
!!$            dcmr3(i) = efac3*topflx*rpdel(i,k-1)
!!$
!!$            q(i,k+1,m) = q(i,k+1,m) + dcmr1(i)
!!$            q(i,k  ,m) = q(i,k  ,m) + dcmr2(i)
!!$            q(i,k-1,m) = q(i,k-1,m) + dcmr3(i)
!!$
!!$40        enddo
!!$50      enddo
!!$
!!$      else
!!$
!!$        do 41 ii=1,len1
!!$!DIR$ IVDEP
!!$           do 51 m=1,pcnst
!!$#endif
!!$! rvk: this comes in both versions:
!!$
!!$            i = indx1(ii)
!!$
!!$! If any of the reported values of the constituent is negative in
!!$! the three adjacent levels, nothing will be done to the profile
!!$
!!$            if ((q(i,k+1,m).lt.0.0) .or. (q(i,k,m).lt.0.0) .or.     &
!!$                (q(i,k-1,m).lt.0.0)) go to 51
!!$
!!$! Specify constituent interface values (linear interpolation)
!!$
!!$            cmrh(i,k  ) = 0.5*(q(i,k-1,m) + q(i,k  ,m))
!!$            cmrh(i,k+1) = 0.5*(q(i,k  ,m) + q(i,k+1,m))
!!$
!!$            cmrc(i) = q(i,k+1,m)
!!$
!!$! Determine fluxes, flux divergence => changes due to convection
!!$! Logic must be included to avoid producing negative values. A bit
!!$! messy since there are no a priori assumptions about profiles.
!!$! Tendency is modified (reduced) when pending disaster detected.
!!$
!!$            botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
!!$            topflx   = beta(i,k)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
!!$            dcmr1(i) = -botflx*rpdel(i,k+1)
!!$            efac1    = 1.0
!!$            efac2    = 1.0
!!$            efac3    = 1.0
!!$
!!$            if (q(i,k+1,m)+dcmr1(i) .lt. 0.0) then
!!$              efac1 = max(small,abs(q(i,k+1,m)/dcmr1(i)) - eps)
!!$            end if
!!$
!!$            if (efac1.eq.small .or. efac1.gt.1.0) efac1 = 0.0
!!$            dcmr1(i) = -efac1*botflx*rpdel(i,k+1)
!!$            dcmr2(i) = (efac1*botflx - topflx)*rpdel(i,k)
!!$!  
!!$            if (q(i,k,m)+dcmr2(i) .lt. 0.0) then
!!$              efac2 = max(small,abs(q(i,k  ,m)/dcmr2(i)) - eps)
!!$            end if
!!$
!!$            if (efac2.eq.small .or. efac2.gt.1.0) efac2 = 0.0
!!$            dcmr2(i) = (efac1*botflx - efac2*topflx)*rpdel(i,k)
!!$            dcmr3(i) = efac2*topflx*rpdel(i,k-1)
!!$
!!$            if (q(i,k-1,m)+dcmr3(i) .lt. 0.0) then
!!$              efac3 = max(small,abs(q(i,k-1,m)/dcmr3(i)) - eps)
!!$            end if
!!$
!!$            if (efac3.eq.small .or. efac3.gt.1.0) efac3 = 0.0
!!$            efac3    = min(efac2,efac3)
!!$            dcmr2(i) = (efac1*botflx - efac3*topflx)*rpdel(i,k)
!!$            dcmr3(i) = efac3*topflx*rpdel(i,k-1)
!!$
!!$            q(i,k+1,m) = q(i,k+1,m) + dcmr1(i)
!!$            q(i,k  ,m) = q(i,k  ,m) + dcmr2(i)
!!$            q(i,k-1,m) = q(i,k-1,m) + dcmr3(i)
!!$51        enddo
!!$41      enddo
!!$
!!$! +rvk  this ends the interactive "if (len1 .ge. pcnst) then"
!!$#if ( PCNST <= 8 )
!!$        endif
!!$#endif
!!$! -rvk
!!$
!!$
!!$! Constituent modifications complete
!!$
!!$70    enddo
!!$      return
!!$    end subroutine hkcmfadj
!===========================================================================================

    subroutine conv_ccm_pjr(plev    ,plevp   ,plon    ,plond            ,  &  !nlev, nlevp1, kproma, kbdim
                            t       ,qh      ,pcpc    ,jctop   ,jcbot   ,  &
                            pblh    ,zm      ,geos    ,zi      ,qtg     ,  &
                            ttg     ,pap     ,paph    ,dpp     ,ts      ,  &
                            delt    ,mcon    ,cme     ,lat     ,           &
                            tpert   ,dlf     ,pflx    ,zdu     ,           &
! +bee
                            cmfdqr  ,mup     ,mdn     ,wpert3d ,capblfc ,  &
                            pblt    ,ireset  ,                             &
                            mu2     ,md2     ,du2     ,eu2     ,ed2     ,  &
                            wdp     ,dsubcld ,jtg     ,maxg    ,ideep   ,  &
                            lengath ,ql      ,cape,    lwc,     iwc,      &
                            rform   ,sform )
! -bee

! ##### MAIN DRIVER FOR ZHANG-MCFARLANE CONVECTION SCHEME #####

!-----------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version 
! of the CCM where it will include a more straightforward formulation 
! and will make use of the standard CCM nomenclature
!-----------------------------------------------------------------------

! same as conv.up except saturation vapor pressure is calculated
! in a different way.

! jul 17/92 - guang jun zhang, m.lazare. calls new buoyan, q1q2
!             and moment (several new work fields added for later).

! nov 21/91 - m.lazare. like previous conv except calls new
!                       clpdprp.
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             previous version conv.
! performs deep convective adjustment based on mass-flux closure
! algorithm.

! ************************ index of variables **********************
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.

!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  wg * betad    downward mass flux at cloud base.
!  wg * betau    upward   mass flux at cloud base.
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpg    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      local sigma half-level thickness (i.e. dshj).
!  ic  * delt     length of model time-step in seconds.
!  wg * wdp      layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * plev     number of model levels.
!  ic  * ilg      lon+2 = size of grid slice.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlat     gaussian latitude index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jtg      top  level index of deep cumulus convection.
!  ic  * kount    current model timestep number.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  wg * long     gathered values of lon.
!  ic  * lev      plev+1.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified 
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  c  * kups     number of points undergoing deep convection.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  i/o * pcp      row of precipitable water in metres.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  i  * pressg   row of surface pressure in pa.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qdb      row of qd at cloud base.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qug      grid slice of mixing ratio in updraft.
!  ic  * rd      dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cpg).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sdb      row of sd at cloud base.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  i  * shbj     grid slice of local bottom interface sigma values.
!  i  * shj      grid slice of local half-level sigma values.
!  i  * shtj     row of local top interfaces of first level.
!  wg * sug      grid slice of dry static energy in updraft.
!  wg * sumde    row of vertically-integrated moist static energy 
!                change.
!  wg * sumdq    row of vertically-integrated scaled mixing ratio 
!                change.
!  wg * sumdt    row of vertically-integrated dry static energy change.
!  wg * sumq     row of vertically-integrated mixing ratio change.
!  i/o * t        grid slice of temperature at mid-layer.
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  o  * jcbot    row of base of cloud indices passed out.
!  w  * tf       grid slice of temperature at interface.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.

!-----------------------------------------------------------------------
     
    IMPLICIT NONE
     integer, INTENT(IN)  ::         &
                 plev,               &    ! number of levels
                 plevp,              &    ! number of levels + 1
                 plon,               &    ! "longitude index", kproma
                 plond                    ! number of longitudes, kbdim

! multi-level i/o fields:

! input/output arguments:

      REAL(dp) :: t(plond,plev)

! jr Added pcnst index to moisture array for consistency

      REAL(dp) :: qh(plond,plev,1) 
      REAL(dp) :: u(plond,plev) 
      REAL(dp) :: v(plond,plev) 
!      REAL(dp) :: utg(plond,plev) 
!      REAL(dp) :: vtg(plond,plev) 
      REAL(dp) :: qtg(plond,plev) 
      REAL(dp) :: ttg(plond,plev)

!       input arguments
      REAL(dp), INTENT(IN) :: pap(plond,plev) 
      REAL(dp), INTENT(IN) :: paph(plond,plev+1) 
      REAL(dp), INTENT(IN) :: dpp(plond,plev) 
      REAL(dp), INTENT(IN) :: zm(plond,plev) 
      REAL(dp), INTENT(IN) :: geos(plond) 
      REAL(dp), INTENT(IN) :: zi(plond,plev+1)
      REAL(dp), INTENT(IN) :: pblh(plond) 
      REAL(dp)             :: zs(plond) 
      REAL(dp), INTENT(IN) :: tpert(plond) 

! output arguments

      REAL(dp) :: pcpck(plond,plev)
      REAL(dp) :: mcon(plond,plev) 
      REAL(dp) :: dlg(plond,plev)     ! gathrd version of the detraining cld h2o tend
      REAL(dp) :: dlf(plond,plev)     ! scattrd version of the detraining cld h2o tend
      REAL(dp) :: pflx(plond,plevp)   ! scattered precip flux at each level
      REAL(dp) :: pflxg(plond,plevp)  ! gather precip flux at each level
      REAL(dp) :: cug(plond,plev)     ! gathered condensation rate 
      REAL(dp) :: evpg(plond,plev)    ! gathered evap rate of rain in downdraft
      REAL(dp) :: mumax(plond) 
      REAL(dp) :: cme(plond,plev)
      REAL(dp) :: zdu(plond,plev)
      REAL(dp) :: mup(plond,plevp)
      REAL(dp) :: mdn(plond,plevp)
      REAL(dp) :: wpert3d(plond,plev)
! +bee
      REAL(dp) :: cmfdqr(plond,plev)
! -bee

! single-level i/o fields:

!       input arguments
      REAL(dp) :: ts(plond) 
      REAL(dp) :: pblt(plond)

!       input/output arguments:

      REAL(dp) :: paprc(plond) 
      REAL(dp) :: paprs(plond) 

!       output arguments:

      REAL(dp) :: jctop(plond) 
      REAL(dp) :: jcbot(plond)
      REAL(dp) :: pcpr(plond) 
      REAL(dp) :: pcps(plond) 
      REAL(dp) :: pcpc(plond)

!-----------------------------------------------------------------------

! general work fields (local variables):

      REAL(dp) :: q(plond,plev) 
      REAL(dp) :: p(plond,plev) 
      REAL(dp) :: z(plond,plev) 
      REAL(dp) :: s(plond,plev) 
      REAL(dp) :: qh0(plond,plev) 
      REAL(dp) :: tp(plond,plev) 
      REAL(dp) :: zf(plond,plev+1) 
      REAL(dp) :: pf(plond,plev+1) 
      REAL(dp) :: qstp(plond,plev) 

      REAL(dp) :: cape(plond) 
      REAL(dp) :: tl(plond) 
      REAL(dp) :: sumq(plond) 
!      REAL(dp) :: sumdt(plond) 
!      REAL(dp) :: sumdq(plond) 
!      REAL(dp) :: sumde(plond)

      INTEGER  :: lcl(plond) 
      INTEGER  :: lfc(plond)
      INTEGER  :: lel(plond) 
      INTEGER  :: maxi(plond) 
      INTEGER  :: ideep(plond) 
      REAL(dp) :: precip

! gathered work fields:

      REAL(dp) :: qg(plond,plev) 
      REAL(dp) :: tg(plond,plev) 
      REAL(dp) :: pg(plond,plev) 
      REAL(dp) :: zg(plond,plev) 
      REAL(dp) :: sg(plond,plev) 
      REAL(dp) :: tpg(plond,plev) 
      REAL(dp) :: zfg(plond,plev+1) 
      REAL(dp) :: qstpg(plond,plev) 
      REAL(dp) :: ug(plond,plev) 
      REAL(dp) :: vg(plond,plev) 
      REAL(dp) :: cmeg(plond,plev)
!!      REAL(dp) :: cu(plond,plev)
!!      REAL(dp) :: evp(plond,plev)

! +bee
      REAL(dp) :: cmfdqrg(plond,plev)
! -bee
      REAL(dp) :: capeg(plond) 
      REAL(dp) :: tlg(plond)

      INTEGER  :: lclg(plond) 
      INTEGER  :: lelg(plond) 
      INTEGER  :: maxg(plond)

! work fields arising from gathered calculations.

      REAL(dp) :: mu(plond,plev) 
      REAL(dp) :: eu(plond,plev) 
      REAL(dp) :: dqdt(plond,plev) 
      REAL(dp) :: dsdt(plond,plev) 
      REAL(dp) :: du(plond,plev) 
      REAL(dp) :: md(plond,plev) 
      REAL(dp) :: ed(plond,plev) 
!!      REAL(dp) :: alpha(plond,plev) 
      REAL(dp) :: sd(plond,plev) 
      REAL(dp) :: qd(plond,plev) 
      REAL(dp) :: mc(plond,plev) 
      REAL(dp) :: qhat(plond,plev) 
      REAL(dp) :: qug(plond,plev) 
      REAL(dp) :: sug(plond,plev) 
      REAL(dp) :: qs(plond,plev) 
      REAL(dp) :: shat(plond,plev) 
      REAL(dp) :: wdp(plond,plev) 
      REAL(dp) :: hmn(plond,plev) 
      REAL(dp) :: hsat(plond,plev) 
      REAL(dp) :: ql(plond,plev) 
! +bee
      REAL(dp) :: qlg(plond,plev) 
! -bee
      REAL(dp) :: dudt(plond,plev) 
      REAL(dp) :: dvdt(plond,plev) 
      REAL(dp) :: ud(plond,plev) 
      REAL(dp) :: vd(plond,plev)

    ! mz_ht_20070918+
      REAL(dp) :: lwc(plond,plev)
      REAL(dp) :: iwc(plond,plev)
      REAL(dp) :: rform(plond,plev)
      REAL(dp) :: sform(plond,plev)
      REAL(dp) :: lwcg(plond,plev)
      REAL(dp) :: iwcg(plond,plev)
      REAL(dp) :: rformg(plond,plev)
      REAL(dp) :: sformg(plond,plev)
      ! mz_ht_20070918-
!      REAL(dp) :: deltat(plond,plev) 
!      REAL(dp) :: deltaq(plond,plev)

      REAL(dp) :: betau(plond) 
      REAL(dp) :: betad(plond) 
!!      REAL(dp) :: qdb(plond) 
!!      REAL(dp) :: sdb(plond) 
      REAL(dp) :: dsubcld(plond) 
      REAL(dp) :: mb(plond) 
      REAL(dp) :: totpcp(plond) 
      REAL(dp) :: totevp(plond)
      REAL(dp) :: mu2(plond,plev) 
      REAL(dp) :: eu2(plond,plev) 
      REAL(dp) :: du2(plond,plev) 
      REAL(dp) :: md2(plond,plev) 
      REAL(dp) :: ed2(plond,plev)

      INTEGER  :: jtg(plond) 
      INTEGER  :: jlcl(plond) 
      INTEGER  :: j0(plond) 
      INTEGER  :: jd(plond)

      REAL(dp) :: capelmt
      REAL(dp) :: delt

      INTEGER  :: i
      INTEGER  :: k
      INTEGER  :: lat
      INTEGER  :: lengath
      REAL(dp) :: psdiss
      REAL(dp) :: psevap
      REAL(dp) :: psheat
      REAL(dp) :: psrain
      REAL(dp) :: qdifr
      REAL(dp) :: qeff
      REAL(dp) :: sdifr
!!      REAL(dp) :: qpmax
      REAL(dp) :: tpmax
      REAL(dp) :: alpha8
      REAL(dp) :: pblhgt
      REAL(dp) :: fac1
      REAL(dp) :: tprime(plond,plev)
      REAL(dp) :: qprime(plond,plev)
      REAL(dp) :: tprg(plond,plev)
      REAL(dp) :: qprg(plond,plev)
!!      REAL(dp) :: tpertg(plond) 
!!      REAL(dp) :: qpertg(plond)
      REAL(dp) :: capblfc(plond)
!!      REAL(dp) :: caped(plond)
      REAL(dp) :: wprime
      INTEGER  :: ireset(plond)
      INTEGER  :: irg(plond)
  



!--------------------------Data statements------------------------------

      logical  :: momentm
      data capelmt/70._dp/
      momentm = .FALSE.

! initialize necessary arrays.
! zero out variables not used in ccm

      do i = 1,plon
         paprc(i) = 0._dp
         paprs(i) = 0._dp
      end do
      psdiss = 0._dp
      psheat = 0._dp
      psevap = 0._dp
      psrain = 0._dp
      lwc    = 0._dp
      iwc    = 0._dp
      rform  = 0._dp
      sform  = 0._dp
! initialize convective tendencies

      do k = 1,plev
         do i = 1,plon
            dqdt(i,k) = 0._dp
            dsdt(i,k) = 0._dp
            dudt(i,k) = 0._dp
            dvdt(i,k) = 0._dp
!            deltaq(i,k) = qh(i,k,1)
!            deltat(i,k) = t(i,k)
            mcon(i,k) = 0._dp
            pcpck(i,k) = 0._dp
            qtg(i,k) = 0._dp
            ttg(i,k) = 0._dp
            pflx(i,k) = 0._dp
            pflxg(i,k) = 0._dp
            cme(i,k) = 0._dp
! +bee
            cmfdqr(i,k) = 0._dp
            ql(i,k) = 0._dp
            qlg(i,k) = 0._dp
! -bee
            zdu(i,k) = 0._dp
            mup(i,k) = 0._dp
            mdn(i,k) = 0._dp
         end do
      end do
      do i = 1,plon
         pcpc(i) = 0._dp
         pflx(i,plevp) = 0._dp
         pflxg(i,plevp) = 0._dp
      end do
      if (.not.momentm) then
         do k = 1,plev
            do i = 1,plon
               u(i,k) = 0._dp
               v(i,k) = 0._dp
            end do
         end do
      end if

      do i = 1,plon
         pblt(i) = real(plev,dp)
         pcpr(i) = 0._dp
         pcps(i) = 0._dp
         dsubcld(i) = 0._dp
         sumq(i) = 0._dp
!         sumdt(i) = 0.
!         sumdq(i) = 0.
         jctop(i) = real(plev,dp)
         jcbot(i) = 1._dp
      end do

! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.

      do i = 1,plon
         zs(i) = geos(i)*rgrav
         pf(i,plev+1) = paph(i,plev+1)*0.01
         zf(i,plev+1) = zi(i,plev+1) + zs(i)
      end do
      do k = 1,plev
         do i = 1,plon
            p(i,k) = pap(i,k)*0.01
            pf(i,k) = paph(i,k)*0.01
            z(i,k) = zm(i,k) + zs(i)
            zf(i,k) = zi(i,k) + zs(i)
         end do
      end do

      do k = plev - 1,msg + 1,-1
         do i = 1,plon
            if (abs(z(i,k)-zs(i)-pblh(i)).lt.           &
                (zf(i,k)-zf(i,k+1))*0.5) pblt(i) =real(k,dp)
         end do
      end do

! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! convert from specific humidity (bounded by qmin) to mixing ratio.
! define dry static energy (normalized by cpg).

      do k = 1,plev
         do i = 1,plon
            qh0(i,k) = qh(i,k,1)
            qeff = max(qh(i,k,1),qmin)
            q(i,k) = qeff
            s(i,k) = t(i,k) + (grav/cpg)*z(i,k)
            tp(i,k)=0.0
            shat(i,k) = s(i,k)
            qhat(i,k) = q(i,k)
            wdp(i,k) = dpp(i,k)*0.01
            qg(i,k) = q(i,k)
            tg(i,k) = t(i,k)
            pg(i,k) = p(i,k)
            zg(i,k) = z(i,k)
            sg(i,k) = s(i,k)
            tpg(i,k) = tp(i,k)
            zfg(i,k) = zf(i,k)
            qstpg(i,k) = q(i,k)
            ug(i,k) = u(i,k)
            vg(i,k) = v(i,k)
            dlg(i,k) = 0.
            dlf(i,k) = 0.
!!            cu(i,k) = 0.
!!            evp(i,k) = 0.
         end do
      end do
      do i = 1,plon
         zfg(i,plev+1) = zf(i,plev+1)
         ireset(i) = -1
         capeg(i) = 0.
         lclg(i) = 1
         lelg(i) = plev
         maxg(i) = 1
         tlg(i) = 400._dp
         dsubcld(i) = 0.
!!         qdb(i) = 0.
!!         sdb(i) = 0.
         betau(i) = 0.
         betad(i) = 0.
      end do

!     Calculate tprime and qprime as in cmfmc **JP/PJR

      tpmax = 1.5_dp
!!      qpmax = 5.e-4_dp
!!      qpmax = 1.5e-3_dp
!!      alpha8 = 2._dp        ! number of bl depths to propogate perturbations
      alpha8 = 1._dp        ! number of bl depths to propogate perturbations

      do k=1,plev
         do i=1,plon
            pblhgt = max(pblh(i),1.0_dp)
            fac1   = max(0.0_dp,1.0_dp-zm(i,k)*rgrav/(alpha8*pblhgt))
            tprime(i,k) = min(tpert(i),tpmax)*fac1
#ifdef QPERT
!!            qprime(i,k) = min(qpert(i),qpmax)*fac1
#else
            qprime(i,k) = 0._dp
#endif
!           write(*,*) tprime(i,k),qprime(i,k),k
         enddo
      enddo

      
! evaluate covective available potential energy (cape).

      call buoyan_pjr(plon,plond   ,plev    ,plevp   ,           &
                  q       ,t       ,p       ,z       ,pf      ,  &
                  tp      ,qstp    ,tl      ,cape    ,           &
                  lcl     ,lel     ,lfc     ,maxi    ,           &
                  lat     ,tprime  ,qprime  ,capblfc )

      if (trigon) then
! +pjr
      do i = 1,plon
         wprime = wpert3d(i,maxi(i))
!!         caped(i) = cape(i)
         if (max(0._dp,wprime)**2 + 2._dp*capblfc(i).lt.0._dp) then
            cape(i) = -.0001_dp
         endif
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         wprime = wpert3d(i,maxi(i))
         write (6,*) 'conv_ccm_pjr: maxi, wprime, capblfc ',           &
              maxi(i), wprime, capblfc(i),                             &
              max(0._dp,wprime)**2 + 2*capblfc(i), wpert3d(i,plev) 
         write (6,*) ' caped, cape, lcl, lel, lfc, maxi ',             &
              caped(i), cape(i), lcl(i), lel(i), lfc(i), maxi(i)
      endif
#endif
! -pjr
      endif

! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).

         call whenfgt(plon,cape,1,capelmt,ideep,lengath)
         if (lengath.eq.0) return
         
#ifdef DEBCONV
         if (lat.eq.latlook) then
            ilookg = 0
            do i = 1,lengath
               if (ideep(i) .eq. ilook) then
                  ilookg = i
               endif
            end do
            write (6,*) ' ilook ', ilook,                          &
                ' corresponds to gather index', ilookg
         endif
#endif
!        jyes = 0
!        jno = plon - 1 + 2
!        do il = 1,plon
!           if (cape(il).gt.capelmt) then
!              jyes = jyes + 1
!              ideep(jyes) = il
!           else
!              jno = jno - 1
!              ideep(jno) = il
!           end if
!        end do
!        lengath = jyes
!        if (lengath.eq.0) return

! obtain gathered arrays necessary for ensuing calculations.

      do k = 1,plev
         do i = 1,lengath
            wdp(i,k) = 0.01*dpp(ideep(i),k)
            qg(i,k) = q(ideep(i),k)
            tg(i,k) = t(ideep(i),k)
            pg(i,k) = p(ideep(i),k)
            zg(i,k) = z(ideep(i),k)
            sg(i,k) = s(ideep(i),k)
            tpg(i,k) = tp(ideep(i),k)
            zfg(i,k) = zf(ideep(i),k)
            qstpg(i,k) = qstp(ideep(i),k)
            ug(i,k) = u(ideep(i),k)
            vg(i,k) = v(ideep(i),k)
            tprg(i,k) = tprime(ideep(i),k)
            qprg(i,k) = qprime(ideep(i),k)
         end do
      end do

      do i = 1,lengath
         zfg(i,plev+1) = zf(ideep(i),plev+1)
      end do
      do i = 1,lengath
         capeg(i) = cape(ideep(i))
         lclg(i) = lcl(ideep(i))
         lelg(i) = lel(ideep(i))
         maxg(i) = maxi(ideep(i))
         tlg(i) = tl(ideep(i))
!!         tpertg(i) = tpert(ideep(i))
!!         qpertg(i) = qpert(ideep(i))
      end do

! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.

      do k = msg + 1,plev
         do i = 1,lengath
            if (k.ge.maxg(i)) then
               dsubcld(i) = dsubcld(i) + wdp(i,k)
            end if
         end do
      end do

! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.

      do k = msg + 2,plev
         do i = 1,lengath
!!            alpha(i,k) = 0.5_dp
            sdifr = 0.
            qdifr = 0.
            if (sg(i,k).gt.0. .or. sg(i,k-1).gt. 0.)                 &
                sdifr = abs((sg(i,k)-sg(i,k-1))/                     &
                max(sg(i,k-1),sg(i,k)))
            if (qg(i,k).gt.0. .or. qg(i,k-1).gt.0.)                  &
                qdifr = abs((qg(i,k)-qg(i,k-1))/                     &
                max(qg(i,k-1),qg(i,k)))
            if (sdifr.gt.1.E-6_dp) then
               shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/ &
                           (sg(i,k-1)-sg(i,k))
            else
               shat(i,k) = 0.5_dp * (sg(i,k)+sg(i,k-1))
            end if
            if (qdifr.gt.1.E-6_dp) then
               qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/ & 
                           (qg(i,k-1)-qg(i,k))
            else
               qhat(i,k) = 0.5_dp * (qg(i,k)+qg(i,k-1))
            end if
         end do
      end do

! obtain cloud properties.

      call cldprp_pjr(plond,   plev  ,plevp   ,                     &
                  qg      ,tg      ,ug      ,vg      ,pg      ,     &
                  zg      ,sg      ,mu      ,eu      ,du      ,     &
                  md      ,ed      ,sd      ,qd      ,ud      ,     &
                  vd      ,mc      ,qug     ,sug     ,zfg     ,     &
                  qs      ,hmn     ,hsat    ,shat    ,              &
                  qlg     ,totpcp  ,totevp  ,cmeg    ,maxg    ,     &
                  lelg    ,jtg     ,jlcl    ,maxg    ,j0      ,     &
                  jd      ,lengath ,                                &
                  pflxg   ,evpg    ,cug     ,mu2     ,eu2     ,     &
                  du2     ,md2     ,ed2     ,cmfdqrg ,              &
                  tprg    ,qprg    ,irg     ,                       &
                  lwcg    ,iwcg    ,rformg  ,sformg  )

! determine cloud base mass flux.

      do i = 1,lengath
!!         qdb(i) = qd(i,maxg(i))
!!         sdb(i) = sd(i,maxg(i))
         betad(i) = md(i,maxg(i))
         betau(i) = mu(i,maxg(i))
      end do

! convert detrainment from units of "1/m" to "1/mb".

      do k = msg + 1,plev
         do i = 1,lengath
            du(i,k) = du(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            eu(i,k) = eu(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            ed(i,k) = ed(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            cug(i,k) = cug(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
! +bee, bug fix from /home/pjr/ccm2/omega0.10.1/fspj01/.
            cmeg(i,k) = cmeg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            cmfdqrg(i,k) = cmfdqrg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
! -bee
            evpg(i,k) = evpg(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            du2(i,k) = du2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            eu2(i,k) = eu2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
            ed2(i,k) = ed2(i,k)* (zfg(i,k)-zfg(i,k+1))/wdp(i,k)
         end do
      end do

      call closure_pjr(plond,   plev  ,                             &
                   qg      ,tg      ,pg      ,sg      ,             &
                   tpg     ,qug     ,sug     ,mc      ,             &
                   du      ,mu      ,md      ,qd      ,sd      ,    &
                   qhat    ,shat    ,wdp     ,qstpg   ,             &
                   zfg     ,qlg     ,dsubcld ,mb      ,capeg   ,    &
                   tlg     ,lclg    ,lelg    ,jtg     ,maxg    ,    &
                   1       ,lengath ,                               &
                   capelmt )

! limit cloud base mass flux to theoretical upper bound.

      do i=1,lengath
        mumax(i) = 0.
      end do
      do k=msg + 2,plev
        do i=1,lengath
          mumax(i) = max(mumax(i), mu(i,k)/wdp(i,k))
        end do
      end do
      do i=1,lengath
        if (mumax(i).gt.0.) then
          mb(i) = min(mb(i),0.5/(delt*mumax(i)))
        else
          mb(i) = 0.
        endif
      end do
      do k=msg+1,plev
        do i=1,lengath
          mu(i,k) = mu(i,k)*mb(i)
          md(i,k) = md(i,k)*mb(i)
          mc(i,k) = mc(i,k)*mb(i)
          du(i,k) = du(i,k)*mb(i)
          eu(i,k) = eu(i,k)*mb(i)
          ed(i,k) = ed(i,k)*mb(i)
          cmeg(i,k) = cmeg(i,k)*mb(i)
! +bee
          cmfdqrg(i,k) = cmfdqrg(i,k)*mb(i)
! -bee
          cug(i,k) = cug(i,k)*mb(i)
          evpg(i,k) = evpg(i,k)*mb(i)
          pflxg(i,k+1) = pflxg(i,k+1)*mb(i)*100./grav
          mu2(i,k) = mu2(i,k)*mb(i)
          md2(i,k) = md2(i,k)*mb(i)
          du2(i,k) = du2(i,k)*mb(i)
          eu2(i,k) = eu2(i,k)*mb(i)
          ed2(i,k) = ed2(i,k)*mb(i)
        end do
      end do
      do i = 1,lengath
         betau(i) = betau(i)*mb(i)
         betad(i) = betad(i)*mb(i)

! totpcp from cldprp has the dimension of kg/kg, here it is 
! converted to kg/(m^2*s), the precipitation rate

         totpcp(i) = totpcp(i)*mb(i)*100./grav
         totevp(i) = totevp(i)*mb(i)*100./grav
      end do

! compute temperature and moisture changes due to convection.

      call q1q2pjr2(plond   ,plev    ,                                &
                    dqdt    ,dsdt    ,                                &
                    qug     ,sug     ,du      ,                       &   
                    qhat    ,shat    ,wdp     ,mu      ,md      ,     &
                    sd      ,qd      ,qlg     ,dsubcld          ,     &
                    jtg     ,maxg    ,1       ,lengath ,              &
                    dlg     ,evpg    ,cug   )

! compute momentum changes due to convection, if desired (i.e
! if logical switch set).

!
!      if(momentm)                                                   
!       then
!        call moment(dudt,dvdt,du,alpha,wdp,ed,eu,mc,md,mu,          &
!                   pg,qd,qu,qhat,sd,su,shat,ud,vd,tg,ug,vg,zg,zfg,  &
!                   dsubcld,maxg,jd,jt,                              &
!                   msg,2.*delt,plev,1,lengath,plond,lat)
!      endif

! +bee, move convective transports outside conv_ccm
!      if (pcnst.gt.1) then
!         call convtran_pjr
!     $                (qh      ,mu2     ,md2     ,du2     ,eu2     ,
!     $                 ed2     ,wdp     ,dsubcld ,jtg     ,maxg    ,
!     $                 ideep   ,1       ,lengath ,nstep   ,lat     ,
!     $                 delt    )
!      endif
! -bee

! gather back temperature and mixing ratio.

      do k = msg + 1,plev
         do i = 1,lengath
            psdiss = psdiss + (dudt(i,k)*u(ideep(i),k)+                  &
                     dvdt(i,k)*v(ideep(i),k))*dpp(ideep(i),k)/grav

! q is updated to compute net precip, and then reset to old value.
! the last line is overwritten. so the input basic variables, i.e.
! q, t, u and v are updated to include convective increments. 
! (5/2/95)

            q(ideep(i),k) = q(ideep(i),k) + 2.*delt*dqdt(i,k)
            t(ideep(i),k) = t(ideep(i),k) + 2.*delt*dsdt(i,k)          
            u(ideep(i),k) = u(ideep(i),k) + 2.*delt*dudt(i,k)
            v(ideep(i),k) = v(ideep(i),k) + 2.*delt*dvdt(i,k)
            cme(ideep(i),k) = cmeg(i,k)
! +bee
            cmfdqr(ideep(i),k) = cmfdqrg(i,k)
            ql(ideep(i),k) = qlg(i,k)
! -bee
            zdu(ideep(i),k) = du2(i,k)
            mup(ideep(i),k) = mu(i,k)
            mdn(ideep(i),k) = md(i,k)
            mcon(ideep(i),k) = mc(i,k)
            qtg(ideep(i),k) = dqdt(i,k)
            ttg(ideep(i),k) = dsdt(i,k)
!            utg(ideep(i),k) = dudt(i,k)
!            vtg(ideep(i),k) = dvdt(i,k)
            dlf(ideep(i),k) = dlg(i,k)
            pflx(ideep(i),k) = pflxg(i,k)
!!            cu(ideep(i),k) = cug(i,k)
!!            evp(ideep(i),k) = evpg(i,k)
            ! mz_ht_20070918+            
            lwc(ideep(i),k) = lwcg(i,k)
            iwc(ideep(i),k) = iwcg(i,k)
            rform(ideep(i),k) = rformg(i,k)
            sform(ideep(i),k) = sformg(i,k)
            ! mz_ht_20070918-
         end do
      end do

      do i = 1,lengath
        ireset(ideep(i)) = irg(i)
         jctop(ideep(i)) = real(jtg(i),dp)
! +bee
         jcbot(ideep(i)) = real(maxg(i),dp)
! -bee
         pflx(ideep(i),plevp) = pflxg(i,plevp)
         psevap = psevap + totevp(i)
         psrain = psrain + totpcp(i)
      end do

! convert back to specific humidity from mixing ratio.
! take into account any moisture added to ensure positiveness
! of specific humidity at start of routine.

      do k = msg + 1,plev
         do i = 1,plon
            qh(i,k,1) = q(i,k)
            qh(i,k,1) = qh(i,k,1) - max((qmin-qh0(i,k)),0._dp)
         end do
      end do
      do k = plev,msg + 1,-1
         do i = 1,plon
            sumq(i) = sumq(i) - dpp(i,k)* (qh(i,k,1)-qh0(i,k))

! account for the detraining cloud water in the precip 

            sumq(i) = sumq(i) - dpp(i,k)*dlf(i,k)*2._dp*delt
            pcpck(i,k) = max(0._dp,sumq(i))
         end do
      end do

! obtain final precipitation rate.

      do i = 1,plon
!         llo1 = ts(i) .ge. tfreez

! here pcpr and pcps are in units of kg/m^2, ie. precip per
! time step

!         pcpr(i) = cvmgt(pcpdh(i)*max(sumq(i),0.),0.,llo1)
!         pcps(i) = cvmgt(0.,pcpdh(i)*max(sumq(i),0.),llo1)
         precip = rgrav*max(sumq(i),0._dp)
         if (ts(i) .ge. tfreez) then
           pcpr(i) = precip
           pcps(i) = 0._dp
         else
           pcpr(i) = 0.
           pcps(i) = precip
         end if
      end do


! accumulate precipitation, the 1000. is the density of water, so
! paprc and paprs are now in units of meters.

      do i = 1,plon
         paprc(i) = paprc(i) + (pcpr(i)+pcps(i))/1000.
         paprs(i) = paprs(i) + pcps(i)/1000.
      end do

! convert precipitation to m/s, ie, precip rate.

      do i = 1,plon
         pcpr(i) = pcpr(i)/ (2.*delt)/1000.
         pcps(i) = pcps(i)/ (2.*delt)/1000.
         pcpc(i) = pcpr(i) + pcps(i)
         psheat = psheat + (pcps(i)+pcpr(i))*rl
      end do
      do k = msg + 1,plev
         do i = 1,plon
            pcpck(i,k) = rgrav*pcpck(i,k)/ (2.*delt)
         end do
      end do

! calculate conservation of quantities.


!       if(lat.eq.jlatpr)then
!        do l=msg+1,plev
!        do i=1,lengath
!          sumdq(i) = sumdq(i) + 2.*delt*(rl/cp)*dpp(ideep(i),l)* dqdt(i,l)
!          sumdt(i) = sumdt(i) + 2.*delt*dpp(ideep(i),l)*dsdt(i,l)
!        end do
!        end do

!        print *,'sumdq,sumdt,sumde in convection subroutine########'
!        do i=1,lengath
!          sumde(i) = sumdt(i) + sumdq(i)
!          write(6, 901) sumdq(i), sumdt(i),sumde(i), i, ideep(i)
!        end do

!        print *,'sumdq,sumdt,sumde ... all points'
!      do i=1,plon
!         sumdq(i) = 0.0
!         sumdt(i) = 0.0
!      end do

!      do l=msg+1,plev
!      do i=1,plon
!        deltaq(i,l) = qh(i,l,1) - deltaq(i,l)
!        deltat(i,l) = t (i,l) - deltat(i,l)
!        sumdq(i) = sumdq(i) + (rl/cp)*dpp(i,l)*deltaq(i,l)
!        sumdt(i) = sumdt(i) + dpp(i,l)*deltat(i,l)
!      end do
!      end do
!      do i=1,plon
!        sumde(i) = sumdt(i) + sumdq(i)
!      end do
!        write(6, 902) (i,sumdq(i),sumdt(i),sumde(i),i=1,plon)

!  901   format(1x,3e20.12, i10, i10)
!  902   format(1x,i10, 3e20.12)

!      endif

      return
    end subroutine conv_ccm_pjr
!============================================================================================
!------------------ subroutines for conv_ccm_pjr start --------------------------------------

    subroutine convtran_pjr                                       &
                         (q     ,mu     ,md      ,du     ,eu     ,&
                          ed    ,wdp    ,jt      ,mx             ,&
                          ideep ,il1g   ,il2g                    ,&
                          delt  ,plond  ,plev    ,pcnst)

    
!-----------------------------------------------------------------------

! Convective transport of trace species

! Note that we are still assuming that the tracers are in a moist mixing ratio
! this will change soon

!-------------------------Code History----------------------------------

! Original version:  P. Rasch, Jan 1996 
! Standardized:      L. Buja,  Feb 1996
! Reviewed:          P. Rasch, Feb 1996      
 
!-----------------------------------------------------------------------
      implicit none

!-----------------------------Arguments---------------------------------
 
! Input
      INTEGER, INTENT(IN) :: plond, &  ! number of columns
                             plev,  &  ! number of levels  
                             pcnst     ! number of tracers
      REAL(dp) :: mu(plond,plev)       ! Mass flux up
      REAL(dp) :: md(plond,plev)       ! Mass flux down
      REAL(dp) :: du(plond,plev)       ! Mass detraining from updraft
      REAL(dp) :: eu(plond,plev)       ! Mass entraining from updraft
      REAL(dp) :: ed(plond,plev)       ! Mass entraining from downdraft
      REAL(dp) :: wdp(plond,plev)       ! Delta pressure between interfaces

      INTEGER  :: jt(plond)         ! Index of cloud top for each column
      INTEGER  :: mx(plond)         ! Index of cloud top for each column
      INTEGER  :: ideep(plond)      ! Gathering array
      INTEGER  :: il1g              ! Gathered min lon indices over which to operate
      INTEGER  :: il2g              ! Gathered max lon indices over which to operate

      REAL(dp) :: delt                 ! Time step

! input/output

      REAL(dp) :: q(plond,plev,pcnst)  ! Tracer array including moisture

!--------------------------Local Variables------------------------------

      INTEGER  :: i                 ! Work index
      INTEGER  :: k                 ! Work index
      INTEGER  :: kbm               ! Highest altitude index of cloud base
      INTEGER  :: kk                ! Work index
      INTEGER  :: kkp1              ! Work index
      INTEGER  :: km1               ! Work index
      INTEGER  :: kp1               ! Work index
      INTEGER  :: ktm               ! Highest altitude index of cloud top
      INTEGER  :: m                 ! Work index

      REAL(dp) :: cabv                 ! Mix ratio of constituent above
      REAL(dp) :: cbel                 ! Mix ratio of constituent below
      REAL(dp) :: cdifr                ! Normalized diff between cabv and cbel
      REAL(dp) :: chat(plond,plev)     ! Mix ratio in env at interfaces
      REAL(dp) :: cond(plond,plev)     ! Mix ratio in downdraft at interfaces
      REAL(dp) :: const(plond,plev)    ! Gathered tracer array 
      REAL(dp) :: conu(plond,plev)     ! Mix ratio in updraft at interfaces
      REAL(dp) :: dcondt(plond,plev)   ! Gathered tend array 
      REAL(dp) :: small                ! A small number
      REAL(dp) :: mbsth                ! Threshold for mass fluxes
      REAL(dp) :: mupdudp              ! A work variable
      REAL(dp) :: minc                 ! A work variable
      REAL(dp) :: maxc                 ! A work variable
      REAL(dp) :: qn                   ! A work variable
      REAL(dp) :: fluxin               ! A work variable
      REAL(dp) :: fluxout              ! A work variable
      REAL(dp) :: netflux              ! A work variable

!-----------------------------------------------------------------------

      small = 1.e-36_dp
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
      mbsth = 1.e-15_dp

! Find the highest level top and bottom levels of convection
      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do

! Loop ever each constituent
      do m = 2,pcnst

! Gather up the constituent and set tend to zero
         do k = 1,plev
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,plev
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc.lt.0._dp) then
                  cdifr = 0.
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr.gt.1.E-6_dp) then
                  cabv = max(const(i,km1),maxc*1.e-12)
                  cbel = max(const(i,k),maxc*1.e-12)
                  chat(i,k) = log(cabv/cbel)               &
                                /(cabv-cbel)               &
                                *cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0.

            end do
         end do

#ifdef CHECKNEG
         do k = 1,plev
            km1 = max(1,k-1)
            do i = il1g, il2g
               if (chat(i,k).lt.0.) then
                  write (6,*) ' negative chat ', i, k, lat, chat(i,k),   &
                       const(i,km1), const(i,k)
!                  stop
               endif
            end do
         end do
#endif


! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = plev 
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*wdp(i,kk)
            if (mupdudp.gt.mbsth) then
               conu(i,kk) = (                                    &
                             +eu(i,kk)*const(i,kk)*wdp(i,kk)      &
                            )/mupdudp
            endif
            if (md(i,k).lt.-mbsth) then
               cond(i,k) =  (                                    &
                             -ed(i,km1)*const(i,km1)*wdp(i,km1)   &
                            )/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = plev-1,1,-1
            kkp1 = min(plev,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*wdp(i,kk)
               if (mupdudp.gt.mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)        &
                                 +eu(i,kk)*const(i,kk)*wdp(i,kk)  &
                               )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,plev
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k).lt.-mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)              &
                                -ed(i,km1)*const(i,km1)*wdp(i,km1)   &
                               )/md(i,k)
               endif
            end do
         end do


         do k = ktm,plev
            km1 = max(1,k-1)
            kp1 = min(plev,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =                                      &
!                       +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))    &
!                         -mu(i,k)*   (conu(i,k)-chat(i,k))        &
!                         +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))    &
!                         -md(i,k)*   (cond(i,k)-chat(i,k))        &
!                        )/wdp(i,k)

! version 2 har to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)      &
!                       -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)  & 
!                       -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)                         &
                       + mu(i,k)*min(chat(i,k),const(i,km1))           &
                       -(md(i,k)  *cond(i,k)                           &
                       + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k)                             &
                        +mu(i,kp1)*min(chat(i,kp1),const(i,k))         &
                       -(md(i,kp1)*cond(i,kp1)                         &
                       + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                  netflux = 0.
               endif
               dcondt(i,k) = netflux/wdp(i,k)
            end do
         end do

         do k = kbm,plev             
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k.eq.mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*                       &
!                   (-mu(i,k)*(conu(i,k)-chat(i,k))                     &
!                    -md(i,k)*(cond(i,k)-chat(i,k))                     &
!                   )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1))        &
                          - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k)                          & 
                          - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux).lt.max(fluxin,fluxout)*1.e-12) then
                     netflux = 0.
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/wdp(i,k)
               else if (k.gt.mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0.
               end if
            end do
         end do

! Update and scatter data back to full arrays

         do k = 1,plev
            kp1 = min(plev,k+1)
            do i = il1g,il2g
               qn = const(i,k)+dcondt(i,k)*2.*delt
               q(ideep(i),k,m) = qn

            end do
         end do
         
      end do                    ! m = 2,pcnst

      return
    end subroutine convtran_pjr
!===============================================================================

      subroutine q1q2pjr2 (  plond  ,plev,                            &
                      dqdt    ,dsdt                               ,   & 
                      qu      ,su      ,du                        ,   &
                      qhat    ,shat    ,wdp     ,mu      ,md      ,   &
                      sd      ,qd      ,ql      ,dsubcld          ,   &
                      jt      ,mx      ,il1g    ,il2g             ,   &
                      dl      ,evp     ,cu     )

      implicit none


! rewritten by phil rasch dec 19 1995
      INTEGER, INTENT(IN) :: plond, plev      ! nmber of columns, levels


      REAL(dp) :: dqdt(plond,plev),    &
                  dsdt(plond,plev)

      REAL(dp) :: fact
      INTEGER  :: il1g
      INTEGER  :: il2g
      INTEGER  :: i
      INTEGER  :: k
      REAL(dp) :: qu(plond,plev),     &
                  su(plond,plev),     &
                  du(plond,plev),     &
                  qhat(plond,plev),   &
                  shat(plond,plev),   &
                  wdp(plond,plev),    &
                  mu(plond,plev),     &
                  md(plond,plev),     &
                  sd(plond,plev),     &
                  qd(plond,plev),     &
                  ql(plond,plev)  
      REAL(dp) :: dl(plond,plev),     &
                  evp(plond,plev),    &
                  cu(plond,plev)
      INTEGER  :: kbm
      INTEGER  :: ktm
      REAL(dp) :: emc
      REAL(dp) :: dsubcld(plond)
                

      INTEGER  :: jt(plond),       &
                  mx(plond)

! work fields:

!-------------------------------------------------------------------
      do k = msg + 1,plev
         do i = il1g,il2g
            dsdt(i,k) = 0.
            dqdt(i,k) = 0.
            dl(i,k) = 0.
         end do
      end do

! find the highest level top and bottom levels of convection

      ktm = plev
      kbm = plev
      do i = il1g, il2g
         ktm = min(ktm,jt(i))
         kbm = min(kbm,mx(i))
      end do

      fact = 0.
!      fact = 1.

      do k = ktm,plev-1
         do i = il1g,il2g
!            fact = 1.
!            if (q(i,k).gt.0.8*qs(i,k) .and. k.lt.plev-3) fact = 0.

            emc = +fact*du(i,k)*ql(i,k+1) & ! evaporating cloud detraining to env
                  -cu(i,k)                & ! condensation in updraft
                  +evp(i,k)                 ! evaporating rain in downdraft
!            emc = 0

            dsdt(i,k) = -rl/cpg*emc                             &
                        + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1))  &
                           -mu(i,k)*   (su(i,k)-shat(i,k))      &
                           +md(i,k+1)* (sd(i,k+1)-shat(i,k+1))  &
                           -md(i,k)*   (sd(i,k)-shat(i,k))      &
                          )/wdp(i,k)
          

            dqdt(i,k) = emc                                     & 
                        +(+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1))   &
                          -mu(i,k)*   (qu(i,k)-qhat(i,k))       &
                          +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1))   &
                          -md(i,k)*   (qd(i,k)-qhat(i,k))       &
                         )/wdp(i,k)
                       

            dl(i,k) = (1._dp - fact)*du(i,k)*ql(i,k+1)

         end do
      end do


      do k = kbm,plev             
         do i = il1g,il2g
            if (k.eq.mx(i)) then
               dsdt(i,k) = (1./dsubcld(i))*               &
                    (-mu(i,k)* (su(i,k)-shat(i,k))        &
                     -md(i,k)* (sd(i,k)-shat(i,k))        &
                    )
               dqdt(i,k) = (1./dsubcld(i))*               &
                    (-mu(i,k)*(qu(i,k)-qhat(i,k))         &
                     -md(i,k)*(qd(i,k)-qhat(i,k))         &
                    )
            else if (k.gt.mx(i)) then
               dsdt(i,k) = dsdt(i,k-1)
               dqdt(i,k) = dqdt(i,k-1)
            end if
         end do
      end do

      return
    end subroutine q1q2pjr2
!==========================================================================

    subroutine buoyan_pjr(plon  ,plond   ,plev    ,plevp   ,          &
                        q       ,t       ,p       ,z       ,pf      , &
                        tp      ,qstp    ,tl      ,cape    ,          & 
                        lcl     ,lel     ,lfc     ,mx      ,          &
                        lat     ,tprime  ,qprime  ,capblfc )
!-----------------------------------------------------------------------
! this subroutine works on UNGATHERED ARRAYS
!-----------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.

! the documentation has been enhanced to the degree that we are able

! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
! Revised            P. Rasch and Jon Petch 1997

! ----------------------------------------------------------------------

! jul 14/92 - guang jun zhang, m.lazare, n.mcfarlane.  as in
!             previous version buoyan except remove pathalogical
!             cases of "zig-zags" in profiles where lel defined
!             far too high (by use of lelten array, which assumes
!             a maximum of five such crossing points).
! feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
!             version buoyan.
      IMPLICIT NONE
! input arguments
      INTEGER, INTENT(IN) :: plon,plond,&   ! number of columns, max columns
                             plev, plevp    ! number of  levels, levels+1

      REAL(dp) :: q(plond,plev)        ! spec. humidity
      REAL(dp) :: t(plond,plev)        ! temperature
      REAL(dp) :: p(plond,plev)        ! pressure
      REAL(dp) :: z(plond,plev)        ! height
      REAL(dp) :: pf(plond,plev+1)     ! pressure at interfaces

! output arguments

      REAL(dp) :: tp(plond,plev)       ! parcel temperature
      REAL(dp) :: qstp(plond,plev)     ! saturation mixing ratio of parcel
      REAL(dp) :: tl(plond)            ! parcel temperature at lcl
      REAL(dp) :: cape(plond)          ! convective aval. pot. energy.
      INTEGER  :: lcl(plond)        ! lifting condensation level
      INTEGER  :: lel(plond)        ! cloud top
!     INTEGER  :: lon(plond)        ! level of onset of deep convection
      INTEGER  :: mx(plond)         ! level of max moist static energy

! -------------------------Local Variables------------------------------

      REAL(dp) :: capeten(plond,5)     ! provisional value of cape
      REAL(dp) :: tv(plond,plev)       ! 
      REAL(dp) :: tpv(plond,plev)      ! 
      REAL(dp) :: buoy(plond,plev)

      REAL(dp) :: a1(plond) 
      REAL(dp) :: a2(plond) 
      REAL(dp) :: estp(plond) 
      REAL(dp) :: pl(plond) 
      REAL(dp) :: plexp(plond) 
!!      REAL(dp) :: hmax(plond) 
      REAL(dp) :: y(plond)

      INTEGER  :: knt(plond) 
      INTEGER  :: lelten(plond,5)

      REAL(dp) :: e

      INTEGER  :: i
      INTEGER  :: k
      INTEGER  :: lat
      INTEGER  :: n
      INTEGER  :: lnb(plond)

      REAL(dp) :: tprime(plond,plev)
      REAL(dp) :: qprime(plond,plev)
      INTEGER  :: lfc(plond)
      REAL(dp) :: capblfc(plond)
      REAL(dp) :: tparc                ! temperature of the parcel at a given level
      REAL(dp) :: qparc                ! spec hum of the parcel at given level
!!      REAL(dp) :: off                  ! value of negative bouyance above which the parcel cannot rise
      REAL(dp) :: bloc
      INTEGER  :: numbllcl
      INTEGER  :: maxlelten

! ----------------------------------------------------------------------
!!      off = -4._dp                ! value of buoy to switch off at
! ----------------------------------------------------------------------


      do n = 1,5
        do i = 1,plon
          lelten(i,n) = plev
          capeten(i,n) = 0.
        end do
      end do

      do k = 1,plev
         do i = 1,plon
            buoy(i,k) = 0.0
         end do
      end do

      do i = 1,plon
!        lon(i) = plev
        knt(i) = 0
        lel(i) = plev
        mx(i) = plev
        cape(i) = 0.
!!        hmax(i) = 0.
      end do

!     find the level (mx) of most buoyant parcel in each column, 
!     and the lcl associated with it, and the convective inhibition 
!     of that parcel
      call findcin (plon, plond, plev, plevp, lat, t, q, tprime, qprime, mx,  &
                    p, pf, z, lcl, capblfc, lfc, lnb, cape)
#ifdef DEBCONV
      if (lat.eq.latlook) then
         write (6,*) ' at findcin ', lat
         i = ilook
         write (6,*) ' findcin returned mx, lcl, capblfc '//      &
              ' lfc, lnb, cape ',                                 &
             mx(i), lcl(i), capblfc(i), lfc(i), lnb(i), cape(i)
      endif
#endif

!      if (.false.) then
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' lat, i, pblt ', lat, i, pblt(i)
         write (6,*) ' buoyan_pjr: t, z, q, qs, rh, h, ha, hpen '
         do k = plev,plev-5,-1
            hmn(i) = cpg*t(i,k) + grav*z(i,k) + rl*q(i,k)
!            e = c1*exp((c2* (t(i,k)-tfreez))/     &
!                      ((t(i,k)-tfreez)+c3))
!            qs = eps1*e/ (p(i,k)-e)
!            esn = estblf(t(i,k))
!            qs2 = epsqs*esn/(p(i,k)*100.-esn)
!            ha = cpg*t(i,k) + grav*z(i,k) + rl*(q(i,k)-qs)
!            ha =  hmn(i)+rl*(q(i,k)-qs)
!            hb = cpg*t(i,k-1) + grav*z(i,k-1) + rl*q(i,k-1)
!            write (6,77) k, t(i,k),           &
!                 z(i,k), q(i,k),              & 
!                 qs, q(i,k)/qs, hmn(i), hb,  ha
!   77       format (i5,2f10.3,1p,6e15.5)
!            if (hb .le. hmn(i) .and. q(i,k).ge.qs*0.99) then
!               write (6,*) ' unconditionally unstable profile '
!            endif
         end do
      endif
#endif

      do i = 1,plon
        knt(i) = 0
        lel(i) = plev
        cape(i) = 0.
!!        hmax(i) = 0.
      end do


! initialize parcel properties in sub-cloud layer below lcl.

      do k = plev,msg + 1,-1
         do i=1,plon
            tv(i,k) = t(i,k)* (1._dp + 1.608_dp*q(i,k))/ (1._dp + q(i,k))
            if (k.ge.lcl(i) .and. k.le.mx(i)) then
               tparc = t(i,mx(i)) + tprime(i,mx(i))
               qparc = q(i,mx(i)) + qprime(i,mx(i))
               qstp(i,k) = qparc
               tp(i,k)   = tparc*(p(i,k)/p(i,mx(i)))**   &
                          (0.2854*(1.0_dp - 0.28*qparc))
               tpv(i,k)  = tp(i,k)*(1._dp + 1.608*qparc)/    &
                          (1._dp + qparc)
               buoy(i,k) = tpv(i,k) - tv(i,k)
            else
               tp(i,k) = t(i,k) + tprime(i,k)
               qparc = q(i,k) + qprime(i,k)
               tpv(i,k)  = tp(i,k)*(1._dp + 1.608*qparc)/    &
                          (1._dp + qparc)
               buoy(i,k) = 0.
            end if
        end do
      end do
!     
!     Define parcel properties at lcl (i.e. level immediately above pl).

      do i=1,plon
         if (lcl(i).ne.0) then
         k = lcl(i)
         e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
         tl(i) = 2840./ (3.5*log(t(i,mx(i)))-log(e)-4.805_dp) + 55._dp
         if (tl(i).lt.t(i,mx(i))) then
            plexp(i) = (1./ (0.2854* (1._dp - 0.28*q(i,mx(i)))))
            pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)
         else
            tl(i) = t(i,mx(i))
            pl(i) = p(i,mx(i))
         end if
!        tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
         qstp(i,k) = q(i,mx(i)) + qprime(i,mx(i))

         tp(i,k) = tl(i)* (p(i,k)/pl(i))**          &
                     (0.2854* (1._dp - 0.28*qstp(i,k)))

         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/   &
                     ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))   

         a1(i) = cpg/rl + qstp(i,k)* (1._dp + qstp(i,k)/eps1)   &
             *rl*eps1/ (rgas*tp(i,k)**2)

         a2(i) = 0.5_dp * (qstp(i,k)* (1._dp + 2._dp/eps1*qstp(i,k))*   &
             (1._dp + qstp(i,k)/eps1)*eps1**2.*rl*rl/                   &
             (rgas**2*tp(i,k)**4)-qstp(i,k)*                            &
             (1._dp + qstp(i,k)/eps1)*2.*eps1*rl/                       &
             (rgas*tp(i,k)**3))
        
         a1(i) = 1./a1(i)
         a2(i) = -a2(i)*a1(i)**3
         y(i)    = q(i,mx(i)) + qprime(i,mx(i)) - qstp(i,k)
         tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2

!     estp(i)  =exp(a-b/tp(i,k))

         estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/          &
             ((tp(i,k)-tfreez)+c3))

         qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))

         tpv(i,k) = (tp(i,k))* (1._dp + 1.608*qstp(i,k))/       &
              (1._dp + q(i,mx(i)) + qprime(i,mx(i)))
         buoy(i,k) = tpv(i,k) - tv(i,k) 
      endif
      end do

! main buoyancy calculation.

      do k = plev - 1,msg + 1,-1
        do i=1,plon
          if (k.lt.lcl(i)) then

!            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**         &
                     (0.2854* (1._dp - 0.28*qstp(i,k)))

!          estp(i) = exp(a-b/tp(i,k))

            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/         &
                     ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))

            a1(i) = cpg/rl + qstp(i,k)* (1._dp + qstp(i,k)/eps1)  &
                   *rl*eps1/ (rgas*tp(i,k)**2)

            a2(i) = 0.5* (qstp(i,k)* (1._dp + 2._dp/eps1*qstp(i,k))*  &
                   (1._dp + qstp(i,k)/eps1)*eps1**2.*rl*rl/           &
                   (rgas**2.*tp(i,k)**4.)-qstp(i,k)*                  &
                   (1._dp + qstp(i,k)/eps1)*2.*eps1*rl/               &
                   (rgas*tp(i,k)**3.))

            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
!          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/         &
                      ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
! jjh          tpv(i,k) =tp(i,k)*(1._dp + 1.608*qstp(i,k))/
! jt            (1._dp + q(i,mx(i)))
            tpv(i,k) = tp(i,k)* (1._dp + 1.608*qstp(i,k))/        &
                        (1._dp + qstp(i,k))
!                       (1._dp + q(i,mx(i)) + qprime(i,mx(i))) 
            buoy(i,k) = tpv(i,k) - tv(i,k)
          end if
        end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         do k=plev,1,-1
            if (k.le.mx(i)) then
               write (6,*) 'buoyan_pjr: k, tpv, tv, buoy ',  &
                   k, tpv(i,k), tv(i,k), buoy(i,k)
            endif
         end do
      endif
#endif

      maxlelten = 5
      do k = msg + 2,plev
         numbllcl = 0
         do i = 1,plon
            if (k.lt.lcl(i) ) then
               if (buoy(i,k+1).gt.0. .and. buoy(i,k).le.0.) then
                  knt(i) = min(maxlelten,knt(i) + 1)
                  lelten(i,knt(i)) = k
               end if
            else
               numbllcl = numbllcl + 1
            endif
         end do
         if (numbllcl.eq.plon) then
            goto 200
         endif
      end do
      write (6,*) 'buoy_pjr: numbllcl problem ', numbllcl
!      stop
  200 continue

!    Calculate convective available potential energy (cape).

      do n = 1,5
        do k = msg + 1,plev
          do i = 1,plon
            if (k.le.mx(i) .and.k.gt.lelten(i,n)) then
              capeten(i,n) = capeten(i,n) +                          &
                            rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
          end do
        end do
      end do

!     find maximum cape from all possible tentative capes from
!     one sounding,
!     and use it as the final cape, april 26, 1995

      do n = 1,maxlelten
        do i = 1,plon
          if (capeten(i,n).gt.cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
          end if
        end do
      end do

      do i = 1,plon
         capblfc(i) = 0.0
         lfc(i) = plev+2
      end do

!     find the first layer above the lcl with positive buoyancy 
!     as the lfc (level of free convection)
      do k = plev,msg + 1,-1
         do i = 1,plon
            if (k.le.lcl(i)           &
               .and.buoy(i,k).gt.0..and.lfc(i).eq.plev+2) then
               lfc(i) = k
            endif
         end do
      end do

!     now calulate cape between cloud base and the lfc
      do k = plev,msg+1,-1
         do i = 1,plon
            if (k.le.mx(i).and.k.gt.lfc(i)) then
               bloc = rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
!               if (bloc .lt. -100. .or. capblfc(i).lt.-1000.) then
!                  write (6,*) ' strange: i, k, lfc, lcl, bloc, capblfc',   &
!                      i, k, lfc(i), lcl(i), bloc, capblfc(i)
!                  stop
!               endif
               capblfc(i) = capblfc(i) + bloc
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
!         do i = 1,plon
            write (6,*) ' buoyan_pjr: capblfc, mx, lcl, lfc, cape ',       &
                i, capblfc(i), mx(i), lcl(i), lfc(i), cape(i)
!         end do
      endif
#endif

!     put lower bound on cape for diagnostic purposes.

      do i = 1,plon
        cape(i) = max(cape(i), 0._dp)
      end do
#ifdef RESETTOP

!     limit the cloud top, recalculate cape

      do i = 1,plon
         btot(i) = 0.0
         lel2(i) = plev+1
      end do
      do k = plev,1,-1
         do i = 1,plon
            bloc = rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            btot(i) = btot(i) + bloc
            if (k.le.lfc(i).and.btot(i).lt.off         &
               .and.lel2(i).eq.plev+1) then
!              write(*,*) 'lel2 found at',k+1,btot(i)
               lel2(i) = k+1
            endif
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' lel reset from ', lel(i),' to ', lel2(i)
      endif
#endif
      do i = 1,plon
         if (lel(i).ne.lel2(i)) then
            lel(i) = lel2(i)
         endif
         cape2(i) = 0.0
      end do
      do k = msg + 1,plev
         do i = 1,plon
            if (k.le.mx(i).and.k.gt.lel(i)) then
               cape2(i) = cape2(i) +                            &
                          rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' new, old cape ', cape2(i), cape(i),  &
             lel(i), mx(i)
         do k = lel(i), mx(i)
            write (6,*)  'cape by level ', k,               &
                  rgas*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
         end do
      endif
#endif
      do i = 1,plon
         cape(i) = max(0.0,cape2(i))
      end do

#endif

      return
    end subroutine buoyan_pjr
!===============================================================================
 
    subroutine cldprp_pjr(plond, plev   ,plevp   ,                  &
                       q       ,t       ,u       ,v       ,p       ,&
                       z       ,s       ,mu      ,eu      ,du      ,&
                       md      ,ed      ,sd      ,qd      ,ud      ,&
                       vd      ,mc      ,qu      ,su      ,zf      ,&
                       qst     ,hmn     ,hsat    ,shat             ,&
                       ql      ,totpcp  ,totevp  ,cmeg    ,jb      ,&
                       lel     ,jt      ,jlcl    ,mx      ,j0      ,&
                       jd      ,il2g    ,                           & 
                       pflx    ,evp     ,cu      ,mu2     ,eu2     ,& 
                       du2     ,md2     ,ed2     ,cmfdqr  ,         & 
                       tprime  ,qprime  ,ireset  ,                  &
                       lwc     ,iwc     ,rform   ,sform )
! -----------------------------------------------------------------------------
! This is contributed code not fully standardized by the CCM core group.

! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release

! the documentation has been enhanced to the degree that we are able

! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996

!**** PLEASE NOTE ****

! we are aware of a specific problem in this code 
! (identified by the string ---> PROBLEM ONE)
! during the calculation of the updraft cloud properties,
! rather than adding a perturbation to the updraft temperature of 
! half a degree, (there was an inadvertant addition of cp*0.5) degrees
! or about 500 degrees. (This problem was in the code prior to its 
! contribution to the NCAR effort)

! Fortunately, the erroneous values
! are overwritten later in the code. The problem is quite subtle.
! The erroneous values would persist between cloud base and the lifting 
! condensation level. The addition of the very high perturbation to the updraft
! temperature causes the saturation mixing ratio to be set to zero, 
! and later the lcl to be set to one level above cloud base.
! There are therefore no levels between cloud base and the lcl. Therefore
! all erroneous values are overwritten.

! The only manifestation we are aware of with respect to this problem
! is that the lifting condensation level is constrained to be one level above
! cloud base.

! We discovered the problem after too much had been invested in
! very long integrations (in terms of computer time)
! to allow for a modification and model retuning. It is our expectation that
! this problem will be fixed in the next release of the model.

! *********** 
! ----------------------------------------------------------------------

! $Id: conv_ccm_pjr.F,v 1.1.6.1 1999/04/03 18:54:29 eaton Exp $
! $Author: eaton $

! ----------------------------------------------------------------------
! nov 20/92 - guang jun zhang,m.lazare. now has deeper (more
!             realistic) downdrafts.
! jul 14/92 - guang jun zhang,m.lazare. add shallow mixing
!             formulation.
! nov 21/91 - m.lazare. like previous cldprop except minimum "f"
!                       now 0.0004 instead of 0.001 (more
!                       realistic with more deep).
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.

! -----------------------------------------------------------------------------
      IMPLICIT NONE
! Input arguments
      INTEGER, INTENT(IN) :: plond,      &! number of max columns
                             plev, plevp  ! number of levels, levels+1


      REAL(dp) :: q(plond,plev)        ! spec. humidity of env
      REAL(dp) :: t(plond,plev)        ! temp of env
      REAL(dp) :: p(plond,plev)        ! pressure of env
      REAL(dp) :: z(plond,plev)        ! height of env
      REAL(dp) :: s(plond,plev)        ! normalized dry static energy of env
      REAL(dp) :: zf(plond,plevp)      ! height of interfaces
      REAL(dp) :: u(plond,plev)        ! zonal velocity of env
      REAL(dp) :: v(plond,plev)        ! merid. velocity of env
      REAL(dp) :: tprime(plond,plev)   ! level dependent tpert
      REAL(dp) :: qprime(plond,plev)   ! level dependent qpert

      INTEGER  :: jb(plond)         ! updraft base level
      INTEGER  :: lel(plond)        ! updraft launch level
      INTEGER  :: jt(plond)         ! updraft plume top
      INTEGER  :: jlcl(plond)       ! updraft lifting cond level
      INTEGER  :: mx(plond)         ! updraft base level (same is jb)
      INTEGER  :: j0(plond)         ! level where updraft begins detraining
      INTEGER  :: jd(plond)         ! level of downdraft

! output

      REAL(dp) :: cmfdqr(plond,plev)   ! rate of production of precip at that layer
      REAL(dp) :: du(plond,plev)       ! detrainement rate of updraft
      REAL(dp) :: ed(plond,plev)       ! entrainment rate of downdraft
      REAL(dp) :: eu(plond,plev)       ! entrainment rate of updraft
      REAL(dp) :: hmn(plond,plev)      ! moist stat energy of env
      REAL(dp) :: hsat(plond,plev)     ! sat moist stat energy of env
      REAL(dp) :: mc(plond,plev)       ! net mass flux
      REAL(dp) :: md(plond,plev)       ! downdraft mass flux
      REAL(dp) :: mu(plond,plev)       ! updraft mass flux
      REAL(dp) :: pflx(plond,plevp)    ! precipitation flux thru layer
      REAL(dp) :: qd(plond,plev)       ! spec humidity of downdraft
      REAL(dp) :: ql(plond,plev)       ! liq water of updraft
      REAL(dp) :: qst(plond,plev)      ! saturation spec humidity of env.
      REAL(dp) :: qu(plond,plev)       ! spec hum of updraft
      REAL(dp) :: sd(plond,plev)       ! normalized dry stat energy of downdraft
      REAL(dp) :: shat(plond,plev)     ! interface values of dry stat energy
      REAL(dp) :: su(plond,plev)       ! normalized dry stat energy of updraft
      REAL(dp) :: ud(plond,plev)       ! downdraft u
      REAL(dp) :: vd(plond,plev)       ! downdraft v

!     these version of the mass fluxes conserve mass (used in tracer transport)

      REAL(dp) :: mu2(plond,plev)      ! updraft mass flux
      REAL(dp) :: eu2(plond,plev)      ! updraft entrainment
      REAL(dp) :: du2(plond,plev)      ! updraft detrainment
      REAL(dp) :: md2(plond,plev)      ! downdraft mass flux
      REAL(dp) :: ed2(plond,plev)      ! downdraft entrainment

      INTEGER  :: il2g              !CORE GROUP REMOVE
! mz_ht_20070918+
      REAL(dp) :: lwc(plond,plev)
      REAL(dp) :: iwc(plond,plev)
      REAL(dp) :: rform(plond,plev)
      REAL(dp) :: sform(plond,plev)
      ! mz_ht_20070918-
! Local workspace

      REAL(dp) :: gamma(plond,plev)  
      REAL(dp) :: dz(plond,plev)  
      REAL(dp) :: iprm(plond,plev)  
      REAL(dp) :: hu(plond,plev)  
      REAL(dp) :: hd(plond,plev)  
      REAL(dp) :: eps(plond,plev)  
      REAL(dp) :: f(plond,plev)  
      REAL(dp) :: k1(plond,plev)  
      REAL(dp) :: i2(plond,plev)  
      REAL(dp) :: ihat(plond,plev)  
      REAL(dp) :: i3(plond,plev)  
      REAL(dp) :: idag(plond,plev)  
      REAL(dp) :: i4(plond,plev)  
      REAL(dp) :: qsthat(plond,plev)  
      REAL(dp) :: hsthat(plond,plev)  
      REAL(dp) :: gamhat(plond,plev)  
      REAL(dp) :: cu(plond,plev)  
      REAL(dp) :: evp(plond,plev)  
      REAL(dp) :: cmeg(plond,plev)  
      REAL(dp) :: qds(plond,plev) 
      REAL(dp) :: hmin(plond)  
      REAL(dp) :: expdif(plond)  
      REAL(dp) :: expnum(plond)  
      REAL(dp) :: ftemp(plond)  
      REAL(dp) :: eps0(plond)  
      REAL(dp) :: rmue(plond)  
      REAL(dp) :: zuef(plond)  
      REAL(dp) :: zdef(plond)  
      REAL(dp) :: epsm(plond)  
      REAL(dp) :: ratmjb(plond)  
      REAL(dp) :: est(plond)  
      REAL(dp) :: totpcp(plond)  
      REAL(dp) :: totevp(plond)  
      REAL(dp) :: alfa(plond) 
      REAL(dp) :: beta
!      REAL(dp) :: c0
      REAL(dp) :: ql1
      REAL(dp) :: weight
      REAL(dp) :: tu
      REAL(dp) :: estu
      REAL(dp) :: qstu

      REAL(dp) :: small
      REAL(dp) :: mdt  

      INTEGER  :: khighest
      INTEGER  :: klowest  
      INTEGER  :: kount 
      INTEGER  :: i,k

      LOGICAL  :: doit(plond)
      LOGICAL  :: done(plond)
      INTEGER  :: ireset(plond)

      REAL(dp) :: tmp

! -----------------------------------------------------------------------------

      do i = 1,il2g
         ftemp(i) = 0.
         expnum(i) = 0.
         expdif(i) = 0.
      end do

! jr Change from msg+1 to 1 to prevent blowup

      do k = 1,plev
        do i = 1,il2g
          dz(i,k) = zf(i,k) - zf(i,k+1)
        end do
      end do


! initialize many output and work variables to zero

      do k = 1,plev
        do i = 1,il2g
          k1(i,k) = 0.
          i2(i,k) = 0.
          i3(i,k) = 0.
          i4(i,k) = 0.
          mu(i,k) = 0.
          f(i,k) = 0.
          eps(i,k) = 0.
          eu(i,k) = 0.
          du(i,k) = 0.
          ql(i,k) = 0.
          cu(i,k) = 0.
          evp(i,k) = 0.
          cmeg(i,k) = 0.
          qds(i,k) = q(i,k)
          md(i,k) = 0.
          ed(i,k) = 0.
          sd(i,k) = s(i,k)
          qd(i,k) = q(i,k)
          ud(i,k) = u(i,k)
          vd(i,k) = v(i,k)
          mc(i,k) = 0.
          qu(i,k)   = q(i,k) + qprime(i,k)
          su(i,k)   = s(i,k) + tprime(i,k)
!        est(i)=exp(a-b/t(i,k))
          est(i)     = c1*exp((c2*(t(i,k) - tfreez))     &
                        /((t(i,k) - tfreez) + c3))
! +bee
          if ( p(i,k)-est(i) .gt. 0. ) then
             qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
          else
             qst(i,k) = 1.0_dp
          end if
! -bee
          gamma(i,k) = qst(i,k)*(1._dp + qst(i,k)/eps1)*eps1*rl/    &
                      (rgas*t(i,k)**2)*rl/cpg
          hmn(i,k) = cpg*t(i,k) + grav*z(i,k) + rl*q(i,k)
          hsat(i,k) = cpg*t(i,k) + grav*z(i,k) + rl*qst(i,k)
! *old          hu(i,k) = hmn(i,k)
          hu(i,k)    = cpg*su(i,k) + rl*qu(i,k)
          hd(i,k) = hmn(i,k)
          mu2(i,k) = 0.
          eu2(i,k) = 0.
          du2(i,k) = 0.
          md2(i,k) = 0.
          ed2(i,k) = 0.
          pflx(i,k) = 0.
          cmfdqr(i,k) = 0.
        end do
      end do


! interpolate the layer values of qst, hsat and gamma to
! layer interfaces

      do i = 1,il2g
        hsthat(i,msg+1) = hsat(i,msg+1)
        qsthat(i,msg+1) = qst(i,msg+1)
        gamhat(i,msg+1) = gamma(i,msg+1)
        totpcp(i) = 0.
        totevp(i) = 0.
        pflx(i,plevp) = 0.
        dz(i,plev) = zf(i,plev) - zf(i,plev+1)
      end do
      do k = msg + 2,plev
        do i = 1,il2g
          if (abs(qst(i,k-1)-qst(i,k)).gt.1.E-6_dp) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*          &
                         qst(i,k)/ (qst(i,k-1)-qst(i,k))
          else
            qsthat(i,k) = qst(i,k)
          end if

          hsthat(i,k) = cpg*shat(i,k) + rl*qsthat(i,k)

           if (abs(gamma(i,k-1)-gamma(i,k)).gt.1.E-6_dp) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*         &
                          gamma(i,k-1)*gamma(i,k)/              &
                         (gamma(i,k-1)-gamma(i,k))
          else
            gamhat(i,k) = gamma(i,k)
          end if
        end do
      end do

!     Initialize cloud top to highest plume top.

      do i = 1,il2g
        jt(i) = max(lel(i),4)
        jd(i) = plev
        jlcl(i) = lel(i)
        hmin(i) = 1.E6_dp
        ireset(i) = 0
      end do
#ifdef DEBCONV      
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' cldprp: jt is initialized to ', jt(i)
         write (6,*) ' lel is ', lel(i)
      endif
#endif

!     Find the level of minimum hsat, where detrainment starts

      do k = msg + 1,plev
        do i = 1,il2g
            if (hsat(i,k).le.hmin(i) .and.          &
               (k.ge.jt(i).and.k.le.jb(i))) then
               hmin(i) = hsat(i,k)
               j0(i) = k
            end if
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 set 1 ', j0(i)
      endif
#endif

      do i = 1,il2g
        j0(i) = min(j0(i),jb(i)-2)
        j0(i) = max(j0(i),jt(i)+2)

! Fix from Guang Zhang to address out of bounds array reference

        j0(i) = min(j0(i),plev)
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 set 2 ', j0(i)
      endif
#endif

! no longer need this since perturbations are applied above

! Initialize certain arrays inside cloud

!      do k = msg + 1,plev
!        do i = 1,il2g
!          if (k.ge.jt(i) .and. k.le.jb(i)) then
!            hu(i,k) = hmn(i,mx(i)) + cpg*0.5
!            su(i,k) = s(i,mx(i)) + 0.5
!          end if
!        end do
!      end do

! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************

      do k = plev - 1,msg + 1,-1
        do i = 1,il2g
          if (k.lt.jb(i) .and. k.ge.jt(i)) then
           tmp = min(hmn(i,k),hsat(i,k))
            k1(i,k) = k1(i,k+1)                      &
                + (hmn(i,mx(i))-tmp)*dz(i,k)
!               + (hmn(i,mx(i))-hmn(i,k))*dz(i,k)
            ihat(i,k) = 0.5* (k1(i,k+1)+k1(i,k))
            i2(i,k) = i2(i,k+1) + ihat(i,k)*dz(i,k)
            idag(i,k) = 0.5* (i2(i,k+1)+i2(i,k))
            i3(i,k) = i3(i,k+1) + idag(i,k)*dz(i,k)
            iprm(i,k) = 0.5* (i3(i,k+1)+i3(i,k))
            i4(i,k) = i4(i,k+1) + iprm(i,k)*dz(i,k)
          end if
        end do
      end do

! re-initialize hmin array for ensuing calculation.

      do i = 1,il2g
        hmin(i) = 1.E6_dp
      end do

      do k = msg + 1,plev
        do i = 1,il2g
!          tmp = hmn(i,k)
!          added to deal with supersaturated layers
           tmp = min(hmn(i,k),hsat(i,k))
          if (k.ge.j0(i).and.k.le.jb(i) .and. tmp.le.hmin(i)) then
            hmin(i) = tmp
            expdif(i) = hmn(i,mx(i)) - hmin(i) 
          end if
        end do
      end do

! *********************************************************
!     compute approximate eps(z) using above taylor series
! *********************************************************

      do k = msg + 2,plev
        do i = 1,il2g
          expnum(i) = 0.
          ftemp(i) = 0.
          if (k.lt.jt(i) .or. k.ge.jb(i)) then
            k1(i,k) = 0.
            expnum(i) = 0.
          else
            expnum(i) = hmn(i,mx(i)) - (hsat(i,k-1)*(zf(i,k)-z(i,k)) +         &
                        hsat(i,k)* (z(i,k-1)-zf(i,k)))/(z(i,k-1)-z(i,k))
          end if
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
             write (6,*) ' expnum, k1 ', k, expnum(i), k1(i,k),       &
                  expnum(i)*dz(i,k), expdif(i),                       &
                  jt(i), jb(i)
          endif
#endif
          if ((expdif(i).gt.100._dp.and.expnum(i).gt.0.) .and.           &
               k1(i,k).gt.expnum(i)*dz(i,k)) then
            ftemp(i) = expnum(i)/k1(i,k)
            f(i,k) = ftemp(i) + i2(i,k)/k1(i,k)*ftemp(i)**2 +         & 
                      (2.*i2(i,k)**2-k1(i,k)*i3(i,k))/k1(i,k)**2*     &
                      ftemp(i)**3 + (-5.*k1(i,k)*i2(i,k)*i3(i,k)+     &
                      5.*i2(i,k)**3+k1(i,k)**2*i4(i,k))/              &
                      k1(i,k)**3*ftemp(i)**4
            f(i,k) = max(f(i,k),0._dp)
            f(i,k) = min(f(i,k),0.0002_dp)
          end if
        end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' f ', f(i,k), ftemp(i)
      endif
#endif
      end do
      do i = 1,il2g
        if (j0(i).lt.jb(i)) then
          if (f(i,j0(i)).lt.1.E-6_dp .and.              &
              f(i,j0(i)+1).gt.f(i,j0(i))) then
             j0(i) = j0(i) + 1
          endif
        end if
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' j0 reset ', j0(i)
      endif
#endif

      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.le.j0(i)) then
            f(i,k) = max(f(i,k),f(i,k-1))
          end if
        end do
      end do
      do i = 1,il2g
        eps0(i) = f(i,j0(i))
        eps(i,jb(i)) = eps0(i)
      end do
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' set eps0 and epsjb ', j0(i), eps0(i),eps(i,jb(i))
      endif
#endif

      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.ge.j0(i) .and. k.le.jb(i)) then
               eps(i,k) = f(i,j0(i))
            end if
         end do
      end do

      do k = plev,msg + 1,-1
         do i = 1,il2g
            if (k.lt.j0(i) .and. k.ge.jt(i)) eps(i,k) = f(i,k)
         end do
      end do

!     specify the updraft mass flux mu, entrainment eu, detrainment du
!     and moist stati! energy hu.
!     here and below mu, eu,du, md and ed are all normalized by mb

      do i = 1,il2g
        if (eps0(i).gt.0.) then
!          mu(i,jb(i)) = 1.
!          eu(i,jb(i)) = eps0(i)/2.
! *pjr NOTE TO CORE GROUP mu and mu2 (and eu and eu2 should have the same vals now)
          mu2(i,jb(i)) = 1._dp
          eu2(i,jb(i)) = mu2(i,jb(i))/dz(i,jb(i))
          mu(i,jb(i)) = mu2(i,jb(i))
          eu(i,jb(i)) = eu2(i,jb(i))
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
             write (6,*) ' jb, mujb eps0 ', jb(i), mu(i,jb(i)), eps0(i)
          endif
#endif
        end if
      end do
      do k = plev,msg + 1,-1
        do i = 1,il2g
           if (eps0(i).gt.0. .and. k.ge.jt(i).and.k.lt.jb(i)) then
            zuef(i) = zf(i,k) - zf(i,jb(i))
            rmue(i) = (1./eps0(i))* (exp(eps(i,k+1)*zuef(i))-1._dp)/zuef(i)
            mu(i,k) = (1./eps0(i))* (exp(eps(i,k)*zuef(i))-1._dp)/zuef(i)
#ifdef DEBCONV
          if (lat.eq.latlook.and.i.eq.ilookg) then
               write (6,*) ' muk, epsk  ', k, mu(i,k), eps(i,k)
          endif
#endif
            eu(i,k) = (rmue(i)-mu(i,k+1))/dz(i,k)
            du(i,k) = (rmue(i)-mu(i,k))/dz(i,k)
            mu2(i,k) = mu(i,k)
            eu2(i,k) = eu(i,k)
            du2(i,k) = du(i,k)
          end if
        end do
      end do

      khighest = plevp
      klowest = 1
      do i=1,il2g
        khighest = min(khighest,lel(i))
        klowest = max(klowest,jb(i))
      end do
      do k = klowest-1,khighest,-1
!dir$ ivdep
        do i = 1,il2g
          if (k.le.jb(i)-1 .and. k.ge.lel(i) .and. eps0(i).gt.0.) then
#ifdef DEBCONV
             if (lat.eq.latlook.and.i.eq.ilookg) then
                write (6,*) ' small muk, hu, hmn ', k, mu(i,k),    &
                    hu(i,k), hmn(i,k), hsat(i,k)
             endif
#endif
            if (mu(i,k).lt.0.01_dp) then
              hu(i,k) = hu(i,jb(i))
              mu(i,k) = 0.
              mu2(i,k) = mu(i,k)
              eu2(i,k) = 0.
              du2(i,k) = mu2(i,k+1)/dz(i,k)
              eu(i,k) = eu2(i,k)
              du(i,k) = du2(i,k)
            else
              hu(i,k) = mu(i,k+1)/mu(i,k)*hu(i,k+1) +           &
                        dz(i,k)/mu(i,k)* (eu(i,k)*hmn(i,k)-     &
                        du(i,k)*hsat(i,k))
            end if
#ifdef DEBCONV
             if (lat.eq.latlook.and.i.eq.ilookg) then
                write (6,*) ' revised muk, hu, hmn ', k, mu(i,k),    &
                    hu(i,k), hmn(i,k)
             endif
#endif
          end if
        end do
      end do
! #ifdef OFFIT

! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset

! *pjr there are diffs here, and I hope to god they dont matter
      do i=1,il2g
         doit(i) = .true.
         ireset(i) = 0
      end do
      do k=klowest-2,khighest-1,-1
         do i=1,il2g
            if (doit(i) .and. k.le.jb(i)-2 .and. k.ge.lel(i)-1) then
#ifdef DEBCONV
               if (lat.eq.latlook.and.i.eq.ilookg) then
                  write (6,*) ' hub, hstb, hua, hsta, mu doit ', k,  &
                      hu(i,k), hsthat(i,k), hu(i,k+1),               &
                      hsthat(i,k+1), mu(i,k), doit(i)
               endif
#endif
               if (hu(i,k  ).le.hsthat(i,k) .and.             &
                   hu(i,k+1).gt.hsthat(i,k+1) .and.           &
!                  mu(i,k).ge.0.02) then
! +mgl set this to 0.01 to be consistent with computation of hu in the cloud column done above
                   mu(i,k).ge.0.01_dp) then
                  if (hu(i,k)-hsthat(i,k).lt.-2000._dp) then
                     jt(i) = k + 1
                     doit(i) = .false.
                     ireset(i) = 1
                  else
                     jt(i) = k
                     ireset(i) = 2
                     if (eps0(i).le.0.) then
                        doit(i) = .false.
                        ireset(i) = 3
                     endif
                  end if
               else if (hu(i,k).gt.hu(i,jb(i)) .or. mu(i,k).lt.0.01_dp        &
! +mgl this condition causes a lot of overturning of only the surface layer!
!                  .or. (hu(i,k)-hsthat(i,k).lt.-2000._dp)) then
                      ) then
                  jt(i) = k + 1
                  doit(i) = .false.
                  ireset(i) = 4
               end if
            end if
#ifdef DEBCONV
            if (lat.eq.latlook.and.i.eq.ilookg) then
               write (6,*) ' jt reset section ', k,        &
                             jt(i), ireset(i), doit(i)
            endif
#endif
        end do
      end do

#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' jt final section ', jt(i), ireset(i)
      endif
#endif
!#endif


      do k = plev,msg + 1,-1
!dir$ ivdep
        do i = 1,il2g
          if (k.ge.lel(i) .and. k.le.jt(i) .and.    &
              eps0(i).gt.0.) then
            mu(i,k) = 0.
            eu(i,k) = 0.
            du(i,k) = 0.
            mu2(i,k) = 0.
            eu2(i,k) = 0.
            du2(i,k) = 0.
!     CPGJR/JP               hu(i,k) = hu(i,jb(i)) !commen out?
          end if
          if (k.eq.jt(i) .and. eps0(i).gt.0.) then
            du(i,k) = mu(i,k+1)/dz(i,k)
            du2(i,k) = mu2(i,k+1)/dz(i,k)
            eu2(i,k) = 0.
            mu2(i,k) = 0.
            eu(i,k) = 0.
            mu(i,k) = 0.
          end if
        end do
      end do

! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.

      do i = 1,il2g

! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1

        alfa(i) = dd
        jt(i) = min(jt(i),jb(i)-1)
        jd(i) = max(j0(i),jt(i)+1)
        jd(i) = min(jd(i),jb(i))
        hd(i,jd(i)) = hmn(i,jd(i)-1)
        ud(i,jd(i)) = u(i,jd(i)-1)
        vd(i,jd(i)) = v(i,jd(i)-1)
        if (jd(i).lt.jb(i) .and. eps0(i).gt.0.) then
          epsm(i) = eps0(i)
!          alfa(i)=2.*epsm(i)*( zf(i,jd(i))-zf(i,jb(i)) )/
!     1         (  exp(2.*epsm(i)*( zf(i,jd(i))-
!               zf(i,jb(i)) ))-1.  )
          md(i,jd(i)) = -alfa(i)*epsm(i)/eps0(i)
          md2(i,jd(i)) = md(i,jd(i))
        end if
      end do
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.gt.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0.) then
            zdef(i) = zf(i,jd(i)) - zf(i,k)
            md(i,k) = -alfa(i)/ (2.*eps0(i))*             &
                      (exp(2.*epsm(i)*zdef(i))-1._dp)/zdef(i)
            md2(i,k) = md(i,k)
          end if
        end do
      end do
      do k = msg + 1,plev
! dir$ ivdep
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.     &
              jd(i).lt.jb(i)) then
            ratmjb(i) = min(abs(mu2(i,jb(i))/md2(i,jb(i))),1._dp)
            md2(i,k) = md2(i,k)*ratmjb(i)
!            ratmjb(i) = min(abs(mu(i,jb(i))/md(i,jb(i))),1.)
!            md(i,k) = md(i,k)*ratmjb(i)
               md(i,k) = md2(i,k)
          end if
        end do
      end do

      small = 1.e-20_dp
      do k = msg + 1,plev
        do i = 1,il2g
          if ((k.ge.jt(i).and.k.le.jb(i)) .and. eps0(i).gt.0._dp) then
            ed2(i,k-1) = (md2(i,k-1)-md2(i,k))/dz(i,k-1)
               ed(i,k-1) = ed2(i,k-1)
            mdt = min(md2(i,k),-small)
            hd(i,k) = (md(i,k-1)*hd(i,k-1) -                        &
                       dz(i,k-1)*ed(i,k-1)*hmn(i,k-1))/mdt
          end if
        end do
      end do

! calculate updraft and downdraft properties.

      do k = msg + 2,plev
        do i = 1,il2g
          if ((k.ge.jd(i).and.k.le.jb(i)) .and. eps0(i).gt.0. .and.      &
              jd(i).lt.jb(i)) then
!         sd(i,k) = shat(i,k)
!    1             +              (hd(i,k)-hsthat(i,k))/
!    2               (cpg    *(1.+gamhat(i,k)))
            qds(i,k) = qsthat(i,k) + gamhat(i,k)*(hd(i,k)-hsthat(i,k))/   &
                       (rl*(1._dp + gamhat(i,k)))
          end if
        end do
      end do

      do i = 1,il2g
         done(i) = .false.
      end do
      kount = 0
      do k = plev,msg + 2,-1
        do i = 1,il2g
          if ((.not.done(i) .and. k.gt.jt(i) .and. k.lt.jb(i)) .and.      &
               eps0(i).gt.0.) then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) +                       &
                      dz(i,k)/mu(i,k)* (eu(i,k)-du(i,k))*s(i,k)
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) +                       &
                      dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)-                   &
                      du(i,k)*qst(i,k))
            tu = su(i,k) - grav/cpg*zf(i,k)
            estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
            qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2.-estu)
            if (qu(i,k).ge.qstu) then
              jlcl(i) = k
              kount = kount + 1
              done(i) = .true.
            end if
          end if
        end do
        if (kount.ge.il2g) goto 690
      end do
 690  continue
#ifdef DEBCONV
      if (lat.eq.latlook.and.ilookg.ne.0) then
         i = ilookg
         write (6,*) ' cldprp: jlcl found to be ', jlcl(i)
         write (6,*) ' jt and jb are currently ', jt(i), jb(i)
      endif
#endif

#ifdef LCLJT
! *PJR*XXX    A Change in SCM JP/PJR 

      do i = il1g, il2g
         if (jlcl(i).le.jt(i)) then
!            write(*,*) ' ***jlcl =',jlcl(1),'jt =', jt(1),nstep
!            write(*,*) 'Convection should be supressed by phil'
            jt(i)   = jb(i)     !jpnew
            jlcl(i) = jb(i)     !jpnew
            ireset(i) = 5
         endif
      end do
#endif

      do k = msg + 2,plev
        do i = 1,il2g
!          if (k.eq.jb(i) .and. eps0(i).gt.0.) then
!            qu(i,k) = q(i,mx(i))
!            su(i,k) = (hu(i,k)-rl*qu(i,k))/cpg
!          end if
          if ((k.gt.jt(i).and.k.le.jlcl(i)) .and. eps0(i).gt.0.) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/              &
                     (cpg* (1._dp + gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*                      &
                     (hu(i,k)-hsthat(i,k))/                           & 
                     (rl* (1._dp + gamhat(i,k)))
          end if
        end do
      end do

      do k = plev,msg + 2,-1
        do i = 1,il2g
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0.) then
            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/          &    
                        dz(i,k)- (eu(i,k)-du(i,k))*s(i,k))/            &
                       (rl/cpg)
            if (k.eq.jt(i)) cu(i,k) = 0.
!               cu(i,k) = max(0.,cu(i,k))
!               cu2     = max(0.,
!               cu2     = max(-1.e99,                                    &
!                         +(eu(i,k)*q(i,k) - du(i,k)*qst(i,k))           &
!                         -(mu(i,k)*qu(i,k)-mu(i,k+1)*qu(i,k+1))/dz(i,k) &
!                         )
!                  
!               if (abs(cu(i,k)-cu2)/(abs(cu(i,k))+abs(cu2)+1.e-50)      &
!                   .gt.0.0000001) then
!                  write (6,*) ' inconsistent condensation rates ',      &
!                       i, k, lat,                                       &
!                       cu(i,k), cu2, jt(i), jb(i), jlcl(i), lel(i)      &
!                       ,mu(i,k)
!               endif
               cu(i,k) = max(0._dp,cu(i,k))
          end if
        end do
      end do

      beta = 0.
!      c0 = 2.E-3_dp
      do k = plev,msg + 2,-1
        do i = 1,il2g
          cmfdqr(i,k) = 0.
! this modification is for test3 run, modified on 6/20/1995
!        if(t(i,jt(i) ).gt.tfreez)    c0=0.
!        if(t(i,jt(i) ).le.tfreez   )    c0=2.e-3
          if (k.ge.jt(i) .and. k.lt.jb(i) .and. eps0(i).gt.0. .and.        &
              mu(i,k).ge.0.0) then
            if (mu(i,k).gt.0.) then
              ql1 = 1./mu(i,k)* (mu(i,k+1)*ql(i,k+1)-                      &
                   dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))
              ql(i,k) = ql1/ (1._dp + dz(i,k)*c0)
            else
              ql(i,k) = 0.
            end if
            totpcp(i) = totpcp(i) + dz(i,k)*(cu(i,k)-du(i,k)*              &
                        (beta*ql(i,k) + (1._dp - beta)*ql(i,k+1)))
            cmfdqr(i,k) = c0*mu(i,k)*ql(i,k)
            ! mz_ht_20070918+
            IF (T(i,k) > Tmelt) THEN
              lwc(i,k)   = ql(i,k)
              rform(i,k) = c0 * ql(i,k)
            ELSE
              iwc(i,k)   = ql(i,k)
              sform(i,k) = c0 * ql(i,k)
            ENDIF
            ! mz_ht_20070918-
          end if
        end do
      end do

      do i = 1,il2g
        qd(i,jd(i)) = qds(i,jd(i))
        sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cpg
      end do

      do k = msg + 2,plev
        do i = 1,il2g
          if (k.ge.jd(i).and.k.lt.jb(i) .and. eps0(i).gt.0.) then
            qd(i,k+1) = qds(i,k+1)
            evp(i,k) = -ed(i,k)*q(i,k) +                               &
                       (md(i,k)*qd(i,k)-md(i,k+1)*qd(i,k+1))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._dp)
            mdt = min(md(i,k+1),-small)
            sd(i,k+1) = ((rl/cpg*evp(i,k)-ed(i,k)*s(i,k))*dz(i,k) +    &
                          md(i,k)*sd(i,k))/mdt
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      end do
      do i = 1,il2g
! *guang         totevp(i) = totevp(i) + md(i,jd(i))*q(i,jd(i)-1) -
        totevp(i) = totevp(i) + md(i,jd(i))*qd(i,jd(i)) -              &
                    md(i,jb(i))*qd(i,jb(i))
      end do
      if (.true.) then
        do i = 1,il2g
          k = jb(i)
          if (eps0(i).gt.0.) then
            evp(i,k) = -ed(i,k)*q(i,k) + (md(i,k)*qd(i,k))/dz(i,k)
            evp(i,k) = max(evp(i,k),0._dp)
            totevp(i) = totevp(i) - dz(i,k)*ed(i,k)*q(i,k)
          end if
        end do
      endif

      do i = 1,il2g
        totpcp(i) = max(totpcp(i),0._dp)
        totevp(i) = max(totevp(i),0._dp)
      end do

      weight = 1.0_dp
      do k = msg + 2,plev
        do i = 1,il2g
          if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
            md2(i,k) = md2(i,k)*min(1._dp,weight*totpcp(i)/         &
                      (totevp(i)+weight*totpcp(i)))
            ed2(i,k) = ed2(i,k)*min(1._dp,weight*totpcp(i)/         & 
                      (totevp(i)+weight*totpcp(i)))
            evp(i,k) = evp(i,k)*min(1._dp,                          &
                       weight*totpcp(i)/ (totevp(i)+             &
                       weight*totpcp(i)))
          else
            md2(i,k) = 0.
            ed2(i,k) = 0.
            evp(i,k) = 0.
          end if
            md(i,k) = md2(i,k)
            ed(i,k) = ed2(i,k)
            cmeg(i,k) = cu(i,k) - evp(i,k)
!           cmeg is the cloud water condensed - rain water evaporated
!           cmfdqr  is the cloud water converted to rain - (rain evaporated)
            cmfdqr(i,k) = cmfdqr(i,k)-evp(i,k)
        end do
      end do
      do k = 2,plevp
        do i = 1,il2g
          pflx(i,k) = pflx(i,k-1) + cmfdqr(i,k-1)*dz(i,k-1)
        end do
      end do
      do i = 1,il2g
        if (totevp(i).gt.0. .and. totpcp(i).gt.0.) then
          totevp(i) = totevp(i)*min(1._dp,                           &   
                      weight*totpcp(i)/(totevp(i) + weight*totpcp(i)))
        else
          totevp(i) = 0.
        end if
      end do

      do k = msg + 1,plev
        do i = 1,il2g
          mc(i,k) = mu(i,k) + md(i,k)
        end do
      end do

      return
    end subroutine cldprp_pjr
!=========================================================================================

    subroutine closure_pjr(plond, plev   ,                           &
                        q       ,t       ,p       ,s       ,         &
                        tp      ,qu      ,su      ,mc      ,         &
                        du      ,mu      ,md      ,qd      ,sd      ,&
                        qhat    ,shat    ,wdp     ,qstp    ,         &
                        zf      ,ql      ,dsubcld ,mb      ,cape    ,&
                        tl      ,lcl     ,lel     ,jt      ,mx      ,&
                        il1g    ,il2g    ,                           &
                        capelmt )
! ----------------------------------------------------------------------

! This is contributed code not fully standardized by the CCM core group.

! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release

! the documentation has been enhanced to the degree that we are able

! Original version:  G. Zhang and collaborators
! Standardized:      Core group staff, 1994 and 195
! Reviewed:          P. Rasch, April 1996
! ----------------------------------------------------------------------

! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.

! ----------------------------Arguments---------------------------------
      IMPLICIT NONE
      INTRINSIC :: TINY
      INTEGER, INTENT(IN) :: plond, plev !number of columns, levels

      REAL(dp) :: q(plond,plev)        ! spec humidity
      REAL(dp) :: t(plond,plev)        ! temperature
      REAL(dp) :: p(plond,plev)        ! pressure (mb)
      REAL(dp) :: s(plond,plev)        ! normalized dry static energy 
      REAL(dp) :: tp(plond,plev)       ! parcel temp
      REAL(dp) :: qu(plond,plev)       ! updraft spec. humidity
      REAL(dp) :: su(plond,plev)       ! normalized dry stat energy of updraft
      REAL(dp) :: mc(plond,plev)       ! net convective mass flux 
      REAL(dp) :: du(plond,plev)       ! detrainment from updraft
      REAL(dp) :: mu(plond,plev)       ! mass flux of updraft
      REAL(dp) :: md(plond,plev)       ! mass flux of downdraft
      REAL(dp) :: qd(plond,plev)       ! spec. humidity of downdraft
      REAL(dp) :: sd(plond,plev)       ! dry static energy of downdraft
      REAL(dp) :: qhat(plond,plev)     ! environment spec humidity at interfaces
      REAL(dp) :: shat(plond,plev)     ! env. normalized dry static energy at intrfcs
      REAL(dp) :: wdp(plond,plev)       ! pressure thickness of layers
      REAL(dp) :: qstp(plond,plev)     ! spec humidity of parcel
      REAL(dp) :: zf(plond,plev+1)     ! height of interface levels
      REAL(dp) :: ql(plond,plev)       ! liquid water mixing ratio

      REAL(dp) :: mb(plond)            ! cloud base mass flux
      REAL(dp) :: cape(plond)          ! available pot. energy of column
      REAL(dp) :: tl(plond)
      REAL(dp) :: dsubcld(plond)       ! thickness of subcloud layer

      INTEGER  :: lcl(plond)        ! index of lcl
      INTEGER  :: lel(plond)        ! index of launch leve
      INTEGER  :: jt(plond)         ! top of updraft
      INTEGER  :: mx(plond)         ! base of updraft

! -------------------------Local variables------------------------------

      REAL(dp) :: dtpdt(plond,plev)
      REAL(dp) :: dqsdtp(plond,plev)
      REAL(dp) :: dtmdt(plond,plev)
      REAL(dp) :: dqmdt(plond,plev)
      REAL(dp) :: dboydt(plond,plev)
      REAL(dp) :: thetavp(plond,plev)
      REAL(dp) :: thetavm(plond,plev)

      REAL(dp) :: dtbdt(plond),dqbdt(plond),dtldt(plond)
      REAL(dp) :: beta
      REAL(dp) :: capelmt
      REAL(dp) :: dadt
      REAL(dp) :: debdt
      REAL(dp) :: dltaa
      REAL(dp) :: eb

      INTEGER  :: i
      INTEGER  :: il1g
      INTEGER  :: il2g
      INTEGER  :: k

! ----------------------------------------------------------------------
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.

      do i = il1g,il2g
         mb(i) = 0.
         eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
         dtbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                 &
                    (shat(i,mx(i))-su(i,mx(i)))+                   &
                     md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
         dqbdt(i) = (1./dsubcld(i))* (mu(i,mx(i))*                 &
                    (qhat(i,mx(i))-qu(i,mx(i)))+                   &
                     md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
         debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
         dtldt(i) = -2840.* (3.5/t(i,mx(i))*dtbdt(i)-debdt/eb)/    &
                    (3.5*log(t(i,mx(i)))-log(eb)-4.805_dp)**2
      end do

!   dtmdt and dqmdt are cumulus heating and drying.

      do k = msg + 1,plev
         do i = il1g,il2g
            dtmdt(i,k) = 0.
            dqmdt(i,k) = 0.
         end do
      end do

      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.eq.jt(i)) then
               dtmdt(i,k) = (1./wdp(i,k))*                              &
                             (mu(i,k+1)* (su(i,k+1)-shat(i,k+1)-        &
                             rl/cpg*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-   &
                             shat(i,k+1)))
               dqmdt(i,k) = (1./wdp(i,k))*(mu(i,k+1)* (qu(i,k+1)-       &
                             qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*          &
                            (qd(i,k+1)-qhat(i,k+1)))
            end if
         end do
      end do

      beta = 0.
      do k = msg + 1,plev - 1
         do i = il1g,il2g
            if (k.gt.jt(i) .and. k.lt.mx(i)) then
               dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+        &
                             mc(i,k+1)* (s(i,k)-shat(i,k+1)))/   &
                             wdp(i,k) - rl/cpg*du(i,k)*          &
                            (beta*ql(i,k)+ (1._dp-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))                &
!                      +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/wdp(i,k) &
!                      +du(i,k)*(qs(i,k)-q(i,k))                 &
!                      +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

               dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+   &
                             cpg/rl* (su(i,k+1)-s(i,k)))-         &
                             mu(i,k)* (qu(i,k)-qhat(i,k)+cpg/rl*  &
                             (su(i,k)-s(i,k)))+md(i,k+1)*         &
                             (qd(i,k+1)-qhat(i,k+1)+cpg/rl*       &
                             (sd(i,k+1)-s(i,k)))-md(i,k)*         &
                             (qd(i,k)-qhat(i,k)+cpg/rl*           &
                             (sd(i,k)-s(i,k))))/wdp(i,k) +        &
                             du(i,k)* (beta*ql(i,k)+              &
                             (1._dp-beta)*ql(i,k+1))
            end if
         end do
      end do

      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.ge.lel(i) .and. k.le.lcl(i)) then
               thetavp(i,k) = tp(i,k)* (1000._dp/p(i,k))** (rgas/cpg)*   &
                             (1._dp+1.608_dp*qstp(i,k)-q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000._dp/p(i,k))** (rgas/cpg)*    &
                             (1._dp+0.608_dp*q(i,k))
               dqsdtp(i,k) = qstp(i,k)* (1._dp+qstp(i,k)/eps1)*eps1*rl/  &
                            (rgas*tp(i,k)**2)

! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.

               dtpdt(i,k) = tp(i,k)/ (1._dp+                             &
                            rl/cpg* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                            (dtbdt(i)/t(i,mx(i))+                     &
                            rl/cpg* (dqbdt(i)/tl(i)-q(i,mx(i))/       &
                            tl(i)**2*dtldt(i)))

! dboydt is the integrand of cape change.

              dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._dp/                 &
                             (1._dp + 1.608_dp*qstp(i,k)-q(i,mx(i)))*   &
                             (1.608_dp * dqsdtp(i,k) * dtpdt(i,k) -     &
                             dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_dp/  &
                             (1._dp+0.608_dp*q(i,k))*dqmdt(i,k)))*grav* &
                             thetavp(i,k)/thetavm(i,k)
            end if
         end do
      end do

      do k = msg + 1,plev
         do i = il1g,il2g
            if (k.gt.lcl(i) .and. k.lt.mx(i)) then
               thetavp(i,k) = tp(i,k)* (1000._dp/p(i,k))** (rgas/cpg)*  &
                              (1._dp+0.608_dp*q(i,mx(i)))
               thetavm(i,k) = t(i,k)* (1000._dp/p(i,k))** (rgas/cpg)*   &
                              (1._dp+0.608_dp*q(i,k))

! dboydt is the integrand of cape change.

               dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+                            &
                              0.608_dp/ (1._dp+0.608_dp*q(i,mx(i)))*dqbdt(i)- &
                              dtmdt(i,k)/t(i,k)-                              &
                              0.608_dp/ (1._dp+0.608_dp*q(i,k))*dqmdt(i,k))*  &
                              grav*thetavp(i,k)/thetavm(i,k)         
            end if
         end do
      end do


! buoyant energy change is set to 2/3*excess cape per 3 hours

      do i = il1g,il2g
         dadt = 0.
         do k = lel(i),mx(i) - 1
            dadt = dadt + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         end do

         dltaa = -1.* (cape(i)-capelmt)
!!         if (dadt.ne.0.) mb(i) = max(dltaa/tau/dadt,0._dp)
         if (ABS(dadt).gt.tiny(0.0_dp)) mb(i) = max(dltaa/tau/dadt,0._dp)
      end do

      return
    end subroutine closure_pjr
!===========================================================================================
    subroutine findcin (plon, plond, plev, plevp, lat, t, q, tprime, qprime,  &
                        mx, pm, pi, zm, lclo, bmax, lfc, lnb, lcape)

!     find the convective inhibition for each parcel
!     define the dry and moist adiabats as lines of constant s and h
!     where s = c_p T + gz, h = s + Lq

!     note that the variations of c_p and R (for hydrostatic equation) with
!     water vapor are ignored

      implicit none

      INTEGER, INTENT(IN) :: plon, plond,  &  ! number of columns, max columns
                             plev, plevp       ! number of levels, levels+1
!     
      INTEGER  :: lat                       ! lat index for this slice
      REAL(dp) :: t(plond,plev)                ! temperature of env
      REAL(dp) :: q(plond, plev)               ! moisture of env
      REAL(dp) :: tprime(plond,plev)           ! tpert for updraft plume
      REAL(dp) :: qprime(plond,plev)           ! qpert for updraft plume
      REAL(dp) :: pm(plond,plev)               ! pressume at layer midpoints
      REAL(dp) :: pi(plond,plevp)              ! pressure at layer interfaces
      REAL(dp), INTENT(IN) :: zm(plond,plev)   ! height at layer midpoints
      INTEGER  :: lclo(plond)               ! lcl of parcel leaving layer mx
      INTEGER  :: mx(plond)                 ! level of parcel with max cin
      REAL(dp) :: bmax(plond)                  ! buoyancy between mx and lcl
      INTEGER  :: klcl(plond,plev)
      INTEGER  :: klfc(plond,plev)
      INTEGER  :: klnb(plond,plev)
      INTEGER  :: lnb(plond)
      INTEGER  :: lfc(plond)
      REAL(dp) :: lcape(plond)
      REAL(dp) :: capen

!     local workspace
      REAL(dp) :: tubase(plond,plev)
      REAL(dp) :: qu(plond,plev)
      REAL(dp) :: bloc2
      REAL(dp) :: tvl
      REAL(dp) :: tev(plond,plev)
      REAL(dp) :: es
      REAL(dp) :: buoy
! rvk:      REAL(dp) :: grav
      REAL(dp) :: c1n
      REAL(dp) :: plcl(plond)
      REAL(dp) :: tparc(plond,plev)
      REAL(dp) :: qparc(plond,plev)
      REAL(dp) :: hold(plond)
      REAL(dp) :: hnew(plond)
      REAL(dp) :: tsat(plond)
      REAL(dp) :: qs(plond)
      REAL(dp) :: dtd
      REAL(dp) :: dtw
      REAL(dp) :: err
      REAL(dp) :: capgam
      REAL(dp) :: dz
      REAL(dp) :: capblfc(plond,plev)
      REAL(dp) :: cape(plond,plev)
      REAL(dp) :: emax
      REAL(dp) :: qold(plond), told(plond)
      REAL(dp) :: ppa(plond), esv(plond)
      REAL(dp) :: gam(plond)
      INTEGER  :: imax

      REAL(dp) :: gmax


      INTEGER  :: kb
      INTEGER  :: i
      INTEGER  :: itmax
      INTEGER  :: k
      INTEGER  :: kt
      INTEGER  :: iter
      INTEGER  :: ktop
      INTEGER  :: nmax

      do k = plev,1,-1
         do i = 1,plon
            qu(i,k) = q(i,k) + qprime(i,k) ! moisture for parcels
            tubase(i,k) = t(i,k) + tprime(i,k) ! temperature for parcels
!            tev(i,k) = t(i,k)*(1.+.608*q(i,k)) ! virtual temperature for env
            tev(i,k) = t(i,k)                           &
                      * (1._dp+1.608_dp*q(i,k))/ (1._dp+q(i,k))
!!            lcl(i,k) = 0
!!          btot(i,k) = -1.e36_dp
            klcl(i,k) = 0
            cape(i,k) = 0.
            capblfc(i,k) = 0.
            klfc(i,k) = 0
            klnb(i,k) = plev+1
         end do
      end do

      do i = 1,plon
         mx(i) = 1
         lclo(i) = 0
         bmax(i) = -1.e36_dp
!!         bpos(i) = 0.
         lcape(i) = -1.e36_dp
      end do

      c1n = grav/cpg
!!      c2n = hlatv/cpg
      capgam = -c1n
!      itmax = 25

      itmax = 500

      ktop = plev+1                     ! highest level to lift a parcel
      kt = plev+1                       ! highest level to start parcels from
      do k = 1,plev
         do i = 1,plon
            if (pm(i,k).gt.100._dp) then
               ktop = min(k-1,ktop)
            endif
            if (pm(i,k).gt.850._dp) then
!            if (pi(i,plev)-pm(i,k).lt.150) then
               kt = min(k-1,kt)
            endif
         end do
      end do
#ifdef DEBCONV
      write (6,*) ' lat ', lat
      write (6,*) ' kt highest possible level of cloud base ', kt
      write (6,*) ' ktop highest possible level of cloud top ', ktop
#endif
            
!     first calculate the properties of parcel at base of lifting
      do kb = plev,kt,-1
#ifdef DEBCONV
         if (lat.eq.latlook) then
            write (6,*) ' kb ', kb
            write (6,*) ' k, lnb, lfc, lcl, tpv, tev, tp, t, p'//       &
                   ' bouy, bloc2, capblfc, cape '
         endif
#endif
         do i = 1,plon
!            write (6,*) ' override ', i
!            tubase(i,kb) = 190.1707794843005
!            qu(i,kb) = 9.833855485759975E-6
!            pm(i,kb) =  647.8308110287944
!            zm(i,kb) = 3117.237132057446
            ppa(i) = pm(i,kb)*100.
         end do
!         tubase(1,kb) = 300.
         call vqsatd(tubase(1,kb), ppa, esv, qs, gam, plon)
!         qu(1,kb) = qs(1)*1.04
         do i = 1,plon

            plcl(i) = 0.

!           start assuming the parcels are unsaturated
            tparc(i,kb) = tubase(i,kb)
            qparc(i,kb) = qu(i,kb)
            hold(i) = cpg*tparc(i,kb)+grav*zm(i,kb)+hlatv*qparc(i,kb)

!           check whether they are saturated
            tsat(i) = tparc(i,kb)

!           if saturated, then first assume water condenses out and
!           releases all latent heat
            if (qs(i) .le. qparc(i,kb)) then
               plcl(i) = max(plcl(i),pm(i,kb))
               klcl(i,kb) = max(klcl(i,kb),kb)
               gmax = max(gam(i),0.1_dp)   ! limit the size of dtw for cold temp
               dtw = hlatv*(qparc(i,kb)-qs(i))/(cpg*gmax)
               tsat(i) = tparc(i,kb) + dtw
!               es = estblf(tsat(i))*0.01
!               qs(i) = eps1*es/(pm(i,kb)-es)
            endif
         end do
#ifdef DEBCONV
         if (lat.eq.latlook) then
            i = ilook
            write (6,*) ' deb 1 ', tsat(i), qs(i), tparc(i,kb), gam(i)
         endif
#endif

!        now correct it so parcel is just saturated and conserves h
         do iter=1,itmax
            do i = 1,plon
               told(i) = tsat(i)
               qold(i) = qs(i)
            end do
            call vqsatd(tsat, ppa, esv, qs, gam, plon)
            nmax = 0
            do i = 1,plon
               if (klcl(i,kb).ne.0) then
                  hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
                  err = hnew(i) - hold(i)
                  tsat(i) = tsat(i) - err/(cpg*(1._dp + gam(i)))
                  es = estblf(tsat(i))*0.01_dp
                  qs(i) = eps1*es/(pm(i,kb)-(1._dp - eps1)*es)

! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere

                  qs(i) = min(1.0_dp,qs(i))

                  if (qs(i) .lt. 0.0) then
                     qs(i) = 0.0
                     es = pm(i,k)
                  end if
                  hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
                  err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
                  if (err.gt.1.e-4_dp) then
                     nmax = nmax+1
                  endif
               endif
            end do
#ifdef DEBCONV
            if (lat.eq.latlook) then
               i = ilook
               write (6,*) ' deb 2 ',                      &
                   iter,hnew(i),hold(i),tsat(i), qs(i)
            endif
#endif


            if (nmax.eq.0) goto 10

          end do     ! iteration loop

!        now check to make sure it converged
         emax = 0.
         imax = 0
         do i = 1,plon
            if (klcl(i,kb).ne.0) then
               hnew(i) = cpg*tsat(i)+grav*zm(i,kb)+hlatv*qs(i)
               tparc(i,kb) = tsat(i)
               qparc(i,kb) = qs(i)
               err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
               if (err.gt.emax) then
                  emax = err
                  imax = i
               endif
            endif
         end do
#ifdef DEBCONV
         if (lat.eq.latlook) then
            i = ilook
            write (6,*) ' deb 3 ',          &
                 iter, hnew(i), hold(i), tsat(i), qs(i)
         endif
#endif

         if (emax.gt.1.e-4_dp) then
            i = imax
            write (6,*) ' problems initializing moist k base ',   &
                 emax, lat, i, klcl(i,kb), hold(i), hnew(i)
            write (6,*) ' told, tnew ',                           &
                 told(i), tsat(i)
            write (6,*) ' qold, qnew ',                           & 
                 qold(i), qs(i)
            write (6,*) ' tubase, qu ', tubase(i,kb), qu(i,kb)
            write (6,*) ' pm, zm, nmax ', pm(i,kb), zm(i,kb), nmax
!            stop
         endif
   10    continue
!         stop 111
         
!        now check it for lfc
         do i = 1,plon
!            tvl = tparc(i,kb)*(1.+.608*qparc(i,kb))
            tvl = tparc(i,kb)                                      &
                    * (1._dp+1.608_dp*qparc(i,kb))/ (1._dp+qparc(i,kb))
            buoy = tvl-tev(i,kb)
            if (klcl(i,kb).gt.0.and.buoy.gt.0.) then
               klfc(i,kb) = max(klfc(i,kb),kb)
            endif
            bloc2 = rgas*buoy*(pi(i,kb+1)-pi(i,kb))/pm(i,kb)
!            bloc2 = rgas*buoy*log(pi(i,kb+1)/pi(i,kb))
            capen = cape(i,kb) + bloc2
            if (klfc(i,kb).eq.0) then
               capblfc(i,kb) = capblfc(i,kb) + bloc2
            endif
            if (buoy.ge.0._dp) then
               klnb(i,kb) = min(klnb(i,kb),kb)
            endif
            if (klfc(i,kb).eq.0.or.buoy.ge.0._dp) then
               cape(i,kb) = capen
            endif
#ifdef DEBCONV
            if (lat.eq.latlook.and.i.eq.ilook) then
               write (6,69) kb, klnb(i,kb), klfc(i,kb),      &
                     klcl(i,kb), tvl, tev(i,kb),             &
                     tparc(i,kb), t(i,kb), pm(i,kb),         &
                     buoy, bloc2,                            &
                     capblfc(i,kb), cape(i,kb)
            endif
#endif
         end do

#ifdef DEBCONV
         write (6,*) ' completed cloud base test for level ',kb
#endif
!        move the parcel up
         do k = kb-1,ktop,-1
            do i = 1,plon
               hold(i) = cpg*tparc(i,k+1)                   &
                    +grav*zm(i,k+1)+hlatv*qparc(i,k+1)

!              first assume it moves up a dry adiabat
               dz = (zm(i,k)-zm(i,k+1))
               dtd = capgam*dz
!!               td = tparc(i,k+1) + dtd
!              es = exp( -6763.6/td - 4.9283*log(td) + 54.2190 )
!               es = estblf(td)*0.01
!               qs(i) = eps1*es/(pm(i,k)-es)
               tparc(i,k) = tparc(i,k+1) + dtd
               qparc(i,k) = qparc(i,k+1)

               hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
               ppa(i) = pm(i,k)*100.
            end do

!           if it hits saturation, note it, and move up a wet adiabat
!           note this is very approximate because
!           1) qs is assumed linear
!           2) the lcl was assumed occurred precisely at level k+1

            call vqsatd(tparc(1,k), ppa, esv, qs, gam, plon)
            do i = 1,plon
               if (qs(i) .lt. qparc(i,k+1)) then
                  plcl(i) = max(plcl(i),pm(i,k))
                  klcl(i,kb) = max(klcl(i,kb),k)
                  dtw = capgam*dz/(1._dp + gam(i))
                  tsat(i) = tparc(i,k+1) + dtw
                  es = estblf(tsat(i))*0.01_dp
                  qs(i) = eps1*es/(pm(i,k)-(1._dp - eps1)*es)

! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere

                  qs(i) = min(1.0_dp,qs(i))

                  if (qs(i) .lt. 0.0) then
                     qs(i) = 0.0
                     es = pm(i,k)
                  end if
                  tparc(i,k) = tsat(i)
               endif
            end do
            
!           now correct any parcels that have hit their lcl
!           so they conserve h, and are saturated.
!           using a newton method iteration
!           this remove the error in assumptions mentioned above
            do iter = 1,itmax
               call vqsatd(tparc(1,k), ppa, esv, qs, gam, plon)
               nmax = 0
               do i = 1,plon
                  if (klcl(i,kb).ne.0) then
                     hnew(i) = cpg*tsat(i)+grav*zm(i,k)+hlatv*qs(i)
                     err = hnew(i) - hold(i)
                     dtw =  - err/(cpg*(1._dp + gam(i)))
                     tsat(i) = tsat(i) + dtw
                     es = estblf(tsat(i))*0.01_dp
                     qs(i) = eps1*es/(pm(i,k)-(1._dp - eps1)*es)

! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere

                     qs(i) = min(1.0_dp,qs(i))

                     if (qs(i) .lt. 0.0) then
                        qs(i) = 0.0
                        es = pm(i,k)
                     end if
#ifdef DEBCONV
                     if (lat.eq.latlook.and.i.eq.ilook) then
                        write (6,*) ' deb 4 ', k,hnew(i), hold(i), dtw,      &
                            tsat(i), pm(i,k),gam(i)
                     endif
#endif
                     tparc(i,k) = tsat(i)
                     qparc(i,k) = qs(i)
                  end if
                  hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
                  err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
                  if (err.gt.1.e-4_dp) then
                     nmax = nmax + 1
                  endif
               end do
               if (nmax.eq.0) goto 20
             end do

!           now check to see whether it converged.
            emax = 0.
            imax = 0
            do i = 1,plon
               hnew(i) = cpg*tparc(i,k)+grav*zm(i,k)+hlatv*qparc(i,k)
               err = abs((hnew(i)-hold(i))/(hnew(i)+hold(i)))
               if (err.gt.emax) then
                  emax = err
                  imax = i
               endif
            end do
            if (emax.gt.1.e-4_dp) then
               i = imax
               write (6,*) ' problems initializing moist k mid',lat,k,       &
                   i, hold(i), hnew(i), tparc(i,k), zm(i,k), gam(i)
!               stop
            end if
   20       continue

!           now check it for lfc
            do i = 1,plon
!               tvl = tparc(i,k)*(1.+.608*qparc(i,k))
               tvl = tparc(i,k)                                     &
                    * (1._dp+1.608_dp*qparc(i,k))/ (1._dp+qparc(i,k))
               buoy = tvl-tev(i,k)
               if (klcl(i,kb).gt.0.and.buoy.gt.0._dp) then
                  klfc(i,kb) = max(klfc(i,kb),k)
               endif
               bloc2 = rgas*buoy*(pi(i,k+1)-pi(i,k))/pm(i,k)
               capen = cape(i,kb) + bloc2
               if (klfc(i,kb).eq.0) then
                  capblfc(i,kb) = capblfc(i,kb) + bloc2
               endif
               if (buoy.ge.0._dp) then
                  klnb(i,kb) = min(klnb(i,kb),k)
               endif
               if (klfc(i,kb).eq.0.or.buoy.ge.0._dp) then
                  cape(i,kb) = capen
               endif
#ifdef DEBCONV
               if (lat.eq.latlook.and.i.eq.ilook) then
                  write (6,*) ' hnew ', hnew(i)
                  write (6,69) k, klnb(i,kb), klfc(i,kb),       &
                       klcl(i,kb), tvl, tev(i,k),               &
                       tparc(i,k), t(i,k), pm(i,k),             &
                       buoy, bloc2,                             &   
                       capblfc(i,kb), cape(i,kb)             
   69              format (4i5,10f11.3)
               endif
#endif
            end do
         end do

      end do

!     now choose the cloud base
      do kb = plev,kt,-1
         do i = 1,plon
!            if (capblfc(i,kb).gt.bmax(i)) then
!            if (capblfc(i,kb).ge.bmax(i)
!     $           .and.cape(i,kb).gt.lcape(i)) then
!           choose a layer with the least inhibition
!           till we find a positive capeblfc, then stop
            if (capblfc(i,kb).ge.bmax(i)        &
                .and.bmax(i).lt.0.) then
               mx(i) = kb
               lclo(i) = klcl(i,kb)
               bmax(i) = capblfc(i,kb)
               lfc(i) = klfc(i,kb)
               lnb(i) = klnb(i,kb)
               lcape(i) = cape(i,kb)
            endif
         end do
      end do
#ifdef DEBCONV
      if (lat.eq.latlook) then
         i = ilook
         write (6,*) ' properties of best parcel ',                     &
             'kb, mx(i), lclo(i), bmax(i), lfc(i), lnb(i), lcape(i)',   &
             kb, mx(i), lclo(i), bmax(i), lfc(i), lnb(i), lcape(i)
      endif
#endif

      return
    end subroutine findcin


!------------------ subroutines for conv_ccm_pjr end ---------------------------
!functions:
!--------------------------------------------------------------_----------------
    ELEMENTAL REAL(dp) FUNCTION  qhalf(sh1, sh2, shbs1, shbs2)  
      
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: sh1, sh2, shbs1, shbs2
      

      qhalf = min(max(sh1,sh2),                         &
              (shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
    END FUNCTION qhalf

!===============================================================================
END MODULE MESSY_CONVECT_ZHANG
