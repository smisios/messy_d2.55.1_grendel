 module messy_vertdiff_camdiffsolver

  !------------------------------------------------------------------------------------ !
  ! Module to solve vertical diffusion equations using a tri-diagonal solver.           !
  ! The module will also apply countergradient fluxes, and apply molecular              ! 
  ! diffusion for constituents.                                                         !
  !                                                                                     !
  ! Public interfaces :                                                                 ! 
  !    compute_vdiff    solves diffusion equations                                      !
  !    vd_lu_solve      tridiagonal solver also used by gwd (should be private)         !
  !    vd_lu_decomp     tridiagonal solver also used by gwd (should be private)         !
  !                                                                                     !
  !------------------------------------ Code History ---------------------------------- !
  ! Initial subroutines :  B. Boville and others, 1991-2004                             !
  ! Modularization      :  J. McCaa, September 2004                                     !
  ! Most Recent Code    :  Sungsu Park, Aug. 2006, Dec. 2008, Jan. 2010.                !
  !------------------------------------------------------------------------------------ !

  use messy_main_constants_mem, ONLY: r8 => dp           &
                                    , cpair  => cp_air   &     ! Specific heat of dry air
                                    , gravit => g        &     ! Acceleration due to gravity
                                    , rair   => rd       &     ! Gas constant for dry air
                                    , zvir   => vtmpc1   &     ! rh2o/rair - 1
                                    , latvap => alv      &     ! Latent heat of vaporization
                                    , karman => c_vKar         ! von Karman constant

  implicit none
  private       
  save


  ! ----------------- !
  ! Public interfaces !
  ! ----------------- !
  public compute_vdiff                                   ! Full routine
  public vd_lu_solve                                     ! Tridiagonal solver also used by gwd ( should be private! )
  public vd_lu_decomp                                    ! Tridiagonal solver also used by gwd ( should be private! )

 contains


  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine compute_vdiff( pcols           , pver               , ncnst         , ncol         , pmid        , &
                            pint            , rpdel              , t             , ztodt        , taux        , &
                            tauy            , shflx              , cflx_q        , cflx_xl      , cflx_xi     , &  
                            cflx            , ntop               , nbot          ,                              &
                            kvh             , kvm                , kvq           , cgs          , cgh         , &
                            ksrftms         , fieldlist,                                                        &
                            u               , v                  , q             , dse          , xl          , &
                            xi              , xt                 ,                                              &
                            tautmsx         , tautmsy            , dtk           ,                              &
                            tauresx         , tauresy            , itaures       , do_iss                       &
                            )

    !-------------------------------------------------------------------------- !
    ! Driver routine to compute vertical diffusion of momentum, moisture, trace !
    ! constituents and dry static energy. The new temperature is computed from  !
    ! the diffused dry static energy.                                           ! 
    ! Turbulent diffusivities and boundary layer nonlocal transport terms are   !
    ! obtained from the turbulence module.                                      !
    !-------------------------------------------------------------------------- !

    ! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
    !                Also, vertical diffusion of cloud droplet number concentration and aerosol number
    !                concentration should be done very carefully in the future version.
    ! --------------- !
    ! Input Arguments !
    ! --------------- !

    integer,  intent(in)    :: pcols
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol            ! Number of atmospheric columns
    integer,  intent(in)    :: ntop            ! Top    interface level to which vertical diffusion is applied ( = 1 ).
    integer,  intent(in)    :: nbot            ! Bottom interface level to which vertical diffusion is applied ( = pver ).
    integer,  intent(in)    :: itaures         ! Indicator determining whether 'tauresx,tauresy'
                                               ! is updated (1) or non-updated (0) in this subroutine.

    real(r8), intent(in)    :: pmid(:,:)       ! Mid-point pressures [ Pa ]
    real(r8), intent(in)    :: pint(:,:)       ! Interface pressures [ Pa ]
    real(r8), intent(in)    :: rpdel(:,:)      ! 1./pdel
    real(r8), intent(in)    :: t(:,:)          ! Temperature [ K ]
    real(r8), intent(in)    :: ztodt           ! 2 delta-t [ s ]
    real(r8), intent(in)    :: taux(:)         ! Surface zonal      stress.
                                               ! Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: tauy(:)         ! Surface meridional stress.
                                               ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: shflx(:)        ! Surface sensible heat flux [ W/m2 ]
    real(r8), intent(in)    :: cflx_q(:)       ! Surface q flux [ kg/m2/s ]
    real(r8), intent(in)    :: cflx_xl(:)      ! Surface xl flux [ kg/m2/s ]
    real(r8), intent(in)    :: cflx_xi(:)      ! Surface xi flux [ kg/m2/s ]
    real(r8), intent(in)    :: cflx(:,:)       ! Surface constituent flux [ kg/m2/s ]
!!$    real(r8), intent(in)    :: zi(:,:)         ! Interface heights [ m ]
    real(r8), intent(in)    :: ksrftms(:)      ! Susrface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
!!$    real(r8), intent(in)    :: qmincg(ncnst)      ! Minimum constituent mixing ratios from cg fluxes
    real(r8), intent(in)    :: kvh(:,:)         ! Eddy diffusivity for heat [ m2/s ]

    logical, intent(in)     :: fieldlist(:)        ! Array of flags selecting which fields to diffuse
    logical, intent(in)     :: do_iss

    ! ---------------------- !
    ! Input-Output Arguments !
    ! ---------------------- !

    real(r8), intent(inout) :: kvm(:,:)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    real(r8), intent(inout) :: kvq(:,:)         ! Eddy diffusivity for constituents
    real(r8), intent(inout) :: cgs(:,:)         ! Counter-gradient star [ cg/flux ]
    real(r8), intent(inout) :: cgh(:,:)         ! Counter-gradient term for heat

    real(r8), intent(inout) :: u(:,:)             ! U wind. This input is the 'raw' input wind to
                                                         ! PBL scheme without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: v(:,:)             ! V wind. This input is the 'raw' input wind to PBL scheme
                                                         ! without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: q(:,:)             ! Moisture [ kg/kg, #/kg ? ]
    real(r8), intent(inout) :: dse(:,:)           ! Dry static energy [ J/kg ]
    real(r8), intent(inout) :: xl(:,:)            ! 
    real(r8), intent(inout) :: xi(:,:)            ! 
    real(r8), intent(inout) :: xt(:,:,:)      ! trace constituents [ kg/kg, #/kg ? ]

    real(r8), intent(inout) :: tauresx(:)            ! Input  : Reserved surface stress at previous time step
    real(r8), intent(inout) :: tauresy(:)            ! Output : Reserved surface stress at current  time step

    ! ---------------- !
    ! Output Arguments !
    ! ---------------- !

    real(r8), intent(out)   :: dtk(:,:)           ! T tendency from KE dissipation
    real(r8), intent(out)   :: tautmsx(:)            ! Implicit zonal      turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: tautmsy(:)            ! Implicit meridional turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]


    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    integer  :: i, k, m, icol                            ! Longitude, level, constituent indices
    integer  :: status                                   ! Status indicator
    integer  :: nbot_molec                               ! Bottom level where molecular diffusivity is applied
    integer  :: ntop_molec                               ! Top level where molecular diffusivity is applied
    logical  :: lqtst(pcols)                             ! Adjust vertical profiles
    logical  :: need_decomp                              ! Whether to compute a new decomposition

    real(r8) :: tmpm(pcols,pver)                         ! Potential temperature, ze term in tri-diag sol'n
    real(r8) :: ca(pcols,pver)                           ! - Upper diag of tri-diag matrix
    real(r8) :: cc(pcols,pver)                           ! - Lower diag of tri-diag matrix
    real(r8) :: dnom(pcols,pver)                         ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))

    real(r8) :: tmp1(pcols)                              ! Temporary storage
    real(r8) :: tmpi1(pcols,pver+1)                      ! Interface KE dissipation
    real(r8) :: tint(pcols,pver+1)                       ! Interface temperature
    real(r8) :: rhoi(pcols,pver+1)                       ! rho at interfaces
    real(r8) :: tmpi2(pcols,pver+1)                      ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: rrho(pcols)                              ! 1./bottom level density 

    real(r8) :: zero(pcols)                              ! Zero array for surface heat exchange coefficients 
    real(r8) :: tautotx(pcols)                           ! Total surface stress ( zonal )
    real(r8) :: tautoty(pcols)                           ! Total surface stress ( meridional )

    real(r8) :: dinp_u(pcols,pver+1)                     ! Vertical difference at interfaces, input u
    real(r8) :: dinp_v(pcols,pver+1)                     ! Vertical difference at interfaces, input v
    real(r8) :: dout_u                                   ! Vertical difference at interfaces, output u
    real(r8) :: dout_v                                   ! Vertical difference at interfaces, output v
    real(r8) :: dse_top(pcols)                           ! dse on top boundary
    real(r8) :: cc_top(pcols)                            ! Lower diagonal at top interface
    real(r8) :: cd_top(pcols)                            ! 
    real(r8) :: rghd(pcols,pver+1)                       ! (1/H_i - 1/H) *(g*rho)^(-1)

    real(r8) :: qtm(pcols,pver)                          ! Temporary copy of q
!!$    real(r8) :: mw_fac(ncnst)                            ! sqrt(1/M_q + 1/M_d) for this constituent
!!$    real(r8) :: cnst_mw(ncnst)                           ! Molecular weight [ kg/kmole ]
!!$    real(r8) :: ubc_mmr(pcols,ncnst)                     ! Upper boundary mixing ratios [ kg/kg ]
!!$    real(r8) :: ubc_flux(ncnst)                          ! Upper boundary flux [ kg/s/m^2 ]
    real(r8) :: ubc_t(pcols)                             ! Upper boundary temperature [ K ]

    real(r8) :: ws(pcols)                                ! Lowest-level wind speed [ m/s ]
    real(r8) :: tau(pcols)                               ! Turbulent surface stress ( not including mountain stress )
    real(r8) :: ksrfturb(pcols)                          ! Surface drag coefficient of 'normal' stress. > 0.
                                                         ! Virtual mass input per unit time per unit area [ kg/s/m2 ]
    real(r8) :: ksrf(pcols)                              ! Surface drag coefficient of 'normal' stress +
                                                         ! Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
    real(r8) :: usum_in(pcols)                           ! Vertical integral of input u-momentum. Total zonal
                                                         ! momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    real(r8) :: vsum_in(pcols)                           ! Vertical integral of input v-momentum. Total meridional
                                                         ! momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    real(r8) :: usum_mid(pcols)                          ! Vertical integral of u-momentum after adding explicit residual stress
    real(r8) :: vsum_mid(pcols)                          ! Vertical integral of v-momentum after adding explicit residual stress
    real(r8) :: usum_out(pcols)                          ! Vertical integral of u-momentum after doing implicit diffusion
    real(r8) :: vsum_out(pcols)                          ! Vertical integral of v-momentum after doing implicit diffusion
    real(r8) :: tauimpx(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: tauimpy(pcols)                           ! Actual net stress added at the current step other than mountain stress
    real(r8) :: wsmin                                    ! Minimum sfc wind speed for estimating frictional
                                                         ! transfer velocity ksrf. [ m/s ]
    real(r8) :: ksrfmin                                  ! Minimum surface drag coefficient [ kg/s/m^2 ]
    real(r8) :: timeres                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]
    real(r8) :: ramda                                    ! dt/timeres [ no unit ]
    real(r8) :: psum
    real(r8) :: u_in, u_res, tauresx_in
    real(r8) :: v_in, v_res, tauresy_in  
    
!!$    real(r8) :: mw_fac_loc(pcols,pver+1,ncnst)           ! Local sqrt(1/M_q + 1/M_d) for this constituent

    !--------------------------------
    ! Variables needed for WACCM-X
    !--------------------------------
    real(r8) :: ttemp(pcols,pver)                   ! temporary temperature array
    real(r8) :: ttemp0(pcols,pver)                   ! temporary temperature array

    ! ------------------------------------------------ !
    ! Parameters for implicit surface stress treatment !
    ! ------------------------------------------------ !

    wsmin    = 1._r8                                     ! Minimum wind speed for ksrfturb computation        [ m/s ]
    ksrfmin  = 1.e-4_r8                                  ! Minimum surface drag coefficient                   [ kg/s/m^2 ]
    timeres  = 7200._r8                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

!!$    errstring = ''
!!$    if( ( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) .and. .not. diffuse(fieldlist,'s') ) then
!!$          errstring = 'diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v'
!!$          return
!!$    end if
    zero(:) = 0._r8

    ! Compute 'rho' and 'dt*(g*rho)^2/dp' at interfaces

    tint(:ncol,1) = t(:ncol,1)
    rhoi(:ncol,1) = pint(:ncol,1) / (rair*tint(:ncol,1))
    do k = 2, pver
       do i = 1, ncol
          tint(i,k)  = 0.5_r8 * ( t(i,k) + t(i,k-1) )
          rhoi(i,k)  = pint(i,k) / (rair*tint(i,k))
          tmpi2(i,k) = ztodt * ( gravit*rhoi(i,k) )**2 / ( pmid(i,k) - pmid(i,k-1) )
       end do
    end do
    tint(:ncol,pver+1) = t(:ncol,pver)
    rhoi(:ncol,pver+1) = pint(:ncol,pver+1) / ( rair*tint(:ncol,pver+1) )

    rrho(:ncol) = rair  * t(:ncol,pver) / pmid(:ncol,pver)
    tmp1(:ncol) = ztodt * gravit * rpdel(:ncol,pver)

    cd_top(:)    = 0._r8
    cc_top(:)    = 0._r8


    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !

!!$    if( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) then

        ! Compute the vertical upward differences of the input u,v for KE dissipation
        ! at each interface.
        ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
        ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

        do i = 1, ncol
           dinp_u(i,1) = 0._r8
           dinp_v(i,1) = 0._r8
           dinp_u(i,pver+1) = -u(i,pver)
           dinp_v(i,pver+1) = -v(i,pver)
        end do
        do k = 2, pver
           do i = 1, ncol
              dinp_u(i,k) = u(i,k) - u(i,k-1)
              dinp_v(i,k) = v(i,k) - v(i,k-1)
           end do
        end do

       ! -------------------------------------------------------------- !
       ! Do 'Implicit Surface Stress' treatment for numerical stability !
       ! in the lowest model layer.                                     !
       ! -------------------------------------------------------------- !

       if( do_iss ) then

         ! Compute surface drag coefficient for implicit diffusion 
         ! including turbulent turbulent mountain stress. 

           do i = 1, ncol
              ws(i)       = max( sqrt( u(i,pver)**2._r8 + v(i,pver)**2._r8 ), wsmin )
              tau(i)      = sqrt( taux(i)**2._r8 + tauy(i)**2._r8 )
              ksrfturb(i) = max( tau(i) / ws(i), ksrfmin )
           end do              
           ksrf(:ncol) = ksrfturb(:ncol) + ksrftms(:ncol)  ! Do all surface stress ( normal + tms ) implicitly

         ! Vertical integration of input momentum. 
         ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
         ! Note (u,v) are the raw input to the PBL scheme, not the
         ! provisionally-marched ones within the iteration loop of the PBL scheme.  

           do i = 1, ncol
              usum_in(i) = 0._r8
              vsum_in(i) = 0._r8
              do k = 1, pver
                 usum_in(i) = usum_in(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_in(i) = vsum_in(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do
           end do              

         ! Add residual stress of previous time step explicitly into the lowest
         ! model layer with a relaxation time scale of 'timeres'.

           ramda         = ztodt / timeres
           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*tauresx(:ncol)*ramda
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauresy(:ncol)*ramda

         ! Vertical integration of momentum after adding explicit residual stress
         ! into the lowest model layer.

           do i = 1, ncol
              usum_mid(i) = 0._r8
              vsum_mid(i) = 0._r8
              do k = 1, pver
                 usum_mid(i) = usum_mid(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_mid(i) = vsum_mid(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do
           end do              

         ! Debug 
         ! icol = phys_debug_col(lchnk) 
         ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
         !      tauresx_in = tauresx(icol)
         !      tauresy_in = tauresy(icol)
         !      u_in  = u(icol,pver) - tmp1(icol) * tauresx(icol) * ramda
         !      v_in  = v(icol,pver) - tmp1(icol) * tauresy(icol) * ramda
         !      u_res = u(icol,pver)
         !      v_res = v(icol,pver)
         ! endif
         ! Debug

       else

         ! In this case, do 'turbulent mountain stress' implicitly, 
         ! but do 'normal turbulent stress' explicitly.
         ! In this case, there is no 'redisual stress' as long as 'tms' is
         ! treated in a fully implicit wway, which is true.

         ! 1. Do 'tms' implicitly

           ksrf(:ncol) = ksrftms(:ncol) 

         ! 2. Do 'normal stress' explicitly

           u(:ncol,pver) = u(:ncol,pver) + tmp1(:ncol)*taux(:ncol)
           v(:ncol,pver) = v(:ncol,pver) + tmp1(:ncol)*tauy(:ncol)

       end if  ! End of 'do iss' ( implicit surface stress )

       ! --------------------------------------------------------------------------------------- !
       ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
       ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
       !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
       !                                   u(pver) : explicitly include 'redisual normal' stress !
       !    For explicit 'normal' stress : ksrf = ksrftms                                        !
       !                                   u(pver) : explicitly include 'normal' stress          !
       ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
       ! --------------------------------------------------------------------------------------- !

       call vd_lu_decomp( pcols , pver , ncol  ,                        &
                          ksrf  , kvm  , tmpi2 , rpdel , ztodt , zero , &
                          ca    , cc   , dnom  , tmpm  , ntop  , nbot )

       call vd_lu_solve(  pcols , pver , ncol  ,                        &
                          u     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       call vd_lu_solve(  pcols , pver , ncol  ,                        &
                          v     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       ! ---------------------------------------------------------------------- !
       ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
       ! have been actually added into the atmosphere at the current time step. ! 
       ! Also, update residual stress, if required.                             !
       ! ---------------------------------------------------------------------- !

       do i = 1, ncol

          ! Compute the implicit 'tms' using the updated winds.
          ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
          ! that has been actually added into the atmosphere both for explicit
          ! and implicit approach. 

          tautmsx(i) = -ksrftms(i)*u(i,pver)
          tautmsy(i) = -ksrftms(i)*v(i,pver)

          if( do_iss ) then

            ! Compute vertical integration of final horizontal momentum

              usum_out(i) = 0._r8
              vsum_out(i) = 0._r8
              do k = 1, pver
                 usum_out(i) = usum_out(i) + (1._r8/gravit)*u(i,k)/rpdel(i,k)
                 vsum_out(i) = vsum_out(i) + (1._r8/gravit)*v(i,k)/rpdel(i,k)
              end do

            ! Compute net stress added into the atmosphere at the current time step.
            ! Note that the difference between 'usum_in' and 'usum_out' are induced
            ! by 'explicit residual stress + implicit total stress' for implicit case, while
            ! by 'explicit normal   stress + implicit tms   stress' for explicit case. 
            ! Here, 'tautotx(i)' is net stress added into the air at the current time step.

              tauimpx(i) = ( usum_out(i) - usum_in(i) ) / ztodt
              tauimpy(i) = ( vsum_out(i) - vsum_in(i) ) / ztodt

              tautotx(i) = tauimpx(i) 
              tautoty(i) = tauimpy(i) 

            ! Compute redisual stress and update if required.
            ! Note that the total stress we should have added at the current step is
            ! the sum of 'taux(i) - ksrftms(i)*u(i,pver) + tauresx(i)'.

              if( itaures .eq. 1 ) then
                  tauresx(i) = taux(i) + tautmsx(i) + tauresx(i) - tauimpx(i)
                  tauresy(i) = tauy(i) + tautmsy(i) + tauresy(i) - tauimpy(i)
              endif

          else

              tautotx(i) = tautmsx(i) + taux(i)
              tautoty(i) = tautmsy(i) + tauy(i)
              tauresx(i) = 0._r8
              tauresy(i) = 0._r8

          end if  ! End of 'do_iss' if

       end do ! End of 'do i = 1, ncol' loop

     ! Debug 
     ! icol = phys_debug_col(lchnk) 
     ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
     !      write(iulog,*)
     !      write(iulog,*)  'diffusion_solver debug'  
     !      write(iulog,*)
     !      write(iulog,*)  'u_in, u_res, u_out'
     !      write(iulog,*)   u_in, u_res, u(icol,pver)
     !      write(iulog,*)  'tauresx_in, tautmsx, tauimpx(actual), tauimpx(derived), tauresx_out, taux'
     !      write(iulog,*)   tauresx_in, tautmsx(icol), tauimpx(icol), -ksrf(icol)*u(icol,pver), tauresx(icol), taux(icol)
     !      write(iulog,*)
     !      write(iulog,*)  'v_in, v_res, v_out'
     !      write(iulog,*)   v_in, v_res, v(icol,pver)
     !      write(iulog,*)  'tauresy_in, tautmsy, tauimpy(actual), tauimpy(derived), tauresy_out, tauy'
     !      write(iulog,*)   tauresy_in, tautmsy(icol), tauimpy(icol), -ksrf(icol)*v(icol,pver), tauresy(icol), tauy(icol)
     !      write(iulog,*)
     !      write(iulog,*)  'itaures, ksrf, ksrfturb, ksrftms'
     !      write(iulog,*)   itaures, ksrf(icol), ksrfturb(icol), ksrftms(icol)
     !      write(iulog,*) 
     ! endif
     ! Debug

       ! ------------------------------------ !
       ! Calculate kinetic energy dissipation !
       ! ------------------------------------ !       

     ! Modification : In future, this should be set exactly same as 
     !                the ones in the convection schemes 

       ! 1. Compute dissipation term at interfaces
       !    Note that 'u,v' are already diffused wind, and 'tautotx,tautoty' are 
       !    implicit stress that has been actually added. On the other hand,
       !    'dinp_u, dinp_v' were computed using non-diffused input wind.

     ! Modification : I should check whether non-consistency between 'u' and 'dinp_u'
     !                is correctly intended approach. I think so.

       k = pver + 1
       do i = 1, ncol
          tmpi1(i,1) = 0._r8
          tmpi1(i,k) = 0.5_r8 * ztodt * gravit * &
                       ( (-u(i,k-1) + dinp_u(i,k))*tautotx(i) + (-v(i,k-1) + dinp_v(i,k))*tautoty(i) )
       end do

       do k = 2, pver
          do i = 1, ncol
             dout_u = u(i,k) - u(i,k-1)
             dout_v = v(i,k) - v(i,k-1)
             tmpi1(i,k) = 0.25_r8 * tmpi2(i,k) * kvm(i,k) * &
                          ( dout_u**2 + dout_v**2 + dout_u*dinp_u(i,k) + dout_v*dinp_v(i,k) )
          end do
       end do

       ! 2. Compute dissipation term at midpoints, add to dry static energy

       do k = 1, pver
          do i = 1, ncol
             dtk(i,k) = ( tmpi1(i,k+1) + tmpi1(i,k) ) * rpdel(i,k)
             dse(i,k) = dse(i,k) + dtk(i,k)
          end do
       end do

!!$    end if ! End of diffuse horizontal momentum, diffuse(fieldlist,'u') routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

  ! Modification : In future, we should diffuse the fully conservative 
  !                moist static energy,not the dry static energy.

!!$    if( diffuse(fieldlist,'s') ) then

      ! Add counter-gradient to input static energy profiles
        do k = 1, pver
           dse(:ncol,k) = dse(:ncol,k) + ztodt * rpdel(:ncol,k) * gravit  *                &
                                       ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgh(:ncol,k+1) &
                                       - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgh(:ncol,k  ) )
       end do

     ! Add the explicit surface fluxes to the lowest layer

       dse(:ncol,pver) = dse(:ncol,pver) + tmp1(:ncol) * shflx(:ncol)

       call vd_lu_decomp( pcols , pver , ncol  ,                         &
                          zero  , kvh  , tmpi2 , rpdel , ztodt , cc_top, &
                          ca    , cc   , dnom  , tmpm  , ntop  , nbot    )

       call vd_lu_solve(  pcols , pver , ncol  ,                         &
                          dse   , ca   , tmpm  , dnom  , ntop  , nbot  , cd_top )



!!$    endif

    !---------------------------- !
    ! Diffuse Water Vapor Tracers !
    !---------------------------- !
       
    q(:ncol,pver) = q(:ncol,pver) + tmp1(:ncol) * cflx_q(:ncol)
    call vd_lu_solve(  pcols , pver , ncol  ,                         &
         q , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )

    xl(:ncol,pver) = xl(:ncol,pver) + tmp1(:ncol) * cflx_xl(:ncol)
    call vd_lu_solve(  pcols , pver , ncol  ,                         &
         xl , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )

    xi(:ncol,pver) = xi(:ncol,pver) + tmp1(:ncol) * cflx_xi(:ncol)
    call vd_lu_solve(  pcols , pver , ncol  ,                         &
         xi , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )

  ! Modification : For aerosols, I need to use separate treatment 
  !                for aerosol mass and aerosol number. 


    ! Loop through constituents
! mz_ab_20130806 decomp not necessary because cc_top = zero
    need_decomp = .true.
!!$    need_decomp = .false.
    do m = 1, ncnst
       if (fieldlist(m)) then

           ! Add the nonlocal transport terms to constituents in the PBL.
           ! Check for neg q's in each constituent and put the original vertical
           ! profile back if a neg value is found. A neg value implies that the
           ! quasi-equilibrium conditions assumed for the countergradient term are
           ! strongly violated.

           qtm(:ncol,:pver) = xt(:ncol,:pver,m)

           !! cflx implemented for co2 only, ignore for now ! mz_ab_20130824
           do k = 1, pver
              xt(:ncol,k,m) = xt(:ncol,k,m) + &
                             ztodt * rpdel(:ncol,k) * gravit  * ( cflx(:ncol,m) * rrho(:ncol) ) * &
                           ( rhoi(:ncol,k+1) * kvh(:ncol,k+1) * cgs(:ncol,k+1)                    &
                           - rhoi(:ncol,k  ) * kvh(:ncol,k  ) * cgs(:ncol,k  ) )
           end do
!!$           lqtst(:ncol) = all(xt(:ncol,1:pver,m) >= qmincg(m), 2)
           lqtst(:ncol) = all(xt(:ncol,1:pver,m) >= 1.e-36_r8, 2)
           do k = 1, pver
              xt(:ncol,k,m) = merge( xt(:ncol,k,m), qtm(:ncol,k), lqtst(:ncol) )
           end do

           ! Add the explicit surface fluxes to the lowest layer
           xt(:ncol,pver,m) = xt(:ncol,pver,m) + tmp1(:ncol) * cflx(:ncol,m)

           ! Diffuse constituents.
           if( need_decomp ) then
               call vd_lu_decomp( pcols , pver , ncol  ,                         &
                                  zero  , kvq  , tmpi2 , rpdel , ztodt , zero  , &
                                  ca    , cc   , dnom  , tmpm  , ntop  , nbot )
               need_decomp =  .false.
           end if

           call vd_lu_solve(  pcols , pver , ncol  ,                         &
                              xt(:,:,m) , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )
       end if
    end do

    return
  end subroutine compute_vdiff

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vd_lu_decomp( pcols, pver, ncol ,                        &
                           ksrf , kv  , tmpi , rpdel, ztodt , cc_top, &
                           ca   , cc  , dnom , ze   , ntop  , nbot, cpairv )
    !---------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the ! 
    ! tridiagonal diffusion matrix.                                         ! 
    ! The diagonal elements (1+ca(k)+cc(k)) are not required by the solver. !
    ! Also determine ze factor and denominator for ze and zf (see solver).  !
    !---------------------------------------------------------------------- !

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    integer,  intent(in)  :: pcols                 ! Number of allocated atmospheric columns
    integer,  intent(in)  :: pver                  ! Number of allocated atmospheric levels 
    integer,  intent(in)  :: ncol                  ! Number of computed atmospheric columns
    integer,  intent(in)  :: ntop                  ! Top level to operate on
    integer,  intent(in)  :: nbot                  ! Bottom level to operate on
    real(r8), intent(in)  :: ksrf(pcols)           ! Surface "drag" coefficient [ kg/s/m2 ]
    real(r8), intent(in)  :: kv(pcols,pver+1)      ! Vertical diffusion coefficients [ m2/s ]
    real(r8), intent(in)  :: tmpi(pcols,pver+1)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in)  :: rpdel(pcols,pver)     ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in)  :: ztodt                 ! 2 delta-t [ s ]
    real(r8), intent(in)  :: cc_top(pcols)         ! Lower diagonal on top interface (for fixed ubc only)

    real(r8), intent(out) :: ca(pcols,pver)        ! Upper diagonal
    real(r8), intent(out) :: cc(pcols,pver)        ! Lower diagonal
    real(r8), intent(out) :: dnom(pcols,pver)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    real(r8), intent(out) :: ze(pcols,pver)        ! Term in tri-diag. matrix system

    real(r8), intent(in), optional  :: cpairv(pcols,pver) ! "Variable" specific heat at constant pressure

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer :: i                                   ! Longitude index
    integer :: k                                   ! Vertical  index

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
    ! a combination of ca and cc; they are not required by the solver.

    !
    !  If switch present and set true, then input kv is kvt for use in diagonal calculations
    !

    if ( present(cpairv) ) then
      do k = nbot-1, ntop, -1
        do i = 1, ncol
           ca(i,k  ) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k  ) / cpairv(i,k)
           cc(i,k+1) = kv(i,k+1)*tmpi(i,k+1)*rpdel(i,k+1) / cpairv(i,k+1)
        end do
      end do
    else
      do k = nbot - 1, ntop, -1
        do i = 1, ncol
          ca(i,k  ) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k  )
          cc(i,k+1) = kv(i,k+1) * tmpi(i,k+1) * rpdel(i,k+1)
        end do
      end do
    endif

    ! The bottom element of the upper diagonal (ca) is zero (element not used).
    ! The subdiagonal (cc) is not needed in the solver.

    do i = 1, ncol
       ca(i,nbot) = 0._r8
    end do

    ! Calculate e(k).  This term is 
    ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

    do i = 1, ncol
       dnom(i,nbot) = 1._r8/(1._r8 + cc(i,nbot) + ksrf(i)*ztodt*gravit*rpdel(i,nbot))
       ze(i,nbot)   = cc(i,nbot)*dnom(i,nbot)
    end do

    do k = nbot - 1, ntop + 1, -1
       do i = 1, ncol
          dnom(i,k) = 1._r8/(1._r8 + ca(i,k) + cc(i,k) - ca(i,k)*ze(i,k+1))
          ze(i,k)   = cc(i,k)*dnom(i,k)
       end do
    end do

    do i = 1, ncol
       dnom(i,ntop) = 1._r8/(1._r8 + ca(i,ntop) + cc_top(i) - ca(i,ntop)*ze(i,ntop+1))
    end do

    return
  end subroutine vd_lu_decomp

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !

  subroutine vd_lu_solve( pcols , pver , ncol , &
                          q     , ca   , ze   , dnom , ntop , nbot , cd_top )
    !----------------------------------------------------------------------------------- !
    ! Solve the implicit vertical diffusion equation with zero flux boundary conditions. !
    ! Procedure for solution of the implicit equation follows Richtmyer and              !
    ! Morton (1967,pp 198-200).                                                          !
    !                                                                                    !
    ! The equation solved is                                                             !
    !                                                                                    !  
    !     -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k),                              !
    !                                                                                    !
    ! where d(k) is the input profile and q(k) is the output profile                     !
    !                                                                                    ! 
    ! The solution has the form                                                          !
    !                                                                                    !
    !     q(k) = ze(k)*q(k-1) + zf(k)                                                    !
    !                                                                                    !
    !     ze(k) = cc(k) * dnom(k)                                                        !
    !                                                                                    !  
    !     zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)                                       !
    !                                                                                    !
    !     dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]   !
    !                                                                                    !
    ! Note that the same routine is used for temperature, momentum and tracers,          !
    ! and that input variables are replaced.                                             !
    ! ---------------------------------------------------------------------------------- ! 

    ! --------------------- !
    ! Input-Output Argument !
    ! --------------------- !

    integer,  intent(in)    :: pcols                  ! Number of allocated atmospheric columns
    integer,  intent(in)    :: pver                   ! Number of allocated atmospheric levels 
    integer,  intent(in)    :: ncol                   ! Number of computed atmospheric columns
    integer,  intent(in)    :: ntop                   ! Top level to operate on
    integer,  intent(in)    :: nbot                   ! Bottom level to operate on
    real(r8), intent(in)    :: ca(pcols,pver)         ! -Upper diag coeff.of tri-diag matrix
    real(r8), intent(in)    :: ze(pcols,pver)         ! Term in tri-diag solution
    real(r8), intent(in)    :: dnom(pcols,pver)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))
    real(r8), intent(in)    :: cd_top(pcols)          ! cc_top * ubc value

    real(r8), intent(inout) :: q(pcols,pver)          ! Constituent field

    ! --------------- !
    ! Local Variables ! 
    ! --------------- !

    real(r8)                :: zf(pcols,pver)         ! Term in tri-diag solution
    integer                    i, k                   ! Longitude, vertical indices

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !

    ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of 
    ! tridiagonal matrix defined by implicit diffusion equation.
    ! Note that only levels ntop through nbot need be solved for.

    do i = 1, ncol
       zf(i,nbot) = q(i,nbot)*dnom(i,nbot)
    end do

    do k = nbot - 1, ntop + 1, -1
       do i = 1, ncol
          zf(i,k) = (q(i,k) + ca(i,k)*zf(i,k+1))*dnom(i,k)
       end do
    end do

    ! Include boundary condition on top element

    k = ntop
    do i = 1, ncol
       zf(i,k) = (q(i,k) + cd_top(i) + ca(i,k)*zf(i,k+1))*dnom(i,k)
    end do

    ! Perform back substitution

    do i = 1, ncol
       q(i,ntop) = zf(i,ntop)
    end do

    do k = ntop + 1, nbot, +1
       do i = 1, ncol
          q(i,k) = zf(i,k) + ze(i,k)*q(i,k-1)
       end do
    end do

    return
  end subroutine vd_lu_solve

  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  

end module messy_vertdiff_camdiffsolver
