!*****************************************************************************
MODULE messy_vertdiff

  USE messy_main_constants_mem, ONLY: DP, HLINE2, STRLEN_SHORT

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vertdiff_read_nml_ctrl           ! read CTRL namelist and initialize
  PUBLIC :: positive_moisture

  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modstr = 'vertdiff' ! name of module
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver = '1.0' ! module version


  ! GLOBAL CTRL-NAMELIST
  INTEGER, PUBLIC      :: eddy_scheme = 0      ! Holtslag and Boville (default) / Holtslag and Boville and Rash
  LOGICAL, PUBLIC      :: do_tms = .TRUE.      ! switch for turbulent mountain stress
  LOGICAL, PUBLIC      :: do_iss = .TRUE.      ! switch for implicit turbulent surface stress
  ! Module data
  REAL(dp), PUBLIC :: rztodt                   ! 1./ztodt [ 1/s ]
  REAL(dp), PUBLIC, POINTER :: dtk(:,:,:)      ! T tendency from KE dissipation
  REAL(dp), PUBLIC, POINTER :: tke(:,:,:)      ! Turbulent kinetic energy [ m2/s2 ]
  INTEGER,  PUBLIC, POINTER :: turbtype(:,:)   ! Turbulent interface types [ no unit ]
  REAL(dp), PUBLIC, POINTER :: smaw(:,:)       ! Normalized Galperin instability function
                                               ! ( 0<= <=4.964 and 1 at neutral )
  REAL(dp), PUBLIC, pointer :: cgs(:,:,:)      ! Counter-gradient star  [ cg/flux ]
  REAL(dp), PUBLIC, pointer :: cgh(:,:)        ! Counter-gradient term for heat
  REAL(dp), PUBLIC, pointer :: ksrftms(:,:)    ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
  REAL(dp), PUBLIC, pointer :: tautmsx(:,:)    ! U component of turbulent mountain stress [ N/m2 ]
  REAL(dp), PUBLIC, pointer :: tautmsy(:,:)    ! V component of turbulent mountain stress [ N/m2 ]
  REAL(dp), PUBLIC, pointer :: tautotx(:,:)    ! U component of total surface stress [ N/m2 ]
  REAL(dp), PUBLIC, pointer :: tautoty(:,:)    ! V component of total surface stress [ N/m2 ]

  REAL(dp), PUBLIC, POINTER :: kvm_in(:,:,:)   ! kvm from previous timestep [ m2/s ]
  REAL(dp), PUBLIC, POINTER :: kvt(:,:,:)      ! Molecular kinematic conductivity for temperature [  ]
  REAL(dp), PUBLIC, POINTER :: kvq(:,:,:)      ! Eddy diffusivity for constituents [ m2/s ]
  REAL(dp), PUBLIC, POINTER :: kvh(:,:,:)      ! Eddy diffusivity for heat [ m2/s ]
  REAL(dp), PUBLIC, POINTER :: kvm(:,:,:)      ! Eddy diffusivity for momentum [ m2/s ]
  REAL(dp), PUBLIC, POINTER :: sfi(:,:)        ! Saturation fraction at interfaces [ fraction ]
  REAL(dp), PUBLIC, POINTER :: sl(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qt(:,:,:)
  REAL(dp), PUBLIC, POINTER :: slv(:,:,:)
  REAL(dp), PUBLIC, POINTER :: slten(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qtten(:,:,:)
  REAL(dp), PUBLIC, POINTER :: slvten(:,:,:)
  REAL(dp), PUBLIC, POINTER :: slflx(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qtflx(:,:,:)
  REAL(dp), PUBLIC, POINTER :: uflx(:,:,:)
  REAL(dp), PUBLIC, POINTER :: vflx(:,:,:)
  REAL(dp), PUBLIC, POINTER :: slflx_cg(:,:,:)
  REAL(dp), PUBLIC, POINTER :: qtflx_cg(:,:,:)
  REAL(dp), PUBLIC, POINTER :: uflx_cg(:,:,:)
  REAL(dp), PUBLIC, POINTER :: vflx_cg(:,:,:)
  REAL(dp), PUBLIC, POINTER :: th(:,:)         ! Potential temperature
  REAL(dp), PUBLIC, POINTER :: wpert(:,:)      ! Turbulent wind gusts

  REAL(dp), PUBLIC, POINTER :: ri(:,:,:)       ! richardson number (HB output)

  REAL(dp), PUBLIC, POINTER :: thvs(:)         ! Virtual potential temperature at surface
  REAL(dp), PUBLIC, POINTER :: rrho(:)         ! Reciprocal of density at surface
  REAL(dp), PUBLIC, POINTER :: khfs(:)         ! sfc kinematic heat flux [mK/s]
  REAL(dp), PUBLIC, POINTER :: kqfs(:)         ! sfc kinematic water vapor flux [m/s]
  REAL(dp), PUBLIC, POINTER :: kbfs(:)         ! sfc kinematic buoyancy flux [m^2/s^3]

  REAL(dp), PUBLIC, POINTER :: ftem(:,:)       ! Saturation vapor pressure before PBL
  REAL(dp), PUBLIC, POINTER :: tem2(:,:)       ! Saturation specific humidity and RH
  REAL(dp), PUBLIC, POINTER :: t_aftPBL(:,:,:) ! Temperature after PBL diffusion
  REAL(dp), PUBLIC, POINTER :: tten(:,:,:)     ! Temperature tendency by PBL diffusion
  REAL(dp), PUBLIC, POINTER :: rhten(:,:,:)    ! RH tendency by PBL diffusion
  REAL(dp), PUBLIC, POINTER :: tauresx(:)      ! Residual stress to be added in vdiff to correct
  REAL(dp), PUBLIC, POINTER :: tauresy(:)      ! for turb stress mismatch between sfc and atm accumulated.
  REAL(dp), PUBLIC, POINTER :: ipbl(:)
  REAL(dp), PUBLIC, POINTER :: kpblh(:)
  REAL(dp), PUBLIC, POINTER :: wstarPBL(:)
  REAL(dp), PUBLIC, POINTER :: rairi(:,:)      ! interface gas constant needed for compute_vdiff

  REAL(dp), PUBLIC, POINTER :: tpert(:,:)
  REAL(dp), PUBLIC, POINTER :: qpert(:)
  REAL(dp), PUBLIC, POINTER :: pblh(:,:)

  REAL(dp), PUBLIC, POINTER :: tmp1(:)         ! Temporary storage

!!$  REAL(dp), PUBLIC, POINTER :: ustar(:,:)      ! Surface friction velocity [ m/s ] --> zust_2d
  REAL(dp), PUBLIC, POINTER :: obklen(:,:)     ! Obukhov length [ m ]
  INTEGER , PUBLIC          :: ntop            ! Top interface level to which vertical diffusion is applied ( = 1 ).
  INTEGER , PUBLIC          :: nbot            ! Bottom interface level to which vertical diffusion is applied ( = pver ).

  REAL(dp), PUBLIC, POINTER :: t_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: s_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: u_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: v_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: q_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: xl_tmp(:,:,:)
  REAL(dp), PUBLIC, POINTER :: xi_tmp(:,:,:)


CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE vertdiff_read_nml_ctrl(status, iou)

    ! READ VERTDIFF NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! LOGICAL I/O unit

    ! LOCAL
    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status
    CHARACTER(LEN=*), PARAMETER :: substr = 'vertdiff_read_nml_ctrl'

    NAMELIST /CTRL/ eddy_scheme, do_tms, do_iss

    ! INITIALIZE
    status = 1 ! error

    ! INPUT NAMELIST
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! no error

  END SUBROUTINE vertdiff_read_nml_ctrl

  ! --------------------------------------------------------------------------



  !***************************************************************************
  subroutine positive_moisture( cp, xlv, xls, ncol, mkx, dt, qvmin, qlmin, qimin, & 
       zdp, qv, ql, qi, t, s, qvten, qlten, qiten, sten )
    ! ------------------------------------------------------------------------------- !
    ! If any 'ql < qlmin, qi < qimin, qv < qvmin' are developed in any layer,         !
    ! force them to be larger than minimum value by (1) condensating water vapor      !
    ! into liquid or ice, and (2) by transporting water vapor from the very lower     !
    ! layer. '2._dp' is multiplied to the minimum values for safety.                  !
    ! Update final state variables and tendencies associated with this correction.    !
    ! If any condensation happens, update (s,t) too.                                  !
    ! Note that (qv,ql,qi,t,s) are final state variables after applying corresponding !
    ! input tendencies.                                                               !
    ! Be careful the order of k : '1': near-surface layer, 'mkx' : top layer          ! 
    ! ------------------------------------------------------------------------------- !
    implicit none
    INTEGER,  intent(in)     :: ncol, mkx
    REAL(dp), intent(in)     :: cp, xlv, xls
    REAL(dp), intent(in)     :: dt, qvmin, qlmin, qimin
    REAL(dp), intent(in)     :: zdp(ncol,mkx)
    REAL(dp), intent(inout)  :: qv(ncol,mkx), ql(ncol,mkx), qi(ncol,mkx), t(ncol,mkx), s(ncol,mkx)
    REAL(dp), intent(inout)  :: qvten(ncol,mkx), qlten(ncol,mkx), qiten(ncol,mkx), sten(ncol,mkx)
    INTEGER ::   i, k
    REAL(dp) :: dql, dqi, dqv, sum, aa, dum 

    ! Modification : I should check whether this is exactly same as the one used in
    !                shallow convection and cloud macrophysics.

    do i = 1, ncol
       do k = mkx, 1, -1    ! From the top to the 1st (lowest) layer from the surface
          dql        = max(0._dp,1._dp*qlmin-ql(i,k))
          dqi        = max(0._dp,1._dp*qimin-qi(i,k))
          qlten(i,k) = qlten(i,k) +  dql/dt
          qiten(i,k) = qiten(i,k) +  dqi/dt
          qvten(i,k) = qvten(i,k) - (dql+dqi)/dt
          sten(i,k)  = sten(i,k)  + xlv * (dql/dt) + xls * (dqi/dt)
          ql(i,k)    = ql(i,k) +  dql
          qi(i,k)    = qi(i,k) +  dqi
          qv(i,k)    = qv(i,k) -  dql - dqi
          s(i,k)     = s(i,k)  +  xlv * dql + xls * dqi
          t(i,k)     = t(i,k)  + (xlv * dql + xls * dqi)/cp
          dqv        = max(0._dp,1._dp*qvmin-qv(i,k))
          qvten(i,k) = qvten(i,k) + dqv/dt
          qv(i,k)    = qv(i,k)    + dqv
          if( k .ne. 1 ) then 
             qv(i,k-1)    = qv(i,k-1)    - dqv*zdp(i,k)/zdp(i,k-1)
             qvten(i,k-1) = qvten(i,k-1) - dqv*zdp(i,k)/zdp(i,k-1)/dt
          endif
          qv(i,k) = max(qv(i,k),qvmin)
          ql(i,k) = max(ql(i,k),qlmin)
          qi(i,k) = max(qi(i,k),qimin)
       end do
       ! Extra moisture used to satisfy 'qv(i,1)=qvmin' is proportionally 
       ! extracted from all the layers that has 'qv > 2*qvmin'. This fully
       ! preserves column moisture. 
       if( dqv .gt. 1.e-20_dp ) then
          sum = 0._dp
          do k = 1, mkx
             if( qv(i,k) .gt. 2._dp*qvmin ) sum = sum + qv(i,k)*zdp(i,k)
          enddo
          aa = dqv*zdp(i,1)/max(1.e-20_dp,sum)
          if( aa .lt. 0.5_dp ) then
             do k = 1, mkx
                if( qv(i,k) .gt. 2._dp*qvmin ) then
                   dum        = aa*qv(i,k)
                   qv(i,k)    = qv(i,k) - dum
                   qvten(i,k) = qvten(i,k) - dum/dt
                endif
             enddo
          else 
             write(*,*) 'Full positive_moisture is impossible in vertical_diffusion'
          endif
       endif
    end do
    return

  end subroutine positive_moisture
  !***************************************************************************

  !*****************************************************************************
END MODULE messy_vertdiff
