MODULE MESSY_cvtrans

!-----------------------------------------------------------------------
!
! Convective transport of trace species
!
! Mass-conserving, monotonic version that allows the possibility for 
! either bulk ensemble transport (in a single central core updraft) or 
! segregated transport in each discrete plume which comprises the ensemble
!
! The code expects all arrays to be passed in with dimensions:
! nlev - levels,
! nlond - longitudes, of which only the first nlon are considered,
! ncnst - tracers (all of which will be transported)
!
! Notes for PORTING:
! 1) closure of the mass fluxes is done internally, and the user is 
!    flagged if the departure from closed mass fluxes exceeds a threshold value
! 2) the downdraft parameters in Zhang-McFarlane (ZM) convection are defined
!    in a somewhat unusual fashion, the closure and conversion to column
!    values will probably need to be adapted to other schemes
!
! "Tunable" parameters include: fle, ???
!
!-------------------------Code History----------------------------------
!
! Original version:        M. Lawrence, Dec 2003
! build in for EMAC:       H. Tost, Dec 2003
! alternative algorithm:   H.G. Ouwersloot, Aug 2014

  USE messy_main_constants_mem, ONLY: DP, i4

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_cvtrans, closure, cv_transport, lcvt_gp, lcvt_lg 
  PUBLIC :: tr_attr, tracer_attr
  ! mz_ho_20140826+
  PUBLIC :: lcvt_sc, trans_fac, hlimit
  PUBLIC :: lmeanconc
  PUBLIC :: lintsteps, maxfrac
  ! mz_ho_20140826-
  PUBLIC :: l_calc_trma, calc_trma
  
! um_ak_20090721+: moved to _si ...
!!$  REAL(dp) , PUBLIC, POINTER :: cumassf(:,:,:)     => NULL()  ! convective upward mass flux after closure
!!$  REAL(dp) , PUBLIC, POINTER :: cdmassf(:,:,:)     => NULL()  ! convective downward mass flux after closure
!!$  REAL(dp) , PUBLIC, POINTER :: cuentr(:,:,:)      => NULL()  ! entrainment of the upward flux after closure
!!$  REAL(dp) , PUBLIC, POINTER :: cudetr(:,:,:)      => NULL()  ! detrainment of the upward flux after closure
!!$  REAL(dp) , PUBLIC, POINTER :: cdentr(:,:,:)      => NULL()  ! entrainment of the downward flux after closure
!!$  REAL(dp) , PUBLIC, POINTER :: cddetr(:,:,:)      => NULL()  ! detrainment of the downward flux after closure
!!$
!!$  REAL(dp) , PUBLIC, POINTER :: updr_velo(:,:,:)   => NULL()  ! corrected updraft velocity after closure
!!$
!!$  REAL(dp) , PUBLIC, POINTER :: kbot(:,:)          => NULL()  ! lowest level of convection in cvtrans
!!$  REAL(dp) , PUBLIC, POINTER :: ktop(:,:)          => NULL()  ! top level of convection in cvtrans
!!$  REAL(dp) , PUBLIC, POINTER :: kraintop(:,:)      => NULL()  ! top level of convective rain in cvtrans
!!$  
!!$  REAL(dp) , PUBLIC, POINTER :: trac_field(:,:,:,:) => NULL()  ! tracer field for coupling to scav
!!$  
! um_ak_20090721-
  CHARACTER(len=*), PUBLIC, PARAMETER :: MODSTR = 'cvtrans'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver = '2.4'

  LOGICAL,     PUBLIC :: segtrans, bulktrans
  INTEGER(I4), PUBLIC :: scav_trans
  TYPE tracer_attr
    INTEGER, DIMENSION(:), POINTER :: amount       => NULL()
  END TYPE tracer_attr
  TYPE(tracer_attr),       POINTER, DIMENSION(:)           :: &
                                      tr_attr      => NULL()

  LOGICAL  :: lcvt_gp = .FALSE.     ! cvtrans for gridpoint tracers
  LOGICAL  :: lcvt_lg = .FALSE.     ! cvtrans for (pseudo)-lagrange tracers

  ! mz_ho_20140826+
  ! Switch whether subcloud processes affect convective transport 
  ! (Ouwersloot et al., 2013); can be set by namelist
  LOGICAL             :: lcvt_sc = .FALSE.
  ! Factor for transported species, equal to -xi_2 in Eq. (13) of 
  ! Ouwersloot et al. (2013); can be set by namelist
  REAL(dp)            :: trans_fac = 1.23_dp
  ! Height (m) below which a cloud base is always assumed to be coupled to a
  ! (convective) boundary layer; can be set by namelist
  REAL(dp)            :: hlimit = 2500._dp
  ! During every time step use the temporal mean concentration at cloud
  ! base to calculate what's entering the plume
  LOGICAL             :: lmeanconc = .FALSE.
  ! Enables (adaptive) intermediate time steps for convective transport
  LOGICAL             :: lintsteps = .FALSE.
  ! Maximum fraction of grid box that is allowed to be emptied by upward 
  ! mass-flux per intermediate time step
  REAL(dp)            :: maxfrac   = 0.1_dp
  ! mz_ho_20140826-
  LOGICAL             :: l_calc_trma = .FALSE.
  LOGICAL             :: calc_trma = .FALSE.

CONTAINS
!===============================================================================

  SUBROUTINE init_cvtrans(iou, fstat)

    ! cvtrans MODULE ROUTINE (CORE)
    !
    ! READ SCAV NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified: Holger Tost, MPICH, 14-11-2003

    IMPLICIT NONE
    
    INTRINSIC :: TRIM
    ! I/O
    INTEGER, INTENT(IN)  :: iou      ! logical I/O unit
    INTEGER, INTENT(out)  :: fstat   ! status
 
    ! (LOCAL) NAMELIST VARIABLES FOR xtsurf MODULE CONTROL

    LOGICAL :: bulk
    LOGICAL :: seg_plume
    INTEGER :: sc_trans
   

    NAMELIST /CTRL/     bulk, seg_plume, sc_trans, lcvt_gp, lcvt_lg &
         , lcvt_sc, trans_fac, hlimit &    ! mz_ho_20140826
         , lmeanconc, lintsteps, maxfrac & ! mz_ho_20140826
         , l_calc_trma

    ! LOCAL
    LOGICAL                           :: lex          ! file exists ?

    ! INITIALIZE GLOBAL CONTROL VARIABLES

    segtrans   = .false.
    bulktrans  = .false.
    scav_trans = 1
    fstat = 0

    ! INITIALIZE NAMELIST VARIABLES
    
    bulk      = .false.
    seg_plume = .false.
    sc_trans  = 1

    ! INPUT NAMELIST
    WRITE(*,*) '*******************************************************'
    WRITE(*,*) 'START CVTRANS MODULE INITIALISATION (init_CVTRANS)'
    WRITE(*,*) '*******************************************************'
    ! CHECK IF FILE EXISTS, YES: SWITCH cvtrans ON
    !                       NO : KEEP cvtrans SWITCHED OFF
    INQUIRE(file=TRIM(modstr)//'.nml', exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) 'WARNING *** FILE '//TRIM(modstr)//'.nml'//'  NOT FOUND !'
       WRITE(*,*) ' CVTRANS not working properly !'
       WRITE(*,*) '******************************************************'
       RETURN
    END IF
    
    ! READ NAMELIST
    OPEN(iou,file=TRIM(modstr)//'.nml')
    WRITE(*,*) 'Reading namelist from '//TRIM(modstr)//'.nml', &
         ' (unit ',iou,') ...'
    READ(iou, NML=CTRL, IOSTAT=fstat)
    IF (fstat /= 0) THEN
       WRITE(*,*) 'ERROR *** READ ERROR in NAMELIST ', &
            TRIM(modstr)//'.nml'//' !'
       WRITE(*,*) '******************************************************'
       RETURN
    END IF
    CLOSE(iou)

    segtrans   = seg_plume
    bulktrans  = bulk
    scav_trans = sc_trans
    if(.not. lcvt_sc) trans_fac = 1.0_dp ! mz_ho_20140826

    if (.not.bulk) then
       fstat=1
       WRITE(*,*) ' ONLY BULK transport available by now!!!!'
       RETURN
    ENDIF
    WRITE(*,*) '******************************************************'
    WRITE(*,*) 'END CVTRANS MODULE INITIALISATION (init_cvtrans)'
    WRITE(*,*) '******************************************************'

  END SUBROUTINE init_cvtrans

!==============================================================================

  subroutine closure(nlev, nlond,                            &
!                    npmax, nlon,                            &
                     dt, kgm2,                               &
!                    zdel, zm,                               &
                     muin, euin, duin, mdin, edin,           &
!                    segtrans,                               & 
                     conv_part,                              &
                     jrow, kb, kt, lwork,                    &
                     mu, md, eu, du, ed, dd,                 & ! um_ak_20090721
                     status)
!-----------------------------------------------------------------------
      implicit none
!--------------------------Commons--------------------------------------
!      include 'pmgrid_sa.h'
!-----------------------------Arguments---------------------------------
      INTRINSIC :: ABS, MIN, TINY
! 
! Input
!
      integer                 &
          nlev                & ! vertical dimension (levels) of arrays
!!          ,npmax              & ! maximum number of discrete plumes
!!                                ! (equal to the max # of detraining levels)
!!          ,nlon               & ! number of columns to consider in loops
          ,nlond              & ! horizontal dimension (longitudes) of arrays
          ,jrow               & ! index of blocks and rows  
          ,lwork(nlond)         ! index field of packed columns

      real(dp) :: dt               ! Time step (*real* time step, not leapfrog!)
      real(dp) :: kgm2(nlond,nlev) ! Column density of air in each cell (kg/m2)
      real(dp) :: muin(nlond,nlev) ! Updraft mass flux, top of layer (kg/m2/s)
      real(dp) :: euin(nlond,nlev) ! Entrainment into updraft (kg/m2/s)
      real(dp) :: duin(nlond,nlev) ! Detrainment from updraft (kg/m2/s)
      real(dp) :: mdin(nlond,nlev) ! Downdraft mass flux, top of layer (kg/m2/s)
      real(dp) :: edin(nlond,nlev) ! Entrainment into downdraft (kg/m2/s)
      REAL(dp) :: conv_part(nlond) ! for convective cloud cover estimate

!!      real(dp) :: zdel(nlond,nlev)   ! Depth of each cell (m)
!!      real(dp) :: zm(nlond,nlev)    ! height above sfc at layer midpoints
   

! Output
      INTEGER, INTENT(INOUT) :: status

! Local Variables

      logical :: segtrans       ! T => transport by segregated plumes,
!                               ! F => bulk transport (one central plume)
   
      integer :: i, k, j        ! Work index

! closed mass fluxes and entrainment/detrainment rates for internal use

      real(dp) :: mu(nlond,nlev)    ! Updraft mass flux, top of layer (kg/m2/s)
      real(dp) :: eu(nlond,nlev)    ! Entrainment into updraft (kg/m2/s)
      real(dp) :: du(nlond,nlev)    ! Detrainment from updraft (kg/m2/s)
      real(dp) :: md(nlond,nlev)    ! Downdraft mass flux, top of layer (kg/m2/s)
      real(dp) :: ed(nlond,nlev)    ! Entrainment into downdraft (kg/m2/s)
      real(dp) :: dd(nlond,nlev)    ! Dntrainment from downdraft (kg/m2/s)


!     segregated plume ensemble member mass fluxes
!     these are all in kg/m2/s 
!!$      real(dp) ::                     &  
!!$          pmu(nlev,npmax)      & ! mass flux in each plume
!!$          ,peu(nlev,npmax)     & ! entrainment in each plume
!!$          ,pdu(nlev,npmax)       ! detrainment from each plume (top level only)

      real(dp), PARAMETER :: minflux = 1.e-7_dp   ! Threshold for mass fluxes
                                     ! typical fluxes are ~1.e-3 - 1.e-2 kg/m2/s

      ! correction factor in case of CFL violation
      REAL(dp) :: corr_fac(nlond,nlev), correc_fac(nlond),r  
      REAL(dp) :: du_new(nlond, nlev)

      integer  :: kt(nlond), kb(nlond), kdd(nlond)

!!      integer  :: np              ! number of discrete (segregated) plumes  
!!      real(dp) :: zdelcol(nlev)
!!      real(dp) :: zmcol(nlev)
   

!!$      real(dp) ::                   &
!!$          pmut(nlev)         & ! total mass flux (sum of plumes)
!!$          ,peut(nlev)        & ! total entrainment flux (sum of plumes)
!!$          ,pdut(nlev)          ! total detrainment flux (sum of plumes)
!-----------------------------------------------------------------------     

      status = 200

      segtrans = .false.         ! bulk transport (one central plume)


!     check for sufficient convection in the column; if none, skip the column
    
      do k=1,nlev
        do i=1,nlond
          if (muin(i,k) < minflux) muin(i,k) = 0.
          if (abs(mdin(i,k)) < minflux) mdin(i,k) = 0.
        enddo
      enddo


!     close the incoming mass fluxes and entrainment/detrainment rates
!     (this ensures mass conservation and column-wise monotonicity)
!     if a small tolerance (minflux) is exceeded, stop and flag user that the 
!     incoming fluxes aren't sufficiently balanced;

      do i=1,nlond
        if (muin(i,1) > 0.0_dp .or. muin(i,1) < 0.0_dp) RETURN
      enddo

!mz_ht_20040106+
!        Incoming mass flux is assumed to be almost correct         
      mu(:,:)=0._dp
      eu(:,:)=0._dp
      du(:,:)=0._dp
      md(:,:)=0._dp
      ed(:,:)=0._dp
      dd(:,:)=0._dp
      corr_fac(:,1:nlev)=1._dp
   
      do k=1,nlev
        do i=1,nlond
!!    limiting mass flux not only by means of mass conservation, but also 
!!    measured values to 100 g/(m^2s)     
!        mu(k)=min(muin(i,k), 0.10_dp/conv_part(i))      
!        mu(k)=min(muin(i,k), 0.03_dp)
          mu(i,k)=muin(i,k)
          du(i,k) = max(duin(i,k), 0.0_dp)
          du_new(i,k) = du(i,k)
        enddo
      enddo
      do k=1,nlev-1
        do i=1,nlond
          du_new(i,k) = min(du_new(i,k), mu(i,k+1))
          ! mz_ho_20140826+
!!$          if (mu(i,k) < minflux) &
!!$            du_new(i,k) = min (mu(i,k+1), (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp))
          if (lintsteps) then 
             ! with this switch mu is allowed to be bigger than kgm2/dt
             if (mu(i,k) < minflux) du_new(i,k) = mu(i,k+1)
          else
             if (mu(i,k) < minflux) &
                  du_new(i,k) = min (mu(i,k+1), &
                  (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp))
          endif
          ! mz_ho_20140826-
          if (mu(i,k+1) < minflux) du_new(i,k) = 0.
          eu(i,k) = mu(i,k) - mu(i,k+1) + du_new(i,k)
          if (eu(i,k) < 0._dp) then
            du_new(i,k) = du_new(i,k) - eu(i,k)
            eu(i,k) = 0.
          endif
        enddo
      enddo
      eu(1:nlond,nlev)= mu(1:nlond,nlev) - du(1:nlond,nlev)

      ! mz_ho_20140826+
      ifnotlintsteps: if (.not. lintsteps) then
         ! In that case, the original CFL criteria need to be applied
      ! mz_ho_20140826-
!     check CFL
      do k=1,nlev
        do i=1,nlond
          if (mu(i,k)*dt > kgm2(i,k)) then
           !    write(*,*) 'mu*dt exceeds cell mass'
           !    write(*,*) i,k,mu(k)*dt, kgm2(i,k)   &
           !         ,mu(k)*dt/kgm2(i,k)-1.
!     if the exceedence isn't large, then it's probably machine prec, "correct"
           !if (abs(1.-mu(k)*dt/kgm2(i,k)).lt.minflux) then
           !      write(*,*) 'correcting', mu(k)
            mu(i,k) = (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp) 
           !      write(*,*) 'after correcting', mu(k)
!     calculate a correction factor of the massflux which is then applied 
!     also to entrainment/detrainment rates
            corr_fac(i,k) = mu(i,k) / muin(i,k)
          end if

          if ((du_new(i,k)*dt > kgm2(i,k)) .AND. (kgm2(i,k) > 0._dp) ) then 
            du_new(i,k) = (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp) 
            corr_fac(i,k) = min(corr_fac(i,k), du(i,k)/du_new(i,k))
          endif
        enddo
      enddo
      do i=1,nlond
        correc_fac(i) = minval(corr_fac(i,1:nlev))
      enddo

      do k=1,nlev
        do i=1,nlond
          mu(i,k) = correc_fac(i) * muin(i,k)
          du(i,k) = max(correc_fac(i) * du_new(i,k),0.0_dp)         
        enddo
      enddo

      do k=1,nlev
        do i=1,nlond
          if ( mu(i,k) > 0._dp .and. (mu(i,k)+du(i,k))*dt >  kgm2(i,k)) then
            r = du(i,k)/mu(i,k)
            mu(i,k) = (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp) / (1 + r)
            corr_fac(i,k) = min(corr_fac(i,k), mu(i,k)/muin(i,k))
          endif
        enddo
      enddo
      ! mz_ho_20140826+
      endif ifnotlintsteps
      ! mz_ho_20140826-
      do i=1,nlond
        correc_fac(i) = minval(corr_fac(i,1:nlev))
      enddo

      do k=1,nlev
        do i=1,nlond
          mu(i,k) = correc_fac(i) * muin(i,k)
          md(i,k) = -correc_fac(i) * mdin(i,k)
          if (md(i,k) > mu(i,k)) md(i,k)=min(md(i,k),0.99*mu(i,k))
          du(i,k) = max(correc_fac(i) * du_new(i,k),0.0_dp)         
          ed(i,k) = min(- correc_fac(i) * edin(i,k),0.0_dp)
          eu(i,k) = correc_fac(i)*eu(i,k)
        enddo
      enddo
      do k=1,nlev-1
        do i=1,nlond
          du_new(i,k) = min(du(i,k), mu(i,k+1))
! mz_ho_20140826+
!!$          if (mu(i,k) < minflux) &
!!$            du_new(i,k) = min (mu(i,k+1), (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp))
          if (lintsteps) then 
             ! with this switch mu is allowed to be bigger than kgm2/dt
             if (mu(i,k) < minflux) du_new(i,k) = mu(i,k+1)
          else
             if (mu(i,k) < minflux) &
                  du_new(i,k) = min (mu(i,k+1), &
                  (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp))
          endif
! mz_ho_20140826-
          if (mu(i,k+1) < minflux) du_new(i,k) = 0.
          eu(i,k) = mu(i,k) - mu(i,k+1) + du_new(i,k)
          if (eu(i,k) < 0._dp) then
            du_new(i,k) = du_new(i,k) - eu(i,k)
            eu(i,k) = 0.
          endif
        enddo
      enddo
      eu(1:nlond,nlev)= mu(1:nlond,nlev) - du(1:nlond,nlev)
! Calculation of cloud top height
      kt(1:nlond) = nlev
 
      do k=nlev-1,1, -1
        do i=1,nlond
          if (mu(i,k+1) > minflux) kt(i) = k
        enddo
      enddo
      
!        do k=1,kt(i)-1
      do k=1,nlev
        do i=1,nlond
          if (k > kt(i)-1) cycle
          du(i,k)=0.
        enddo
      enddo
! Assuming values for entrainment/detrainment at cloud top
      do i=1,nlond
        eu(i,kt(i)) = 0.
      enddo
! Maximizing detrainment with massflux
      do k=1,nlev-1
        do i=1,nlond
          du(i,k) = min ( du_new(i,k), mu(i,k+1) )
        enddo
      enddo
! Calculate entrainment rates assuming detrainment and massflux are correct   
!!$      do k=1,nlev-1
!!$        if (mu(k).lt.minflux) du(k)=min(mu(k+1), (kgm2(i,k)/dt)*(1._dp - 1.e-12_dp))
!!$        if (mu(k+1).lt.minflux) du(k) = 0.
!!$        eu(k) = mu(k) - mu(k+1) + du(k)
!!$        if (eu(k).lt.0._dp) then
!!$          du(k) = du(k) - eu(k)
!!$          eu(k) = 0.
!!$        endif
!!$      enddo
!!$      eu(nlev)= mu(nlev) - du(nlev)
!      du(kt-1) = min(mu(kt), (kgm2(i,kt-1)/dt)*(1._dp - 1.e-12_dp))         

!     also close the downdraft mass fluxes in the same way as the updrafts;
!     in the case of ZM convection, only the downdraft entrainment 
!     is specified, for convenience also specify the detrainment locally;
!     also, note that the downdraft mass fluxes for ZM are negative,
!     while the entrainment rates are positive, except sometimes they are 
!     negative (but not always), in which case they represent the detrainment;
!     make both of these positive for use in the transport algorithm

!     for ease of formulation, TRANSFORM DOWNDRAFT MASS FLUXES AND ENTRAINMENT
!     TO POSITIVE VALUES, and determine a downdraft detrainment profile

!     for the ZM scheme, ed *sometimes* also appears to include dd!  
!     whenever ed<0, it is detrainment; ignore these terms and diagnose
!     downdraft detrainment from the mass flux

      kdd(1:nlond) = 1
      do k = nlev,1,-1
        do i=1,nlond
          if (md(i,k) > minflux) kdd(i) = k
          ed(i,k) = -min(0._dp,edin(i,k))  !mz_ht_20040106
        end do
      enddo
      do k=1,nlev-1
         ed(1:nlond,k) = min( ed(1:nlond,k), md(1:nlond,k+1))
      enddo

      do i=1,nlond
        if  (kdd(i) /= 1) ed(i,kdd(i)-1) = md(i,kdd(i))
        ed(i,nlev) = 0.
      enddo
          
      do k = 1,nlev-1
        do i=1,nlond
          dd(i,k) = md(i,k) + ed(i,k) - md(i,k+1)
          if (dd(i,k) < 0._dp)  then
            ed(i,k) = ed(i,k) - dd(i,k)
            dd(i,k) = 0.
          endif
!         if (md(k).lt.minflux.and.ed(k).eq.dd(k)) then     
          if ( (md(i,k) < minflux) .and. &
            (ABS(ed(i,k)-dd(i,k)) < tiny(0.0_dp)) ) then
            dd(i,k) = 0.
            ed(i,k) = 0.
          endif
        end do
      enddo
        
      dd(1:nlond,nlev) = md(1:nlond,nlev)


!     determine whether to transport in segregated plumes (segtrans=T) or bulk

      if (segtrans) then
            
!     first segregate the bulk convective mass flux and entrainment profiles
!     into a set of discrete plumes for each column;
!     this is done in a subroutine for development purposes, later for
!     efficiency these operations should be done directly here
!     (note that this is only done for updrafts, downdrafts all detrain at
!     the same level)

!     currently calling with npmax=nlev, could constrain this based on the
!     convection code later

!     create a zdelcol (g77, can use arguments in f90)

!!         do k = 1,nlev
!!            zdelcol(k) = zdel(i,k)
!!            zmcol(k) = zm(i,k)
!!         end do



      else                   ! DO BULK TRANSPORT (SEGPLUMES=F)
            

!     set cloud base and cloud top

         kb(1:nlond) = 1
         do k = 1,nlev          ! could start from ktop, top conv level allowed
           do i=1,nlond
             if (mu(i,k) > 0._dp) kb(i) = k
             if (md(i,k) > 0._dp) kb(i) = k 
           end do
         enddo
! um_ak_20090721+: moved to SMIL
!!$         do k=1,nlev
!!$           do i=1,nlond
!!$             j = lwork(i)
!!$             cumassf(j,k,jrow) = mu(i,k)
!!$             cdmassf(j,k,jrow) = md(i,k)
!!$             cuentr (j,k,jrow) = eu(i,k)
!!$             cudetr (j,k,jrow) = du(i,k)
!!$             cdentr (j,k,jrow) = ed(i,k)
!!$             cddetr (j,k,jrow) = dd(i,k)
!!$           enddo
!!$         enddo
! um_ak_20090721-
       end if                 ! on segtrans

       status = 0

     END subroutine closure

!===============================================================================

    SUBROUTINE cv_transport(nlev, nlon, ncnst,                       &
!                           ntrans, npmax,                           &
                            dt, kgm2,                                &
!                           zdel,zm, segtrans,                       & 
                            fracis, q, jrow, kb, kt, status,  &
                            lwork, option, js,                       &
                            mu, md, eu, du, ed, dd & ! um_ak_20090721
                            , belowh               & ! mz_ho_20140826
                            )

!-----------------------------------------------------------------------
      implicit none
!--------------------------Commons--------------------------------------
!      include 'pmgrid_sa.h'
!-----------------------------Arguments---------------------------------
      INTRINSIC :: ABS, MAX, MIN
! 
! Input
!
      integer  ::             &
          nlev                & ! vertical dimension (levels) of arrays
!          ,npmax              & ! maximum number of discrete plumes
!                                ! (equal to the max # of detraining levels)
          ,nlon               & ! number of columns to consider in loops
          ,ncnst              & ! number of tracers (dimension)
!          ,ntrans             & ! number of tracers to transport (start w/ #1)
          ,jrow               & ! index of blocks and rows  
          ,lwork(nlon)          ! index field of packed columns

      real(dp) :: dt                 ! Time step (*real* time step, not leapfrog!)
      real(dp) :: kgm2(nlon,nlev)   ! Column density of air in each cell (kg/m2)

!      real(dp) :: zdel(nlon,nlev)   ! Depth of each cell (m)
!      real(dp) :: zm(nlon,nlev)    ! height above sfc at layer midpoints
!      real(dp) :: zdelcol(nlev)
!      real(dp) :: zmcol(nlev)
!      logical  :: docol
     
      real(dp) :: fracis(nlon,nlev,ncnst) ! fraction of tracer that is insoluble

      LOGICAL, INTENT(IN), OPTIONAL :: belowh(nlon) ! mz_ho_20140826

! input/output

      real(dp) :: q(nlon,nlev,ncnst)  ! Tracer array 
      INTEGER, INTENT(INOUT) :: status
      INTEGER, INTENT(IN)    :: option, js
! Local Variables

      integer  :: i, j, k, m       ! Work index
      LOGICAL  :: bh(nlon) ! mz_ho_20140826

! closed mass fluxes and entrainment/detrainment rates for a column

      real(dp) :: mu(nlon,nlev)    ! Updraft mass flux, top of layer (kg/m2/s)
      real(dp) :: eu(nlon,nlev)    ! Entrainment into updraft (kg/m2/s)
      real(dp) :: du(nlon,nlev)    ! Detrainment from updraft (kg/m2/s)
      real(dp) :: md(nlon,nlev)    ! Downdraft mass flux, top of layer (kg/m2/s)
      real(dp) :: ed(nlon,nlev)    ! Entrainment into downdraft (kg/m2/s)
      real(dp) :: dd(nlon,nlev)    ! Dntrainment from downdraft (kg/m2/s)

!     segregated plume ensemble member mass fluxes
!     these are all in kg/m2/s 
!!$      real(dp) ::                     &  
!!$          pmu(nlev,npmax)    & ! mass flux in each plume
!!$          ,peu(nlev,npmax)   & ! entrainment in each plume
!!$          ,pdu(nlev,npmax)     ! detrainment from each plume (top level only)

      real(dp) :: const(nlon,nlev)    ! Gathered tracer array 
      real(dp) :: fism(nlon,nlev)     ! insoluble fraction of tracer m

      real(dp), parameter :: small   = 1.0e-15_dp   ! A small number

      integer  ::  kt(nlon), kb(nlon)

      real(dp) :: qu(nlon,nlev),qud(nlon,nlev)
      real(dp) :: qd(nlon,nlev),qdd(nlon,nlev)
      real(dp) :: fle
      real(dp) :: new_const(nlon,nlev)  ! Gathered tracer array, after transport
      ! mz_ho_20140826+
      real(dp) :: qave
      real(dp) :: fact
      ! mz_ho_20140826-

      INTEGER, PARAMETER :: min_lev = 2

!c++debug

!!$      real(dp) ::                   &
!!$          pmut(nlev)         & ! total mass flux (sum of plumes)
!!$          ,peut(nlev)        & ! total entrainment flux (sum of plumes)
!!$          ,pdut(nlev)          ! total detrainment flux (sum of plumes)
            
!-----------------------------------------------------------------------
      ! mz_ho_20140826+
      if (present(belowh)) then
        bh(:) = belowh(:)
      else
        bh(:) = .FALSE.
      endif
      ! mz_ho_20140826-
  
      status = 100
      
! um_ak_20090721+: moved to SMIL
!!$      do k=1,nlev
!!$        do i=1,nlon
!!$          j = lwork(i)
!!$          mu(i,k) = cumassf(j,k,jrow) 
!!$          md(i,k) = cdmassf(j,k,jrow) 
!!$          eu(i,k) = cuentr (j,k,jrow) 
!!$          du(i,k) = cudetr (j,k,jrow) 
!!$          ed(i,k) = cdentr (j,k,jrow) 
!!$          dd(i,k) = cddetr (j,k,jrow) 
!!$        enddo
!!$      enddo
! um_ak_20090721-

! Loop ever each constituent
      do m = 1, ncnst

! Gather up the constituent
        do k = 1,nlev
          do i=1,nlon
            const(i,k) = q(i,k,m)
            new_const(i,k)=q(i,k,m)
!!!$       if (const(i,k) < -1.e-20_dp .or. const(i,k) > 1._dp)  &
!!!$       print*, "bad input parameter, negative at input", q(i,k,m),i,k,jrow,m

  
            fism(i,k) = fracis(i,k,m)
! not set up to deal with soluble gases yet!
            if (abs(fism(i,k)-1._dp) > small) status=200
          end do
        enddo
        if (status == 200) RETURN

!     set updraft mixing ratio at cloud base layer
        
        do i=1,nlon
          qu(i,:)     = 0._dp
          ! mz_ho_20140826+
!!$       qu(i,kb(i)) = const(i,kb(i))
          if (lmeanconc .and. (mu(i,kb(i)) .gt. 0._dp) ) then
             if (bh(i)) then
                fact = trans_fac * mu(i,kb(i)) * dt / kgm2(i,kb(i))
             else  !below h
                fact = mu(i,kb(i)) * dt / kgm2(i,kb(i))
             endif !below h
             qave = const(i,kb(i)-1) + (const(i,kb(i)) &
                  - const(i,kb(i)-1)) * (1.0_dp-exp(-fact))/fact
             ! from eq. d(kgmt * q(kb))/dt = 
             !   - mu (q(kb) + (trans_fac-1) (qu(kb)-qu(kb-1))) + mu q(kb-1)
             qu(i,kb(i)) = qave
             if (bh(i)) qu(i,kb(i)) = qu(i,kb(i)) + &
                  (trans_fac - 1.0_dp) * (qave - const(i,kb(i)-1))
             !prevent unphysical situations
             if (qu(i,kb(i)) .lt. 0.0_dp) qu(i,kb(i)) = 0._dp
          else    !lmeanconc & mu > 0
             qu(i,kb(i)) = const(i,kb(i))
             if (bh(i)) qu(i,kb(i)) = qu(i,kb(i)) + &
                  (trans_fac - 1.0_dp) * (const(i,kb(i)) - const(i,kb(i)-1))
             if (qu(i,kb(i)) .lt. 0.0_dp) &
                  qu(i,kb(i)) = 0._dp ! prevent unphysical situations
          endif
          ! mz_ho_20140826-
!     initialize the updraft detraining mixing ratio                 
          qud(i,:) = 0._dp
        enddo
!     determine the updraft and mixing ratios and detrainment mixing ratios

!                do k = kb-1,kt,-1
        do k=nlev-1,min_lev,-1
!                do k = kb-1,kt-1,-1
          do i=1,nlon
            if ( (k > kb(i)-1) .or. (k < kt(i)-1) ) cycle 

!     detrainment mixing ratio is equal to the updraft mixing ratio upon
!     entering the layer, plus the contribution from air which has entrained
!     in this layer and also detrains in this layer (which is not negligible
!     in ZM, which assumes entrainment in the plumes while they are passing
!     through their detrainment layer); this is accounted for with the 
!     parameter "fle", which is the fraction of local entrainment which 
!     also detrains in this layer (only non-zero when du>0 and eu>0)

            if (du(i,k) > 0.0_dp .and. eu(i,k) > 0.0_dp) then

!     extreme assumption: all the detraining air comes from below and all
!     the air entraining in this layer heads out the top to the next layer:

!                     fle = 0.0

!     other extreme: all the air which entrains in this layer directly 
!     detrains in this layer, so that the plume passing through does not
!     take on any of the tracer from this layer:

!                     fle = 1.0


!     based on the output from segplumes, a "typical" value of fle can 
!     be estimated; from a 10-day global run, this turned out to average
!     fle = 0.45 +/- 0.21 (mean +/- 1 sigma), close to the average for a 
!     single column over Africa for one month of fle = 0.49 +/- 0.19,
!     so use the value of 0.5 as a standard for this code

              fle = 0.5_dp

! there are some limitations on fle (max and min):

! 1) it cannot be so large that the amount entrained locally which should 
!    detrain locally exceeds the amount which actually does detrain locally:

              fle = min(fle,du(i,k)/eu(i,k)-small) 
                        ! fle = min(fle,du(k)/eu(k)) 

! 2) it cannot be so small that the amount which detrains locally exceeds the
!    sum of the amount which enters and the amount which is entrained locally 
!    and detrains locally (k+1 is okay since loop starts at kb-1):

              if (du(i,k) > mu(i,k+1)) then
                fle = max(fle,((du(i,k)-mu(i,k+1))/eu(i,k))*(1._dp+small))
#if !(defined(__SX__))
                write(*,*) 'du(k) > mu(k+1), adjusting fle' 
                write(*,*) i,jrow,k,fle,du(i,k),mu(i,k+1), &
                  eu(i,k),(du(i,k)-mu(i,k+1))/eu(i,k), kb(i),kt
#endif
              end if

            else
! set to zero where there is no detrainment so the equations below can be
! formulated more generally
              fle = 0._dp
            end if

            if ( du(i,k) < TINY(0.0_dp) ) then
               qud(i,k) = 0.0_dp
            ELSE
               qud(i,k) = ( (du(i,k)-fle*eu(i,k))*qu(i,k+1)  &
!                      + fle*eu(i,k)*fism(i,k)*const(i,k) ) &
! will have to think through carefully how to do interstitial phase, ignore for now
                    + fle*eu(i,k)*const(i,k) )         &
                       / du(i,k)
            END if
               
!!$                      if (qud(k).lt.0._dp) print*, "QUD",qud(k),&
!!$                        (du(k)-fle*eu(k))*qu(k+1) , &
!!$                        (du(k)-fle*eu(k)),  du(k), fle, &
!!$                        min(0.5_dp,du(k)/eu(k)), spacing(fle), &
!!$                        eu(k), qu(k+1), const(i,k), i,k, option

! update the updraft
            if ( mu(i,k) < tiny(0.0_dp) ) then
               qu(i,k) = 0.0_dp
            else
               qu(i,k) = ( (mu(i,k+1)-(du(i,k)-fle*eu(i,k)))   &
                    *qu(i,k+1)                           &
!                   + (1.-fle)*eu(k)*fism(i,k)*const(i,k) ) &
! will have to think through carefully how to do interstitial phase, 
! ignore for now
                    + (1._dp-fle)*eu(i,k)*const(i,k) )      &
                    / mu(i,k)
            END if
!!$         if (qu(k).lt.0._dp) print*, "QU",qu(k),&
!!$                     (mu(k+1)-(du(k)-fle*eu(k)))*qu(k+1) , &
!!$                      du(k), fle, eu(k), qu(k+1), const(i,k), i,k, kb-1, kt-1
          end do
        enddo





!     now do the same for the downdraft (only entrains)
        do i=1,nlon
          qd(i,:) = 0. 
     !!!     qd(i,kt(i)-1) = 0. 
!     initialize the downdraft detraining mixing ratio
          qdd(i,:) = 0. 
        enddo
!               do k = kt,kb 
        do k = min_lev,nlev
!            do k = kt-1,kb 
          do i=1,nlon
            if ( k < kt(i)-1 .or. k > kb(i) ) cycle
                    

!     detrainment mixing ratio is equal to the downdraft mixing ratio upon
!     entering the layer (assume entrained air from layer does not detrain
!     in the same layer - okay for ZM)

            qdd(i,k) = qd(i,k-1)
               
!     update the downdraft
! mz_pj_20060316 : to 'if' statements instead of one ...
            if (k < nlev) then
               if (md(i,k+1) > 0._dp) then
                  qd(i,k) = ( qd(i,k-1)*(md(i,k)-dd(i,k)) &
                       + const(i,k)*ed(i,k) ) / md(i,k+1)                     
               else
                  qd(i,k) = 0. ! no flux out the bottom of the model
               end if
            else
               qd(i,k) = 0. ! no flux out the bottom of the model
            end if
          end do
       enddo
 
!     now update the environment (includes the redistribution due to 
!     mass balance subsidence)

        do k=nlev,min_lev,-1
!              do k = kb,kt,-1
          do i=1,nlon
!              do k = kb,kt-1,-1
            if ( k > kb(i) .or. k < kt(i)-1 ) cycle

!     updraft only

!     mass leftover from original local mass (after subsidence and detrainment)
!     times local mixing ratio plus subsidence flux plus detrainment flux
!     divided by total local air mass (which should not change!)

!                  newconst(i,k) = 
!     $                 ( (kgm2(i,k) - (mu(k)+du(k))*dt) * const(i,k) 
!     $                 + mu(k)*dt*const(i,k-1)  
!     $                 + du(k)*dt*qud(k) )
!     $                 / kgm2(i,k)

!     including downdraft
!     compute the amount of subsidence across each interface as the 
!     residual of the updraft minus downdraft flux at that interface

!     debug: check that the downdraft flux does not exceed the updraft

!!$                  if (md(k).gt.mu(k)) then

!!$                     write(*,*) 'downdraft > updraft, ',i,k,md(k),mu(k)
!!$                     do kk = 1,nlev
!!$                        write(*,'(i3,12(1pe12.3))') kk &
!!$                            , muin(i,kk),mu(kk),mdin(i,kk),md(kk)
!!$                     end do
!!$                  end if
                     
               
!               newconst(i,k) = 
!     since only considering subsidence, can still overwrite const
!     (as long as loop goes upwards so you don't use updated const(k-1))

!!$                  if (mu(k).gt.md(k)) then

!!$            if (kgm2(i,k) == 0._dp) print*, i,k,kb(i), kt(i), &
!!$              kgm2(i,k), kgm2(i,kb(i)),kgm2(i,kt(i)), kgm2(i,kt(i)-1), &
!!$              "2.Term", (mu(i,k)-md(i,k))*dt*const(i,k-1), &
!!$              "3. Term ",du(i,k)*dt*qud(i,k),&
!!$              "4. Term ",dd(i,k)*dt*qdd(i,k),&
!!$              "teile k:",mu(i,k),md(i,k), du(i,k), dd(i,k), &
!!$              "alles mu:",mu(i,:), &
!!$              "alles kgm2: ", kgm2(i,:)

            new_const(i,k) = ( (kgm2(i,k) - (mu(i,k)-md(i,k)          &
                             + du(i,k) + dd(i,k)) * dt) * const(i,k)  &
                             + (mu(i,k)-md(i,k))*dt*const(i,k-1)      &   ! subsidence
                             + du(i,k)*dt*qud(i,k)                    &   ! detrained from updraft
                             + dd(i,k)*dt*qdd(i,k) )                  &   ! detrained from downdraft
                             / kgm2(i,k)
            ! mz_ho_20140826
            if ((k==kb(i)) .and. (bh(i) .or. lmeanconc)) then
               new_const(i,k) = ( new_const(i,k)*kgm2(i,k) - &
                    mu(i,k)*dt*(qu(i,k)-const(i,k)) ) / kgm2(i,k)
               ! Pumping effect of air leaving with different 
               ! properties below cloud
            endif
            ! mz_ho_20140826

            IF (.NOT. calc_trma) THEN
#if !(defined(__SX__))
! for DEBUG purposes
            if (tr_attr(js)%amount(m) == 1) then
              if ( (new_const(i,k) < -1.e-20_dp .or. new_const(i,k) > 1._dp) &
                .and. (const(i,k) > 0._dp) )  &
                print*, "bad transport", &
                option, i,k,jrow,m,kt(i), &
                new_const(i,k)*kgm2(i,k), kgm2(i,k),   &
                "first",&
                (kgm2(i,k)-(mu(i,k)-md(i,k)+du(i,k)+dd(i,k))*dt)*const(i,k), &
                mu(i,k), md(i,k), du(i,k), dd(i,k), &
                mu(i,k)-md(i,k) + du(i,k) + dd(i,k),&
                (mu(i,k)-md(i,k) + du(i,k) + dd(i,k)) * dt,&
                "second",(mu(i,k)-md(i,k))*dt*const(i,k-1),(mu(i,k)-md(i,k)), &
                const(i,k-1), &
                "third", du(i,k)*dt*qud(i,k), du(i,k), qud(i,k),eu(i,k),      &
                const(i,k),           &
                "fourth", dd(i,k)*dt*qdd(i,k), qdd(i,k), md(i,k), dd(i,k),    &
                ed(i,k), md(i,k-1), dd(i,k-1),ed(i,k-1) 
            else if (tr_attr(js)%amount(m) == 2) then            
              if ((new_const(i,k) < -1.e-10_dp .or. new_const(i,k) > 1.e20_dp) &
                .and. (const(i,k) > 0._dp) ) &
                print*, "bad transport",     &
                option, i,k,jrow,m,kt(i),       &
                new_const(i,k)*kgm2(i,k), kgm2(i,k), &
                "first", &
                (kgm2(i,k)-(mu(i,k)-md(i,k)+du(i,k)+dd(i,k))*dt)*const(i,k), &
                mu(i,k), md(i,k), du(i,k), dd(i,k), &
                mu(i,k)-md(i,k) + du(i,k) + dd(i,k),&
                (mu(i,k)-md(i,k) + du(i,k) + dd(i,k)) * dt,&
                "second",(mu(i,k)-md(i,k))*dt*const(i,k-1),(mu(i,k)-md(i,k)), &
                const(i,k-1), &
                "third",du(i,k)*dt*qud(i,k), du(i,k), qud(i,k),eu(i,k),       &
                const(i,k),           &
                "fourth", dd(i,k)*dt*qdd(i,k), qdd(i,k),md(i,k), dd(i,k),     &
                ed(i,k), md(i,k-1), dd(i,k-1),ed(i,k-1) 
            endif
#endif
            ENDIF
!!$                  else
!!$                     if(kb > nlev) then
!!$                        const(i,k) =                              &
!!$                             ( (kgm2(i,k) - (md(k)-mu(k)          &
!!$                             + du(k) + dd(k)) * dt) * const(i,k)  &
!!$                             + (md(k)-mu(k))*dt*const(i,k+1)      &
!!$                             + du(k)*dt*qud(k)                    &
!!$                             + dd(k)*dt*qdd(k) )                  &
!!$                             / kgm2(i,k)
!!$                     else
!!$                        const(i,k) =                              &
!!$                             ( (kgm2(i,k) - (md(k)-mu(k)          &
!!$                             + du(k) + dd(k)) * dt) * const(i,k)  &
!!$                        !     + (md(k)-mu(k))*dt*const(i,k+1)      &
!!$                             + du(k)*dt*qud(k)                    &
!!$                             + dd(k)*dt*qdd(k) )                  &
!!$                             / kgm2(i,k)
!!$                     endif
!!$                  endif
          end do
        enddo

! Return tracer back to full arrays

! the code sometimes produces *extremely* small negatives (order 1.e-30), 
! can set these to zero without a significant violation of mass conservation
! (note that this also results in a small non-monotonicity, but again, this is
!  near the machine precision level)

        do k = 1,nlev
          do i=1,nlon
            q(i,k,m) = new_const(i,k)
!                  q(i,k,m) = max(0.,const(i,k))
          end do
        enddo

      end do              ! on m
      
      status = 0

      return

    END SUBROUTINE cv_transport

!===============================================================================

END MODULE MESSY_cvtrans
