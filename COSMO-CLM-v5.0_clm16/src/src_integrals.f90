!+ Source module for computing volume- and surface integrals
!------------------------------------------------------------------------------

MODULE src_integrals

!------------------------------------------------------------------------------
!
! Description:
!   Contained are tools to calculate 
!   - volume integrals of arbitrary fields over an arbitrary 
!     cuboid (= a rectangular parallelepiped, a 'Quader' in German) 
!     defined in the numerical (!) (i.e. terrain-following) grid and 
!   - surface integrals of arbitrary vector fields ('fluxes')
!     over the surface of this cuboid  
!   for inspecting conservation properties.
!
! Contained are the following:
!
!   - init_integral_3d
!     Subroutine to perform some initializations
!
!   - integral_3d_total
!     Function to calculate the integral of a field in the total domai
!
!   - integral_3d
!     Function to calculate the integral of a field in the domain of the 
!     processor
!
!   - integral_3d_3tl
!     Function to calculate the integral of a field in the domain of the 
!     processor
!
!   - integral_3d_cond
!     Function to calculate the integral of a field in the subdomain of the 
!     processor, which is described by 'integ_cuboid % qmy'.
!
!   - surface_integral_total
!     Subroutine to calculate the surface integrals over fluxes
!
!   - surface_integral_cond
!     Subroutine to  calculate surface integrals over a processor area
!
!   - organize_integrals
!     the main driver routine (this probably has TO BE CHANGED BY THE USER, see below!)
!
!   - init_check_conservation / final_check_conservation
!     initialization / deallocation  routines
!
!   - check_rho_conservation
!     the main routine (if other variables than total density has to be inspected
!     this has TO BE CHANGED BY THE USER, see below!)
!
!   - adv_flux_upwind / adv_flux_laxwen
!     these routine calculate advective fluxes
!
! Quick users guide:
! ==================
!
!   - make a copy of the example subr. 'check_rho_conservation' and give it
!     a new name (e.g. 'check_xy_conservation')
!   - in Subr. 'organize_integrals' (this is the main driver routine):
!     call this new routine 'check_xy_conservation' additionally to/instead of
!     'check_rho_conservation'
!   - in 'init_check_conservation':
!     - define a new variable like 'rho_xy' and allocate it
!       (analogous to variable 'rho' or 'rho_q_tot')
!     - define additional storage variables like 'integ_old_xy' and 'deltyQ_xy'
!   - the main work has to be done in the new routine 'check_xy_conservation':
!     - calculate the new scalar field 'rho_xy'
!       from which you want to calculate volume integrals
!     - if you have only advective fluxes of this scalar field (through the
!       boundaries) then things are easy. Otherwise you have to add additional
!       flux components (from diffusion, sedimentation, ...) to the fields flux_x,
!       flux_y, flux_z (an example of such a more elaborate subroutine is 'check_qx_conservation'=
!   - in the NAMELIST: set l_integrals=.TRUE., and if needed, define your cubic edges
!   - ... after the model run: have a look in your model ascii-output  :-)
!   - more information can be found in: Baldauf (2008), COSMO-Newsletter Nr. 7
!
! Current Code Owner: DWD, Michael Baldauf
!  phone:  +49 69 8062 2733
!  fax :   +49 69 8062 3721
!  email:  Michael.Baldauf@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V3_23        2007/03/30 Michael Baldauf
!  Initial release
! V4_4         2008/07/16 Ulrich Schaettler
!  Eliminated duplicated variables
! V4_5         2008/09/10 Ronny Petrik, Michael Baldauf
!  Additional flux calculation method (Lax-Wendroff)
! V4_8         2009/02/16 Oliver Fuhrer
!  Use global values only if num_compute greater 1
! V4_9         2009/07/16 Ulrich Schaettler
!  Eliminated SR calc_sqrtg_r and moved it to numeric_utilities
! V4_10        2009/09/11 Michael Baldauf
!  Editorial changes
! V4_13        2010/05/11 Michael Gertz
!  Adaptions to SVN
! V4_18        2011/05/26 Michael Baldauf
!  Introduced new SRs check_qx_conservation and calc_cuboid_geometry
! V4_23        2012/05/10 Ulrich Schaettler
!  Use field sqrtg_r_* from new module grid_metrics_utilities
! V4_25        2012/09/28 Anne Roches, Oliver Fuhrer
!  Replaced qx-variables by using them from the tracer module
! V4_27        2013/03/19 Astrid Kerkweg, Ulrich Schaettler
!  MESSy interface introduced
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE data_parameters , ONLY :   &
  ireals,     & ! KIND-type parameter for "normal" integer variables
  iintegers     ! KIND-type parameters for real variables

USE data_constants,   ONLY:    &
  r_earth, pi, g, r_d, r_v

USE data_fields,      ONLY:    &
  crlat, rho0, hhl, u, v, w, pp, p0, T

USE data_modelconfig, ONLY :   &
  ie,         & ! number of grid points in zonal direction
  je,         & ! number of grid points in meridional direction
  ke,         & ! number of grid points in vertical direction
  istart,     & ! start index for the forecast of w, t, qd, qw and pp
  iend,       & ! end index for the forecast of w, t, qd, qw and pp
  jstart,     & ! start index for the forecast of w, t, qd, qw and pp
  jend,       & ! end index for the forecast of w, t, qd, qw and pp
  dt,         & ! long time-step
  dlon,       & ! grid point distance in zonal  direction (in degrees)
  dlat,       & ! grid point distance in merid. direction (in degrees)
  lalloc_prr_gsp, lalloc_prr_con, lalloc_prs_gsp, lalloc_prs_con, & ! 
  lalloc_prg_gsp,                                                 & !
  idt_qv, idt_qc, idt_qs, idt_qg, idt_qr, idt_qi

USE data_parallel,    ONLY:    &
  my_cart_id,      & ! rank of this subdomain in the cartesian communicator
  isubpos,         & ! positions of the subdomains in the total domain
  nboundlines,     & ! number of boundary lines of the domain for which
                     ! no forecast is computed
  icomm_cart,      & ! communicator for the virtual cartesian topology
  num_compute,     & ! number of compute PEs
  imp_reals          ! determines the correct REAL type used in the model

USE data_runcontrol,    ONLY:  ntstep, nstart, nnow

USE grid_metrics_utilities, ONLY:  sqrtg_r_s, sqrtg_r_u, sqrtg_r_v
USE parallel_utilities, ONLY:  global_values
USE environment,        ONLY:  model_abort

USE src_tracer,         ONLY:  trcr_get, trcr_errorstr

USE data_tracer,        ONLY:  T_ERR_NOTFOUND

!------------------------------------------------------------------------------

IMPLICIT NONE

!------------------------------------------------------------------------------

! Local Variables
! corner coordinates for the cuboid for integration
  INTEGER   (KIND=iintegers)       ::           &
    imin_integ, imax_integ,  &
    jmin_integ, jmax_integ,  &
    kmin_integ, kmax_integ

  LOGICAL :: l_integrals    ! main switch to (de)activate integrals

! Type declarations

  TYPE cuboid_type
    INTEGER (KIND=iintegers) :: imin
    INTEGER (KIND=iintegers) :: imax
    INTEGER (KIND=iintegers) :: jmin
    INTEGER (KIND=iintegers) :: jmax
    INTEGER (KIND=iintegers) :: kmin
    INTEGER (KIND=iintegers) :: kmax
  END TYPE cuboid_type


  TYPE cuboid_multiproz_type
    TYPE (cuboid_type) :: qfd      ! full domain of the cuboid (global indizes)
    TYPE (cuboid_type) :: qmy      ! part of the cuboid, which covers the 
            ! prozessor area; this is of course itself a cuboid (local indizes)
    LOGICAL :: is_in_integ_domain  ! is the prozessor domain part of the cuboid?
      ! are boundaries contained in this prozessor domain? :
    LOGICAL :: has_north_boundary
    LOGICAL :: has_south_boundary
    LOGICAL :: has_west_boundary
    LOGICAL :: has_east_boundary

    REAL (KIND=ireals) :: volume
    REAL (KIND=ireals) :: area_top
    REAL (KIND=ireals) :: area_bottom
    REAL (KIND=ireals) :: area_N
    REAL (KIND=ireals) :: area_S
    REAL (KIND=ireals) :: area_W
    REAL (KIND=ireals) :: area_E

  END TYPE cuboid_multiproz_type

! Variable declarations

  TYPE( cuboid_multiproz_type ) :: integ_cuboid

! And for applying this module
REAL (KIND=ireals), ALLOCATABLE, DIMENSION(:,:,:) :: &
   flux_x, flux_y, flux_z   ! moisture mass flux components

! declare storage variables:

! for the balance of total mass:
REAL (KIND=ireals), ALLOCATABLE :: &
   rho(:,:,:)   ! total density (i.e. dry air + all water constituents)
REAL (KIND=ireals) :: integ_old_rho
REAL (KIND=ireals) :: sum_Res_rho

! for the balance of total water mass:
REAL (KIND=ireals), ALLOCATABLE :: &
   rho_q_tot(:,:,:) ! total density of all water constituents
REAL (KIND=ireals) :: integ_old_q_tot
REAL (KIND=ireals) :: sum_Res_q_tot
 
!------------------------------------------------------------------------------

CONTAINS

!==============================================================================
!==============================================================================

SUBROUTINE init_integral_3d

!------------------------------------------------------------------------------
!
! Description:
!   calculate the struct 'integ_cuboid' from the global indizes
!   imin_integ, ..., kmax_integ of the cuboid 
!   (read as NAMELIST-Parameters in src_meanvalues) 
! 
!   Input:  imin_integ, imax_integ,
!           jmin_integ, jmax_integ,
!           kmin_integ, kmax_integ 
!   Output: integ_cuboid
!
! Method:
!
!------------------------------------------------------------------------------

! global indizes for corners of prozessor domain:
INTEGER(KIND=iintegers) :: my_imin, my_imax, my_jmin, my_jmax  

! global indizes for corners of the
! part of the cuboid lying in the prozessor domain :
INTEGER(KIND=iintegers) :: imin_glob, imax_glob, jmin_glob, jmax_glob 

!------------------------------------------------------------------------------

  ! --- define global corner-points of the cuboid from NAMELIST-variables ---

  integ_cuboid % qfd % imin = imin_integ
  integ_cuboid % qfd % imax = imax_integ

  integ_cuboid % qfd % jmin = jmin_integ
  integ_cuboid % qfd % jmax = jmax_integ

  integ_cuboid % qfd % kmin = kmin_integ
  integ_cuboid % qfd % kmax = kmax_integ

  ! --- global indizes of corners for 'my' prozessor -------

  integ_cuboid % is_in_integ_domain = .TRUE.

  my_imin = isubpos(my_cart_id,1)
  my_jmin = isubpos(my_cart_id,2)
  my_imax = isubpos(my_cart_id,3)
  my_jmax = isubpos(my_cart_id,4)

  imin_glob = MAX( integ_cuboid % qfd % imin, my_imin )
  IF ( imin_glob > my_imax ) THEN
    integ_cuboid % is_in_integ_domain = .FALSE.
  END IF

  imax_glob = MIN( integ_cuboid % qfd % imax, my_imax )
  IF ( imax_glob < my_imin ) THEN
    integ_cuboid % is_in_integ_domain = .FALSE.
  END IF

  jmin_glob = MAX( integ_cuboid % qfd % jmin, my_jmin )
  IF ( jmin_glob > my_jmax ) THEN
    integ_cuboid % is_in_integ_domain = .FALSE.
  END IF

  jmax_glob = MIN( integ_cuboid % qfd % jmax, my_jmax )
  IF ( jmax_glob < my_jmin ) THEN
    integ_cuboid % is_in_integ_domain = .FALSE.
  END IF

  ! -- local indizes for corners of the cuboid at the local prozessor domain 
  !     (are calculated here, even if proc. is not in the cuboid)

  integ_cuboid % qmy % imin = imin_glob - my_imin + istart
  integ_cuboid % qmy % imax = imax_glob - my_imin + istart
  integ_cuboid % qmy % jmin = jmin_glob - my_jmin + jstart
  integ_cuboid % qmy % jmax = jmax_glob - my_jmin + jstart
  integ_cuboid % qmy % kmin = integ_cuboid % qfd % kmin
  integ_cuboid % qmy % kmax = integ_cuboid % qfd % kmax

  ! --- does my processor contain any boundaries of the cuboid? ---

  IF ( integ_cuboid % is_in_integ_domain ) THEN

    IF ( ( my_imin <= integ_cuboid % qfd % imin ) .AND.  &
         (            integ_cuboid % qfd % imin <= my_imax ) )  THEN
      integ_cuboid % has_west_boundary = .TRUE.
    ELSE        
      integ_cuboid % has_west_boundary = .FALSE.
    END IF

    IF ( ( my_imin <= integ_cuboid % qfd % imax ) .AND.  &
         (            integ_cuboid % qfd % imax <= my_imax ) )  THEN
      integ_cuboid % has_east_boundary = .TRUE.
    ELSE        
      integ_cuboid % has_east_boundary = .FALSE.
    END IF

    IF ( ( my_jmin <= integ_cuboid % qfd % jmin ) .AND.   &
         (            integ_cuboid % qfd % jmin <= my_jmax ) ) THEN
      integ_cuboid % has_south_boundary = .TRUE.
    ELSE        
      integ_cuboid % has_south_boundary = .FALSE.
    END IF

    IF ( ( my_jmin <= integ_cuboid % qfd % jmax ) .AND.  &
         (            integ_cuboid % qfd % jmax <= my_jmax ) ) THEN
      integ_cuboid % has_north_boundary = .TRUE.
    ELSE        
      integ_cuboid % has_north_boundary = .FALSE.
    END IF

  END IF

  ! --- control output -------

  !WRITE(*,'(A,I3,L2,6(A,I4))') "integ_cuboid: proz=", my_cart_id, &
  !  &    integ_cuboid%is_in_integ_domain, &
  !  &   ", imin=", integ_cuboid % qmy % imin, &
  !  &   ", imax=", integ_cuboid % qmy % imax, &
  !  &   ", jmin=", integ_cuboid % qmy % jmin, &
  !  &   ", jmax=", integ_cuboid % qmy % jmax, &
  !  &   ", kmin=", integ_cuboid % qmy % kmin, &
  !  &   ", kmax=", integ_cuboid % qmy % kmax

END SUBROUTINE init_integral_3d

!==============================================================================
!==============================================================================

REAL (KIND=ireals)  FUNCTION integral_3d_total( field, integ_type )

!------------------------------------------------------------------------------
!
! Description:
!   calculate the integral of 'field' in the total domain
!   This is the driver routine for Subr. 
!   'integral_3d', 'integral_3d_3tl' or 'integral_3d_cond'.
!
! Method:
!
!------------------------------------------------------------------------------

REAL    (KIND=ireals),    INTENT(in) ::  &
  field(1:ie, 1:je, 1:ke)

INTEGER (KIND=iintegers), INTENT(in) ::  &
  integ_type

! Local Variables
REAL (KIND=ireals)       :: z_int 
CHARACTER (LEN=80)       :: yzerrmsg
INTEGER (KIND=iintegers) :: izerror

!------------------------------------------------------------------------------

  SELECT CASE ( integ_type ) 
  CASE (1)
    ! integral over the whole domain, sqrtg_r_s is already calculated
    z_int = integral_3d     ( field, nboundlines+1, ie-nboundlines,   &
                                     nboundlines+1, je-nboundlines,   &
                                                 3,           ke-2    ) 
  CASE (2)
    ! integral over the whole domain, sqrtg_r_s is 
    !   not available (3-timelevel scheme)
    z_int = integral_3d_3tl ( field, nboundlines+1, ie-nboundlines,   &
                                    nboundlines+1, je-nboundlines,    &
                                                3,           ke-2    ) 
  CASE (3)
    ! integral over a cuboid domain defined by 'integ_cuboid'
    z_int = integral_3d_cond( field, integ_cuboid )
  !CASE (4)
  !  ! integral over a cuboid domain defined by 'integ_cuboid',
  !  ! quadratic interpolation
  !  z_int = integral_3d_quad_int_cond( field, integ_cuboid )
  CASE default
    PRINT*, "ERROR in integral_3d_total: false value in integ_type"
  END SELECT

  IF (num_compute > 1) THEN
    CALL global_values( z_int, 1, 'SUM', imp_reals, icomm_cart ,-1,     &
                        yzerrmsg, izerror )
  ENDIF

  integral_3d_total = z_int

END FUNCTION integral_3d_total

!==============================================================================
!==============================================================================

REAL (KIND=ireals) FUNCTION integral_3d( field, istart, iend, jstart, jend, &
                                         kstart, kend )

!------------------------------------------------------------------------------
!
! Description:
!   calculate the integral of 'field' in the domain of the processor
!   here: sqrt(g) = measure of volume of the grid box
!   is already calculated
!
! Method:
!
!------------------------------------------------------------------------------

REAL (KIND=ireals),       INTENT(in) :: field(:,:,:)
INTEGER (KIND=iintegers), INTENT(in) ::  &
            istart, iend, jstart, jend, kstart, kend
REAL (KIND=ireals), PARAMETER        :: dzeta = 1.0_ireals
REAL (KIND=ireals)                   :: h1, h2
INTEGER (KIND=iintegers)             :: i, j, k

!------------------------------------------------------------------------------

  integral_3d  = 0.0_ireals

  DO j = jstart, jend
    h1 = 0.0_ireals
    DO k = kstart, kend
      h2 = 0.0_ireals
      DO i = istart, iend
        h2 = h2 + field(i,j,k) / sqrtg_r_s(i,j,k)
      END DO
      h1 = h1 + h2
    END DO
    integral_3d = integral_3d  + h1 * crlat(j,1)
  END DO

  integral_3d = integral_3d * r_earth**2 * (pi/180.0_ireals)**2   &
                     * dlon * dlat * dzeta

END FUNCTION integral_3d

!==============================================================================
!==============================================================================

REAL (KIND=ireals)  FUNCTION integral_3d_3tl( field, istart, iend,    &
                                                jstart, jend, kstart, kend )
!------------------------------------------------------------------------------
!
! Description:
!   calculate the integral of 'field' in the domain of the processor
!
! Method:
!
!------------------------------------------------------------------------------

REAL (KIND=ireals),       INTENT(in) :: field(1:ie, 1:je, 1:ke)
INTEGER (KIND=iintegers), INTENT(in) :: &
               istart, iend, jstart, jend, kstart, kend
INTEGER (KIND=iintegers)  :: i,j,k
REAL (KIND=ireals)        :: fdet ! absolute value of functional determinant
REAL (KIND=ireals)        :: dxdlam, dydphi, dzdzeta
REAL (KIND=ireals), PARAMETER ::  dzeta = 1.0_ireals

!------------------------------------------------------------------------------

  integral_3d_3tl  = 0.0_ireals

  dxdlam  = r_earth * pi / 180.0_ireals        ! (lambda in Grad)

  DO j = jstart, jend
    dydphi  = r_earth * pi / 180.0_ireals    ! * cos(phi) vernachlaess.

    DO i = istart, iend

      dzdzeta = 1.0_ireals/g/rho0(i,j,1) * ( p0(i,j,2)-p0(i,j,1) )/dzeta
      fdet = ABS( dxdlam * dydphi * dzdzeta )
      integral_3d_3tl = integral_3d_3tl + fdet * field(i,j,1)

      DO k=kstart, kend-1
        dzdzeta = 0.5_ireals/g/rho0(i,j,k) * ( p0(i,j,k+1)-p0(i,j,k-1) )/dzeta
        fdet = ABS( dxdlam * dydphi * dzdzeta )
        integral_3d_3tl = integral_3d_3tl + fdet * field(i,j,k)
      END DO

      dzdzeta = 1.0_ireals/g/rho0(i,j,kend) * ( p0(i,j,kend)-p0(i,j,kend-1) )/dzeta
      fdet = ABS( dxdlam * dydphi * dzdzeta )
      integral_3d_3tl = integral_3d_3tl + fdet * field(i,j,kend)

    END DO
  END DO

  integral_3d_3tl = integral_3d_3tl * dlon * dlat * dzeta

END FUNCTION integral_3d_3tl

!=============================================================================
!=============================================================================

REAL (KIND=ireals) FUNCTION integral_3d_cond( field, integ_cuboid )

!-----------------------------------------------------------------------------
!
! Description:
!   calculate the integral of 'field' in the subdomain of the processor
!   which is described by 'integ_cuboid % qmy'.
!   It is assumed, that sqrt(g) (= measure of volume of the grid box)
!   is already calculated!
!
! Method:
!
!-----------------------------------------------------------------------------

REAL (KIND=ireals),           INTENT(IN) :: field(1:ie, 1:je, 1:ke)
TYPE (cuboid_multiproz_type), INTENT(IN) :: integ_cuboid

REAL (KIND=ireals), PARAMETER       :: dzeta = 1.0_ireals
REAL (KIND=ireals)                  :: h1, h2
INTEGER (KIND=iintegers)            :: i, j, k

!------------------------------------------------------------------------------

  integral_3d_cond  = 0.0_ireals

  IF ( integ_cuboid % is_in_integ_domain ) THEN

    DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
      h1 = 0.0_ireals
      DO k = integ_cuboid % qmy % kmin, integ_cuboid % qmy % kmax
        h2 = 0.0_ireals
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
          h2 = h2 + field(i,j,k) / sqrtg_r_s(i,j,k)
        END DO
        h1 = h1 + h2
      END DO
      integral_3d_cond = integral_3d_cond  + h1 * crlat(j,1)
    END DO

    integral_3d_cond = integral_3d_cond * r_earth**2 * (pi/180.0_ireals)**2 &
                           * dlon * dlat * dzeta

  END IF

  !WRITE(*,'(A,E21.14)') "integral= ", integral_3d_cond

END FUNCTION integral_3d_cond

!==============================================================================
!==============================================================================

SUBROUTINE surface_integral_total( flux_x, flux_y, flux_z,   &
             int_top, int_bottom, int_n, int_s, int_w, int_e)

!------------------------------------------------------------------------------
!
! Description:
!   calculate the surface integrals over fluxes.
!   The fluxes flux_x, flux_y, flux_z are assumed to be defined at
!   the positions of u, v and w, respectively. 
!   The components are cartesian components (more precisely: 
!   physical components, i.e. for normalized base vectors, 
!   in the earth spherical coordinates).
!
!   This is the driver routine for SUBROUTINE surface_integral_cond
!
! Method:
!
!------------------------------------------------------------------------------

! TYPE( cuboid_multiproz_type ), intent(in) :: integ_cuboid
! fluxes:
REAL (KIND=ireals), INTENT(IN) :: flux_x(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(IN) :: flux_y(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(IN) :: flux_z(1:ie, 1:je, 1:ke+1)

! flux integrals over the 6 faces of the cuboid:
REAL (KIND=ireals), INTENT(OUT) :: &
        int_top, int_bottom, int_n, int_s, int_w, int_e

CHARACTER (LEN=80)  :: yzerrmsg
INTEGER             :: izerror

!------------------------------------------------------------------------------

  CALL surface_integral_cond( flux_x, flux_y, flux_z,   &
             int_top, int_bottom, int_n, int_s, int_w, int_e, integ_cuboid)

  IF (num_compute > 1) THEN
    CALL global_values( int_top,    1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
    CALL global_values( int_bottom, 1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
    CALL global_values( int_n,      1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
    CALL global_values( int_s,      1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
    CALL global_values( int_w,      1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
    CALL global_values( int_e,      1, 'SUM', imp_reals, icomm_cart ,-1,     &
                                              yzerrmsg, izerror )
  ENDIF

END SUBROUTINE surface_integral_total

!==============================================================================
!==============================================================================

SUBROUTINE surface_integral_cond( flux_x, flux_y, flux_z,    &
        int_top, int_bottom, int_n, int_s, int_w, int_e, integ_cuboid)

!------------------------------------------------------------------------------
!
! Description:
!   calculate surface integral over on eprocessor area
! Method:
!
!------------------------------------------------------------------------------

! fluxes:
REAL (KIND=ireals), INTENT(IN) :: flux_x(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(IN) :: flux_y(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(IN) :: flux_z(1:ie, 1:je, 1:ke+1)

! flux integrals over the 6 faces of the cuboid:
REAL (KIND=ireals), INTENT(INOUT) ::   &
        int_top, int_bottom, int_n, int_s, int_w, int_e

TYPE( cuboid_multiproz_type ), INTENT(IN)  :: integ_cuboid

INTEGER :: i, j, k
REAL (KIND=ireals) :: dx_aequ, dx, dy, dzdx, dzdy, h
REAL (KIND=ireals) :: flux
REAL (KIND=ireals) :: da_x, da_y, da_z

!------------------------------------------------------------------------------

  int_top    = 0.0_ireals
  int_bottom = 0.0_ireals
  int_n      = 0.0_ireals
  int_s      = 0.0_ireals
  int_w      = 0.0_ireals
  int_e      = 0.0_ireals

  IF ( integ_cuboid % is_in_integ_domain ) THEN

    IF ( integ_cuboid % has_north_boundary ) THEN
      j = integ_cuboid % qmy % jmax
      DO k = integ_cuboid % qmy % kmin, integ_cuboid % qmy % kmax
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
         !int_n = int_n + flux_y(i,j,k) * ( hhl(i,j,k) - hhl(i,j,k+1) )
          int_n = int_n + flux_y(i,j,k) * 1.0_ireals/sqrtg_r_v(i,j,k)
        END DO
      END DO
      int_n = int_n * r_earth * (pi/180.0_ireals) * dlon * crlat(j, 2)
    END IF

    IF ( integ_cuboid % has_south_boundary ) THEN
      j = integ_cuboid % qmy % jmin
      DO k = integ_cuboid % qmy % kmin, integ_cuboid % qmy % kmax
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
          !int_s = int_s - flux_y(i,j-1,k) * ( hhl(i,j,k) - hhl(i,j,k+1) )
           int_s = int_s - flux_y(i,j-1,k) * 1.0_ireals/sqrtg_r_v(i,j-1,k)
        END DO
      END DO
      int_s = int_s * r_earth * (pi/180.0_ireals) * dlon * crlat(j-1, 2)
    END IF

    IF ( integ_cuboid % has_west_boundary ) THEN
      i = integ_cuboid % qmy % imin
      DO k = integ_cuboid % qmy % kmin, integ_cuboid % qmy % kmax
        DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
          !int_w = int_w - flux_x(i-1,j,k) * ( hhl(i,j,k) - hhl(i,j,k+1) )
           int_w = int_w - flux_x(i-1,j,k) * 1.0_ireals/sqrtg_r_u(i-1,j,k)
        END DO
      END DO
      int_w = int_w * r_earth * (pi/180.0_ireals) * dlat
    END IF

    IF ( integ_cuboid % has_east_boundary ) THEN
      i = integ_cuboid % qmy % imax
      DO k = integ_cuboid % qmy % kmin, integ_cuboid % qmy % kmax
        DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
          !int_e = int_e + flux_x(i,j,k) * ( hhl(i,j,k) - hhl(i,j,k+1) )
           int_e = int_e + flux_x(i,j,k) * 1.0_ireals/sqrtg_r_u(i,j,k)
        END DO
      END DO
      int_e = int_e * r_earth * (pi/180.0_ireals) * dlat
    END IF

    dx_aequ = r_earth * (pi/180.0_ireals) * dlon
    dy      = r_earth * (pi/180.0_ireals) * dlat 

    k = integ_cuboid % qmy % kmin
    IF ( k > 1 ) THEN
      DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
        dx = dx_aequ * crlat(j,1) 
        h  = dx * dy
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
          !dzdx = ( hhl(i+1,j,k) - hhl(i,j,k) ) / dx
          !dzdy = ( hhl(i,j+1,k) - hhl(i,j,k) ) / dy
          ! slope in 4th order accuracy:
          dzdx = ( 2.0_ireals/ 3.0_ireals * ( hhl(i+1,j,k)-hhl(i-1,j,k) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i+2,j,k)-hhl(i-2,j,k) ) ) / dx
          dzdy = ( 2.0_ireals/ 3.0_ireals * ( hhl(i,j+1,k)-hhl(i,j-1,k) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i,j+2,k)-hhl(i,j-2,k) ) ) / dy
          ! normal vector of the surface element (1. order approx.)
          da_x = - dzdx * h
          da_y = - dzdy * h
          da_z = h

          int_top = int_top                                                  &
            &    + 0.25 * ( flux_x(i,j,k  ) + flux_x(i-1,j,k  )            &
            &             + flux_x(i,j,k-1) + flux_x(i-1,j,k-1) ) * da_x   &
            &    + 0.25 * ( flux_y(i,j,k  ) + flux_y(i,j-1,k  )            &
            &             + flux_y(i,j,k-1) + flux_y(i,j-1,k-1) ) * da_y   &
            &    + flux_z(i,j,k) * da_z

        END DO
      END DO
    ELSE
      ! k=1:
      DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
        dx = dx_aequ * crlat(j,1)
        h  = dx * dy
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
          !dzdx = ( hhl(i+1,j,k) - hhl(i,j,k) ) / dx
          !dzdy = ( hhl(i,j+1,k) - hhl(i,j,k) ) / dy
          ! slope in 4th order accuracy:
          dzdx = ( 2.0_ireals/ 3.0_ireals * ( hhl(i+1,j,k)-hhl(i-1,j,k) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i+2,j,k)-hhl(i-2,j,k) ) ) / dx
          dzdy = ( 2.0_ireals/ 3.0_ireals * ( hhl(i,j+1,k)-hhl(i,j-1,k) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i,j+2,k)-hhl(i,j-2,k) ) ) / dy

          ! normal vector of the surface element (1. order approx.)
          da_x = - dzdx * h
          da_y = - dzdy * h
          da_z = h

          int_top = int_top                 &
            + 0.5 * ( flux_x(i,j,k) + flux_x(i-1,j,k) ) * da_x   &
            + 0.5 * ( flux_y(i,j,k) + flux_y(i,j-1,k) ) * da_y   &
               + flux_z(i,j,k) * da_z

        END DO
      END DO
    END IF

    k = integ_cuboid % qmy % kmax
    IF ( k < ke ) THEN
      DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
        dx = dx_aequ * crlat(j,1)
        h  = dx * dy
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax
          !dzdx = ( hhl(i+1,j,k+1) - hhl(i,j,k+1) ) / dx
          !dzdy = ( hhl(i,j+1,k+1) - hhl(i,j,k+1) ) / dy
          ! slope in 4th order accuracy:
          dzdx = ( 2.0_ireals/ 3.0_ireals * ( hhl(i+1,j,k+1)-hhl(i-1,j,k+1) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i+2,j,k+1)-hhl(i-2,j,k+1) ) ) / dx
          dzdy = ( 2.0_ireals/ 3.0_ireals * ( hhl(i,j+1,k+1)-hhl(i,j-1,k+1) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i,j+2,k+1)-hhl(i,j-2,k+1) ) ) / dy

          ! normal vector of the surface element (1. order approx.)
          da_x = dzdx * h
          da_y = dzdy * h
          da_z = -h

          int_bottom = int_bottom               &
            &      + 0.25 * ( flux_x(i,j,k  ) + flux_x(i-1,j,k  )            &
            &               + flux_x(i,j,k+1) + flux_x(i-1,j,k+1) ) * da_x   &
            &      + 0.25 * ( flux_y(i,j,k  ) + flux_y(i,j-1,k  )            &
            &               + flux_y(i,j,k+1) + flux_y(i,j-1,k+1) ) * da_y   &
            &       + flux_z(i,j,k+1) * da_z
        END DO
      END DO
    ELSE
      ! k=ke:
      DO j = integ_cuboid % qmy % jmin, integ_cuboid % qmy % jmax
        dx = dx_aequ * crlat(j,1)
        h  = dx * dy
        DO i = integ_cuboid % qmy % imin, integ_cuboid % qmy % imax

          !dzdx = ( hhl(i+1,j,k+1) - hhl(i,j,k+1) ) / dx
          !dzdy = ( hhl(i,j+1,k+1) - hhl(i,j,k+1) ) / dy
          ! slope in 2nd order accuracy:
          !dzdx = 0.5_ireals * ( hhl(i+1,j,k+1)-hhl(i-1,j,k+1) ) / dx
          !dzdy = 0.5_ireals * ( hhl(i,j+1,k+1)-hhl(i,j-1,k+1) ) / dy
          ! slope in 4th order accuracy:
          dzdx = ( 2.0_ireals/ 3.0_ireals * ( hhl(i+1,j,k+1)-hhl(i-1,j,k+1) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i+2,j,k+1)-hhl(i-2,j,k+1) ) ) / dx
          dzdy = ( 2.0_ireals/ 3.0_ireals * ( hhl(i,j+1,k+1)-hhl(i,j-1,k+1) )  &
                 - 1.0_ireals/12.0_ireals * ( hhl(i,j+2,k+1)-hhl(i,j-2,k+1) ) ) / dy

          ! normal vector of the surface element (1. order approx.)
          da_x = dzdx * h
          da_y = dzdy * h
          da_z = -h

          flux = 0.5 * ( flux_x(i,j,k) + flux_x(i-1,j,k) ) * da_x   &
            &  + 0.5 * ( flux_y(i,j,k) + flux_y(i,j-1,k) ) * da_y   &
            &  + flux_z(i,j,k+1) * da_z

          int_bottom = int_bottom + flux

        END DO
      END DO

    END IF

  END IF

END SUBROUTINE surface_integral_cond

!==============================================================================
!==============================================================================

SUBROUTINE organize_integrals( yaction )

!-----------------------------------------------------------------------------
! Description:
!   Driver routine
!-----------------------------------------------------------------------------
  
CHARACTER(*) :: yaction
!-----------------------------------------------------------------------------

  IF ( yaction=="init" ) THEN
    CALL  init_integral_3d
    CALL  init_check_conservation

  ELSE IF ( yaction=="compute" ) THEN

    CALL check_qx_conservation  ( 1 )

    CALL check_rho_conservation ( 1 )
    ! flux calculation: upwind = 1, Lax-Wendroff = 2

    ! hints:
    ! 1.) call 'check_rho_conservation' as the last of these 'check_...'-subr.,
    !     because rho(:,:,:) should be updated here
    ! 2.) if you don't call 'check_rho_conservation' then call explicitely
    ! call calc_total_density( nnow)

  ELSE IF ( yaction=="finalize" ) THEN
    CALL final_check_conservation

  END IF

END SUBROUTINE organize_integrals

!==============================================================================
!==============================================================================

SUBROUTINE init_check_conservation
 
!-----------------------------------------------------------------------------
!
! Description:
!   Initialise the testing tool
! Method:
!
!-----------------------------------------------------------------------------
 
  INTEGER (KIND=iintegers) :: istat 
  CHARACTER(100) :: yzerrmsg
  INTEGER (KIND=iintegers) :: izerror

  izerror=0   ! correct value ??
  IF (my_cart_id == 0) THEN
    WRITE(*,*) " INITIALIZE INTEGRATION TOOL:"
  END IF

  ! allocate one storage field for each variable that will be inspected about conservation:
  ALLOCATE( rho       ( 1:ie, 1:je, 1:ke   ), STAT=istat )   ! total density of air
  IF ( istat /= 0 ) THEN
    yzerrmsg="allocation of rho"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'init_check_conservation')
  END IF

  ALLOCATE( rho_q_tot ( 1:ie, 1:je, 1:ke   ), STAT=istat  )  ! total density of moisture
  IF ( istat /= 0 ) THEN
    yzerrmsg="allocation of rho_q_tot"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'init_check_conservation')
  END IF
 
  ! initialize these fields with zero and define 2 storage variables:
  rho(:,:,:)       = 0.0_ireals
  integ_old_rho    = 0.0_ireals
  sum_Res_rho      = 0.0_ireals

  rho_q_tot(:,:,:) = 0.0_ireals
  integ_old_q_tot  = 0.0_ireals
  sum_Res_q_tot    = 0.0_ireals
  
  ! check geometrical properties of the cuboid:
  CALL calc_cuboid_geometry

END SUBROUTINE init_check_conservation
 
!==============================================================================
!==============================================================================
 
SUBROUTINE final_check_conservation
 
!------------------------------------------------------------------------------
! Description:
!   Deallocations
!------------------------------------------------------------------------------
 
  DEALLOCATE( rho )
  DEALLOCATE( rho_q_tot )
 
END SUBROUTINE final_check_conservation
 
!==============================================================================
!==============================================================================
 
SUBROUTINE calc_cuboid_geometry

!------------------------------------------------------------------------------
!
! Description:
!   Calculate the volume of the cuboid and the areas of the 6 sides
!
!------------------------------------------------------------------------------
 
!REAL (kind=ireals) :: integ
!REAL (kind=ireals) :: int_top, int_bottom, int_n, int_s, int_w, int_e
REAL (kind=ireals) :: dummy

REAL (kind=ireals), ALLOCATABLE, DIMENSION(:,:,:) :: field_0, field_1, field_1_
INTEGER (kind=iintegers) :: i, j, k
INTEGER (kind=iintegers) :: istat
INTEGER (KIND=iintegers) :: izerror

!------------------------------------------------------------------------------

  ! allocate fields containing 1 or 0, for a simple application 
  ! of the integration routines:
  ALLOCATE( field_0  ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( field_1  ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( field_1_ ( 1:ie, 1:je, 1:ke+1 ) )

  ! areas by integrating over a flux = 1 

  field_0 (:,:,:) = 0.0_ireals
  field_1 (:,:,:) = 1.0_ireals
  field_1_(:,:,:) = 1.0_ireals

  ! volume of the cuboid = volume integral over 1 :
  integ_cuboid%volume = integral_3d_total( field_1, 3)
 
  ! areas of the 4 side faces = integral over 1 :
  CALL surface_integral_total( field_1, field_1, field_1_,  &
           dummy, dummy, integ_cuboid%area_N, integ_cuboid%area_S, &
           integ_cuboid%area_W, integ_cuboid%area_E )

  integ_cuboid%area_S      = -integ_cuboid%area_S
  integ_cuboid%area_W      = -integ_cuboid%area_W
 
  ! areas of the projections of the top and bottom surfaces to the flat earth:
  CALL surface_integral_total( field_0, field_0, field_1_,  &
           integ_cuboid%area_top, integ_cuboid%area_bottom, &
           dummy, dummy, dummy, dummy )

  integ_cuboid%area_bottom = -integ_cuboid%area_bottom

  IF (my_cart_id == 0) THEN
    WRITE(*,'(A,E15.7,A)') "   cuboid volume:", integ_cuboid%volume, " m^3"
    WRITE(*,'(A,6(A,E14.7),A)') "   cuboid face areas:",                     &
        & " T:", integ_cuboid%area_top, " B:", integ_cuboid%area_bottom,   &
        & " N:", integ_cuboid%area_N,   " S:", integ_cuboid%area_S,        &
        & " W:", integ_cuboid%area_W,   " E:", integ_cuboid%area_E,        &
        & " m^2"
    WRITE(*,*)
  END IF

  DEALLOCATE( field_0, field_1, field_1_ )
 
END SUBROUTINE calc_cuboid_geometry

!==============================================================================
!==============================================================================

SUBROUTINE calc_total_density( n )

!-----------------------------------------------------------------------------
! Description:
!   calculate the total density rho(:,:,:) at timelevel n
!
! Input:
!   fields T, p0, pp, qv, and qc, qi, qr, qs, qg (if they exist)
!   at timelevel n
!
! Output:
!   rho(:,:,:)
!
!-----------------------------------------------------------------------------

  INTEGER, INTENT(in) :: n
  CHARACTER(100) :: yzerrmsg
  INTEGER (KIND=iintegers) :: izerror

  REAL (kind=ireals), ALLOCATABLE :: &
      q_cond(:,:,:) ! specific mass of condensed (liquid+frozen) water 
  INTEGER :: i, j, k, istat

  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:) => NULL(),          &       ! QV at n
    qc  (:,:,:) => NULL(),          &       ! QC at n
    qi  (:,:,:) => NULL(),          &       ! QI at n
    qg  (:,:,:) => NULL(),          &       ! QG at n
    qr  (:,:,:) => NULL(),          &       ! QR at n
    qs  (:,:,:) => NULL()                   ! QS at n

  CHARACTER(LEN=25) :: yzroutine = 'calc_total_density'

  ALLOCATE( q_cond(1:ie, 1:je, 1:ke), STAT=istat)
  IF ( istat /= 0 ) THEN
    yzerrmsg="allocation of q_cond"
    izerror=0
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'calc_total_density')
  END IF

  ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = n, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = n, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = n, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = n, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = n, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = n, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! add up moisture specific mass:

  q_cond(:,:,:) = 0.0_ireals
  q_cond(:,:,:) = q_cond(:,:,:) + qc(:,:,:)
  IF ( ASSOCIATED(qi) )  q_cond(:,:,:) = q_cond(:,:,:) + qi(:,:,:)
  IF ( ASSOCIATED(qr) )  q_cond(:,:,:) = q_cond(:,:,:) + qr(:,:,:)
  IF ( ASSOCIATED(qs) )  q_cond(:,:,:) = q_cond(:,:,:) + qs(:,:,:)
  IF ( ASSOCIATED(qg) )  q_cond(:,:,:) = q_cond(:,:,:) + qg(:,:,:)

  ! calculate density:

  DO i=1, ie
    DO j=1, je
      DO k=1, ke
 
        rho(i,j,k) = ( p0(i,j,k) + pp(i,j,k,n) ) / R_d / T(i,j,k,n)          &
            / (1.0_ireals+(R_v/R_d-1.0_ireals)*qv(i,j,k) - q_cond(i,j,k) )
 
      END DO
    END DO
  END DO

  DEALLOCATE( q_cond, STAT=istat )

END SUBROUTINE calc_total_density

!==============================================================================
!==============================================================================
 
SUBROUTINE check_rho_conservation (fluxmethod)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine is an example of how to apply the 
!   conservation inspection tool 
!   here: check validity of the continuity equation (for density)
!
! Input:
!   fluxmethod = 1: upwind calculation of fluxes 
!   fluxmethod = 2: Lax-Wendroff calculation of fluxes 
!
!   it uses the field 'rho' from the previous timestep
!
! Output:
!   field 'rho' is updated
!   storage variables 'integ_old_rho' and 'sum_Res_rho' are updated
!
! Method:
!
!------------------------------------------------------------------------------
 
INTEGER (KIND=iintegers), INTENT(in) :: fluxmethod

REAL (kind=ireals) :: integ
REAL (kind=ireals) :: int_top, int_bottom, int_n, int_s, int_w, int_e
REAL (kind=ireals) :: dQdt, residuum, divQv
REAL (kind=ireals) :: norm  ! normalization factor for fluxes/volume integrals

 
INTEGER (kind=iintegers) :: i, j, k
INTEGER (kind=iintegers) :: istat
INTEGER (KIND=iintegers) :: izerror

!------------------------------------------------------------------------------
 
  !norm = integ_cuboid%area_bottom ! all integrals and fluxes are
                                   ! normalized by the bottom area of the cuboid
  norm = 1.0_ireals             ! take fluxes/integrals as they are

  ! allocate the flux vector field:
  ALLOCATE( flux_x ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( flux_y ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( flux_z ( 1:ie, 1:je, 1:ke+1 ) )

  ! --- calculate the density field and the appropriate density fluxes ----

  ! -- advective flux (upwind; with rho at the old timestep n-1) --

  SELECT CASE (fluxmethod)

  CASE(1) !!Upwind method with Euler in time

    CALL adv_flux_upwind( rho, u(:,:,:,nnow), v(:,:,:,nnow), w(:,:,:,nnow), &
                          flux_x, flux_y, flux_z )

  CASE(2) !!Lax Wendroff flux measurement with Euler in time

    CALL adv_flux_laxwen( rho, u(:,:,:,nnow), v(:,:,:,nnow), w(:,:,:,nnow), &
                          flux_x, flux_y, flux_z )

  END SELECT

  ! -- calculate rho(:,:,:) at the current timestep n -- 
  CALL calc_total_density( nnow )

  ! ---- apply integration tool -------

  ! -- volume integral:
  integ = integral_3d_total( rho, 3)

  dQdt = ( integ - integ_old_rho ) / dt
 
  integ_old_rho = integ   ! this value is stored for the next timestep!
 
  ! -- surface integral:
  CALL surface_integral_total( flux_x, flux_y, flux_z,  &
           int_top, int_bottom, int_n, int_s, int_w, int_e )
 
  divQv = int_top + int_bottom + int_n + int_s + int_w + int_e
 
  residuum = dQdt + divQv
 
  IF ( ntstep > nstart ) THEN
    ! first timestep always gives a false value for dQdt
    sum_Res_rho = sum_Res_rho + residuum * dt

    IF (my_cart_id == 0) THEN
      WRITE(*,'(A11,6(A,E15.7))') "rho-flx:", &
        &    " T:", int_top/norm, " B:", int_bottom/norm, & 
        &    " N:", int_n  /norm, " S:", int_s     /norm, &
        &    " W:", int_w  /norm, " E:", int_e     /norm
      WRITE(*,'(A11,5(A,E14.7))') "rho_cons:",  &
        &  " Q = ", integ/norm,  ", dQdt = ", dQdt/norm,   &
        &  ", divQv = ", divQv/norm,                           &
        &  ", Res = ", residuum/norm, ", sum_Res = ", sum_Res_rho/norm
      WRITE(*,*)
    ENDIF

  END IF

  DEALLOCATE( flux_x, flux_y, flux_z )
 
END SUBROUTINE check_rho_conservation
 
!==============================================================================
!==============================================================================
 
SUBROUTINE check_qx_conservation (fluxmethod)

  !------------------------------------------------------------------------------
  !
  ! Description:
  !   This subroutine is an example of how to apply the 
  !   conservation inspection tool 
  !   here: check validity of the continuity equation (for moisture density)
  !
  ! Input:
  !   fluxmethod = 1: upwind calculation of fluxes 
  !   fluxmethod = 2: Lax-Wendroff calculation of fluxes  
  ! 
  !   it uses the fields 'rho_q_tot' and 'rho' from the previous timestep
  !
  ! Output:
  !   field 'rho_q_tot' is updated
  !   storage variables 'integ_old_q_tot' and 'sum_Res_q_tot' are updated
  !
  ! Method:
  !
  !------------------------------------------------------------------------------

  USE data_fields     , ONLY :   &
    prr_gsp    ,    & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp    ,    & ! precipitation rate of snow, grid-scale        (kg/m2*s)
    prg_gsp    ,    & ! precipitation rate of graupel, grid-scale     (kg/m2*s)
    prr_con    ,    & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con    ,    & ! precipitation rate of snow, convective        (kg/m2*s)
    qvsflx            ! surface flux of water vapour        (1/m2s) Einheit ???? !MB<

  INTEGER (KIND=iintegers), INTENT(in) :: fluxmethod

  REAL (kind=ireals) :: integ
  REAL (kind=ireals) :: int_top, int_bottom, int_n, int_s, int_w, int_e
  REAL (kind=ireals) :: dQdt, residuum, divQv

   ! specific mass of condensed (liquid+frozen) moisture:
  REAL (kind=ireals), ALLOCATABLE :: q_cond(:,:,:)

  INTEGER (kind=iintegers) :: i, j, k
  INTEGER (kind=iintegers) :: istat

  REAL (kind=ireals), DIMENSION(:,:,:), ALLOCATABLE :: uh, vh, wh
  REAL (kind=ireals), DIMENSION(:,:,:), ALLOCATABLE :: &
    rho_qx  ! temporal field for the density of any kind of water variable
  REAL (kind=ireals) :: rho_qr 
  REAL (kind=ireals) :: rho_tot  ! temporal variable for the total density

  REAL (kind=ireals), DIMENSION(:,:,:), ALLOCATABLE ::    &
    nullflux, sedim_flux_z, diffus_flux_z

  CHARACTER(100) :: yzerrmsg
  INTEGER (KIND=iintegers) :: izerror

  REAL (kind=ireals) :: norm

  REAL (KIND=ireals), POINTER :: &
    qv  (:,:,:) => NULL(),          &       ! QV at nnow
    qc  (:,:,:) => NULL(),          &       ! QC at nnow
    qi  (:,:,:) => NULL(),          &       ! QI at nnow
    qg  (:,:,:) => NULL(),          &       ! QG at nnow
    qr  (:,:,:) => NULL(),          &       ! QR at nnow
    qs  (:,:,:) => NULL()                   ! QS at nnow

  CHARACTER (LEN=25) :: yzroutine = 'check_qx_conservation'

  !------------------------------------------------------------------------------

  !norm = integ_cuboid%area_bottom ! all volumes and fluxes are
                                   ! normalized by the bottom area of the cuboid
  norm = 1.0_ireals             ! take fluxes/integrals as they are

  izerror=0   ! correct value ??

  ! --- Allocations --------------

  ALLOCATE( uh(1:ie, 1:je, 1:ke  ) )
  ALLOCATE( vh(1:ie, 1:je, 1:ke  ) )
  ALLOCATE( wh(1:ie, 1:je, 1:ke+1) )

  ! allocate the flux vector field:
  ALLOCATE( flux_x ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( flux_y ( 1:ie, 1:je, 1:ke   ) )
  ALLOCATE( flux_z ( 1:ie, 1:je, 1:ke+1 ) )

  ALLOCATE( nullflux(1:ie, 1:je, 1:ke  ) )

  ALLOCATE( rho_qx(1:ie, 1:je, 1:ke) ) 

  ! --- calculate the density field and the appropriate density fluxes ----

  nullflux(:,:,:) = 0.0_ireals

  ! Retrieve the required microphysics tracers
  CALL trcr_get(izerror, idt_qv, ptr_tlev = nnow, ptr = qv)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qc, ptr_tlev = nnow, ptr = qc)
  IF (izerror /= 0) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qi, ptr_tlev = nnow, ptr = qi)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qg, ptr_tlev = nnow, ptr = qg)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qr, ptr_tlev = nnow, ptr = qr)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF
  CALL trcr_get(izerror, idt_qs, ptr_tlev = nnow, ptr = qs)
  IF (izerror /= 0 .AND. izerror /= T_ERR_NOTFOUND) THEN
    yzerrmsg = trcr_errorstr(izerror)
    CALL model_abort(my_cart_id, izerror, yzerrmsg, yzroutine)
  ENDIF

  ! === fluxes ==================

  ! --- advective flux (with rho_q_tot at the old timestep n-1) ---

  SELECT CASE (fluxmethod)

  CASE(1) !!Upwind method with Euler in time

    CALL adv_flux_upwind( rho_q_tot, u(:,:,:,nnow), v(:,:,:,nnow), w(:,:,:,nnow), &
      flux_x, flux_y, flux_z )

  CASE(2) !!Lax Wendroff flux measurement with Euler in time

    CALL adv_flux_laxwen( rho_q_tot, u(:,:,:,nnow), v(:,:,:,nnow), w(:,:,:,nnow), &
      flux_x, flux_y, flux_z )

  END SELECT
  
  !CALL flux_bottom_boundary_corr( flux_x, flux_y, flux_z )

  CALL surface_integral_total( flux_x, flux_y, flux_z,  &
    int_top, int_bottom, int_n, int_s, int_w, int_e )
  IF (my_cart_id == 0) THEN
    WRITE(*,'(A11,6(A,E15.7))') "q:adv-flux:", &
        &    " T:", int_top/norm, " B:", int_bottom/norm, & 
        &    " N:", int_n  /norm, " S:", int_s     /norm, &
        &    " W:", int_w  /norm, " E:", int_e     /norm
  END IF

  ! --- sedimentation flux for rain (upwind; qr, rho at old timestep n-1) ---
  !     (gridscale and convective)

  ALLOCATE( sedim_flux_z (1:ie, 1:je, 1:ke+1) )

  sedim_flux_z(:,:,:) = 0.0_ireals

  IF ( ASSOCIATED(qr) ) THEN

    ! use the ground flux of sedimentation of qr
    ! --> this works only, if the lower box boundary is adjacent to the ground!
    IF ( kmax_integ == ke ) THEN
      IF ( lalloc_prr_gsp ) THEN
        !WRITE(*,'(A11,3(A6,E14.7))') "prr_gsp:", "Min=", MINVAL(prr_gsp), "Max=", MAXVAL(prr_gsp), "Mean=", SUM(prr_gsp)/ie/je
        ! Vorzeichen-Definition: prr_gsp > 0 <==> Fluss nach unten, d.h. in den Boden gerichtet
        sedim_flux_z(:,:,ke+1) = sedim_flux_z(:,:,ke+1) - prr_gsp(:,:)
      END IF
      IF ( lalloc_prr_con ) THEN
        !WRITE(*,'(A11,3(A6,E14.7))') "prr_con:", "Min=", MINVAL(prr_con), "Max=", MAXVAL(prr_con), "Mean=", SUM(prr_con)/ie/je
        ! Vorzeichen-Definition: prr_con > 0 <==> Fluss nach unten, d.h. in den Boden gerichtet
        sedim_flux_z(:,:,ke+1) = sedim_flux_z(:,:,ke+1) - prr_con(:,:)
      END IF
    ELSE
      yzerrmsg="box must be adjacent to the ground (k=ke!)"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'check_qx_conservation')
    END IF

  END IF

  ! --- sedimentation flux for snow (upwind; qs at old timestep n-1) ---
  !     (gridscale and convective)
  IF ( ASSOCIATED(qs) ) THEN

    ! use the ground flux of sedimentation of qs
    ! --> this works only, if the lower box boundary is adjacent to the ground!

    IF ( kmax_integ == ke ) THEN
      IF ( lalloc_prs_gsp ) THEN
        !WRITE(*,'(A11,3(A6,E14.7))') "prs_gsp:", "Min=", MINVAL(prs_gsp), "Max=", MAXVAL(prs_gsp), "Mean=", SUM(prs_gsp)/ie/je
        sedim_flux_z(:,:,ke+1) = sedim_flux_z(:,:,ke+1) - prs_gsp(:,:)
      END IF
      IF ( lalloc_prs_con ) THEN
        !WRITE(*,'(A11,3(A6,E14.7))') "prs_con:", "Min=", MINVAL(prs_con), "Max=", MAXVAL(prs_con), "Mean=", SUM(prs_con)/ie/je
        sedim_flux_z(:,:,ke+1) = sedim_flux_z(:,:,ke+1) - prs_con(:,:)
      END IF
    ELSE
      yzerrmsg="box must be adjacent to the ground (k=ke!)"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'check_qx_conservation')
    END IF

  END IF

  ! --- sedimentation flux for graupel (upwind; qg at old timestep n-1) ---
  !     (gridscale)
  IF ( ASSOCIATED(qg) ) THEN

    ! use the ground flux of sedimentation of qg
    ! --> this works only, if the lower box boundary is adjacent to the ground!

    IF ( kmax_integ == ke ) THEN
      IF ( lalloc_prg_gsp ) THEN
        !WRITE(*,'(A11,3(A6,E14.7))') "prg_gsp:", "Min=", MINVAL(prg_gsp), "Max=", MAXVAL(prg_gsp), "Mean=", SUM(prg_gsp)/ie/je
        sedim_flux_z(:,:,ke+1) = sedim_flux_z(:,:,ke+1) - prg_gsp(:,:)
      END IF
    ELSE
      yzerrmsg="box must be adjacent to the ground (k=ke!)"
      CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'check_qx_conservation')
  END IF

  END IF

  CALL surface_integral_total( nullflux, nullflux, sedim_flux_z,  &
    int_top, int_bottom, int_n, int_s, int_w, int_e )
  IF (my_cart_id == 0) THEN
    WRITE(*,'(A11,6(A,E15.7))') "sedim-flx:",              &
        &    " T:", int_top/norm, " B:", int_bottom/norm,  &
        &    " N:", int_n  /norm, " S:", int_s     /norm,  &
        &    " W:", int_w  /norm, " E:", int_e     /norm
  END IF

  flux_z(:,:,:) = flux_z(:,:,:) + sedim_flux_z(:,:,:)

  DEALLOCATE( sedim_flux_z )

  ! --- diffusion fluxes for qv, qc, qi ---

  ! da ich keine Diffusionsfluesse habe, nehme ich einfach den
  ! gesamten Feuchtefluss des Transferschemas
  ! das funktioniert also nur, wenn die Box am Boden aufliegt!!

  ALLOCATE( diffus_flux_z(1:ie, 1:je, 1:ke+1) )

  IF ( kmax_integ == ke ) THEN
    !WRITE(*,'(A11,3(A6,E14.7))') "qvsflx:", "Min=", MINVAL(qvsflx), "Max=", MAXVAL(qvsflx), "Mean=", SUM(qvsflx)/ie/je
    ! Vorzeichen-Definition: qvsflx > 0  <==> Fluss nach unten, d.h. in den Boden gerichtet
    flux_z (:,:,ke+1) = flux_z(:,:,ke+1) - qvsflx(:,:)
    diffus_flux_z(:,:,ke+1) = - qvsflx(:,:)
  ELSE
    yzerrmsg="box must be adjacent to the ground (k=ke!)"
    CALL model_abort (my_cart_id, 100+izerror, yzerrmsg, 'check_qx_conservation')
  END IF

  CALL surface_integral_total( nullflux, nullflux, diffus_flux_z,  &
    int_top, int_bottom, int_n, int_s, int_w, int_e )
  IF (my_cart_id == 0) THEN
    WRITE(*,'(A11,6(A,E15.7))') "qvsflx:",  &
        &    " T:", int_top/norm, " B:", int_bottom/norm, & 
        &    " N:", int_n  /norm, " S:", int_s     /norm, &
        &    " W:", int_w  /norm, " E:", int_e     /norm
  END IF

  DEALLOCATE( diffus_flux_z )
  DEALLOCATE( rho_qx )

  ! === variable rho_q_tot at the current timestep n =========

  ALLOCATE( q_cond(1:ie, 1:je, 1:ke), STAT=istat)

  q_cond(:,:,:) = 0.0_ireals
  q_cond(:,:,:) = q_cond(:,:,:) + qc(:,:,:)
  IF ( ASSOCIATED(qi) )  q_cond(:,:,:) = q_cond(:,:,:) + qi(:,:,:)
  IF ( ASSOCIATED(qr) )  q_cond(:,:,:) = q_cond(:,:,:) + qr(:,:,:)
  IF ( ASSOCIATED(qs) )  q_cond(:,:,:) = q_cond(:,:,:) + qs(:,:,:)
  IF ( ASSOCIATED(qg) )  q_cond(:,:,:) = q_cond(:,:,:) + qg(:,:,:)

      DO k=1, ke
    DO j=1, je
      DO i=1, ie
 
        rho_tot = ( p0(i,j,k) + pp(i,j,k,nnow) ) / R_d / T(i,j,k,nnow)   &
          / ( 1.0_ireals + (R_v/R_d-1.0_ireals)*qv(i,j,k) - q_cond(i,j,k) )

        rho_q_tot(i,j,k) = rho_tot * ( qv(i,j,k) + q_cond(i,j,k) )
 
      END DO
    END DO
  END DO

  DEALLOCATE( q_cond, STAT=istat )
 
  ! ======= apply integration tool ========
 
  ! -- volume integral:
  integ = integral_3d_total( rho_q_tot, 3)
 
  dQdt = ( integ - integ_old_q_tot ) / dt
 
  integ_old_q_tot = integ
 
  ! -- surface integral:
  CALL surface_integral_total( flux_x, flux_y, flux_z,  &
    int_top, int_bottom, int_n, int_s, int_w, int_e )
 
  divQv = int_top + int_bottom + int_n + int_s + int_w + int_e
 
  residuum = dQdt + divQv
 
  IF ( ntstep > nstart ) THEN
    ! first timestep always gives a false value for dQdt
    sum_Res_q_tot = sum_Res_q_tot + residuum * dt

  IF (my_cart_id == 0) THEN
      WRITE(*,'(A11,6(A,E15.7))') "rho_q-flx:", & 
        &    " T:", int_top/norm, " B:", int_bottom/norm, & 
        &    " N:", int_n  /norm, " S:", int_s     /norm, &
        &    " W:", int_w  /norm, " E:", int_e     /norm

      WRITE(*,'(A11,5(A,E14.7))') "rho_q_cons:",  &
        &  " Q = ", integ/norm,  ", dQdt = ", dQdt/norm,   &
        &  ", divQv = ", divQv/norm,                           &
        &  ", Res = ", residuum/norm, ", sum_Res = ", sum_Res_q_tot/norm
    WRITE(*,*)
  ENDIF
 
  END IF

  DEALLOCATE( nullflux )

  DEALLOCATE( flux_x )
  DEALLOCATE( flux_y )
  DEALLOCATE( flux_z )

  DEALLOCATE( uh )
  DEALLOCATE( vh )
  DEALLOCATE( wh )

END SUBROUTINE check_qx_conservation
 
!==============================================================================
!==============================================================================

SUBROUTINE adv_flux_upwind( rho, uh, vh, wh, flux_x, flux_y, flux_z)

!------------------------------------------------------------------------------
!
! Description:
!   calculate advective fluxes (i.e. 'flux = v * rho') upwind 
!
! Method:
!
!------------------------------------------------------------------------------

REAL (KIND=ireals), INTENT(in)  :: rho   (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: uh    (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: vh    (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: wh    (1:ie, 1:je, 1:ke+1)
REAL (KIND=ireals), INTENT(out) :: flux_x(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(out) :: flux_y(1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(out) :: flux_z(1:ie, 1:je, 1:ke+1)

INTEGER :: i,j,k

!------------------------------------------------------------------------------

  DO k=1, ke
    DO j = 1, je
      DO i = 1, ie-1

        IF ( uh(i,j,k) > 0 ) THEN
          flux_x(i,j,k) = uh(i,j,k) * rho(i,j,k)
        ELSE 
          flux_x(i,j,k) = uh(i,j,k) * rho(i+1,j,k)
        END IF

      ENDDO
    ENDDO
  ENDDO

  DO k=1, ke
    DO j = 1, je-1
      DO i = 1, ie

        IF ( vh(i,j,k) > 0 ) THEN
          flux_y(i,j,k) = vh(i,j,k) * rho(i,j,k)
        ELSE 
          flux_y(i,j,k) = vh(i,j,k) * rho(i,j+1,k)
        END IF

      ENDDO
    ENDDO
  ENDDO

  DO k=1+1, ke
    DO j = 1, je
      DO i = 1, ie

        IF ( wh(i,j,k) > 0 ) THEN
          flux_z(i,j,k) = wh(i,j,k) * rho(i,j,k)
        ELSE 
          flux_z(i,j,k) = wh(i,j,k) * rho(i,j,k-1)
        END IF

      ENDDO
    ENDDO
  ENDDO

  DO j = 1, je
    DO i = 1, ie

      flux_z(i,j,1)    = wh(i,j,1)    * rho(i,j,1)

      ! Am Unterrand etwas 'genauer':
      IF ( wh(i,j,ke+1) < 0 ) THEN
        flux_z(i,j,ke+1) = wh(i,j,ke+1) * rho(i,j,ke)
      ELSE
        flux_z(i,j,ke+1) = wh(i,j,ke+1) * rho(i,j,ke) ! d.h. rho wird 'gespiegelt'
        !flux_z(i,j,ke+1) = wh(i,j,ke+1) * ( 0.5_ireals * rho(i,j,ke) ) ! zentr. diff.
        !flux_z(i,j,ke+1) = 0.0_ireals
      END IF

    END DO
  END DO

END SUBROUTINE adv_flux_upwind
 
!==============================================================================

SUBROUTINE adv_flux_laxwen( rho_skal, uh, vh, wh, flux_x, flux_y, flux_z)

!------------------------------------------------------------------------------
!
! Description:
!   calculate advective fluxes (i.e. 'flux = v * rho * skalar') with
!   Lax-Wendroff-scheme.
!
! Method:
!
!------------------------------------------------------------------------------

REAL (KIND=ireals), INTENT(in)  :: rho_skal (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: uh       (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: vh       (1:ie, 1:je, 1:ke)
REAL (KIND=ireals), INTENT(in)  :: wh       (1:ie, 1:je, 1:ke+1)
REAL (KIND=ireals), INTENT(out) :: flux_x   (1:ie, 1:je, 1:ke  )
REAL (KIND=ireals), INTENT(out) :: flux_y   (1:ie, 1:je, 1:ke  )
REAL (KIND=ireals), INTENT(out) :: flux_z   (1:ie, 1:je, 1:ke+1)

INTEGER :: i,j,k

!------------------------------------------------------------------------------

  DO k=1, ke
    DO j = 1, je
      DO i = 1, ie-1

        flux_x(i,j,k) = 0.5_ireals*( uh(i,j,k)*rho_skal(i+1,j,k)      &
                                    +uh(i,j,k)*rho_skal(i  ,j,k) )    &
                       -0.5_ireals* uh(i,j,k)**2.0_ireals*dt          &
                          /(r_earth*crlat(j,1)*dlon*pi/180.0_ireals)  &
                       *( rho_skal(i+1,j,k) - rho_skal(i,j,k) )

      ENDDO
    ENDDO
  ENDDO


  DO k=1, ke
    DO j = 1, je-1
      DO i = 1, ie

        flux_y(i,j,k) = 0.5_ireals*( vh(i,j,k)*rho_skal(i,j+1,k)    &
                                    +vh(i,j,k)*rho_skal(i,  j,k) )  &
                       -0.5_ireals* vh(i,j,k)**2.0_ireals*dt        &
                          /(r_earth*dlat*pi/180.0_ireals)           &
                       *( rho_skal(i,j+1,k) - rho_skal(i,j,k) )

      ENDDO
    ENDDO
  ENDDO


  DO k=2, ke
    DO j = 1, je
      DO i = 1, ie

        flux_z(i,j,k) = 0.5_ireals*( wh(i,j,k)*rho_skal(i,j,k-1)    &
                                    +wh(i,j,k)*rho_skal(i,j,k  ) )  &
                       -0.5_ireals* wh(i,j,k)**2.0_ireals*dt        &
                          /(0.5_ireals*(hhl(i,j,k-1)-hhl(i,j,k+1)) )&
                       *( rho_skal(i,j,k-1) - rho_skal(i,j,k) )

      ENDDO
    ENDDO
  ENDDO
  flux_z(:,:,1)    = wh(:,:,1)    * rho_skal(:,:,1)
  flux_z(:,:,ke+1) = wh(:,:,ke+1) * rho_skal(:,:,ke)

END SUBROUTINE adv_flux_laxwen

!==============================================================================

END MODULE src_integrals
