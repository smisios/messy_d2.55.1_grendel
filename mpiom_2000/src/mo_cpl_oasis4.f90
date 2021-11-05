MODULE cpl_oasis4
#ifndef MESSY
   !!======================================================================
   !!                    ***  MODULE cpl_oasis4  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS4
   !!               special case: OPA/LIM coupled to ECHAM5
   !!=====================================================================
   !! History :   
   !!   9.0  !  04-06  (R. Redler, NEC CCRLE, Germany) Original code
   !!   " "  !  04-11  (R. Redler, N. Keenlyside) revision
   !!   " "  !  04-11  (V. Gayler, MPI M&D) Grid writing
   !!   " "  !  05-08  (R. Redler, W. Park) frld initialization, paral(2) revision
   !!   " "  !  05-09  (R. Redler) extended to allow for communication over root only
   !!   " "  !  05-09  (R. Redler) extended to allow for communication over root only
   !!   " "  !  05-12  (R. Hill, Met. Office) Tweaks and hacks to get NEMO/O4 working
   !!   " "  !  06-02  (R. Redler, W. Park) Bug fixes and updates according to the OASIS3 interface
   !!   " "  !  06-02  (R. Redler) app/grid/grid_name from namelist
   !!   " "  !  08-06  (H. Haak) adjustments for MPIOM/ECHAM5
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_oasis4'                    coupled Ocean/Atmosphere via OASIS4
   !!----------------------------------------------------------------------
   !!   cpl_prism_init     : initialization of coupled mode communication
   !!   cpl_prism_define   : definition of grid and fields
   !!   cpl_prism_send     : send out fields in coupled mode
   !!   cpl_prism_recv     : receive fields in coupled mode
   !!   cpl_prism_finalize : finalize the coupled mode communication
   !!----------------------------------------------------------------------
   !! * Modules used
!##################### WARNING coupled mode ###############################
!##################### WARNING coupled mode ###############################
!   Following line must be enabled if coupling with OASIS
!   USE prism                        ! prism module
!##################### WARNING coupled mode ###############################
!##################### WARNING coupled mode ###############################

   USE mo_parallel, ONLY : nprocxy, p_pe   ! message passing
   USE mo_parallel, ONLY : gather            ! message passing
   USE mo_parallel, ONLY : scatter            ! message passing
   USE mo_units, ONLY : io_stdout
   USE mo_commo1, ONLY : dt

!   USE daymod                       ! date and time info
!   USE dom_oce                      ! ocean space and time domain
!   USE in_out_manager               ! I/O manager
!   USE par_oce                      !
!   USE phycst, only : rt0           ! freezing point of sea water
!   USE oasis4_date                  ! OASIS4 date declarations in
                                    ! PRISM compatible format
   IMPLICIT NONE
!
! Exchange parameters for coupling ORCA-LIM with ECHAM5
!
#ifndef __cpl_co2
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nrecv = 16 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nrecv = 15 ! no. of fields passed from atmosphere to ocean
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nsend = 8  ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nsend = 6  ! no. of fields passed from ocean to atmosphere
#endif
#else /*__cpl_co2*/
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nrecv = 18 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nrecv = 17 ! no. of fields passed from atmosphere to ocean
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nsend = 10 ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nsend = 8  ! no. of fields passed from ocean to atmosphere
#endif
#endif /*__cpl_co2*/


   INTEGER, DIMENSION(nsend)  :: send_id
   INTEGER, DIMENSION(nrecv)  :: recv_id

   CHARACTER(len=32)          :: cpl_send (nsend)
   CHARACTER(len=32)          :: cpl_recv (nrecv)

   CHARACTER(len=16)          :: app_name       ! application name for OASIS use
   CHARACTER(len=16)          :: comp_name      ! name of this PRISM component
   CHARACTER(len=16)          :: grid_name(3)      ! name of the grid

! The following now come in via new module oasis4_date
   TYPE(PRISM_Time_struct), PUBLIC    :: dates          ! date info for send operation
   TYPE(PRISM_Time_struct), PUBLIC    :: dates_bound(2) ! date info for send operation
   TYPE(PRISM_Time_struct), PUBLIC    :: dater          ! date info for receive operation
   TYPE(PRISM_Time_struct), PUBLIC    :: dater_bound(2) ! date info for receive operation
   TYPE(PRISM_Time_struct), PUBLIC    :: tmpdate

   PRIVATE

   INTEGER, PARAMETER         :: localRoot  = 0

   INTEGER                    :: localRank      ! local MPI rank
   INTEGER                    :: localSize      ! local MPI size
   INTEGER                    :: localComm      ! local MPI size
   LOGICAL                    :: commRank       ! true for ranks doing OASIS communication
   INTEGER                    :: comp_id        ! id returned by prism_init_comp

   INTEGER                    :: RANGE(5)

   LOGICAL, SAVE              :: prism_was_initialized
   LOGICAL, SAVE              :: prism_was_terminated
   LOGICAL, PARAMETER         :: lwp=.TRUE.

   INTEGER                    :: ierror         ! return error code

   LOGICAL                    :: rootexchg =.FALSE.    ! logical switch


   REAL(wp), DIMENSION(:,:), ALLOCATABLE :: exfld  ! Temporary buffer for exchange
   REAL(wp), DIMENSION(:),   ALLOCATABLE :: buffer ! Temporary buffer for exchange
   INTEGER, DIMENSION(:,:),  ALLOCATABLE :: ranges ! Temporary buffer for exchange

   DOUBLE PRECISION           :: date_incr

   !! Routine accessibility
   PUBLIC cpl_prism_init
   PUBLIC cpl_prism_define
   PUBLIC cpl_prism_send
   PUBLIC cpl_prism_recv
   PUBLIC cpl_prism_finalize

   PUBLIC send_id, recv_id


CONTAINS

   SUBROUTINE cpl_prism_init( localCommunicator )

      IMPLICIT NONE

      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_prism_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS4 software)
      !!
      !! ** Method  :   OASIS4 MPI communication 
      !!--------------------------------------------------------------------
      !! * Arguments
      !!
      INTEGER, INTENT(OUT)       :: localCommunicator
      !!
      !! * Local declarations
      !!

!      NAMELIST/nam_mpp/ app_name, comp_name, c_mpi_send, grid_name

       app_name='mpiom'
       comp_name='ocean'
       grid_name(1)='scalar'
       grid_name(2)='upoint'
       grid_name(3)='vpoint'

      !!
      !!--------------------------------------------------------------------
      !!
      IF(lwp) WRITE(io_stdout,*)
      IF(lwp) WRITE(io_stdout,*) 'cpl_prism_init : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(io_stdout,*) '~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(io_stdout,*)
     


      !------------------------------------------------------------------
      ! 1st Initialize the PRISM system for the application
      !------------------------------------------------------------------

      CALL prism_initialized (prism_was_initialized, ierror)
      IF ( ierror /= PRISM_Success ) &
        CALL prism_abort( comp_id, 'MPIOM', 'cpl_prism_init: Failure in prism_initialized' )

      IF ( .NOT. prism_was_initialized ) THEN
         CALL prism_init( app_name, ierror )
         IF ( ierror /= PRISM_Success ) &
            CALL prism_abort(comp_id, 'MPIOM', 'cpl_prism_init: Failure in prism_init')
         prism_was_initialized = .TRUE.
      ELSE
         CALL prism_abort(comp_id, 'MPIOM', 'cpl_prism_init: Do not initialize prism twice!')
      ENDIF
      !
      ! Obtain the actual dates and date bounds
      !
      ! date is determined by adding days since beginning of
      !   the run to the corresponding initial date. Note that
      !   model internal info about the start date of the experiment
      !   is bypassed. Instead we rely sololy on the info provided
      !   by the SCC.xml file. 
      !
      dates   = PRISM_Jobstart_date

      WRITE(io_stdout,*) "PRISM JOB START DATE IS", dates

      !
      ! upper bound is determined by adding half a time step
      !
      tmpdate = dates
      date_incr = dt/2.0
      CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
      dates_bound(2) = tmpdate
      !
      ! lower bound is determined by half distance to date from previous run
      !
      tmpdate   = dates
      date_incr = -dt/2.0
      CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
      dates_bound(1) = tmpdate

      dater = dates
      dater_bound(1) = dates_bound(1) 
      dater_bound(2) = dates_bound(2) 

      WRITE(io_stdout,*) "DATE send and rec BOUNDS",dater_bound
      WRITE(io_stdout,*) "OTHER BITS FOR DATE",dt


      !------------------------------------------------------------------
      ! 2nd Initialize the PRISM system for the component
      !------------------------------------------------------------------

      CALL prism_init_comp ( comp_id, comp_name, ierror )
      IF ( ierror /= PRISM_Success ) &
         CALL prism_abort (comp_id, 'MPIOM', 'cpl_prism_init: Failure in prism_init_comp')

      WRITE(io_stdout,*) "COMPLETED INIT_COMP",comp_name,comp_id

      !------------------------------------------------------------------
      ! 3rd Get an MPI communicator for MPIOM local communication
      !------------------------------------------------------------------

      CALL prism_get_localcomm ( comp_id, localComm, ierror )
      IF ( ierror /= PRISM_Success ) &
         CALL prism_abort (comp_id, 'MPIOM', 'cpl_prism_init: Failure in prism_get_localcomm' )

      localCommunicator = localComm

       WRITE(io_stdout,*) "COMPLETED GET_LOCALCOMM",comp_name,comp_id


   END SUBROUTINE cpl_prism_init


   SUBROUTINE cpl_prism_define ()

      IMPLICIT NONE

      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_prism_define  ***
      !!
      !! ** Purpose :   Define grid and field information for ocean
      !!    exchange between AGCM, OGCM and COUPLER. (OASIS4 software)
      !!
      !! ** Method  :   OASIS4 MPI communication 
      !!--------------------------------------------------------------------
      !! * Arguments
      !!
      !! * Local declarations

      INTEGER                    :: grid_id(3)     ! id returned by prism_def_grid

      INTEGER                    :: upoint_id(1), &
                                    vpoint_id(1), &
                                    tpoint_id(1)   ! ids returned by prism_set_points

      INTEGER                    :: umask_id(1), &
                                    vmask_id(1), &
                                    tmask_id(1)    ! ids returned by prism_set_mask

      INTEGER                    :: grid_type      ! PRISM grid type

      INTEGER                    :: SHAPE(2,3)     ! shape of arrays passed to PSMILe
      INTEGER                    :: nodim(2)
      INTEGER                    :: data_type      ! data type of transients

      INTEGER                    :: nbr_corners

      LOGICAL                    :: new_points
      LOGICAL                    :: new_mask
      LOGICAL                    :: mask(ie,je,1)

      INTEGER                    :: ji, jj, jk     ! local loop indicees

      CHARACTER(len=32)          :: cpl_send (nsend)
      CHARACTER(len=32)          :: cpl_recv (nrecv)

      CHARACTER(len=32)          :: grid_name      ! name of the grid
      CHARACTER(len=32)          :: point_name     ! name of the grid points

      REAL(kind=wp), ALLOCATABLE :: rclon(:,:,:)
      REAL(kind=wp), ALLOCATABLE :: rclat(:,:,:)
      REAL(kind=wp), ALLOCATABLE :: rcz  (:,:)

      !!--------------------------------------------------------------------
     
      IF(lwp) WRITE(io_stdout,*)
      IF(lwp) WRITE(io_stdout,*) 'cpl_prism_define : initialization in coupled ocean/atmosphere case'
      IF(lwp) WRITE(io_stdout,*) '~~~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(io_stdout,*)
     

      ! -----------------------------------------------------------------
      ! ... Some initialisation
      ! -----------------------------------------------------------------

      send_id = 0
      recv_id = 0

#ifndef NOMPI

      ! -----------------------------------------------------------------
      ! ... Some MPI stuff relevant for optional exchange via root only
      ! -----------------------------------------------------------------

      commRank = .FALSE.

      localRank = p_pe ! from lib_mpp
      localSize = nprocxy ! from lib_mpp

      IF(lwp) WRITE(io_stdout,*) "CALLING DEFINE"

      IF ( rootexchg ) THEN
         IF ( localRank == localRoot ) commRank = .TRUE.
      ELSE
         commRank = .TRUE.
      ENDIF

#else
      !
      ! For non-parallel configurations the one and only process ("localRoot")
      ! takes part in the communication
      ! 
      localRank = localRoot
      commRank = .TRUE.

#endif

      ! -----------------------------------------------------------------
      ! ... Allocate memory for data exchange
      ! -----------------------------------------------------------------


      IF(lwp) WRITE(io_stdout,*) "Abbout to allocate exfld",jpi,jpj

      ALLOCATE(exfld(1:ie,1:je), stat = ierror)
      IF (ierror > 0) THEN
         CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in allocating Reals')
         RETURN
      ENDIF

      IF ( rootexchg .AND. localRank == localRoot ) THEN
         ALLOCATE(ranges(5,0:localSize-1), stat = ierror)
         IF (ierror > 0) THEN
            CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in allocating Integer')
            RETURN
         ENDIF
      ENDIF

      !------------------------------------------------------------------
      ! 1st Declare the local grid characteristics for
      !     surface coupling. The halo regions must be excluded. For
      !     surface coupling it is sufficient to specify only one
      !     vertical z-level.
      !------------------------------------------------------------------

      grid_type = PRISM_irrlonlat_regvrt

      IF(lwp) WRITE(io_stdout,*) "Set grid type"


      ! -----------------------------------------------------------------
      ! ... Define the shape of the valid region without the halo.
      !     For serial configuration (-DNOMPI not being active)
      !     nl* is set to the global values 1 and jp*glo.
      ! -----------------------------------------------------------------

      IF ( rootexchg ) THEN
         SHAPE(1,1) = 2
         SHAPE(2,1) = ie_g-1
         SHAPE(1,2) = 2
         SHAPE(2,2) = je_g-1
         SHAPE(1,3) = 1
         SHAPE(2,3) = 1
      ELSE
         SHAPE(1,1) = 2
         SHAPE(2,1) = ie-1
         SHAPE(1,2) = 2
         SHAPE(2,2) = je-1
         SHAPE(1,3) = 1
         SHAPE(2,3) = 1
      ENDIF

      IF(lwp) WRITE(io_stdout,*) "commrank is", commRank

      IF ( commRank ) THEN

         IF(lwp) WRITE(io_stdout,*) "CALLING DEF_GRID"

         IF(lwp) WRITE(io_stdout,*) "grid name",grid_name
         IF(lwp) WRITE(io_stdout,*) " shape",shape
         IF(lwp) WRITE(io_stdout,*) "grid type",grid_type

         DO ji=1,3        
         CALL prism_def_grid ( grid_id(ji), grid_name(ji), comp_id, shape, &
              grid_type, ierror )
         IF ( ierror /= PRISM_Success ) THEN
            PRINT *, 'cpl_prism_define: Failure in prism_def_grid ',i
            CALL prism_abort (comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_def_grid')
         ENDIF
         ENNDO

         !------------------------------------------------------------------
         ! 2nd Declare the geometic information for this grid.
         !------------------------------------------------------------------

         ! -----------------------------------------------------------------
         ! ... Redefine shape which may now include the halo region as well.
         ! -----------------------------------------------------------------

      IF ( rootexchg ) THEN
         SHAPE(1,1) = 1
         SHAPE(2,1) = ie_g
         SHAPE(1,2) = 1
         SHAPE(2,2) = je_g
         SHAPE(1,3) = 1
         SHAPE(2,3) = 1

! missing gather der corners weiter unten       

      ELSE
         SHAPE(1,1) = 1
         SHAPE(2,1) = ie
         SHAPE(1,2) = 1
         SHAPE(2,2) = je
         SHAPE(1,3) = 1
         SHAPE(2,3) = 1
      ENDIF

         IF(lwp) WRITE(io_stdout,*) "redefined shape",shape

         ! -----------------------------------------------------------------
         ! ... Define the elements, i.e. specify the corner points for each
         !     volume element. In case OPA runs on level coordinates (regular
         !     in the vertical) we only need to give the 4 horizontal corners
         !     for a volume element plus the vertical position of the upper
         !     and lower face. Nevertheless the volume element has 8 corners.
         ! -----------------------------------------------------------------

         !
         ! ... Treat corners in the horizontal plane
         !
         ALLOCATE(rclon(SHAPE(1,1):SHAPE(2,1),SHAPE(1,2):SHAPE(2,2),4), &
              STAT=ierror)
         IF ( ierror /= 0 ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: error in allocate for rclon')

         ALLOCATE(rclat(SHAPE(1,1):SHAPE(2,1),SHAPE(1,2):SHAPE(2,2),4), &
              STAT=ierror)
         IF ( ierror /= 0 ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: error in allocate for rclon')

         nbr_corners = 8

         !define scalar grid  
         !
         ! ... Set longitudes and latitudes
         !          
         DO jj = SHAPE(1,2), SHAPE(2,2)
            DO ji = SHAPE(1,1), SHAPE(2,1)
               rclon(ji,jj,1) = gila(2*ji+1,2*jj+1)*GRARAD
               rclon(ji,jj,2) = gila(2*ji+1,2*jj-1)*GRARAD
               rclon(ji,jj,3) = gila(2*ji-1,2*jj-1)*GRARAD
               rclon(ji,jj,4) = gila(2*ji-1,2*jj+1)*GRARAD

               rclat(ji,jj,1) = giph(2*ji+1,2*jj+1)*GRARAD
               rclat(ji,jj,2) = giph(2*ji+1,2*jj-1)*GRARAD
               rclat(ji,jj,3) = giph(2*ji-1,2*jj-1)*GRARAD
               rclat(ji,jj,4) = giph(2*ji-1,2*jj+1)*GRARAD

            ENDDO
         ENDDO
         !
         ! ... Treat corners along the vertical axis
         !
         ALLOCATE(rcz(SHAPE(1,3):SHAPE(2,3),2), STAT=ierror)
         IF ( ierror /= 0 ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: error in allocate for rcz')

         DO jk = SHAPE(1,3), SHAPE(2,3)
            rcz(jk,1) = tiestw(jk)
            rcz(jk,2) = tiestw(jk+1)
         ENDDO

         IF(lwp) WRITE(io_stdout,*) "ABOUT TO CALL SET CORNERS",shape 

         CALL prism_set_corners ( grid_id(1), nbr_corners, shape, rclon, rclat, &
              rcz, ierror)
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_corners')

         ! upoints

         ! ... Set longitudes and upper latitudes
         !          
         DO jj = SHAPE(1,2), SHAPE(2,2)
            DO ji = SHAPE(1,1), SHAPE(2,1)
               rclon(ji,jj,1) = gila(2*ji+2,2*jj+1)*GRARAD
               rclon(ji,jj,2) = gila(2*ji+2,2*jj-1)*GRARAD
               rclon(ji,jj,3) = gila(2*ji,2*jj-1)*GRARAD
               rclon(ji,jj,4) = gila(2*ji,2*jj+1)*GRARAD

               rclat(ji,jj,1) = giph(2*ji+2,2*jj+1)*GRARAD
               rclat(ji,jj,2) = giph(2*ji+2,2*jj-1)*GRARAD
               rclat(ji,jj,3) = giph(2*ji,2*jj-1)*GRARAD
               rclat(ji,jj,4) = giph(2*ji,2*jj+1)*GRARAD

            ENDDO
         ENDDO
         !
         ! ... Treat corners along the vertical axis
         !
         ALLOCATE(rcz(SHAPE(1,3):SHAPE(2,3),2), STAT=ierror)
         IF ( ierror /= 0 ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: error in allocate for rcz')

         DO jk = SHAPE(1,3), SHAPE(2,3)
            rcz(jk,1) = tiestw(jk)
            rcz(jk,2) = tiestw(jk+1)
         ENDDO

         IF(lwp) WRITE(io_stdout,*) "ABOUT TO CALL SET CORNERS",shape 

         CALL prism_set_corners ( grid_id(1), nbr_corners, shape, rclon, rclat, &
              rcz, ierror)
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_corners')


         ! vpoints

         !
         ! ... Set longitudes and upper latitudes
         !          
         DO jj = SHAPE(1,2), SHAPE(2,2)
            DO ji = SHAPE(1,1), SHAPE(2,1)
               rclon(ji,jj,1) = gila(2*ji+1,2*jj+1)*GRARAD
               rclon(ji,jj,2) = gila(2*ji,2*jj-1)*GRARAD
               rclon(ji,jj,3) = gila(2*ji,2*jj-1)*GRARAD
               rclon(ji,jj,4) = gila(2*ji,2*jj+1)*GRARAD

               rclat(ji,jj,1) = giph(2*ji+2,2*jj+1)*GRARAD
               rclat(ji,jj,2) = giph(2*ji+2,2*jj-1)*GRARAD
               rclat(ji,jj,3) = giph(2*ji,2*jj-1)*GRARAD
               rclat(ji,jj,4) = giph(2*ji,2*jj+1)*GRARAD
            ENDDO
         ENDDO
         !
         ! ... Treat corners along the vertical axis
         !
         ALLOCATE(rcz(SHAPE(1,3):SHAPE(2,3),2), STAT=ierror)
         IF ( ierror /= 0 ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: error in allocate for rcz')

         DO jk = SHAPE(1,3), SHAPE(2,3)
            rcz(jk,1) = tiestw(jk)
            rcz(jk,2) = tiestw(jk+1)
         ENDDO

         IF(lwp) WRITE(io_stdout,*) "ABOUT TO CALL SET CORNERS",shape 

         CALL prism_set_corners ( grid_id(2), nbr_corners, shape, rclon, rclat, &
              rcz, ierror)
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_corners')


         DEALLOCATE(rclon, rclat, rcz)

         ! -----------------------------------------------------------------
         ! ... Define the gridpoints  
         ! -----------------------------------------------------------------

         new_points = .TRUE.

         IF(lwp) WRITE(io_stdout,*) "CALLING SET_POINTS"

         !
         ! ... the u-points
         !
         point_name = 'upoint'
         CALL prism_set_points ( upoint_id(1), point_name, grid_id(2), shape,      &
              glamu, gphiu, gdept(SHAPE(1,3):SHAPE(2,3)), new_points, ierror )
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_points upoint_id')
         !
         ! ... the v-points
         !

         IF(lwp) WRITE(io_stdout,*) "CALLING SET_POINTS done u doing v"

         point_name = 'vpoint'
         CALL prism_set_points ( vpoint_id(1), point_name, grid_id(3), shape,      &
              glamv, gphiv, gdept(SHAPE(1,3):SHAPE(2,3)), new_points, ierror )      
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_points vpoint_id')
         !
         ! ... the t-points

         point_name = 'scalar'
         CALL prism_set_points ( tpoint_id(1), point_name, grid_id(1), shape,   &
              glamt, gphit, gdept(SHAPE(1,3):SHAPE(2,3)), new_points, ierror )
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_points tpoint_id')
         !
         ! ... the f-points
         !
!         point_name = 'f-points'
!         CALL prism_set_points ( fpoint_id(1), point_name, grid_id(1), shape,   &
!              glamf, gphif, gdept(shape(1,3):shape(2,3)), new_points, ierror )
!         IF ( ierror /= PRISM_Success ) &
!              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_points fpoint_id')


         IF(lwp) WRITE(io_stdout,*) "CALLING SET_POINTS done f"

         ! -----------------------------------------------------------------
         ! ... Convert OPA masks to logicals and define the masks
         ! -----------------------------------------------------------------

         new_mask = .TRUE.

         IF ( rootexchg ) THEN
            mask = ( amsuo_g(:,:,1) == 1)
         ELSE
            mask = ( amsuo(:,:,1) == 1)
         ENDIF

         CALL prism_set_mask (umask_id(1), grid_id(2), shape, &
                 mask(SHAPE(1,1):SHAPE(2,1),                  &
                      SHAPE(1,2):SHAPE(2,2),                  &
                      SHAPE(1,3):SHAPE(2,3)),                 &
              new_mask, ierror )
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_mask umask_id')



         IF ( rootexchg ) THEN
            mask = ( amsue_g(:,:,1) == 1)
         ELSE
            mask = ( amsue(:,:,1) == 1)
         ENDIF
        
         CALL prism_set_mask (vmask_id(1), grid_id(3), shape, &
                 mask(SHAPE(1,1):SHAPE(2,1),                  &
                      SHAPE(1,2):SHAPE(2,2),                  &
                      SHAPE(1,3):SHAPE(2,3)),                 &
              new_mask, ierror )
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_mask umask_id')

         IF ( rootexchg ) THEN
            mask = ( weto_g(:,:,1) == 1)
         ELSE
            mask = ( weto(:,:,1) == 1)
         ENDIF

         CALL prism_set_mask (tmask_id(1), grid_id(1), shape, &
                 mask(SHAPE(1,1):SHAPE(2,1),                  &
                      SHAPE(1,2):SHAPE(2,2),                  &
                      SHAPE(1,3):SHAPE(2,3)),                 &
              new_mask, ierror )
         IF ( ierror /= PRISM_Success ) &
              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_mask umask_id')

!         mask = (fmask == 1)
!         CALL prism_set_mask (fmask_id(1), grid_id(1), shape, &
!                 mask(shape(1,1):shape(2,1),                  &
!                      shape(1,2):shape(2,2),                  &
!                      shape(1,3):shape(2,3)),                 &
!              new_mask, ierror )
!         IF ( ierror /= PRISM_Success ) &
!              CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_set_mask umask_id')

         IF(lwp) WRITE(io_stdout,*) "DONE ALL THE SET MASKS"

         

         ! -----------------------------------------------------------------
         ! ... Define the angles
         !   This is needed if zonal tau is not oriented E-W and meridional
         !   tau is not oriented along N-S but rather along local coordinate
         !   axis. Please check!!!!
         ! -----------------------------------------------------------------

!rr      cal prism_set_angles ( ..., ierror ) ! not yet supported by OASIS4

         ! -----------------------------------------------------------------
         ! ... Define the partition 
         ! -----------------------------------------------------------------
         
!         IF ( rootexchg ) THEN
!
!            range(1) = nimpp-1+nldi   ! global start in i
!            range(2) = nlei-nldi+1    ! local size in i of valid region
!            range(3) = njmpp-1+nldj   ! global start in j
!            range(4) = nlej-nldj+1    ! local size in j of valid region
!            range(5) = range(2) &
!                     * range(4)       ! local horizontal size
!            !
!            ! Collect ranges from all NEMO procs on the local root process
!            !
!            CALL mpi_gather(range,  5, MPI_INTEGER, &
!                            ranges, 5, MPI_INTEGER, localRoot, localComm, ierror)
!
!            IF ( localRank == localRoot ) THEN
!
!               maxlen = maxval(ranges(5,:))
!
!               ALLOCATE(buffer(1:maxlen), stat = ierror)
!               IF (ierror > 0) THEN
!                  CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in allocating buffer')
!                  RETURN
!               ENDIF
!
!            ENDIF
!
!         ENDIF

         ! -----------------------------------------------------------------
         ! ... Define the scalefactors 
         ! -----------------------------------------------------------------

!rr      WRITE(io_stdout,*) "CALLING SCALEFACTOR"
!rr      call prism_set_scalefactor ( grid_id(1), shape, e1t, e2t, e3t, ierror )  ! not yet supported by OASIS4
!rr      WRITE(io_stdout,*) "ABOUT TO DEFINE THE TRANSIENTS"

         !------------------------------------------------------------------
         ! 3rd Declare the transient variables
         !------------------------------------------------------------------
         !
         ! ... Define symbolic names for the transient fields send by the ocean
         !     These must be identical to the names specified in the SMIOC file.
         !
         cpl_send( 1)='SSTOCEAN' ! sea surface temperature [K]              
         cpl_send( 2)='SITOCEAN' ! sea ice thickness       [m]             
         cpl_send( 3)='SICOCEAN' ! sea ice concentration   [frac.] 
         cpl_send( 4)='SNTOCEAN' ! surface snow thickness over sea ice [m]  
         cpl_send( 5)='OCUOCEAN' ! ocean eastward velocity  [m/s]                         
         cpl_send( 6)='OCVOCEAN' ! ocean northward velocity  [m/s]                         
#ifdef __cpl_co2
         cpl_send( 7)='CO2TRAOC' ! CO2 transfer coefficient [??]
         cpl_send( 8)='CO"OCEAN' ! pCO2 in the uppermost ocean layer [ppm CO2]             
#endif


         !
         ! ...  Define symbolic names for transient fields received by the ocean.
         !      These must be identical to the names specified in the SMIOC file.
         !
         ! ...  a) U-Grid fields
         !
         cpl_recv( 1)='TXWOCEAU' ! surface downward eastward stress where open sea [pa]
         cpl_recv( 2)='TYWOCEAU' ! surface downward northward stress where open sea [pa]
         cpl_recv( 3)='TXIOCEAU' ! surface downward eastward stress where sea ice [pa]
         cpl_recv( 4)='TYIOCEAU' ! surface downward northward stress where sea ice [pa]
         !
         ! ...  a) V-Grid fields
         !
         cpl_recv( 5)='TXWOCEAV' ! surface downward eastward stress where open sea [pa]
         cpl_recv( 6)='TYWOCEAV' ! surface downward northward stress where open sea [pa]
         cpl_recv( 7)='TXIOCEAV' ! surface downward eastward stress where sea ice [pa]
         cpl_recv( 8)='TYIOCEAV' ! surface downward northward stress where sea ice [pa]
         !
         ! ...  a) T-Grid fields
         !
         cpl_recv( 9)='FRIOCEAN' ! P-E surface downward solid [m/s]
         cpl_recv(10)='FRWOCEAN' ! P-E surface downward liquid [m/s]
         cpl_recv(11)='RHIOCEAN' ! surface downward residual heat flux where sea ice [W/m2]
         cpl_recv(12)='CHIOCEAN' ! surface downward heat flux where sea ice [W/m2]
         cpl_recv(13)='NHWOCEAN' ! surface downward heat flux where open sea [W/m2]
         cpl_recv(14)='SHWOCEAN' ! surface net downward shortwave radiation [W/m2]
         cpl_recv(15)='WSVOCEAN' ! wind speed at 10m [m/s]

#ifdef __cpl_co2
         cpl_recv(16)='CO2CONOC' ! atm. CO2 concentration [ppm]
         cpl_recv(17)='CO2FLXOC' ! CO2 surface upward flux [?]
#endif
         
         data_type = PRISM_DOUBLE_PRECISION

         nodim(1) = 3 ! check
         nodim(2) = 0
         !
         ! ... Announce send variables, all on T points. 
         !
         DO ji = 1, nsend
            CALL prism_def_var (send_id(ji), cpl_send(ji), grid_id(1), &
                 tpoint_id(1), tmask_id(1), nodim, shape, data_type, ierror)
            IF ( ierror /= PRISM_Success ) THEN
               PRINT *, 'Failed to define transient ', ji, TRIM(cpl_send(ji))
               CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_def_var')
            ENDIF
         ENDDO
         !
         nodim(1) = 3 ! check
         nodim(2) = 0
         !
         ! ... Announce recv variables. 
         !
         ! ... a) on U points
         !
         DO ji = 1, 4
            CALL prism_def_var (recv_id(ji), cpl_recv(ji), grid_id(2), &
                 upoint_id(1), umask_id(1), nodim, shape, data_type, ierror)
            IF ( ierror /= PRISM_Success ) THEN
               PRINT *, 'Failed to define transient ', ji, TRIM(cpl_recv(ji))
               CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_def_var')
            ENDIF
         ENDDO
         !
         ! ... b) on V points
         !
         DO ji = 5, 8 
            CALL prism_def_var (recv_id(ji), cpl_recv(ji), grid_id(3), &
                 vpoint_id(1), vmask_id(1), nodim, shape, data_type, ierror)
            IF ( ierror /= PRISM_Success ) THEN
               PRINT *, 'Failed to define transient ', TRIM(cpl_recv(ji))
               CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_def_var')
            ENDIF
         ENDDO
         !
         ! ... c) on T points
         !
         DO ji = 9, nrecv
            CALL prism_def_var (recv_id(ji), cpl_recv(ji), grid_id(1), &
                 tpoint_id(1), tmask_id(1), nodim, shape, data_type, ierror)
            IF ( ierror /= PRISM_Success ) THEN
               PRINT *, 'Failed to define transient ', TRIM(cpl_recv(ji))
               CALL prism_abort ( comp_id, 'MPIOM', 'OPA cpl_prism_define: Failure in prism_def_var')
            ENDIF
         ENDDO

      ENDIF ! commRank

      !------------------------------------------------------------------
      ! 4th End of definition phase
      !------------------------------------------------------------------

      IF(lwp) WRITE(io_stdout,*) "ABOUT TO CALL PRISM_ENDDEF" 

      CALL prism_enddef(ierror)

      IF(lwp) WRITE(io_stdout,*) "DONE ENDDEF",ierror

      IF ( ierror /= PRISM_Success ) &
         CALL prism_abort ( comp_id, 'MPIOM', 'cpl_prism_define: Failure in prism_enddef')
        
      IF(lwp) WRITE(io_stdout,*) "ALL DONE, EXITING PRISM SET UP PHASE"
 
   END SUBROUTINE cpl_prism_define



   SUBROUTINE cpl_prism_send( var_id, date, data_array, info )

      IMPLICIT NONE

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_send  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine sends fields
      !!      like sst or ice cover to the coupler or remote application.
      !!
      !! ** Method  : OASIS4
      !!----------------------------------------------------------------------
      !! * Arguments
      !!
      INTEGER, INTENT( IN )  :: var_id    ! variable Id
      INTEGER, INTENT( OUT ) :: info      ! variable Id
      INTEGER, INTENT( IN )  :: date      ! ocean time-step in seconds
      REAL(wp)               :: data_array(:,:)
      !!
      !! * Local declarations
      !!

      REAL(dp),POINTER               :: global_array(:,:)  nur auf localroot alloaten

      !
      INTEGER                :: request    ! MPI isend request
      INTEGER                :: ji, jj, jn ! local loop indicees

      !!
      INTEGER, SAVE          :: ncount = 0
      !!
      !!--------------------------------------------------------------------
      !!
      ncount = ncount + 1

#ifndef NOMPI

      request = 0

      IF ( rootexchg ) THEN

         ! collect data on the local root process

         IF ( localRank == localRoot ) THEN
            ALLOCATE(global_array(ie_g,je_g))
         ELSE
            global_array => NULL()
         ENDIF

         call gather(data_array,global_array,localRoot)

         ! send data from local root to OASIS4
         !
         CALL prism_put ( var_id, dates, dates_bound, global_array, info, ierror )      

      ELSE
         !
         ! send local data from every process to OASIS4
         !
         CALL prism_put ( var_id, dates, dates_bound, data_array, info, ierror )      

      ENDIF !rootexchg

#else

      !
      ! send local data from every process to OASIS4
      !
      IF ( commRank ) &
      CALL prism_put ( var_id, dates, dates_bound, data_array, info, ierror )

#endif

      IF ( commRank ) THEN

         IF (lwp) THEN

            IF ( info==PRISM_Cpl ) THEN
               WRITE(io_stdout,*) '****************'
               DO ji = 1, nsend
                  IF (var_id == send_id(ji) ) THEN
                     WRITE(io_stdout,*) 'prism_put_proto: Outgoing ', cpl_send(ji)
                     EXIT
                  ENDIF
               ENDDO
               WRITE(io_stdout,*) 'prism_put: var_id       ', var_id
               WRITE(io_stdout,*) 'prism_put:   date       ', date
               WRITE(io_stdout,*) 'prism_put:   info       ', info
               WRITE(io_stdout,*) '     - Minimum value is ', MINVAL(data_array)
               WRITE(io_stdout,*) '     - Maximum value is ', MAXVAL(data_array)
               WRITE(io_stdout,*) '     -     Sum value is ', SUM(data_array)
               WRITE(io_stdout,*) '****************'
            ENDIF

         ENDIF

         IF ( ncount == nsend ) THEN
            !
            !  3. Update dates and dates_bound for next step. We assume that cpl_prism_send
            !  is called for all send fields at each time step. Therefore we update
            !  the date argument to prism_put only every nsend call to cpl_prism_send.
            !
            dates_bound(1) = dates_bound(2)

            tmpdate    = dates_bound(2)
            date_incr  = dt/2.0

            CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
            dates = tmpdate
            CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
            dates_bound(2) = tmpdate

            ncount = 0

         ENDIF

      ENDIF ! commRank

   END SUBROUTINE cpl_prism_send



   SUBROUTINE cpl_prism_recv(  var_id, date, data_array, info )

      IMPLICIT NONE

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_recv  ***
      !!
      !! ** Purpose : - At each coupling time-step,this routine receives fields
      !!      like stresses and fluxes from the coupler or remote application.
      !!
      !! ** Method  : OASIS4
      !!----------------------------------------------------------------------
      !! * Arguments
      !!
      INTEGER, INTENT( IN )  :: var_id    ! variable Id
      INTEGER, INTENT( OUT ) :: info      ! variable Id
      INTEGER, INTENT( IN )  :: date      ! ocean time-step in seconds
      REAL(wp),INTENT( INOUT ) :: data_array(:,:)
      !!
      !! * Local declarations
      !!
      REAL(dp),ALLOCATABLE               :: global_array(:,:)  
      !
      LOGICAL                :: action = .FALSE.
      INTEGER                :: request    ! MPI isend request
      INTEGER                :: ji, jj, jn ! local loop indicees

      INTEGER, SAVE          :: ncount = 0
      !!
      !!--------------------------------------------------------------------
      !!
      ncount  = ncount + 1

#ifndef NOMPI

      request = 0

      IF ( rootexchg ) THEN

         !
         ! receive data from OASIS4 on local root
         !
         IF ( localRank == localRoot ) THEN
            ALLOCATE(global_array(ie_g,je_g))
         ELSE
            ALLOCATE(global_array(0,0))
         ENDIF
         

         IF ( commRank ) &
         CALL prism_get (var_id, dater, dater_bound, global_array, info, ierror)

      ELSE
         !
         ! receive local data from OASIS4 on every process
         !
         CALL prism_get (var_id, dater, dater_bound, exfld, info, ierror)

      ENDIF

      CALL MPI_BCAST ( info, 1, MPI_INTEGER, localRoot, localComm, ierror )

      action = (info==PRISM_CplIO)

      IF ( action ) THEN

      IF ( rootexchg ) THEN

         CALL scatter(data_array,global_array,localroot)

      ELSE

         data_array=exfld
         
      ENDIF


         IF (lwp) THEN        
            WRITE(io_stdout,*) '****************'
            DO ji = 1, nrecv
               IF (var_id == recv_id(ji) ) THEN
                  WRITE(io_stdout,*) 'prism_get: Incoming ', cpl_recv(ji)
                  EXIT
               ENDIF
            ENDDO
            WRITE(io_stdout,*) 'prism_get: var_id       ', var_id
            WRITE(io_stdout,*) 'prism_get:   date       ', date
            WRITE(io_stdout,*) 'prism_get:   info       ', info
            WRITE(io_stdout,*) '     - Minimum value is ', MINVAL(data_array)
            WRITE(io_stdout,*) '     - Maximum value is ', MAXVAL(data_array)
            WRITE(io_stdout,*) '     -     Sum value is ', SUM(data_array)
            WRITE(io_stdout,*) '****************'
         ENDIF

      ENDIF

#else

      CALL prism_get (var_id, dater, dater_bound, exfld, info, ierror)

      IF ( info==PRISM_CplIO ) THEN
               data_array=exfld

         IF (lwp) THEN        
            WRITE(io_stdout,*) '****************'
            DO ji = 1, nrecv
               IF (var_id == recv_id(ji) ) THEN
                  WRITE(io_stdout,*) 'prism_get: Incoming ', cpl_recv(ji)
                  EXIT
               ENDIF
            ENDDO
            WRITE(io_stdout,*) 'prism_get: var_id       ', var_id
            WRITE(io_stdout,*) 'prism_get:   date       ', date
            WRITE(io_stdout,*) 'prism_get:   info       ', info
            WRITE(io_stdout,*) '     - Minimum value is ', MINVAL(data_array)
            WRITE(io_stdout,*) '     - Maximum value is ', MAXVAL(data_array)
            WRITE(io_stdout,*) '     -     Sum value is ', SUM(data_array)
            WRITE(io_stdout,*) '****************'
         ENDIF

      ENDIF

#endif

      IF ( ncount == nrecv ) THEN
         !
         !  3. Update dater and dater_bound for next step. We assume that cpl_prism_recv
         !  is called for all recv fields at each time step. Therefore we update
         !  the date argument to prism_get only every nrecv call to cpl_prism_recv.
         !
         dater_bound(1) = dater_bound(2)

         tmpdate    = dater_bound(2)
         date_incr  = dt/2.0

         CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
         dater = tmpdate
         CALL PRISM_calc_newdate ( tmpdate, date_incr, ierror )
         dater_bound(2) = tmpdate

         ncount = 0

      ENDIF

   END SUBROUTINE cpl_prism_recv



   SUBROUTINE cpl_prism_finalize

      IMPLICIT NONE

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_prism_finalize  ***
      !!
      !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
      !!      called explicitly before cpl_prism_init it will also close
      !!      MPI communication.
      !!
      !! ** Method  : OASIS4
      !!----------------------------------------------------------------------

!      DEALLOCATE(exfld)

      IF ( prism_was_initialized ) THEN

         CALL prism_terminated ( prism_was_terminated, ierror )
         
         IF ( prism_was_terminated ) THEN
            PRINT *, 'prism has already been terminated.'
         ELSE
            CALL prism_terminate ( ierror )
            prism_was_terminated = .TRUE.
         ENDIF

      ELSE

         PRINT *, 'Initialize prism before terminating it.'

      ENDIF


   END SUBROUTINE cpl_prism_finalize


   !!----------------------------------------------------------------------
   !!   Default case           Dummy module         forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
! Removed in elimination of conditional compilation flags
! CONTAINS
!    SUBROUTINE cpl_prism_init             ! Dummy routine
!    END SUBROUTINE cpl_prism_init
!    SUBROUTINE cpl_prism_define           ! Dummy routine
!    END SUBROUTINE cpl_prism_define
!    SUBROUTINE cpl_prism_send             ! Dummy routine
!    END SUBROUTINE cpl_prism_send
!    SUBROUTINE cpl_prism_recv             ! Dummy routine
!    END SUBROUTINE cpl_prism_recv
!    SUBROUTINE cpl_prism_finalize         ! Dummy routine
!    END SUBROUTINE cpl_prism_finalize

#endif
END MODULE cpl_oasis4
