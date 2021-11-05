      MODULE mo_control_bgc
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/mo_control_bgc.f90,v $\\
!$Revision: 1.3.22.1 $\\
!$Date: 2006/04/03 11:27:49 $\\

!***********************************************************************
!
!**** *MODULE mo_control_bgc* - control variables for bgc modules.
!
!     S.Legutke,        *MPI-MaD, HH*    28.02.02
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     - declaration
!
!
!**********************************************************************
      USE mo_kind, ONLY: wp
      implicit none

      LOGICAL :: lspinbgc = .false.    !  Switch for BGC spinup (alkalinity correction)
      REAL(wp)    :: fspinbgc = 1._wp  !  Factor for manual adjustment of alkalinity

      LOGICAL :: invreset = .false.    !  Switch for resetting of carbon inventory

      REAL(wp)    :: deltacalc = 0._wp     !  Correction for calcite pool [kMol/s]
      REAL(wp)    :: deltaorg  = 0._wp     !  Correction for organic carbon pool [kMol/s]
      REAL(wp)    :: deltasil  = 0._wp     !  Correction for silicate pool [kMol/s]

! Logical unit number for I/O.

      INTEGER :: io_stdo_bgc = 6       !  standard out.
      INTEGER :: io_stdi_bgc = 5       !  standard in.

      INTEGER :: io_rsti_bgc = 27      !  restart in.
      INTEGER :: io_rsto_bgc = 27      !  restart out.

! Control variables

      REAL(wp)    :: dtbgc            !  time step length [sec].
      REAL(wp)    :: dtb              !  time step length [days].
      INTEGER :: ndtdaybgc        !  time steps per day.

      INTEGER :: ldtbgc           !  time step number from bgc restart file
      INTEGER :: ldtrunbgc        !  actual time steps of run.


      INTEGER :: icyclibgc        !  switch for cyclicity.
      INTEGER :: ndtrunbgc        !  total no. of time steps of run.

      INTEGER :: bgcstartyear            !  year of ocean restart file
      INTEGER :: bgcstartmonth           !  month of ocean restart file
      INTEGER :: bgcstartday             !  day of ocean restart file

! MPIOM is using variable lyear already !
!      INTEGER :: ldtoce           !  time step number from bgc ocean file

      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL(wp)    :: rmasks = 0.0_wp    !  value at wet cells in sediment.
!      REAL(wp)    :: rmasko = 99999.00  !  value at wet cells in ocean.    ! will overwrite value from namelist!
      REAL(wp)    :: rmasko = -9.e33_wp !  value at wet cells in ocean.    ! will overwrite value from namelist!

      INTEGER :: kchck = 0          !  switch for extended print control (0/1).
      END MODULE mo_control_bgc
