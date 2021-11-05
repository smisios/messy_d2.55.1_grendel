      SUBROUTINE GET_DUST(kpie,kpje,kpke,pddpo)
!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/get_dust.f90,v $\\
!$Revision: 1.2.10.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!****************************************************************
!
!**** *GET_DUST* -
!
!     Iris Kriest,    *MPI-Met, HH*    18.10.02
!
!     Modified
!     --------
!     Patrick Wetzel : read data in NetCDF format.
!
!     Purpose
!     -------
!     - read monthly dust fluxes from nudged data set by C. Timmreck.
!
!     Method
!     -------
!     - read monthly dust fluxes [kg/m2/month] data into array dustin(i,j,k).
!     - if there is a wet cell write value into array dusty.
!     - if there is a dry cell assign 0.
!
!**   Interface.
!     ----------
!
!     *CALL*       *GET_DUST*
!
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      use mo_param1_bgc

      USE MO_PARALLEL
      USE mo_commo1, ONLY: dlxp,dlyp
      USE mo_mpi

implicit none
      INTEGER :: kpie,kpje,kpke,i,j,l

      REAL(wp) ::pddpo(kpie,kpje,kpke)

! define the fields

      REAL(wp) :: dustin(kpie,kpje,12)

! annual iron input
      REAL(wp) :: iron_sum,dust_sum

#ifdef PNETCDF
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat

!
! Open netCDF data file
!
       IF(p_pe==p_io) THEN
        ncstat = NF_OPEN('INPDUST.nc',NF_NOWRITE, ncid)
        IF (ncstat.NE.NF_NOERR ) CALL STOP_ALL('get_dust: Problem with netCDF1')
       END IF
!
! Read  data
       call read_netcdf_var(ncid,'DUST',dustin(1,1,1),12)

!
! Close file
       IF(p_pe==p_io) THEN
        ncstat = NF_CLOSE(ncid)
        IF ( ncstat .NE. NF_NOERR ) CALL STOP_ALL('get_dust: Problem with netCDF200')
       END IF
#else

! define logical io unit (is this safe, so that it doesn't override the
! io definitions in other parts of the program?)
      INTEGER :: io_inbgc5

      io_inbgc5=175

      IF(p_pe == p_io) THEN
       open(io_inbgc5,file='INPDUST',status='unknown',    &
     &        access='sequential',form='unformatted')
      ENDIF

      do l=1,12
       call read_slice(io_inbgc5,dustin(:,:,l))
      enddo

      IF(p_pe == p_io) close(io_inbgc5)

#endif /*PNETCDF*/

! set to missing value (0.) over land
      do l=1,12
        do j=1,kpje
          do i=1,kpie
            IF (pddpo(i, j, 1) .GT. 0.5_wp) THEN
              dusty(i,j,l) = dustin(i,j,l)
            else
              dusty(i,j,l) = 0._wp
            endif
          enddo
         enddo
      enddo

! sum up annual iron input from dust
! units of dustin : kg dust /m**2/yr
! converted to input of bioavailable iron

      dust_sum = 0._wp

      do l=1,12
        do j=1,kpje
          do i=1,kpie
            IF (pddpo(i, j, 1) .GT. 0.5_wp) THEN
              dust_sum = dust_sum &
                   + dustin(i, j, l) * dlxp(i, j) * dlyp(i, j)/12._wp
            endif
          enddo
         enddo
      enddo

      CALL global_sum(dust_sum)

! convert to bioavailable iron
     iron_sum = dust_sum * 0.035_wp * 0.01_wp * 1.e-7_wp
     if (p_pe == p_io) write(io_stdo_bgc,*) 'total annual bioavailable &
     &iron input [1E7 kg Fe(bio)]',iron_sum

      IF(65-p_ioff>=1 .AND. 65-p_ioff<=kpie .AND. &
     &   35-p_joff>=1 .AND. 35-p_joff<=kpje) THEN
        write(io_stdo_bgc,*) 'DUST input at NABE'
        write(io_stdo_bgc,*) 'i= 65 (47.34N) j=35 (19.29N)'
        DO l=1,12
          write(io_stdo_bgc,*) dustin(65-p_ioff,35-p_joff,l)
        ENDDO
      ENDIF

      RETURN
      END
