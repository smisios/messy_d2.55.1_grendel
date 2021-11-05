! ************************************************************************
MODULE messy_qbo
! ************************************************************************

  !     ------------------------------------------------------------
  !     assimilation of QBO zonal wind observations extended 
  !     over 19 levels between 90 and 3 hPa,
  !     data currently
  !     from January 1954 to October 2003 => 598 months
  !      -> can be easily extended
  !     ------------------------------------------------------------
  !
  ! REFERENCES:
  ! * Giorgetta, M. A., and L. Bengtsson (1999), The potential role of the 
  !   quasi-biennial oscillation in the stratosphere-troposphere exchange as 
  !   found in water vapour in general circulation model experiments, J. 
  !   Geophys. Res., 104, 6003-6019.
  ! * Naujokat, B. (1986), An update of the observed quasi-biennial 
  !   oscillation of the stratospheric winds over the tropics.  J. Atmos. 
  !   Sci., 43, 1873-1877.
  ! * Marquardt, C., and B. Naujokat (1997) An update of the equatorial QBO 
  !   and its variability. 1st SPARC General Assembly, Melbourne Australia, 
  !   WMO/TD-No. 814, Vol. 1, 87-90.
  !
  ! AUTHORS:
  !   Marco Giorgetta, MPI for Meteorology, Hamburg, October 1999
  !     - original code for ECHAM4
  !   Maarten van Aalst, MPI for Chemistry, Mainz, 2003/2004
  !     - implementation in ECHAM5 (as preliminary MESSy submodel)
  !   Patrick Joeckel, MPI for Chemistry, Mainz, September 2004
  !     - revision for MESSy 0.9
  !

  USE messy_main_constants_mem,  ONLY: DP, SP, I4, I8

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: EXP, LOG, MIN

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr='qbo'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver='2.3'

  PUBLIC :: DP, SP, I4, I8

  ! GLOBAL PARAMETER
  ! MAXIMUM NUMBER OF LEVELS IN NUDGING DATA
  INTEGER, PARAMETER :: MAXLEV = 25

  ! number of levels
  INTEGER, PUBLIC :: nqbolev = 0
  ! levels [Pa]
  REAL(DP) , DIMENSION(:),   POINTER, PUBLIC :: r_p => NULL()
  ! data
  !     --------------------------------------------------------------
  !              monthly mean zonal wind               [m/s]
  !              Observations at Canton Island, Gan/Maledives and
  !              Singapore, observed at at 10-70/90 hPa and extended
  !              to 3-90hPa, for 598 months from Jan.1954 to Jan.2001.
  !     --------------------------------------------------------------
  ! data-level dependent weight of nudging
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: r_nw => NULL()
  ! data-level dependent half width of QBO [deg]
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: r_hw => NULL()

  ! GLOBAL CTRL NAMELIST VARIABLES
  ! full nudging for |latitude| <= r_lat1 [deg]
  REAL(DP), PUBLIC :: r_lat1 = 10.0_DP
  ! no   nudging for |latitude| >  r_lat2 [deg]
  REAL(DP), PUBLIC :: r_lat2 = 20.0_DP
  ! basic amplitude of the nudging field, [1/s]
  REAL(DP), PUBLIC :: r_nudg0 = 0.0_DP
  ! data-level dependent weight of nudging
  REAL(DP), DIMENSION(MAXLEV), PUBLIC :: r_nweight = -1.0_DP
  ! data-level dependent half width of QBO [deg]
  REAL(DP), DIMENSION(MAXLEV), PUBLIC :: r_hwidth  = -1.0_DP

  ! QBO DATA LINEARLY INTERPOLATED IN TIME
  REAL(DP), DIMENSION(:), POINTER, PUBLIC :: qbodd => NULL()

  PUBLIC :: qbo_read_nml_ctrl
  PUBLIC :: qbo_add_mem
  PUBLIC :: qbo_interpolate_index

  PUBLIC :: qbo_1
  PUBLIC :: qbo_2
  PUBLIC :: qbo_3
  PUBLIC :: qbo_clean

CONTAINS

!--------------------------------------------------------------------------
  SUBROUTINE qbo_read_nml_ctrl(status, iou)

    ! QBO MODULE ROUTINE
    !
    ! read namelist

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'qbo_read_nml_ctrl'

    NAMELIST /CTRL/ r_lat1, r_lat2, r_nudg0, r_nweight, r_hwidth

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <qbo>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CHECK NAMELIST
    WRITE(*,*) 'FULL NUDGING AT   : |latitude| <= ',r_lat1,' deg'
    WRITE(*,*) 'NO   NUDGING AT   : |latitude| >  ',r_lat2,' deg'
    WRITE(*,*) 'NUDGING AMPLITUDE : ',r_nudg0
    WRITE(*,*) 'NUDGING WEIGHT    : ',r_nweight(:)
    WRITE(*,*) 'HALF WIDTH        : ',r_hwidth(:)

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE qbo_read_nml_ctrl
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  SUBROUTINE qbo_add_mem

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'qbo_add_mem'

    ! ALLOCATE SPACE FOR NAMELIST DATA
    ALLOCATE(r_nw(nqbolev))
    ALLOCATE(r_hw(nqbolev))

    r_nw(:) = r_nweight(1:nqbolev)

    r_hw(:) = r_hwidth(1:nqbolev)

  END SUBROUTINE qbo_add_mem
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  SUBROUTINE qbo_interpolate_index(ndim1,ndim2,g1,g2,k2first,k2last,k2i)

    !       This subroutine initializes the module variables which
    !       are used by the interpolation subroutines.
    !
    !      Input:
    !       ------
    !       ndim1   : dimension of grid 1
    !       g1      : coordinates of grid 1
    !       ndim2   : dimension of grid 2
    !       g2      : coordinates of grid 2
    !
    !       Output:
    !       -------
    !       k2first : index of first coordinate of grid 2
    !                 in the range of grid 1
    !       k2last  : index of last coordinate of grid 2
    !                 in the range of grid 1
    !       k2i     : contains for each coordinate of grid2
    !                 in the range of grid1 the index of the
    !                 nearest bigger coordinate of grid1
    !
    !       Data on grid 1 may be interpolated to the index
    !       range [k2first,k2last] of grid 2.
    !
    !       Author: Marco Giorgetta, MPI for Meteorology, Hamburg
    !               October 1999
    
    IMPLICIT NONE
    
    ! I/O
    INTEGER  :: ndim1, ndim2
    REAL(DP) :: g1(ndim1), g2(ndim2)
    INTEGER  :: k2first, k2last, k2i(ndim2)
    INTEGER  :: k1, k2

    !       Find index of first element of grid 2
    !       which is in the range of grid 1
    !       -------------------------------------
    k2first=0 !***CHECKloopstructure***
    k2=1
    loop1: DO
       IF (g2(k2).ge.g1(1)) THEN
          k2first=k2
       ELSE
          k2=k2+1
       ENDIF
       IF (k2first.ne.0) EXIT loop1 
       !CYCLE loop1
    ENDDO loop1
    
    !       Find index of last element of grid 2
    !       which is in the range of grid 1
    !       ------------------------------------
    k2last=0
    k2=ndim2
    loop2: DO
       IF (g2(k2).le.g1(ndim1)) THEN
          k2last=k2
       ELSE
          k2=k2-1
       ENDIF
       IF (k2last.ne.0) EXIT loop2
    ENDDO loop2
    
    !       Find indices k2i for elements k2first
    !       to k2last of grid 2
    !       -------------------------------------
    k2i=0
    k2=k2first
    k1=2
    loop3: DO
       IF (g2(k2).le.g1(k1)) THEN
          k2i(k2)=k1
          k2=k2+1
       ELSE
          k1=k1+1
       ENDIF
       IF (k2.gt.k2last) EXIT loop3
    ENDDO loop3
    
  END SUBROUTINE qbo_interpolate_index
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
#ifndef MBM_QBO
  PURE SUBROUTINE qbo_1(uqbo, iind0, iind1, pfp, lati)
#else
       SUBROUTINE qbo_1(uqbo, iind0, iind1, pfp, lati)
#endif

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: uqbo     ! nudging data (interpolated)
    INTEGER,  INTENT(IN)  :: iind0    ! vertical index
    INTEGER,  INTENT(IN)  :: iind1    ! vertical index
    REAL(DP), INTENT(IN)  :: pfp      ! full level pressure at time t [Pa]
    REAL(DP), INTENT(IN)  :: lati     ! current latitude [deg]

    ! LOCAL
    REAL(DP) :: ulat

    ! compute QBO at equator at level pfp
    uqbo =(qbodd(iind0)*(r_p(iind1)-pfp)+   &
         qbodd(iind1)*(pfp-r_p(iind0)))/  &
         (r_p(iind1)-r_p(iind0))

    ! compute half width of QBO at level pfp
    ulat =(r_hw(iind0)*(r_p(iind1)-pfp)+   &
         r_hw(iind1)*(pfp-r_p(iind0)))/   &
         (r_p(iind1)-r_p(iind0))

    ! compute QBO at latitude lati and level pfp
    uqbo = uqbo*EXP(LOG(.5_DP)*(lati/ulat)**2)

#ifdef MBM_QBO
!for testing in boxmodel
    write(*,"(A,f6.0,A)") 'Pressure = ',pfp,' Pa'
    write(*,"(A,f6.0,A,f7.2,A,f5.2)") 'QBO at ',r_p(iind0),' Pa: ' &
         , qbodd(iind0), ', weight ', (r_p(iind1)-pfp)/ (r_p(iind1)-r_p(iind0))
    write(*,"(A,f6.0,A,f7.2,A,f5.2)") 'QBO at ',r_p(iind1),' Pa: ' &
         , qbodd(iind1), ', weight ', (pfp-r_p(iind0))/ (r_p(iind1)-r_p(iind0))
    write(*,"(A,f6.2)") '--> QBO at this level: ',uqbo
    write(*,*) '------------------------------------------------------'
    write(*,"(A,f6.0,A)") 'Pressure = ',pfp,' Pa'
    write(*,"(A,f6.0,A,f7.2,A,f5.2)") 'QBO at ',r_p(iind0),' Pa: ' &
         , qbodd(iind0), ', weight ' &
         , (log10(r_p(iind1))-log10(pfp))/ (log10(r_p(iind1))-log10(r_p(iind0)))
    write(*,"(A,f6.0,A,f7.2,A,f5.2)") 'QBO at ',r_p(iind1),' Pa: ' &
         , qbodd(iind1), ', weight ' &
         , (log10(pfp)-log10(r_p(iind0)))/ (log10(r_p(iind1))-log10(r_p(iind0)))
    write(*,"(A,f6.2)") '--> QBO at this level: ',uqbo
    write(*,*) '------------------------------------------------------'
    write(*,*) '------------------------------------------------------'
#endif

  END SUBROUTINE qbo_1
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  PURE SUBROUTINE qbo_2(unudg, iind0, iind1, pfp, lati)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: unudg    ! nudging amplitude
    INTEGER,  INTENT(IN)  :: iind0    ! vertical index
    INTEGER,  INTENT(IN)  :: iind1    ! vertical index
    REAL(DP), INTENT(IN)  :: pfp      ! full level pressure at time t [Pa]
    REAL(DP), INTENT(IN)  :: lati     ! current latitude [deg]

    unudg = r_nudg0
    unudg= unudg*(r_nw(iind0)*(r_p(iind1)-pfp)+      &
         r_nw(iind1)*(pfp-r_p(iind0)))/     &
         (r_p(iind1)-r_p(iind0))
    unudg= unudg * MIN((r_lat2-lati)/(r_lat2-r_lat1), 1.0_DP)

  END SUBROUTINE qbo_2
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  PURE SUBROUTINE qbo_3(dudtn, uqbo, unudg, pum1, pvom, twodt)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: dudtn  ! qbo nudging tendency du/dt
    REAL(DP), INTENT(IN)  :: uqbo   ! nudging data (interpolated)
    REAL(DP), INTENT(IN)  :: unudg  ! nudging amplitude    
    REAl(DP), INTENT(IN)  :: pum1   ! zonal wind at t+dt       [m/s]
    REAl(DP), INTENT(IN)  :: pvom   ! zonal wind tendency at t [m/s^2]
    REAl(DP), INTENT(IN)  :: twodt  ! 2 x time step

    dudtn = (uqbo - (pum1 +twodt*pvom) ) / (1.0_DP + twodt*unudg ) * unudg

  END SUBROUTINE qbo_3
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
  SUBROUTINE qbo_clean

    IMPLICIT NONE
    
    ! NAMELIST
    IF (ASSOCIATED(r_nw))      DEALLOCATE(r_nw)
    NULLIFY(r_nw)
    IF (ASSOCIATED(r_hw))      DEALLOCATE(r_hw)
    NULLIFY(r_hw)

  END SUBROUTINE qbo_clean
!--------------------------------------------------------------------------

! ************************************************************************
END MODULE messy_qbo
! ************************************************************************

