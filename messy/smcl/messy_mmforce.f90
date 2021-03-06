! **********************************************************************
MODULE messy_mmforce
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mmforce'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.0'

  ! CTRL NAMELIST SWITCH(ES)
  INTEGER,  PUBLIC :: i_month = 1
  REAL(DP), PUBLIC :: r_dudt  = 2.0_dp

  REAL(DP), PUBLIC, SAVE :: dayno ! number of day in year (set in global_start)

  PUBLIC :: mmforce_read_nml_ctrl
  PUBLIC :: mmforce_physc_int
  !PRIVATE :: mforc_physct
  !PRIVATE :: mforc_physc1

CONTAINS

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_physc_int(vom, press, philat, forc, hw)

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: vom
    REAL(DP), DIMENSION(:,:), INTENT(IN)    :: press
    REAL(DP), DIMENSION(:),   INTENT(IN)    :: philat
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: forc
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: hw

    SELECT CASE(i_month)
    CASE(0)
       CALL mforc_physct(vom, press, philat, forc, hw)
    CASE(1)
       CALL mforc_physc1(vom, press, philat, forc, hw)
    CASE(7)
       CALL mforc_physc1(vom, press, philat, forc, hw)
    CASE DEFAULT
       ! CANNOT BE REACHED
    END SELECT

  END SUBROUTINE mmforce_physc_int
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
    INTEGER, INTENT(OUT) :: status ! error status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'mmforce_read_nml_ctrl'
    LOGICAL                      :: lex          ! file exists?
    INTEGER                      :: fstat        ! file status

    NAMELIST /CTRL/ i_month, r_dudt

    ! initialize
    status = 1 ! error

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML = CTRL, IOSTAT = fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! diagnose namelist and set global switches
    ! check namelist entries
    SELECT CASE(i_month)
    CASE(0)
       WRITE(*,*) substr,': ... TRANSIENT'
    CASE(1)
       WRITE(*,*) substr,': ... PERPETUAL JANUARY'
    CASE(7)
       WRITE(*,*) substr,': ... PERPETUAL JULY'
    CASE DEFAULT
       WRITE(*,*) substr,': ERROR: UNKNOWN SELECTION FOR i_month'
       RETURN
    END SELECT

    WRITE(*,*) substr,': FORCING TENDENCY = ',r_dudt,' m/s/day'

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! no error    

  END SUBROUTINE mmforce_read_nml_ctrl
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mforc_physc1(vom, press, philat, forc, hw)

    IMPLICIT NONE

    INTRINSIC :: ASIN, SIN, SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: vom
    REAL(DP), DIMENSION(:,:), INTENT(IN)    :: press
    REAL(DP), DIMENSION(:),   INTENT(IN)    :: philat
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: forc
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: hw

    ! LOCAL
    REAL(DP) :: lati, latib, pi, dudtf
    INTEGER ::  jk, jp
    INTEGER :: nproma, nlev

    nproma = SIZE(vom,1)
    nlev   = SIZE(vom,2)

    pi=2.0_dp*asin(1.0_dp)

    dudtf=r_dudt/86400.0_dp   ! => 2.3148*10^-5 m/s^2 , 
    !                                         corresponds to 2 m/s /day 

    DO  jp = 1, nproma    ! vector loop

       lati = philat(jp) ! degrees
       latib=lati*pi/180.0_dp    ! radian
       
       DO jk = 1, nlev    ! level loop
                 
!
!     hw(p) = 1.0                    for p < 10 hPa        
!             (100 hPa -p)/90 hPa    for 10 hPa < p < 100 hPa
!             0.                     for p > 100 hPa


          IF (press(jp,jk) > 10000._dp) THEN
             
             hw(jp,jk)=0.0_dp  
             
          ELSE IF (press(jp,jk) <= 10000._dp .and. press(jp,jk) > 1000._dp) THEN
             
             hw(jp,jk)=(10000.0_dp - press(jp,jk))/9000.0_dp 
             
          ELSE

             hw(jp,jk)=1.0_dp
                  
          ENDIF


          ! calculate forcing
          ! for perpetual july runs set IF (lati >= 0.0_dp) THEN ...

          month: IF (i_month == 1) THEN

             IF (lati <= 0.0_dp) THEN   ! zero forcing on SH -> perp. january

                forc(jp,jk)=0.0_dp

             ELSE

                forc(jp,jk)=(sin(2.0_dp*latib))**2 * hw(jp,jk) * dudtf

             ENDIF

          ELSE

             IF (lati >= 0.0_dp) THEN   ! zero forcing on NH -> perp. july

                forc(jp,jk)=0.0_dp

             ELSE

                forc(jp,jk)=(sin(2.0_dp*latib))**2 * hw(jp,jk) * dudtf

             ENDIF

          END IF month

          vom(jp,jk) = vom(jp,jk) + forc(jp,jk)

                 
       ENDDO   ! level loop


    ENDDO   ! vector loop    

END SUBROUTINE mforc_physc1
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mforc_physct(vom, press, philat, forc, hw)

    ! this subroutine contains a time-dependent zonal wind forcing 
    ! and should NOT be used in perpetual january runs

    IMPLICIT NONE

    INTRINSIC :: ASIN, SIN, SIZE

    ! I/O
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: vom
    REAL(DP), DIMENSION(:,:), INTENT(IN)    :: press
    REAL(DP), DIMENSION(:),   INTENT(IN)    :: philat
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: forc
    REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: hw

    ! LOCAL
    REAL(DP) :: lati, latib, pi, dudtf
    REAL(DP) :: twnh, twsh, count
    INTEGER  :: jk, jp
    INTEGER  :: nproma, nlev

    nproma = SIZE(vom,1)
    nlev   = SIZE(vom,2)

    pi=2.0_dp*asin(1.0_dp)
           
    dudtf=r_dudt/86400.0_dp   ! => 2.3148*10^-5 m/s^2 , 
!                                         corresponds to 2 m/s /day 

!      temporal weighting factors
!      forcing on northern hemisphere (twnh) starts on 15th october, 
!      max on 15th january, end of forcing on 16th april
!      forcing on southern hemisphere (twsh) starts on 16th april,
!      max on 16th july, end of forcing on 15th october


    count = dayno - 288.75_dp

    IF (dayno .ge. 106.25_dp .and. dayno .le. 288.75_dp) THEN
           
       twnh = 0.0_dp
       twsh = sin(pi*(count-182.5_dp)/182.5_dp)

    ELSE
                  
       twnh = sin(pi*count/182.5_dp)
       twsh = 0.0_dp

    ENDIF

    DO  jp = 1, nproma    ! vector loop

       lati = philat(jp) ! degrees
       latib=lati*pi/180.0_dp    ! radian
       
       DO jk = 1, nlev    ! level loop

!
!     hw(p) = 1.0                    for p < 10 hPa        
!             (100 hPa -p)/90 hPa    for 10 hPa < p < 100 hPa
!             0.                     for p > 100 hPa


          IF (press(jp,jk) > 10000._dp) THEN
                    
             hw(jp,jk)=0.0_dp  
             
          ELSE IF (press(jp,jk) <= 10000._dp .and. press(jp,jk) > 1000._dp) THEN
             
             hw(jp,jk)=(10000.0_dp - press(jp,jk))/9000.0_dp 
             
          ELSE

             hw(jp,jk)=1.0_dp
             
          ENDIF


          ! calculate forcing
          
          IF (lati <= 0.0_dp) THEN   

             forc(jp,jk)= twsh*(sin(2.0_dp*latib))**2 * hw(jp,jk) * dudtf

          ELSE


             forc(jp,jk)= twnh*(sin(2.0_dp*latib))**2 * hw(jp,jk) * dudtf

          ENDIF


          vom(jp,jk) = vom(jp,jk) + forc(jp,jk)
                 
       ENDDO   ! level loop

    ENDDO   ! vector loop    

  END SUBROUTINE mforc_physct
! ------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_mmforce
! **********************************************************************
