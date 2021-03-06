MODULE messy_ubcnox

  USE messy_main_constants_mem, ONLY: &
    dp, &                                 ! kind parameter for real
    Pi, STRLEN_ULONG

  IMPLICIT NONE

!-----------------------------------------------------------------
! Everything is PRIVATE, except when explicitely stated otherwise:
!-----------------------------------------------------------------
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modver = '1.5.3'
  CHARACTER(LEN=*), PARAMETER :: modstr = 'ubcnox'
 
  INTEGER, SAVE :: nox_switch, top_levels, ese_switch
  REAL(dp), SAVE :: tempthresh
  CHARACTER (LEN=STRLEN_ULONG), PUBLIC        ::  hn_file=''

  PUBLIC :: ubcnox_read_nml_ctrl, ubcnox_read_noxlut, ubcnox_provide_data &
       , ubcnox_calc_current_ubc, ubcnox_linear_interp

  PUBLIC :: dp, modver, modstr, nox_switch, top_levels, ese_switch, tempthresh



CONTAINS

  SUBROUTINE ubcnox_read_nml_ctrl(status, iou)

  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  USE messy_main_blather, ONLY: start_message, end_message
  
  IMPLICIT NONE
  
  !-----
  ! input/output variables
  !-----
  INTEGER, INTENT(out) :: status   ! error status
  INTEGER, INTENT(in) :: iou       ! logical I/O unit
  
  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'ubcnox_read_nml_ctrl'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status

  NAMELIST /CTRL/ nox_switch, top_levels, ese_switch, tempthresh, hn_file

  CALL start_message(TRIM(modstr),'INITIALISATION', substr)
  !-----
  ! initialisation
  !-----
  status = 1

  CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
  IF (.not.lex) RETURN   ! ubcnox.nml does not exist

  READ(iou, NML=CTRL, IOSTAT=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  IF (fstat /= 0) RETURN   ! error while reading namelist

  write(*,*) ' /----------------------------------\ '
  write(*,*) ' | UPPER BOUNDARY CONDITION FOR NOX | '
  write(*,*) ' \----------------------------------/ '
  write(*,*) ' '
  write(*,*) ' TYPE OF CALCULATION '
  SELECT CASE (nox_switch)
    CASE (1)
      write(*,*) ' You have chosen nox_switch = 1 in ubcnox.nml.'
      write(*,*) ' Number density in Molec/cm?? expected.'
    CASE (2)
      write(*,*) ' You have chosen nox_switch = 2 in ubcnox.nml.'
      write(*,*) ' Flux in Molec/cm??/s expected.'
    CASE (3)
      write(*,*) ' You have chosen nox_switch = 3 in ubcnox.nml.'
      write(*,*) ' Production rate based on Nieder will be used.'
      write(*,*) ' File is used: ',hn_file
      write(*,*) ' WARNING: Do not use for production runs!'
    CASE (4)
      write(*,*) ' You have chosen nox_switch = 4 in ubcnox.nml.'
      write(*,*) ' Number density for Top4 Level in L90 in Molec/cm?? expected.'
    CASE (5)
      write(*,*) ' You have chosen nox_switch = 5 in ubcnox.nml.'
      write(*,*) ' Online calculation of NOy upper boundary condition according to Funke et al 2016'
      write(*,*) ' Level 1 to ', top_levels, ' are used in the parameterization'
      
      SELECT CASE (ese_switch)
        CASE (0)
          write(*,*) ' You have chosen ese_switch = 0 in ubcnox.nml.'
          write(*,*) ' No increase of NOx during Elevated stratopause events'
          write(*,*) ' This setting is not recommended for production runs!'
        CASE (1)
          write(*,*) ' You have chosen ese_switch = 1 in ubcnox.nml.'
          write(*,*) ' Temperature difference between high and low latitudes used to detect ESE.'
          write(*,*) ' You have chosen a temperature threshold of ',tempthresh
        CASE (2)
          write(*,*) ' You have chosen ese_switch = 2 in ubcnox.nml.'
          write(*,*) ' Not implemented yet'
        CASE DEFAULT
          write(*,*) ' Warning: You have chosen an unsupported number for ese_switch'
      END SELECT   ! ese_switch
      
    CASE DEFAULT
      write(*,*) ' Warning: You have chosen an unsupported number for nox_switch'
  END SELECT  ! nox_switch
  
  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error

  CALL end_message(TRIM(modstr),'INITIALISATION', substr)

END SUBROUTINE ubcnox_read_nml_ctrl


! ----------------------------------------------------------------------
  SUBROUTINE ubcnox_read_noxlut(noxlut)
  use netcdf
  USE messy_main_blather, ONLY: start_message, end_message
  
  IMPLICIT NONE
  
  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'ubcnox_read_noxlut'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status
  !character (len = *), parameter :: FILE_NAME = "/home/kit/imk-asf/hj6804/input/messy_ini/raw/ubcnox/noxdbtestnew2-1l_n3.nc"
  !character (len = *) :: FILE_NAME
  ! We are reading 2D data, a 6 x 12 grid. 
  integer, parameter :: nkp = 14, nf107 = 3, nmonth=12, nlat=18
  REAL(dp) :: noxlut(nkp, nmonth, nf107, nlat)

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid
  
  CALL start_message(TRIM(modstr),'READ NOX Look Up Table', substr)
  !-----
  ! initialisation
  !-----

  write(*,*) 'Reading file ',hn_file
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  !call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
  call check( nf90_open(hn_file, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, "pNO", varid) )

  ! Read the data.
  call check( nf90_get_var(ncid, varid, noxlut) )


  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  !print *,"*** SUCCESS reading example file ", FILE_NAME, "! "
  print *,"*** SUCCESS reading example file ", hn_file, "! "



  CALL end_message(TRIM(modstr),'READ NOX Look Up Table', substr)
  
  contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

END SUBROUTINE ubcnox_read_noxlut


! ----------------------------------------------------------------------
SUBROUTINE ubcnox_provide_data(ubc_lat_bin,ubc_kp,ubc_f107)

 REAL(dp) :: ubc_lat_bin(19)
 REAL(dp) :: ubc_kp(14)
 REAL(dp) :: ubc_f107(3)
 
 INTEGER :: i
 
 DO i=1,19
   ubc_lat_bin(i)=-i*10._dp+100._dp
   ! 90,80,70,...,-80,-90
 END DO
 ubc_kp(1)=-1._dp
 ubc_kp(2)=0._dp
 ubc_kp(3)=0.5_dp
 DO i=4,14
   ubc_kp(i)=i*0.5_dp-1._dp
 END DO
 ubc_f107(1)=70._dp
 ubc_f107(2)=130._dp
 ubc_f107(3)=220._dp
 

END SUBROUTINE ubcnox_provide_data


! ----------------------------------------------------------------------
SUBROUTINE ubcnox_calc_current_ubc(ubc,ep,bg,ese,top_levels_in,ap,press,lat,gboxarea,nlon,fac,ese_switch, daysofese, face, dn , ds)

    USE messy_main_constants_mem, ONLY: N_A
    INTEGER, INTENT(in) :: top_levels_in   ! number of levels from top where UBC should be prescribed
    REAL(dp), INTENT(out), DIMENSION(top_levels_in) :: ubc        ! EPP-noy field (complete)
    REAL(dp), INTENT(out), DIMENSION(top_levels_in) :: ep         ! EPP-noy field (without background)
    REAL(dp), INTENT(out), DIMENSION(top_levels_in) :: bg         ! EPP-noy field (only background)
    REAL(dp), INTENT(out), DIMENSION(top_levels_in) :: ese        ! EPP-noy field (only due to ESE)
    REAL(dp), DIMENSION(250), INTENT(in) :: ap          ! Ap Index read in import_ts
    REAL(dp), INTENT(in), DIMENSION(top_levels_in) :: press   ! Pressure in hPascal
    REAL(dp), INTENT(in), DIMENSION(top_levels_in) :: fac   ! Factor calculated in global_start (saved there)
    REAL(dp), INTENT(in), DIMENSION(top_levels_in) :: face  ! Factor for ESE calculated in global_start (saved there)
    REAL(dp), INTENT(in) :: dn,ds          ! days since Jan, 1st (ds) and Jul, 1st (dn)

    REAL(dp), INTENT(in) :: lat             ! latitude in degree
    REAL(dp), INTENT(in) :: gboxarea        ! gridboxarea in cm??
    INTEGER, INTENT(in) :: ese_switch      ! 1: ESE should be detected and additional NOx will be calculated
    INTEGER, INTENT(in) :: daysofese       ! number of day in the ESE
    INTEGER             :: i, j
    REAL(dp), PARAMETER :: tol = 1.e-5
    REAL(dp), PARAMETER :: dw = 250.
    REAL(dp), DIMENSION(250) :: dl
    REAL(dp), DIMENSION(12) :: pref, txnr, txsr, tmnr, tmsr, tfnr, tfsr, &
                               nmnr, nmsr, nfnr, nfsr, wmnr, wmsr, wfnr, wfsr
    REAL(dp), DIMENSION(18) :: lref
    REAL(dp), DIMENSION(9) :: lrefs, lrefn, hs, hn
    REAL(dp), DIMENSION(12,9) :: lds, ldn, lde
    REAL(dp), DIMENSION(12,18) :: amr, a1r, a2r, a3r, p1r, p2r, p3r
    
    REAL(dp) :: tn, ts, nn, ns, wn, ws, seasn, seass, txn, txs, xn, xs, xldn, xlds, &
                facs, facn
    REAL(dp), DIMENSION(250) :: filtern, filters, filtern_flip, filters_flip
    REAL(dp) :: lp, tm, xu, wu, fm, wm, xb, nne, we, rfac, tl, sease, xe           ! variables for ESE
    
    ! local temporary variables
    REAL(dp) :: tmnr1, tmnr2, pref1, pref2, lat1, lat2
    LOGICAL, DIMENSION(12) :: mask
    INTEGER :: ind1, ind2, indlat1, indlat2, nlon
    REAL(dp), DIMENSION(18) :: a1,a2,a3,p1,p2,p3,am,bgh
! ka_sv_20170523+
    REAL(dp) :: startdateofese
! ka_sv_20170523-
    

! start of definition for several variables; values are from Bernd Funke
        DO i=2,250
          dl(i)=i-1
        END DO
        dl(1) = 0.5
        pref=(/1.00, 0.70, 0.50, 0.30, 0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02, 0.01/)
        lref=(/-85.,-75.,-65.,-55.,-45.,-35.,-25.,-15., -5., 5., 15., 25., 35., 45., 55., 65., 75., 85./)
        !%; references latitude centers per hemisphere
        lrefs=(/-85.,-75.,-65.,-55.,-45.,-35.,-25.,-15., -5./)
        lrefn=(/  5., 15., 25., 35., 45., 55., 65., 75., 85./)

        !; transport times (NH and SH)
        txnr=(/40.3,  33.3,  29.8,  38.0,  41.4,  41.1,  36.8,  29.9,  22.7,  13.7,  12.2,  11.7/)
        txsr=(/40.3,  33.3,  28.1,  22.7,  20.2,  19.2,  18.0,  17.0,  15.8,  13.3,  10.6,  10.4/)

        !; maximum density occurence times (NH and SH)
        tmnr=(/179.8, 183.1, 187.4, 193.6, 196.1, 195.9, 192.7, 187.5, 182.1, 175.3, 174.2, 173.8/)
        tmsr=(/195.3, 190.0, 186.1, 182.1, 180.2, 179.4, 178.5, 177.8, 176.9, 175.0, 173.0, 172.8/)

        !; maximum flux occurence times (NH and SH)
        tfnr=(/158.7, 155.4, 152.8, 150.0, 148.5, 147.8, 147.1, 146.8, 146.6, 146.4, 146.0, 145.8/)
        tfsr=(/185.6, 180.6, 176.6, 171.8, 168.9, 167.2, 165.4, 164.0, 163.0, 161.5, 160.1, 159.8/)

        !%;maximum densities (NH and SH)
        nmnr=(/0.00128, 0.00110, 0.00091, 0.00067, 0.00058, 0.00057, 0.00066, 0.00079, 0.00092, 0.00104, 0.00100, 0.00115/)
        nmsr=(/0.00413, 0.00346, 0.00288, 0.00216, 0.00174, 0.00150, 0.00128, 0.00115, 0.00111, 0.00114, 0.00122, 0.00135/)

        !;maximum fluxes (NH and SH)
        nfnr=(/0.00048, 0.00054, 0.00057, 0.00058, 0.00059, 0.00060, 0.00061, 0.00062, 0.00062, 0.00063, 0.00064, 0.00070/)
        nfsr=(/0.00156, 0.00163, 0.00169, 0.00178, 0.00186, 0.00191, 0.00200, 0.00207, 0.00214, 0.00223, 0.00227, 0.00240/)

        !; width parameter densities (NH and SH)
        wmnr=(/0.0558, 0.0469, 0.0383, 0.0280, 0.0242, 0.0244, 0.0284, 0.0345, 0.0409, 0.0479, 0.0475, 0.0474/)
        wmsr=(/0.0631, 0.0598, 0.0562, 0.0508, 0.0480, 0.0470, 0.0469, 0.0477, 0.0483, 0.0479, 0.0473, 0.0472/)

        !; width parameter fluxes (NH and SH)
        wfnr=(/0.1032, 0.0994, 0.0939, 0.0838, 0.0759, 0.0711, 0.0660, 0.0633, 0.0622, 0.0622, 0.0625, 0.0627/)
        wfsr=(/0.0754, 0.0751, 0.0751, 0.0752, 0.0753, 0.0754, 0.0754, 0.0753, 0.0752, 0.0751, 0.0752, 0.0754/)
        
        lds(:, 1)=(/ 0.290, 0.285, 0.288, 0.299, 0.311, 0.319, 0.327, 0.340, 0.357, 0.374, 0.393, 0.412/)
        lds(:, 2)=(/ 0.259, 0.248, 0.241, 0.240, 0.246, 0.260, 0.281, 0.302, 0.309, 0.306, 0.297, 0.294/)
        lds(:, 3)=(/ 0.210, 0.201, 0.195, 0.193, 0.196, 0.202, 0.207, 0.208, 0.206, 0.202, 0.197, 0.187/)
        lds(:, 4)=(/ 0.153, 0.160, 0.161, 0.155, 0.147, 0.136, 0.121, 0.104, 0.092, 0.086, 0.083, 0.078/)
        lds(:, 5)=(/ 0.070, 0.084, 0.092, 0.091, 0.082, 0.069, 0.053, 0.038, 0.028, 0.024, 0.022, 0.021/)
        lds(:, 6)=(/ 0.016, 0.020, 0.021, 0.020, 0.017, 0.013, 0.010, 0.008, 0.007, 0.006, 0.006, 0.006/)
        lds(:, 7)=(/ 0.002, 0.002, 0.002, 0.001, 0.001, 0.000, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001/)
        lds(:, 8)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lds(:, 9)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        WHERE (lds<1.e-20) lds = 1.e-20
        
        ldn(:, 1)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        ldn(:, 2)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        ldn(:, 3)=(/ 0.003, 0.004, 0.004, 0.004, 0.003, 0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002/)
        ldn(:, 4)=(/ 0.011, 0.016, 0.021, 0.023, 0.023, 0.020, 0.016, 0.011, 0.008, 0.006, 0.005, 0.005/)
        ldn(:, 5)=(/ 0.036, 0.049, 0.060, 0.066, 0.067, 0.057, 0.044, 0.030, 0.023, 0.019, 0.018, 0.018/)
        ldn(:, 6)=(/ 0.096, 0.106, 0.120, 0.130, 0.139, 0.131, 0.115, 0.099, 0.088, 0.083, 0.081, 0.080/)
        ldn(:, 7)=(/ 0.185, 0.184, 0.191, 0.206, 0.229, 0.245, 0.250, 0.246, 0.238, 0.230, 0.228, 0.226/)
        ldn(:, 8)=(/ 0.307, 0.291, 0.273, 0.263, 0.261, 0.278, 0.299, 0.318, 0.325, 0.326, 0.325, 0.325/)
        ldn(:, 9)=(/ 0.362, 0.351, 0.332, 0.308, 0.278, 0.266, 0.273, 0.293, 0.316, 0.333, 0.342, 0.344/)
        WHERE (ldn<1.e-20) ldn = 1.e-20

        lde(:, 1)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 2)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 3)=(/ 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        lde(:, 4)=(/ 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.001, 0.001, 0.001, 0.001/)
        lde(:, 5)=(/ 0.009, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.007, 0.007, 0.006, 0.006/)
        lde(:, 6)=(/ 0.044, 0.042, 0.039, 0.036, 0.034, 0.035, 0.037, 0.040, 0.042, 0.042, 0.042, 0.041/)
        lde(:, 7)=(/ 0.146, 0.147, 0.148, 0.139, 0.135, 0.132, 0.137, 0.145, 0.150, 0.153, 0.153, 0.150/)
        lde(:, 8)=(/ 0.333, 0.331, 0.330, 0.334, 0.339, 0.341, 0.344, 0.343, 0.341, 0.336, 0.332, 0.330/)
        lde(:, 9)=(/ 0.467, 0.471, 0.472, 0.481, 0.481, 0.482, 0.471, 0.462, 0.459, 0.462, 0.467, 0.473/)
        WHERE (lde<1.e-20) lde = 1.e-20
        
        !; background densities
        amr(:, 1)=(/1.53e+08, 8.81e+07, 5.44e+07, 2.73e+07, 1.61e+07, 1.09e+07, 6.21e+06, 3.89e+06, 1.57e+06, 1.27e+06, 1.06e+06, 9.00e+05/)
        amr(:, 2)=(/1.59e+08, 9.16e+07, 5.70e+07, 2.84e+07, 1.64e+07, 1.09e+07, 5.99e+06, 3.65e+06, 1.87e+06, 1.45e+06, 9.96e+05, 8.53e+05/)
        amr(:, 3)=(/1.70e+08, 9.86e+07, 6.15e+07, 2.97e+07, 1.69e+07, 1.11e+07, 5.90e+06, 3.59e+06, 2.44e+06, 1.65e+06, 9.38e+05, 8.13e+05/)
        amr(:, 4)=(/2.01e+08, 1.17e+08, 7.23e+07, 3.40e+07, 1.91e+07, 1.25e+07, 6.50e+06, 3.97e+06, 2.70e+06, 1.63e+06, 9.35e+05, 8.10e+05/)
        amr(:, 5)=(/2.41e+08, 1.44e+08, 8.92e+07, 4.11e+07, 2.23e+07, 1.41e+07, 6.95e+06, 4.11e+06, 2.72e+06, 1.67e+06, 9.58e+05, 8.30e+05/)
        amr(:, 6)=(/2.68e+08, 1.60e+08, 9.80e+07, 4.46e+07, 2.39e+07, 1.49e+07, 6.99e+06, 3.95e+06, 2.77e+06, 1.59e+06, 9.19e+05, 7.97e+05/)
        amr(:, 7)=(/2.81e+08, 1.63e+08, 9.65e+07, 4.28e+07, 2.22e+07, 1.37e+07, 6.35e+06, 3.47e+06, 2.37e+06, 1.40e+06, 8.08e+05, 7.00e+05/)
        amr(:, 8)=(/2.81e+08, 1.57e+08, 8.94e+07, 3.83e+07, 1.93e+07, 1.17e+07, 5.33e+06, 2.87e+06, 1.95e+06, 1.17e+06, 6.73e+05, 5.83e+05/)
        amr(:, 9)=(/2.74e+08, 1.46e+08, 8.09e+07, 3.45e+07, 1.74e+07, 1.05e+07, 4.66e+06, 2.46e+06, 1.64e+06, 9.87e+05, 5.69e+05, 4.93e+05/)
        amr(:,10)=(/2.79e+08, 1.48e+08, 8.15e+07, 3.44e+07, 1.70e+07, 1.01e+07, 4.32e+06, 2.22e+06, 1.48e+06, 8.93e+05, 5.15e+05, 4.47e+05/)
        amr(:,11)=(/2.96e+08, 1.62e+08, 9.04e+07, 3.70e+07, 1.75e+07, 1.02e+07, 4.32e+06, 2.24e+06, 1.46e+06, 8.53e+05, 4.92e+05, 4.27e+05/)
        amr(:,12)=(/3.04e+08, 1.73e+08, 9.97e+07, 4.16e+07, 2.02e+07, 1.19e+07, 5.23e+06, 2.71e+06, 1.73e+06, 9.93e+05, 5.73e+05, 4.97e+05/)
        amr(:,13)=(/2.81e+08, 1.65e+08, 9.82e+07, 4.28e+07, 2.16e+07, 1.28e+07, 5.76e+06, 3.04e+06, 1.95e+06, 1.09e+06, 6.27e+05, 5.43e+05/)
        amr(:,14)=(/2.36e+08, 1.41e+08, 8.52e+07, 3.82e+07, 1.97e+07, 1.20e+07, 5.63e+06, 3.08e+06, 1.97e+06, 1.07e+06, 6.15e+05, 5.33e+05/)
        amr(:,15)=(/1.99e+08, 1.18e+08, 7.16e+07, 3.28e+07, 1.76e+07, 1.10e+07, 5.29e+06, 2.89e+06, 1.92e+06, 1.05e+06, 6.04e+05, 5.23e+05/)
        amr(:,16)=(/1.82e+08, 1.04e+08, 6.26e+07, 2.97e+07, 1.64e+07, 1.05e+07, 5.18e+06, 2.78e+06, 1.92e+06, 1.01e+06, 5.73e+05, 4.97e+05/)
        amr(:,17)=(/1.75e+08, 9.65e+07, 5.76e+07, 2.81e+07, 1.64e+07, 1.10e+07, 5.97e+06, 3.42e+06, 2.66e+06, 1.47e+06, 8.38e+05, 7.20e+05/)
        amr(:,18)=(/1.71e+08, 9.43e+07, 5.66e+07, 2.77e+07, 1.76e+07, 1.30e+07, 8.19e+06, 5.14e+06, 3.04e+06, 1.65e+06, 1.06e+06, 8.87e+05/)

        a1r(:, 1)=(/0.548, 0.592, 0.635, 0.691, 0.714, 0.697, 0.678, 0.648, 0.815, 1.035, 1.098, 1.158/)
        a1r(:, 2)=(/0.483, 0.532, 0.570, 0.621, 0.636, 0.615, 0.596, 0.576, 0.742, 0.962, 1.004, 1.036/)
        a1r(:, 3)=(/0.333, 0.381, 0.412, 0.463, 0.481, 0.472, 0.458, 0.440, 0.570, 0.813, 0.838, 0.843/)
        a1r(:, 4)=(/0.167, 0.235, 0.271, 0.315, 0.333, 0.332, 0.331, 0.356, 0.505, 0.805, 0.823, 0.823/)
        a1r(:, 5)=(/0.164, 0.187, 0.216, 0.227, 0.209, 0.185, 0.177, 0.228, 0.345, 0.583, 0.592, 0.592/)
        a1r(:, 6)=(/0.217, 0.211, 0.202, 0.147, 0.107, 0.103, 0.149, 0.206, 0.254, 0.375, 0.376, 0.376/)
        a1r(:, 7)=(/0.245, 0.221, 0.201, 0.159, 0.158, 0.172, 0.215, 0.258, 0.256, 0.285, 0.285, 0.285/)
        a1r(:, 8)=(/0.213, 0.172, 0.151, 0.139, 0.151, 0.177, 0.238, 0.280, 0.270, 0.268, 0.268, 0.268/)
        a1r(:, 9)=(/0.117, 0.099, 0.099, 0.114, 0.142, 0.174, 0.232, 0.273, 0.272, 0.264, 0.264, 0.264/)
        a1r(:,10)=(/0.097, 0.172, 0.220, 0.220, 0.186, 0.193, 0.213, 0.230, 0.242, 0.243, 0.243, 0.243/)
        a1r(:,11)=(/0.170, 0.244, 0.303, 0.299, 0.229, 0.201, 0.169, 0.165, 0.152, 0.148, 0.148, 0.148/)
        a1r(:,12)=(/0.192, 0.246, 0.297, 0.278, 0.216, 0.187, 0.137, 0.095, 0.097, 0.106, 0.106, 0.106/)
        a1r(:,13)=(/0.132, 0.190, 0.239, 0.227, 0.190, 0.205, 0.207, 0.128, 0.147, 0.112, 0.112, 0.112/)
        a1r(:,14)=(/0.112, 0.166, 0.231, 0.246, 0.237, 0.248, 0.238, 0.173, 0.193, 0.139, 0.143, 0.143/)
        a1r(:,15)=(/0.084, 0.136, 0.229, 0.283, 0.290, 0.301, 0.278, 0.200, 0.161, 0.220, 0.231, 0.231/)
        a1r(:,16)=(/0.053, 0.126, 0.227, 0.291, 0.327, 0.351, 0.329, 0.259, 0.294, 0.623, 0.648, 0.649/)
        a1r(:,17)=(/0.191, 0.234, 0.324, 0.425, 0.433, 0.407, 0.348, 0.185, 0.265, 0.552, 0.580, 0.591/)
        a1r(:,18)=(/0.278, 0.342, 0.446, 0.602, 0.443, 0.251, 0.178, 0.324, 0.296, 0.538, 0.646, 0.724/)

        a2r(:, 1)=(/0.039, 0.115, 0.175, 0.228, 0.254, 0.270, 0.337, 0.359, 0.477, 0.540, 0.514, 0.484/)
        a2r(:, 2)=(/0.020, 0.090, 0.136, 0.180, 0.215, 0.243, 0.318, 0.333, 0.427, 0.496, 0.491, 0.484/)
        a2r(:, 3)=(/0.097, 0.139, 0.164, 0.166, 0.185, 0.206, 0.251, 0.278, 0.369, 0.416, 0.413, 0.411/)
        a2r(:, 4)=(/0.155, 0.194, 0.222, 0.200, 0.182, 0.184, 0.229, 0.266, 0.321, 0.349, 0.346, 0.346/)
        a2r(:, 5)=(/0.150, 0.192, 0.210, 0.192, 0.155, 0.157, 0.225, 0.268, 0.268, 0.262, 0.259, 0.259/)
        a2r(:, 6)=(/0.162, 0.192, 0.201, 0.187, 0.169, 0.172, 0.211, 0.250, 0.239, 0.220, 0.219, 0.219/)
        a2r(:, 7)=(/0.162, 0.168, 0.168, 0.172, 0.170, 0.175, 0.203, 0.213, 0.213, 0.205, 0.205, 0.205/)
        a2r(:, 8)=(/0.128, 0.129, 0.133, 0.149, 0.156, 0.167, 0.196, 0.203, 0.207, 0.202, 0.202, 0.202/)
        a2r(:, 9)=(/0.106, 0.106, 0.114, 0.142, 0.150, 0.161, 0.188, 0.200, 0.225, 0.229, 0.229, 0.229/)
        a2r(:,10)=(/0.113, 0.096, 0.099, 0.144, 0.160, 0.176, 0.221, 0.235, 0.253, 0.250, 0.250, 0.250/)
        a2r(:,11)=(/0.095, 0.097, 0.112, 0.159, 0.176, 0.189, 0.225, 0.225, 0.239, 0.237, 0.237, 0.237/)
        a2r(:,12)=(/0.095, 0.096, 0.115, 0.159, 0.192, 0.211, 0.271, 0.293, 0.295, 0.296, 0.296, 0.296/)
        a2r(:,13)=(/0.116, 0.125, 0.142, 0.173, 0.215, 0.237, 0.312, 0.344, 0.363, 0.368, 0.368, 0.368/)
        a2r(:,14)=(/0.126, 0.144, 0.166, 0.196, 0.231, 0.243, 0.303, 0.332, 0.341, 0.374, 0.370, 0.370/)
        a2r(:,15)=(/0.077, 0.116, 0.169, 0.240, 0.270, 0.270, 0.288, 0.284, 0.323, 0.357, 0.353, 0.353/)
        a2r(:,16)=(/0.107, 0.144, 0.197, 0.295, 0.327, 0.329, 0.312, 0.215, 0.260, 0.380, 0.383, 0.383/)
        a2r(:,17)=(/0.154, 0.214, 0.269, 0.371, 0.421, 0.449, 0.476, 0.430, 0.557, 0.674, 0.680, 0.680/)
        a2r(:,18)=(/0.185, 0.249, 0.301, 0.399, 0.502, 0.603, 0.728, 0.785, 0.820, 0.720, 0.689, 0.665/)

        a3r(:, 1)=(/0.096, 0.097, 0.094, 0.105, 0.084, 0.082, 0.053, 0.027, 0.117, 0.174, 0.210, 0.250/)
        a3r(:, 2)=(/0.082, 0.056, 0.038, 0.046, 0.040, 0.041, 0.021, 0.026, 0.106, 0.131, 0.137, 0.143/)
        a3r(:, 3)=(/0.082, 0.048, 0.022, 0.026, 0.034, 0.037, 0.040, 0.062, 0.119, 0.092, 0.094, 0.097/)
        a3r(:, 4)=(/0.102, 0.098, 0.077, 0.058, 0.070, 0.080, 0.072, 0.072, 0.090, 0.076, 0.077, 0.077/)
        a3r(:, 5)=(/0.084, 0.087, 0.081, 0.074, 0.074, 0.076, 0.058, 0.035, 0.043, 0.027, 0.025, 0.025/)
        a3r(:, 6)=(/0.044, 0.045, 0.044, 0.034, 0.022, 0.027, 0.048, 0.048, 0.067, 0.073, 0.073, 0.073/)
        a3r(:, 7)=(/0.028, 0.032, 0.027, 0.011, 0.018, 0.021, 0.023, 0.019, 0.029, 0.048, 0.048, 0.048/)
        a3r(:, 8)=(/0.006, 0.011, 0.017, 0.003, 0.015, 0.025, 0.040, 0.051, 0.048, 0.052, 0.052, 0.052/)
        a3r(:, 9)=(/0.052, 0.025, 0.019, 0.005, 0.014, 0.019, 0.023, 0.037, 0.047, 0.053, 0.053, 0.053/)
        a3r(:,10)=(/0.067, 0.035, 0.009, 0.003, 0.016, 0.026, 0.041, 0.050, 0.033, 0.026, 0.026, 0.026/)
        a3r(:,11)=(/0.054, 0.043, 0.027, 0.021, 0.020, 0.023, 0.037, 0.033, 0.014, 0.020, 0.020, 0.020/)
        a3r(:,12)=(/0.036, 0.027, 0.023, 0.031, 0.034, 0.035, 0.055, 0.082, 0.078, 0.068, 0.068, 0.068/)
        a3r(:,13)=(/0.035, 0.028, 0.026, 0.029, 0.033, 0.044, 0.083, 0.103, 0.070, 0.060, 0.060, 0.060/)
        a3r(:,14)=(/0.051, 0.051, 0.063, 0.060, 0.052, 0.055, 0.076, 0.096, 0.076, 0.090, 0.091, 0.091/)
        a3r(:,15)=(/0.076, 0.084, 0.103, 0.091, 0.069, 0.056, 0.056, 0.056, 0.028, 0.050, 0.051, 0.051/)
        a3r(:,16)=(/0.081, 0.094, 0.125, 0.112, 0.079, 0.043, 0.050, 0.104, 0.185, 0.171, 0.168, 0.167/)
        a3r(:,17)=(/0.084, 0.089, 0.106, 0.093, 0.087, 0.095, 0.186, 0.230, 0.326, 0.272, 0.260, 0.252/)
        a3r(:,18)=(/0.071, 0.098, 0.144, 0.156, 0.125, 0.257, 0.461, 0.570, 0.666, 0.572, 0.545, 0.526/)

        p1r(:, 1)=(/ 1.602,  1.679,  1.741,  1.801,  1.818,  1.809,  1.755,  1.747,  1.790,  1.785,  1.773,  1.762/)
        p1r(:, 2)=(/ 1.651,  1.723,  1.748,  1.806,  1.828,  1.814,  1.764,  1.777,  1.845,  1.846,  1.846,  1.846/)
        p1r(:, 3)=(/ 1.694,  1.749,  1.748,  1.836,  1.850,  1.842,  1.824,  1.902,  1.939,  1.879,  1.877,  1.877/)
        p1r(:, 4)=(/ 1.219,  1.443,  1.499,  1.702,  1.804,  1.845,  1.919,  2.041,  1.996,  1.895,  1.893,  1.893/)
        p1r(:, 5)=(/ 0.432,  0.846,  1.110,  1.452,  1.695,  1.857,  2.231,  2.427,  2.249,  2.031,  2.028,  2.028/)
        p1r(:, 6)=(/ 0.630,  0.843,  1.082,  1.499,  2.313,  2.828, -2.911, -2.987,  2.863,  2.421,  2.419,  2.419/)
        p1r(:, 7)=(/ 1.038,  1.199,  1.444,  1.983,  2.627,  2.939, -2.945, -2.925,  3.090,  2.774,  2.774,  2.774/)
        p1r(:, 8)=(/ 1.477,  1.656,  1.879,  2.346,  2.854,  3.101, -2.988, -2.942, -3.014,  3.133,  3.133,  3.133/)
        p1r(:, 9)=(/ 2.157,  2.820, -2.956, -2.741, -2.698, -2.716, -2.815, -2.912, -2.936, -2.981, -2.981, -2.981/)
        p1r(:,10)=(/-2.317, -2.115, -1.991, -2.006, -2.281, -2.494, -2.779, -2.903, -2.961, -3.033, -3.033, -3.033/)
        p1r(:,11)=(/-1.631, -1.698, -1.693, -1.765, -1.971, -2.168, -2.546, -2.668, -2.933,  3.120,  3.120,  3.120/)
        p1r(:,12)=(/-1.416, -1.487, -1.500, -1.557, -1.596, -1.693, -1.977, -2.292, -2.981, -2.993, -2.993, -2.993/)
        p1r(:,13)=(/-1.726, -1.725, -1.704, -1.801, -1.948, -2.124, -2.327, -2.928,  2.714, -2.950, -2.949, -2.949/)
        p1r(:,14)=(/-2.752, -2.334, -2.184, -2.241, -2.358, -2.446, -2.586, -3.115,  2.718, -2.271, -2.233, -2.233/)
        p1r(:,15)=(/ 3.018, -2.473, -2.344, -2.417, -2.477, -2.476, -2.567, -3.038,  3.025, -1.889, -1.855, -1.855/)
        p1r(:,16)=(/-1.407, -1.861, -2.034, -2.075, -2.125, -2.135, -2.178, -2.254, -1.874, -1.433, -1.422, -1.423/)
        p1r(:,17)=(/-1.115, -1.459, -1.661, -1.662, -1.660, -1.654, -1.594, -1.550, -1.308, -1.283, -1.294, -1.307/)
        p1r(:,18)=(/-1.160, -1.413, -1.562, -1.595, -1.492, -1.167,  0.099,  0.825,  0.323, -1.031, -1.156, -1.227/)

        p2r(:, 1)=(/ 1.221,  0.838,  0.889,  1.131,  1.276,  1.233,  1.254,  1.334,  1.657,  1.800,  1.824,  1.845/)
        p2r(:, 2)=(/-0.509,  0.295,  0.541,  0.951,  1.152,  1.157,  1.229,  1.271,  1.629,  1.752,  1.743,  1.729/)
        p2r(:, 3)=(/-1.334, -0.756, -0.422,  0.216,  0.520,  0.643,  0.873,  0.991,  1.469,  1.646,  1.644,  1.637/)
        p2r(:, 4)=(/-1.420, -1.015, -0.756, -0.228,  0.208,  0.539,  1.001,  1.257,  1.574,  1.730,  1.731,  1.731/)
        p2r(:, 5)=(/-1.209, -0.827, -0.579, -0.177,  0.253,  0.667,  1.172,  1.400,  1.564,  1.675,  1.671,  1.671/)
        p2r(:, 6)=(/-0.682, -0.389, -0.094,  0.371,  0.794,  1.033,  1.401,  1.556,  1.661,  1.757,  1.756,  1.756/)
        p2r(:, 7)=(/-0.226,  0.024,  0.325,  0.776,  1.054,  1.170,  1.406,  1.531,  1.638,  1.672,  1.672,  1.672/)
        p2r(:, 8)=(/ 0.332,  0.554,  0.771,  1.138,  1.290,  1.347,  1.433,  1.496,  1.526,  1.487,  1.487,  1.487/)
        p2r(:, 9)=(/ 1.278,  1.177,  1.096,  1.287,  1.349,  1.368,  1.412,  1.458,  1.519,  1.504,  1.504,  1.504/)
        p2r(:,10)=(/ 1.611,  1.572,  1.471,  1.584,  1.614,  1.622,  1.721,  1.816,  1.862,  1.854,  1.854,  1.854/)
        p2r(:,11)=(/ 1.050,  1.384,  1.604,  1.703,  1.701,  1.731,  1.844,  1.993,  2.003,  1.965,  1.965,  1.965/)
        p2r(:,12)=(/ 0.388,  0.798,  1.191,  1.530,  1.687,  1.767,  1.886,  1.987,  2.122,  2.186,  2.186,  2.186/)
        p2r(:,13)=(/-0.033,  0.302,  0.607,  1.055,  1.386,  1.542,  1.814,  1.966,  2.128,  2.175,  2.175,  2.175/)
        p2r(:,14)=(/-0.132,  0.136,  0.390,  0.850,  1.219,  1.406,  1.719,  1.938,  2.087,  2.089,  2.086,  2.086/)
        p2r(:,15)=(/ 0.322,  0.513,  0.635,  0.973,  1.237,  1.353,  1.604,  1.868,  2.040,  1.987,  1.978,  1.978/)
        p2r(:,16)=(/ 1.113,  1.019,  1.040,  1.258,  1.398,  1.474,  1.604,  1.686,  1.882,  1.766,  1.756,  1.759/)
        p2r(:,17)=(/ 1.123,  1.206,  1.337,  1.524,  1.573,  1.558,  1.550,  1.407,  1.567,  1.660,  1.673,  1.687/)
        p2r(:,18)=(/ 0.943,  1.215,  1.490,  1.677,  1.459,  1.291,  1.161,  1.023,  1.110,  1.306,  1.370,  1.425/)

        p3r(:, 1)=(/ 0.448,  0.187, -0.035, -0.485, -0.565, -0.414,  0.045,  0.360,  2.115,  1.859,  1.852,  1.861/)
        p3r(:, 2)=(/ 0.737,  0.413,  0.061, -0.658, -0.809, -0.647,  0.268,  0.527,  1.914,  1.613,  1.633,  1.676/)
        p3r(:, 3)=(/ 1.398,  1.610,  2.634, -2.644, -2.495, -2.476,  3.027,  3.027,  2.738,  2.250,  2.216,  2.238/)
        p3r(:, 4)=(/ 1.887,  2.285,  2.792, -2.928, -2.646, -2.635, -3.067,  2.843,  2.545,  2.115,  2.089,  2.089/)
        p3r(:, 5)=(/ 2.027,  2.352,  2.725, -3.046, -2.612, -2.314, -2.126, -2.121, -2.048, -1.658, -1.688, -1.688/)
        p3r(:, 6)=(/ 2.861, -2.965, -2.795, -2.678, -2.540, -1.907, -1.385, -1.339, -1.487, -1.397, -1.403, -1.403/)
        p3r(:, 7)=(/-2.067, -1.755, -1.456, -0.213,  0.761,  0.865,  0.864,  0.363, -1.469, -1.543, -1.543, -1.543/)
        p3r(:, 8)=(/ 2.714, -1.745, -1.472, -0.564,  1.315,  1.389,  1.420,  1.417,  1.636,  1.689,  1.689,  1.689/)
        p3r(:, 9)=(/ 2.188, -3.141, -1.970, -2.523,  2.400,  2.127,  2.204,  2.330,  2.112,  1.834,  1.834,  1.834/)
        p3r(:,10)=(/ 1.937,  2.421, -3.100, -2.837,  2.314,  2.349,  2.700,  2.701,  2.348,  1.823,  1.823,  1.823/)
        p3r(:,11)=(/ 1.582,  2.197,  2.809, -2.635, -2.687, -2.845, -2.902, -2.912, -1.696, -0.804, -0.804, -0.804/)
        p3r(:,12)=(/ 0.977,  1.943,  2.775, -2.766, -2.609, -2.459, -1.773, -1.823, -1.656, -1.346, -1.346, -1.346/)
        p3r(:,13)=(/ 0.678,  1.592,  2.017,  2.682, -3.079, -2.848, -2.254, -2.171, -1.930, -1.349, -1.349, -1.349/)
        p3r(:,14)=(/ 0.591,  1.260,  1.521,  1.800,  2.296,  2.699, -2.759, -2.599, -2.212, -1.564, -1.534, -1.534/)
        p3r(:,15)=(/ 0.356,  0.909,  1.232,  1.488,  1.817,  2.057,  3.141, -2.852, -1.689, -1.104, -1.054, -1.054/)
        p3r(:,16)=(/ 0.306,  0.897,  1.089,  1.094,  1.102,  0.988, -0.897, -0.704, -0.600, -0.877, -0.888, -0.895/)
        p3r(:,17)=(/ 0.568,  1.026,  1.309,  1.313,  0.804,  0.214, -0.353, -0.250, -0.277, -0.392, -0.400, -0.412/)
        p3r(:,18)=(/ 0.989,  1.509,  1.872,  2.185,  0.844,  0.276,  0.081,  0.081, -0.090, -0.413, -0.502, -0.572/)

! end of definition for several variables

        ! in the original program by bernd funke here ap index is read
        ! Here we get Ap from import_ts
        
        ! start calculation here
        ! Loop top_levels
        DO i=1,top_levels_in
          mask=(/ .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE., .TRUE./)
          ind1=minloc(abs(pref-press(i)),DIM=1,MASK=mask)
          mask(ind1)=.FALSE.
          pref1=pref(ind1)
          ind2=MAX(ind1-1,1)
          IF (pref1 < press(i)) THEN
            ind2=MIN(ind1+1,12)
          END IF
          pref2=pref(ind2)
          tmnr1=tmnr(ind1)
          tmnr2=tmnr(ind2)
          CALL ubcnox_linear_interp(tn,press(i),pref1,pref2,tmnr1,tmnr2)
          CALL ubcnox_linear_interp(ts,press(i),pref1,pref2,tmsr(ind1),tmsr(ind2))
          CALL ubcnox_linear_interp(nn,press(i),pref1,pref2,nmnr(ind1),nmnr(ind2))
          CALL ubcnox_linear_interp(ns,press(i),pref1,pref2,nmsr(ind1),nmsr(ind2))
          CALL ubcnox_linear_interp(wn,press(i),pref1,pref2,wmnr(ind1),wmnr(ind2))
          CALL ubcnox_linear_interp(ws,press(i),pref1,pref2,wmsr(ind1),wmsr(ind2))
          seasn=4._dp*nn*exp(-wn*(dn-tn))/(1._dp+exp(-wn*(dn-tn)))**2._dp
          seass=4._dp*ns*exp(-ws*(ds-ts))/(1._dp+exp(-ws*(ds-ts)))**2._dp

          
          ! filter function for calculation of weighted Ap 
          CALL ubcnox_linear_interp(txn,press(i),pref1,pref2,txnr(ind1),txnr(ind2))
          CALL ubcnox_linear_interp(txs,press(i),pref1,pref2,txsr(ind1),txsr(ind2))
          filtern(1:250)=sqrt(1./(dl(1:250)*dl(1:250)*dl(1:250)))*exp(-(dl(1:250)-txn)*(dl(1:250)-txn)/(2*(sqrt(0.7*txn)+6.)**2.*dl(1:250)/txn))
          filters(1:250)=sqrt(1./(dl(1:250)*dl(1:250)*dl(1:250)))*exp(-(dl(1:250)-txs)*(dl(1:250)-txs)/(2*(sqrt(0.7*txs)+6.)**2.*dl(1:250)/txs))
          filtern=filtern/sum(filtern(:))
          filters=filters/sum(filters(:))
          ! flip indices for filtern and filters
          filtern_flip(1:250)=filtern(250:1:-1)
          filters_flip(1:250)=filters(250:1:-1)
          xn=sum(ap(1:250)*filtern_flip(1:250))*seasn
          xs=sum(ap(1:250)*filters_flip(1:250))*seass
          
          IF (lat<0._dp) THEN
            ep(i)=xs*fac(i)
          END IF
          IF (lat>=0._dp) THEN
            ep(i)=xn*fac(i)
          END IF
          
          ! calculate background NOy
          indlat1=minloc(abs(lref-lat),DIM=1)
          IF (lat<lref(indlat1)) THEN
            indlat2=MAX(indlat1-1,1)
          END IF
          IF (lat>=lref(indlat1)) THEN
            indlat2=MIN(indlat1+1,18)
          END IF
          lat1=lref(indlat1)
          lat2=lref(indlat2)
          
          DO j=1,18
            CALL ubcnox_linear_interp(am(j),press(i),pref1,pref2,amr(ind1,j),amr(ind2,j))
            CALL ubcnox_linear_interp(a1(j),press(i),pref1,pref2,a1r(ind1,j),a1r(ind2,j))
            CALL ubcnox_linear_interp(a2(j),press(i),pref1,pref2,a2r(ind1,j),a2r(ind2,j))
            CALL ubcnox_linear_interp(a3(j),press(i),pref1,pref2,a3r(ind1,j),a3r(ind2,j))
            CALL ubcnox_linear_interp(p1(j),press(i),pref1,pref2,p1r(ind1,j),p1r(ind2,j))
            CALL ubcnox_linear_interp(p2(j),press(i),pref1,pref2,p2r(ind1,j),p2r(ind2,j))
            CALL ubcnox_linear_interp(p3(j),press(i),pref1,pref2,p3r(ind1,j),p3r(ind2,j))
            bgh(j)=am(j)*(1._dp+a1(j)*sin(ds/365._dp*2._dp*pi+p1(j))+a2(j)*sin(ds/365._dp*4._dp*pi+p2(j))+a3(j)*sin(ds/365._dp*6._dp*pi+p3(j)))
          END DO
          CALL ubcnox_linear_interp(bg(i),lat,lat1,lat2,bgh(indlat1),bgh(indlat2))
          IF (lat>lref(18)) THEN
            bg(i)=bgh(18)
          END IF
          IF (lat<lref(1)) THEN
            bg(i)=bgh(1)
          END IF
          IF ((lat<=5._dp) .and. (lat>=0._dp)) THEN
            bg(i)=bgh(10)
          END IF
          IF ((lat>=-5._dp) .and. (lat<0._dp)) THEN
            bg(i)=bgh(9)
          END IF
          
          
          
        END DO ! top_levels
        
        IF (ese_switch==1) THEN
          DO i=1,top_levels_in
                ese(i)=0._dp;
! ka_sv_20170523+
                startdateofese=dn-daysofese+1
! ka_sv_20170523-
                lp=log(press(i))                                                                            ! log pressure levels
                tm = 62.7637+23.3374*lp+3.34175*lp*lp+0.2589*lp*lp*lp+0.0106088*lp*lp*lp*lp               !% vertical time lag variation
! ka_sv_20170523+
!bugfix: dn -> startdateofese
                tm = MIN(tm + exp((tm+startdateofese-279.)/4.),270.)                                           !% variation at equinox transition
                xu = 4._dp*exp(-0.046*(startdateofese-173.))/(1.+exp(-0.046*(startdateofese-173.)))**2._dp*0.0075                !% seasonal dependence of amount at source region
                wu = 4._dp*exp(-0.043*(startdateofese-173.))/(1.+exp(-0.043*(startdateofese-173.)))**2._dp*1.25                  !% seasonal dependence of ESE wbar at source region
                fm = 0.357087-0.239236*lp+0.00420932*lp*lp+0.0105685*lp*lp*lp+0.00107633*lp*lp*lp*lp      !% vertical flux variation
                fm = MAX(fm/(1.+exp((tm+startdateofese-273.)/8.))*xu*wu,0.)                                    !% scale with source region amount*wbar, consider equinox transition
                wm = exp(-1.69674-0.493714*lp+0.151089*lp*lp+0.00082302*lp*lp*lp-0.0139315*lp*lp*lp*lp-0.000871843*lp*lp*lp*lp*lp+0.000161791*lp*lp*lp*lp*lp*lp)  !% vertical wbar variation
                wm = wm/(1.+exp((tm+startdateofese-280.)/9.))*wu !% scale with source region wbar, consider equinox transition
                !xb=4._dp*nn*exp(-wn*(dn+int(tm)-tn))/(1._dp+exp(-wn*(dn+int(tm)-tn)))**2._dp
                xb=4._dp*nn*exp(-wn*(startdateofese+int(tm)-tn))/(1._dp+exp(-wn*(startdateofese+int(tm)-tn)))**2._dp
! ka_sv_20170523-
                nne = fm/wm-xb
                we = 0.15

                if (daysofese < tm) THEN
                  rfac = ((daysofese)/tm)**0.3
                else
                  rfac = 1.    !%fade in after ESE onset
                end if
                if (dn > 304. .AND. dn < 324.) THEN
                  rfac = rfac*((324.-dn)/20.)**0.5         !%fade out after 1st May
                end if
! ka_sv_20170523+
!bugfix: dn -> startdateofese
                tl = startdateofese+tm
! ka_sv_20170523-
                sease = 4.*nne*exp(-we*(dn-tl))/(1.+exp(-we*(dn-tl)))**2._dp*rfac
                xe=sum(ap(1:250)*filtern_flip(1:250))*sease
                ese(i)=xe*face(i)
           end do   ! top_levels_in
        END IF   ! ese_switch==1
        
        IF (ese_switch==1) THEN
! ka_sv_20170523+
          IF (daysofese>0._dp) THEN
            ubc=ep+bg+ese
          ELSE
            ubc=ep+bg
          END IF
! ka_sv_20170523-
        ELSE
          ubc=ep+bg
        END IF

END SUBROUTINE ubcnox_calc_current_ubc

! ----------------------------------------------------------------------


SUBROUTINE ubcnox_linear_interp(interp,xq,x1,x2,y1,y2)
! linear interpolation between two points
   REAL(dp), INTENT(out) :: interp         ! interpolated value
   REAL(dp), INTENT(in) :: xq              ! x value for interpolated value
   REAL(dp), INTENT(in) :: x1,x2,y1,y2     ! two points for interpolation
   
   REAL(dp) :: dx                          ! total distance for x
   REAL(dp) :: w1, w2                      ! weights for point1 and point2
   
   dx=x2-x1
   IF (abs(dx)<1.e-20) THEN
     interp=y1
   ELSE
     w1=1._dp-(xq-x1)/dx
     w2=1._dp-(x2-xq)/dx
     interp=w1*y1+w2*y2
   END IF
   
END SUBROUTINE


END MODULE messy_ubcnox
