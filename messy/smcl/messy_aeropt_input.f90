MODULE MESSY_AEROPT_INPUT

  USE MESSY_AEROPT_MEM

  IMPLICIT NONE

  INTEGER :: i

CONTAINS
! subroutines to build the lookup tables including a NETCDF reading routine, 
! memory allocation etc.

! can and will be replaced if the files can be accessed using the 
! generic input module


!===============================================================================
 SUBROUTINE lookup_initialize_rad0(L_IO, idx)

   !---------------------------------------------------------------------
   ! *** lookup_INITIALIZE_RAD0:
   !
   !     Initialize all variables and look-up-tables required by
   !     the ECHAM5 radiation module (sw+lw) to calculate aerosol forcing:
   !
   !     Open 2 (sw/lw) NetCDF-files containing look-up-tables and determine
   !     size of fields for memory allocation.
   ! ----------------------------------------------------------------------

   USE netcdf

   IMPLICIT NONE

   LOGICAL, INTENT(IN) :: L_IO
   INTEGER, INTENT(IN) :: idx

   ! IDs of dimensions
   INTEGER :: nrDimID(nmodrad), niDimID(nmodrad), mspDimID(nmodrad)  
   INTEGER :: lambdaDimID, aerosolDimID, e5bandDimID
   INTEGER :: n             ! loop counter
   INTEGER :: status        ! status of NetCDF-subroutines
   INTEGER :: num_aerosols2 ! number of aerosol species in arrays (lw)
   INTEGER :: num_e5_bands  ! number of spectral bands (from NetCDF-data)
   INTEGER :: band

   CHARACTER (len=16)            :: varname ! name of NetCDF-variable
   CHARACTER (len=NF90_MAX_NAME) :: dimname ! name of NetCDF-dimension

   CHARACTER (LEN=*), PARAMETER  :: subroutine_name = 'lookup_initialize_rad0'

   ! --- code starts ---
   
   IF(l_io) THEN
     ! print filename of shortwave look-up table
    WRITE(*,*) "rad_sw_filename =", trim(rad_sw_filename)
     ! print filename of longwave look-up table
     WRITE(*,*) "rad_lw_filename =", trim(rad_lw_filename)
   END IF ! l_io

     ! * * *   s h o r t w a v e   * * *
     
     band = 1
     
     ! 1. Read parameter of sw look-up tables from NetCDF-file.
     
     ! 1.1 Open NetCDF-file containing look-up tables of aerosol optical
     !     properties and inquire dimensions of fields.
     
     !     open NetCDF-file (read only)
     
     status = NF90_OPEN(trim(rad_sw_filename), NF90_NoWrite, ncID_sw)
     IF (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
     
     ! 1.2 Inquire dimensions.
     
     !     number of spectral sub-intervals
     
     status = NF90_INQ_DIMID(ncID_sw, 'lambda', lambdaDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_sw, lambdaDimID, dimname, &
       num_sw_intervals)
     
     !     number of aerosol species
     
     status = NF90_INQ_DIMID(ncID_sw, 'aerosol', aerosolDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_sw, aerosolDimID, dimname, num_aerosols)
     
     
     !     number of ECHAM5 spectral bands (shortwave)
     
     status = NF90_INQ_DIMID(ncID_sw, 'e5band', e5bandDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_sw, e5bandDimID, dimname, num_e5_bands)
     
     ! 1.3 Check for consistency.
     
     IF (num_sw_intervals > max_sw_intervals) THEN
       IF(l_io) &
         WRITE(*,*) trim(subroutine_name),' error: number of spectral bands',&
         ' in NetCDF is greater than maximum', &
         ' number currently allowed (shortwave)'
       STATUS = 1
       RETURN
     END IF
     
     IF (num_aerosols > max_aerosols) THEN
       IF(l_io) WRITE(*,*) trim(subroutine_name),&
         ' error: number of aerosol species',&
         ' in NetCDF is greater than maximum', &
         ' number currently allowed (shortwave)'
       STATUS = 1
       RETURN
     END IF
     
     IF (num_e5_bands /= nsw) THEN
       IF(l_io) &
         WRITE(*,*) trim(subroutine_name),' error: number of spectral bands',&
         ' in NetCDF file does not match',  &
         ' number of ECHAM5 spectral bands',&
         ' (shortwave)'
       STATUS = 1
       RETURN
     END IF
     
     DO n = 1, nmodrad
       
       !  dimension "nr"
       
       WRITE (varname,'(A,I1)') 'nr', n
       status = NF90_INQ_DIMID(ncID_sw, trim(varname), nrDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_sw, nrDimID(n), dimname, &
         idim_nr(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
       !  dimension "ni"
       
       WRITE (varname,'(A,I1)') 'ni', n
       status = NF90_INQ_DIMID(ncID_sw, trim(varname), niDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_sw, niDimID(n), dimname, &
         idim_ni(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
       !  dimension "msp"
       
       WRITE (varname,'(A,I1)') 'msp', n
       status = NF90_INQ_DIMID(ncID_sw, trim(varname), mspDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_sw, mspDimID(n), dimname, &
         idim_msp(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
     END DO
     
     ! Do not close NetCDF-file. This will be done by lookup_initialize_rad2.
     
     ! * * *   l o n g w a v e   * * *
     
     band = 2
     
     ! 2. Read parameters of lw look-up tables from NetCDF-file.
     
     ! 2.1 Open NetCDF-file containing look-up tables of aerosol optical
     !     properties and inquire dimensions of fields.
     
     !     open NetCDF-file (read only)
     
     status = NF90_OPEN(trim(rad_lw_filename), NF90_NoWrite, ncID_lw)
     IF (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
     
     ! 2.2 Inquire dimensions.
     
     !     number of spectral sub-intervals
     
     status = NF90_INQ_DIMID(ncID_lw, 'lambda', lambdaDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_lw, lambdaDimID, dimname, &
       num_lw_intervals)
     
     !     number of aerosol species
     
     status = NF90_INQ_DIMID(ncID_lw, 'aerosol', aerosolDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_lw, aerosolDimID, dimname, &
       num_aerosols2)
     
     !     number of ECHAM5 spectral bands (longwave)
     
     status = NF90_INQ_DIMID(ncID_lw, 'e5band', e5bandDimID)
     IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
     status = NF90_INQUIRE_DIMENSION(ncID_lw, e5bandDimID, dimname, &
       num_e5_bands)
     
     ! 2.3 Check for consistency.
     
     IF (num_lw_intervals > max_lw_intervals) THEN
       IF(l_io) &
         WRITE(*,*) trim(subroutine_name),' error: number of spectral bands',&
         ' in NetCDF is greater than maximum',  &
         ' number currently allowed (longwave)'
       STATUS = 1
       RETURN
     END IF
     
     IF (num_aerosols > num_aerosols2) THEN
       IF(l_io) &
         WRITE(*,*) trim(subroutine_name),' error: number of aerosol types', &
         ' in (longwave) NetCDF file does not', &
         ' match shortwave NetCDF file'
       STATUS = 1
       RETURN
     END IF
     
     IF (num_e5_bands /= jpband) THEN
       IF(l_io) &
         WRITE(*,*) trim(subroutine_name),' error: number of spectral bands',&
         ' in NetCDF file does not match',  &
         ' number of ECHAM5 spectral bands', &
         '(longwave)'
       STATUS = 1
       RETURN
     END IF
     
     DO n = 1, nmodrad
       
       !  dimension "nr"
       
       WRITE (varname,'(A,I1)') 'nr', n
       status = NF90_INQ_DIMID(ncID_lw, trim(varname), nrDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_lw, nrDimID(n), dimname, &
         idim_nr(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
       !  dimension "ni"
       
       WRITE (varname,'(A,I1)') 'ni', n
       status = NF90_INQ_DIMID(ncID_lw, trim(varname), niDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_lw, niDimID(n), dimname, &
         idim_ni(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
       !  dimension "msp"
       
       WRITE (varname,'(A,I1)') 'msp', n
       status = NF90_INQ_DIMID(ncID_lw, trim(varname), mspDimID(n))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       status = NF90_INQUIRE_DIMENSION(ncID_lw, mspDimID(n), dimname, &
         idim_msp(band))
       IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
       
     END DO
     
     ! Do no close NetCDF-file. This will be done by lookup_initialize_rad2
     
   END SUBROUTINE lookup_initialize_rad0

!============================================================================

   SUBROUTINE lookup_initialize_rad1(idx)


   ! *** lookup_INITIALIZE_RAD1:
   !
   !     Allocate memory for look-up-tables required by the ECHAM5 radiation
   !     module (sw+lw) to calculate aerosol forcing.
   ! ----------------------------------------------------------------------

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: idx
     INTEGER             :: v1, v2, v3, i

     ! Allocate memory for look-up table of aerosol optical properties.

     ! *** shortwave ***

     DO i=1,2
       v1 = idim_nr(i)
       v2 = idim_ni(i)
       v3 = idim_msp(i)
       SELECT CASE (i)
       CASE (1)
         ALLOCATE(tab_set(idx)%lut_sw_sigma(v1,v2,v3,nmodrad))
         ALLOCATE(tab_set(idx)%lut_sw_omega(v1,v2,v3,nmodrad))
         ALLOCATE(tab_set(idx)%lut_sw_gamma(v1,v2,v3,nmodrad))

         ALLOCATE(tab_set(idx)%ref_re_sw(num_aerosols,num_sw_intervals))
         ALLOCATE(tab_set(idx)%ref_im_sw(num_aerosols,num_sw_intervals))
         
       CASE(2)
         ALLOCATE(tab_set(idx)%lut_lw_sigma(v1,v2,v3,nmodrad))

         ALLOCATE(tab_set(idx)%ref_re_lw(num_aerosols,num_lw_intervals))
         ALLOCATE(tab_set(idx)%ref_im_lw(num_aerosols,num_lw_intervals))

       END SELECT
     END DO

   END SUBROUTINE lookup_initialize_rad1

!=============================================================================

   SUBROUTINE lookup_initialize_rad2(L_IO)

   !---------------------------------------------------------------------
   ! *** lookup_INITIALIZE_RAD2:
   !
   !     Read the look-up tables for the aerosol optical properties
   !     from 2 (sw/lw) NetCDF-files.
   ! ----------------------------------------------------------------------

   USE netcdf

   IMPLICIT NONE

   INTRINSIC ASSOCIATED

   LOGICAL, INTENT(IN) :: L_IO
   ! IDs of vars
   INTEGER :: nrVarID(nmodrad), niVarID(nmodrad), mspVarID(nmodrad)  
   INTEGER :: sigmaVarID(nmodrad), omegaVarID(nmodrad), gammaVarID(nmodrad)
   INTEGER :: lambdaVarID=0, weightVarID=0, ref_reVarID=0, ref_imVarID=0
   INTEGER :: numVarID=0
   INTEGER :: n      ! loop counter
   INTEGER :: status ! status of NetCDF-subroutines
   INTEGER :: band

   REAL(dp), DIMENSION(nsw)        :: tmp_sw
   REAL(dp), DIMENSION(jpband)     :: tmp_lw

   CHARACTER (len=16)            :: varname ! name of NetCDF-variable

   CHARACTER (LEN=*), PARAMETER  :: subroutine_name = 'lookup_initialize_rad2'

   REAL(dp), DIMENSION(1) :: aux

   ! look-up table (temporary)
   REAL(dp), DIMENSION(:,:,:), POINTER   :: lut_sw_tmp   => NULL()
   REAL(dp), DIMENSION(:,:,:), POINTER   :: lut_lw_tmp   => NULL()

   ! --- code starts ---

   ! * * *   s h o r t w a v e   * * *

   band = 1

   ! 1. Read sw look-up tables from NetCDF-file. The file has already been
   !    opened by lookup_initialize_rad0.

   ! 1.1 Inquire variables.

   DO n = 1, nmodrad

      !  dimension "nr"

     WRITE (varname,'(A,I1)') 'nr', n

     status = NF90_INQ_VARID(ncID_sw, trim(varname), nrVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_sw, nrVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      nr_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_sw, nrVarID(n), aux, start = (/idim_nr(band)/),&
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      nr_max(band,n) = aux(1)
      
      !  dimension "ni"
      
      WRITE (varname,'(A,I1)') 'ni', n
      
      status = NF90_INQ_VARID(ncID_sw, trim(varname), niVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_sw, niVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ni_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_sw, niVarID(n), aux, start = (/idim_ni(band)/),&
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ni_max(band,n) = aux(1)
      
      !  dimension "msp"
      
      WRITE (varname,'(A,I1)') 'msp', n

      status = NF90_INQ_VARID(ncID_sw, trim(varname), mspVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_sw, mspVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      msp_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_sw, mspVarID(n), aux, start = &
        (/idim_msp(band)/), count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      msp_max(band,n) = aux(1)
      
      !  variable "sigma"
      
      WRITE (varname,'(A,I1)') 'sigma_', n
      status = NF90_INQ_VARID(ncID_sw, TRIM(varname), sigmaVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ! get attribute 'sigma_g' (standard deviation of log-normal distribution)
      status = NF90_GET_ATT(ncID_sw, sigmaVarID(n), 'sigma_g', sigma_g(band,n))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      
      !  variable "omega"
      
      WRITE (varname,'(A,I1)') 'omega_', n
      status = NF90_INQ_VARID(ncID_sw, TRIM(varname), omegaVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      !  variable "gamma"
      
      WRITE (varname,'(A,I1)') 'gamma_', n
      status = NF90_INQ_VARID(ncID_sw, TRIM(varname), gammaVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
    END DO
    
    !     variable "lambda"
    
    status = NF90_INQ_VARID(ncID_sw, "lambda", lambdaVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "weight"
    
    status = NF90_INQ_VARID(ncID_sw, "weight", weightVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "ri_r"
    
    status = NF90_INQ_VARID(ncID_sw, "ri_r", ref_reVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "ri_i"
    
    status = NF90_INQ_VARID(ncID_sw, "ri_i", ref_imVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "num_intervals"
    
    status = NF90_INQ_VARID(ncID_sw, "num_intervals", numVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    ! 1.2 Allocate memory for temporary variables.
    
    ! lut (temporary)
    ALLOCATE(lut_sw_tmp(idim_nr(band),idim_ni(band),idim_msp(band)))
    
    ! 1.3 Read aerosol optical properties from NetCDF-file for each mode.
    
    DO n = 1, nmodrad
      
      !  "sigma"
      status = NF90_GET_VAR(ncID_sw, sigmaVarID(n), lut_sw_tmp(:,:,:))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      lut_sw_sigma(:,:,:,n) = lut_sw_tmp(:,:,:)
      
      !  "omega"
      status = NF90_GET_VAR(ncID_sw, omegaVarID(n), lut_sw_tmp(:,:,:))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      lut_sw_omega(:,:,:,n) = lut_sw_tmp(:,:,:)
      
      !  "gamma"
      status = NF90_GET_VAR(ncID_sw, gammaVarID(n), lut_sw_tmp(:,:,:))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      lut_sw_gamma(:,:,:,n) = lut_sw_tmp(:,:,:)
      
    END DO
    
    ! 1.4 Read other variables.
    
    !     "lambda"
    status = NF90_GET_VAR(ncID_sw, lambdaVarID, lambda_sw, start = (/1/), &
      count = (/num_sw_intervals/))
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)

    !     "weight"
    status = NF90_GET_VAR(ncID_sw, weightVarID, weight_sw, start = (/1/), &
      count = (/num_sw_intervals/))
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "ri_r"
    status = NF90_GET_VAR(ncID_sw, ref_reVarID, ref_re_sw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "ri_i"
    status = NF90_GET_VAR(ncID_sw, ref_imVarID, ref_im_sw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "num_intervals"
    status = NF90_GET_VAR(ncID_sw, numVarID, tmp_sw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    DO n = 1, nsw
      sw_intervals(n) = NINT(tmp_sw(n))
    END DO

    ! 1.5 Close NetCDF-file which has been opened by lookup_initialize_rad0.

    status = NF90_CLOSE(ncID_sw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)

    ! * * *   l o n g w a v e   * * *

    band = 2

    ! 2. Read lw look-up tables from NetCDF-file. The file has already been
    !    opened by lookup_initialize_rad0.
    
    ! 2.1 Inquire variables.
    
    DO n = 1, nmodrad
     
      !  dimension "nr"

      WRITE (varname,'(A,I1)') 'nr', n
      
     
      status = NF90_INQ_VARID(ncID_lw, trim(varname), nrVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_lw, nrVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      nr_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_lw, nrVarID(n), aux, start = (/idim_nr(band)/),&
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      nr_max(band,n) = aux(1)
      
      !  dimension "ni"
      
      WRITE (varname,'(A,I1)') 'ni', n
      
      status = NF90_INQ_VARID(ncID_lw, trim(varname), niVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_lw, niVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ni_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_lw, niVarID(n), aux, start = (/idim_ni(band)/),&
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ni_max(band,n) = aux(1)
      
      !  dimension "msp"
      
      WRITE (varname,'(A,I1)') 'msp', n
      
      status = NF90_INQ_VARID(ncID_lw, trim(varname), mspVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      
      status = NF90_GET_VAR(ncID_lw, mspVarID(n), aux, start = (/1/), &
        count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      msp_min(band,n) = aux(1)
      
      status = NF90_GET_VAR(ncID_lw, mspVarID(n), aux, start = &
        (/idim_msp(band)/), count = (/1/))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      msp_max(band,n) = aux(1)
      
      !  variable "sigma"
      
      WRITE (varname,'(A,I1)') 'sigma_', n
      status = NF90_INQ_VARID(ncID_lw, TRIM(varname), sigmaVarID(n))
      IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
      ! get attribute 'sigma_g' (standard deviation of log-normal distribution)
      status = NF90_GET_ATT(ncId_lw, sigmaVarID(n), 'sigma_g', sigma_g(band,n))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      
    END DO
    
    ! Check for consistency of pre-calculated table with chosen standard
    ! deviation of aerosol modes.
    
    !     variable "lambda"
    
    status = NF90_INQ_VARID(ncID_lw, "lambda", lambdaVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "weight"
    
    status = NF90_INQ_VARID(ncID_lw, "weight", weightVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "ri_r"
    
    status = NF90_INQ_VARID(ncID_lw, "ri_r", ref_reVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "ri_i"
    
    status = NF90_INQ_VARID(ncID_lw, "ri_i", ref_imVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    !     variable "num_intervals"
    
    status = NF90_INQ_VARID(ncID_lw, "num_intervals", numVarID)
    IF (status /= nf90_NoErr) CALL handle_err(subroutine_name, status, l_io)
    
    ! 2.2 Allocate memory for temporary variables.
    
    ! lut (temporary)
    ALLOCATE(lut_lw_tmp(idim_nr(band),idim_ni(band),idim_msp(band)))
    
    ! 2.3 Read aerosol optical properties from NetCDF-file for each mode.
    
    DO n = 1, nmodrad
      
      !  "sigma"
      status = NF90_GET_VAR(ncID_lw, sigmaVarID(n), lut_lw_tmp(:,:,:))
      if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
      lut_lw_sigma(:,:,:,n) = lut_lw_tmp(:,:,:)
      
    END DO
    
    ! 2.4 Read other variables.
    
    !     "lambda"
    status = NF90_GET_VAR(ncID_lw, lambdaVarID, lambda_lw, start = (/1/), &
      count = (/num_lw_intervals/))
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "weight"
    status = NF90_GET_VAR(ncID_lw, weightVarID, weight_lw, start = (/1/), &
      count = (/num_lw_intervals/))
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "ri_r"
    !   status = NF90_GET_VAR(ncID_lw, ref_reVarID, ref_re_lw, start = (/1,1/), &
    !                         count = (/num_aerosols,num_lw_intervals/))
    status = NF90_GET_VAR(ncID_lw, ref_reVarID, ref_re_lw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    !     "ri_i"
    !   status = NF90_GET_VAR(ncID_lw, ref_imVarID, ref_im_lw, start = (/1,1/), &
    !                         count = (/num_aerosols,num_lw_intervals/))
    status = NF90_GET_VAR(ncID_lw, ref_imVarID, ref_im_lw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    

    !     "num_intervals"
    status = NF90_GET_VAR(ncID_lw, numVarID, tmp_lw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    DO n = 1, jpband
      lw_intervals(n) = NINT(tmp_lw(n))
    END DO

    ! 2.5 Close NetCDF-file, which has been opened by lookup_initialize_rad0.

    status = NF90_CLOSE(ncID_lw)
    if (status /= nf90_NoErr) call handle_err(subroutine_name, status, l_io)
    
    ! 3. Clean up temporary arrays.

    IF (ASSOCIATED(lut_sw_tmp)) DEALLOCATE(lut_sw_tmp)
    IF (ASSOCIATED(lut_lw_tmp)) DEALLOCATE(lut_lw_tmp)
    
  END SUBROUTINE lookup_initialize_rad2
  
!=============================================================================

  SUBROUTINE lookup_initialize_rad3(L_IO)

   !---------------------------------------------------------------------
   ! *** lookup_INITIALIZE_RAD3:
   !
   !     Finish initialization of look-up tables for the aerosol optical
   !     properties.
   ! ----------------------------------------------------------------------

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: L_IO

    INTEGER :: n
    INTEGER :: band

    CHARACTER (LEN=*), PARAMETER  :: subroutine_name = 'lookup_initialize_rad3'

    
    ! 2. Calculate variables used to speed up table look-up.
    
    DO band = 1, 2
      DO n = 1, nmodrad
        nr_step(band,n)     = (nr_max(band,n) - nr_min(band,n)) &
          / REAL(idim_nr(band) - 1)
        log_ni_min(band,n)  = LOG(ni_min(band,n))
        ni_step(band,n)     = (LOG(ni_max(band,n)) - log_ni_min(band,n)) &
          / REAL(idim_ni(band) - 1)
        log_msp_min(band,n) = LOG(msp_min(band,n))
        msp_step(band,n)    = (LOG(msp_max(band,n)) - log_msp_min(band,n)) &
          / REAL(idim_msp(band) - 1)
      END DO
    END DO
   
  END SUBROUTINE lookup_initialize_rad3

!===============================================================================

  SUBROUTINE handle_err ( name, status, L_IO)

    USE netcdf

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(in) :: name
    INTEGER, INTENT(in) :: status
    LOGICAL, INTENT(in) :: L_IO
    
    IF(l_io) &
      WRITE(*,*) trim(name), 'NetCDF error: ', trim(NF90_STRERROR(status))
    
  END SUBROUTINE handle_err
  

  ! -------------------------------------------------------------------------

  SUBROUTINE handle_table_err (name, filename, sigma_g, i_mod, L_IO, i)
    
    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER (LEN=*),          INTENT(in) :: name, filename
    REAL(dp), INTENT(in) :: sigma_g
    INTEGER,  INTENT(in) :: i_mod
    LOGICAL, INTENT(in) :: L_IO
    INTEGER, INTENT(IN) :: I

    CHARACTER(LEN=12) :: size

    IF(l_io) THEN
      size = 'ERROR'
      SELECT CASE (radmod(i_mod))
      CASE (1)
        size = 'Aitken'
      CASE (2)
        size = 'accumulation'
      CASE (3)
        size = 'coarse'
      END SELECT
      WRITE(*,*)
      WRITE(*,*) trim(name),' - CAUTION: CRITICAL WARNING:'
      WRITE(*,*)
      WRITE(*,*) 'The '//TRIM(size)//' mode standard deviation in the'
      WRITE(*,*) 'look-up table for the calculation of the aerosol'
      WRITE(*,*) 'optical properties and the geometric standard'
      WRITE(*,'(A,I1,A)') 'deviation of mode', i_mod, 'do not match.'
      WRITE(*,*)
      WRITE(*,*) 'Look-up table (',trim(filename),'):'
      WRITE(*,'(1X,A,I1,A,F4.2)')     ' sigma_g(', radmod(i_mod), ') = ', &
        sigma_g
      WRITE(*,*)
      WRITE(*,*) 'aerosol submodel sigma definition:'
      WRITE(*,'(1X,A,I1,F4.2)') ' sigma(', i_mod, ') = ', sigma(i_mod)
      WRITE(*,*)
    END IF
    
  END SUBROUTINE handle_table_err

!==============================================================================

  SUBROUTINE aeropt_check_lut(opts, l_io)

    USE MESSY_AEROPT_SETS,          ONLY: MAP_OPTSET_STRUCT, MAP_LUT_STRUCT

    IMPLICIT NONE

    TYPE(t_aero_set), POINTER, INTENT(IN) :: opts
    LOGICAL, INTENT(in) :: L_IO

    CHARACTER(LEN=256)  :: name
    INTEGER             :: band, n
    CHARACTER (LEN=*), PARAMETER  :: subroutine_name = 'check_lookuptables'
    

    CALL MAP_OPTSET_STRUCT(opts)
    CALL MAP_LUT_STRUCT(lut_number)
    
    ! Check for consistency of pre-calculated table with chosen standard
    ! deviation of aerosol modes.
    
    DO band=1,2
      SELECT CASE (BAND)
      CASE(1)
        name = rad_sw_filename
      CASE(2)
        name = rad_lw_filename
      END SELECT

      DO n = ns+1, nmod
        IF (ABS(sigma_g(band,radmod(n)) - sigma(n)) .GT. 1.0e-6_dp) THEN
          CALL handle_table_err(&
               subroutine_name, name, sigma_g(band,radmod(n)), n, l_io, i)
        END IF
      END DO
    END DO
  END SUBROUTINE aeropt_check_lut

!==============================================================================
  SUBROUTINE aeropt_finalise_lut(opts,l_io, status)

    USE MESSY_AEROPT_SETS,          ONLY: MAP_OPTSET_STRUCT, MAP_LUT_STRUCT

    IMPLICIT NONE

    TYPE (t_aero_set), POINTER, INTENT(IN) :: opts
    INTEGER, INTENT(INOUT) :: status
    LOGICAL, INTENT(IN)    :: L_IO

    INTEGER  :: num_wave, band
    INTEGER  :: i, j, k, l, n
    CHARACTER (LEN=*), PARAMETER  :: subroutine_name = 'finalise_lut'
    CHARACTER (LEN=1024)          :: error_str

    ! refractive index, real part
    REAL(dp), DIMENSION(max_aerosols,max_diag_wavelens) :: ref_re_opt
    ! refractive index, imaginary part
    REAL(dp), DIMENSION(max_aerosols,max_diag_wavelens) :: ref_im_opt

    ! refractive index, real part
    REAL(dp), DIMENSION(max_aerosols,n_jv_calc) :: ref_re_jval
    ! refractive index, imaginary part
    REAL(dp), DIMENSION(max_aerosols,n_jv_calc) :: ref_im_jval

    num_wave = max_sw_intervals+max_lw_intervals+max_diag_wavelens+n_jv_calc
    ALLOCATE(opts%lambda(num_wave))
    ALLOCATE(opts%lambda_squared(num_wave))
    ALLOCATE(opts%weight(num_wave))
    ALLOCATE(opts%ref_re(max_aerosols,num_wave))
    ALLOCATE(opts%ref_im(max_aerosols,num_wave))
    ALLOCATE(opts%int2band(num_wave))

    CALL MAP_OPTSET_STRUCT(opts)
    CALL MAP_LUT_STRUCT(opts%lut_number)

    status = 1

    IF (opts%num_opt_wavelens > 0) THEN

      ! linear interpolation

      DO i = 1, opts%num_opt_wavelens

        IF (rad_diag_wavelen(i) < lambda_sw(num_sw_intervals)) THEN
          DO j = 1, num_sw_intervals - 1
            IF ((lambda_sw(j) <= rad_diag_wavelen(i)).and. &
              (lambda_sw(j+1) >= rad_diag_wavelen(i))) EXIT
          END DO
          
          DO n = 1, num_aerosols
            ref_re_opt(n,i) = ref_re_sw(n,j)                         &
              + (rad_diag_wavelen(i) - lambda_sw(j)) &
              * (ref_re_sw(n,j+1) - ref_re_sw(n,j))  &
              / (lambda_sw(j+1) - lambda_sw(j))
            ref_im_opt(n,i) = ref_im_sw(n,j)                         &
              + (rad_diag_wavelen(i) - lambda_sw(j)) &
              * (ref_im_sw(n,j+1) - ref_im_sw(n,j))  &
              / (lambda_sw(j+1) - lambda_sw(j))
          END DO
        ELSE
          DO j = 1, num_lw_intervals - 1
            IF ((lambda_lw(j) <= rad_diag_wavelen(i)).and. &
              (lambda_lw(j+1) >= rad_diag_wavelen(i))) EXIT
          END DO
          
          DO n = 1, num_aerosols
            ref_re_opt(n,i) = ref_re_lw(n,j)                         &
              + (rad_diag_wavelen(i) - lambda_lw(j)) &
              * (ref_re_lw(n,j+1) - ref_re_lw(n,j))  &
              / (lambda_lw(j+1) - lambda_lw(j))
            ref_im_opt(n,i) = ref_im_lw(n,j)                         &
              + (rad_diag_wavelen(i) - lambda_lw(j)) &
              * (ref_im_lw(n,j+1) - ref_im_lw(n,j))  &
              / (lambda_lw(j+1) - lambda_lw(j))
          END DO
        END IF
      END DO

    END IF

    IF (n_jv_bands > 0) THEN

      ! linear interpolation

      DO i = 1, n_jv_calc

        IF ((jval_wavelen(i) < lambda_sw(1)).or. &
          (jval_wavelen(i) > lambda_sw(num_sw_intervals))) THEN
          IF(l_io) &
            WRITE (error_str, '(A,F5.3,A,F5.3,A)')                   &
            'Warning: wavelengths between ', lambda_sw(1),       &
            ' and ', lambda_sw(num_sw_intervals),              &
            ' micro-m for diagnostics supported only; choosing lowest'//&
           &' available wavelength for lower values (jval_wavelen).'
          IF(l_io) &
            WRITE(*,*) trim(subroutine_name), trim(error_str)
        END IF
        IF (jval_wavelen(i) < lambda_sw(1)) THEN
          ref_re_jval(1:num_aerosols,i) = ref_re_sw(1:num_aerosols,1)
          ref_im_jval(1:num_aerosols,i) = ref_im_sw(1:num_aerosols,1)     
        ELSE IF (jval_wavelen(i) > lambda_sw(num_sw_intervals)) THEN
          ref_re_jval(1:num_aerosols,i) = ref_re_sw(1:num_aerosols,num_sw_intervals)
          ref_im_jval(1:num_aerosols,i) = ref_im_sw(1:num_aerosols,num_sw_intervals)     
        ELSE
          DO j = 1, num_sw_intervals - 1
            IF ((lambda_sw(j) <= jval_wavelen(i)).and. &
              (lambda_sw(j+1) >= jval_wavelen(i))) EXIT
          END DO
          
          DO n = 1, num_aerosols
            ref_re_jval(n,i) = ref_re_sw(n,j)                         &
              + (jval_wavelen(i) - lambda_sw(j)) &
              * (ref_re_sw(n,j+1) - ref_re_sw(n,j))  &
              / (lambda_sw(j+1) - lambda_sw(j))
            ref_im_jval(n,i) = ref_im_sw(n,j)                         &
              + (jval_wavelen(i) - lambda_sw(j)) &
              * (ref_im_sw(n,j+1) - ref_im_sw(n,j))  &
              / (lambda_sw(j+1) - lambda_sw(j))
          END DO
        END IF

      END DO
      
    END IF

    ! 2.1 wavelength (center) of sub-intervals [m]:
    !     convert from [micro-m] to [m]
    
    DO n = 1, num_sw_intervals
      lambda(n) = lambda_sw(n) * 1.0e-6
    END DO
    DO n = 1, num_opt_wavelens  ! optional wavelengths (sw) for diagnostics
      lambda(n+num_sw_intervals) = rad_diag_wavelen(n) * 1.0e-6
    END DO
    DO n = 1, n_jv_bands  ! optional wavelengths (sw) for JVAL coupling
      lambda(n+num_sw_intervals+num_opt_wavelens) = jval_wavelen(n) * 1.0e-6
    END DO

    DO n = 1, num_lw_intervals
      lambda(n+num_sw_intervals+num_opt_wavelens+n_jv_bands) = lambda_lw(n) * 1.0e-6
    END DO
    
    DO n = 1, num_sw_intervals+num_opt_wavelens+n_jv_bands+num_lw_intervals
      lambda_squared(n) = lambda(n) * lambda(n)  ! [m2]
    END DO
    
    ! 2.2 indices for mapping sub-intervals to ECHAM5-bands
    
    l = 1
    DO n = 1, nsw
      DO k = 1, sw_intervals(n)
        int2band(l) = n
        l = l + 1
      END DO
    END DO
    
    DO n = 1, num_opt_wavelens
      int2band(l) = 0
      l = l + 1
    END DO

    DO n = 1, n_jv_bands
      int2band(l) = 0
      l = l + 1
    END DO
    
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ECHAM5 longwave bands are sorted by ascending wavenumber,
    ! _NOT_ as the shortwave bands by ascending wavelength!!!
    ! The wavelength intervals in the NetCDF-file are sorted
    ! by ascending wavelength
    !
    ! ---> swap order of assignment
    !      wavelength interval <---> ECHAM5 band!!!
    !      set assignment variable "int2band" accordingly
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    DO n = 1, jpband
      DO k = 1, lw_intervals(n)
        int2band(l) = jpband + 1 - n  ! reverse oder!!! see comment above!!!
        l = l + 1
      END DO
    END DO
    
    ! 2.3 weight of sub-intervals and sum for each ECHAM5 band
    
    do n = 1, nsw
      sumweight_sw(n) = 0.0_dp
    end do
    do n = 1, jpband
      sumweight_lw(n) = 0.0_dp
    end do
    
    do n = 1, num_sw_intervals+num_opt_wavelens+n_jv_bands+num_lw_intervals
      band = int2band(n)              ! ECHAM band index for current interval
      weight(n) = 0.0_dp              ! diagnostic wavelength
      if (n <= num_sw_intervals) then ! shortwave interval
        weight(n) = weight_sw(n)
        sumweight_sw(band) = sumweight_sw(band) + weight(n)
      else if (n > num_sw_intervals+num_opt_wavelens+n_jv_bands) then ! longwave interval
        weight(n) = weight_lw(n-num_sw_intervals-num_opt_wavelens-n_jv_bands)
        sumweight_lw(band) = sumweight_lw(band) + weight(n)
      else                            ! diagnostic or JVAL wavelength
        weight(n) = 0.0_dp
      end if
    end do
    

    ! 2.4 refractive indices

    DO n = 1, num_aerosols
      DO k = 1, num_sw_intervals
        opts%ref_re(n,k) = ref_re_sw(n,k)
        opts%ref_im(n,k) = ref_im_sw(n,k)
      END DO
      DO k = 1, num_opt_wavelens
        opts%ref_re(n,k+num_sw_intervals) = ref_re_opt(n,k)
        opts%ref_im(n,k+num_sw_intervals) = ref_im_opt(n,k)
      END DO
      DO k = 1, n_jv_bands
        opts%ref_re(n,k+num_sw_intervals+num_opt_wavelens) = ref_re_jval(n,k)
        opts%ref_im(n,k+num_sw_intervals+num_opt_wavelens) = ref_im_jval(n,k)
      END DO
      DO k = 1, num_lw_intervals
        opts%ref_re(n,k+num_opt_wavelens+num_sw_intervals+n_jv_bands) = ref_re_lw(n,k)
        opts%ref_im(n,k+num_opt_wavelens+num_sw_intervals+n_jv_bands) = ref_im_lw(n,k)
      END DO
    END DO
    

   ! 3. Release temporary arrays.

    IF (ASSOCIATED(ref_re_sw))      DEALLOCATE(ref_re_sw)
    IF (ASSOCIATED(ref_im_sw))      DEALLOCATE(ref_im_sw)
    IF (ASSOCIATED(ref_re_lw))      DEALLOCATE(ref_re_lw)
    IF (ASSOCIATED(ref_im_lw))      DEALLOCATE(ref_im_lw)

    znwave = num_sw_intervals+num_opt_wavelens+num_lw_intervals+n_jv_bands
    status = 0

  END SUBROUTINE aeropt_finalise_lut

!==============================================================================
END MODULE MESSY_AEROPT_INPUT
