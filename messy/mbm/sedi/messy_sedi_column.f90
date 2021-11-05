! Time-stamp: <2014-01-28 17:18:39 joec_pa> -*- f90 -*-
! Authors:
! Astrid Kerkweg, MPICH, 2006

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE messy_sedi_column

  USE messy_main_constants_mem, ONLY: DP
  USE sedi_column_netcdf,       ONLY: open_nc_file, write_nc_file, close_nc_file
  USE messy_sedi,               ONLY: modstr, order

  IMPLICIT NONE

  PRIVATE

  INTEGER :: ncid_messy, ncid_aerosol, ncid_loss

  ! INIT DIMENSIONS
  INTEGER, PUBLIC :: DT = 200 ! DEFAULT MODEL TIMESTEP
  INTEGER, PUBLIC :: NT = 100 ! DEFAULT NUMBER OF TIMESTEPs  
  INTEGER, PUBLIC            :: profile = -1
  CHARACTER(LEN=120), PUBLIC :: prof_meteo =''

  INTEGER, PARAMETER :: nproma = 1  ! number of columns
  INTEGER :: nlev  ! number of vertical levels
  REAL(dp), POINTER, DIMENSION(:,:) :: temp 
  REAL(dp), POINTER, DIMENSION(:,:) :: press 
  REAL(dp), POINTER, DIMENSION(:,:) :: pressi 
  REAL(dp), POINTER, DIMENSION(:,:) :: sphum
  REAL(dp), POINTER, DIMENSION(:,:) :: xt    ! aerosol number concentration
                                               ! in 1/mol
  REAL(dp), POINTER, DIMENSION(:,:) :: densaer ! aerosol density
  REAL(dp), POINTER, DIMENSION(:,:) :: radaer  ! aerosol radius

  REAL(dp), POINTER, DIMENSION(:,:) :: vt_term ! terminal velocity
  REAL(dp), POINTER, DIMENSION(:,:) :: sfloss  ! loss at surface [mol/mol/s]


  PUBLIC :: sedi_column_initialize ! initialize particles
                                   ! CORE:     sedi_read_nml_ctrl
                                   ! PRIVATE   sedi_read_nml_cpl_column
  PUBLIC :: sedi_column_init_memory
  PUBLIC :: sedi_column_init_column
       ! PRIAVTE :: sedi_init_aerosol_profiles
       ! PRIVATE :: sedi_init_profiles_meteo
  PUBLIC :: sedi_column_physc  ! calculate  terminal velocities
       ! CORE:     calc_vt
       !           sedi_0
       !           sedi_1 
  PUBLIC :: sedi_column_result      ! print results
  PUBLIC :: sedi_column_free_memory ! free memory

CONTAINS

  !*****************************************************************************

  SUBROUTINE sedi_column_initialize
    
    USE messy_sedi,   ONLY: sedi_read_nml_ctrl
    
    IMPLICIT NONE

    INTEGER :: status =-1
    INTEGER :: IOU

    ! read required order of sedimentation scheme
    IOU =101
    CALL sedi_read_nml_ctrl(status,IOU)
    
    IF (status /= 0) THEN
       write(*,*) "READ ERROR in SEDI CTRL"
       STOP
    ELSE
       write (*,*) " SEDIMENTATION SCHEME: ", order
    ENDIF
    
   ! determine aerosol profile
    IOU =102
    CALL sedi_read_nml_cpl_column(status,IOU)

    IF (status /= 0) THEN
       write(*,*) "READ ERROR in SEDI CTRL_COLUMN"
       STOP 
    ELSE
       write(*,*) '---------------------------------------'
       write(*,*) 'TIMESTEP LENGTH: ',DT
       write(*,*) '---------------------------------------'
       write(*,*) 'NUMBER OF TIMESTEPS: ',NT
       write(*,*) '---------------------------------------'

       write(*,*) " PROFILE:               ", profile
       IF (profile > 5) THEN
          WRITE(*,*) "PROFILE NOT AVAILABLE"
          STOP
       ENDIF

       write(*,*) "COLUMN PARAMETERS ARE READ FROM ",prof_meteo
       write(*,*) '---------------------------------------'
    ENDIF

  END SUBROUTINE sedi_column_initialize

! --------------------------------------------------------------------------
  SUBROUTINE sedi_column_init_column
    
    IMPLICIT NONE

    CALL sedi_init_profiles_meteo

    CALL sedi_init_aerosol_profiles

  END SUBROUTINE sedi_column_init_column
! --------------------------------------------------------------------------

  SUBROUTINE sedi_column_init_memory

    IMPLICIT NONE

    ! open netcdf files and write headers

    write(*,*) 'opening file sedi_messy'
    ! write meteorology 
    CALL open_nc_file(ncid_messy, 'sedi_messy', &
      (/ 'press ','temp  ','sphum ' /), &
      (/ 'Pa   ','K    ','kg/kg'/),klev=nlev )

    ! terminal velocity
    ALLOCATE(vt_term(nproma,nlev))
    vt_term(:,:) = 0._dp
    ! sedimantational loss
    ALLOCATE(sfloss(nproma,nlev))
    sfloss(:,nlev) = 0._dp

    write(*,*) 'opening file sedi_aerosol'
    CALL open_nc_file(ncid_aerosol, 'sedi_aerosol', &
         (/'vt_term','densaer','radaer ','xt     '/),   &
         (/'m/s  ','kg/m3','m    ','1/mol'/),klev=nlev )

    write(*,*) 'opening file sedi_loss'

    CALL open_nc_file(ncid_loss, 'sedi_loss', (/'sfloss'/), &
         (/'1/(mol s)'/),klev=nlev)

  END SUBROUTINE sedi_column_init_memory

!*****************************************************************************

  SUBROUTINE sedi_column_physc

    USE messy_sedi,              ONLY:  mean_free_path &!!$, density &
                                      , air_viscosity           &  
                                      , calc_vt, sedi_0, sedi_1
    USE messy_main_tools,        ONLY: mass_density ! op_pj_20140128

    ! MESSy
    USE messy_main_constants_mem,  ONLY: g, M_air

    IMPLICIT NONE
    
    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='sedi_column_physc'
    REAL(dp) :: rho_air(nproma,nlev)    ! density of air
    REAL(dp) :: lambda_air(nproma,nlev) ! mean free path of air
    REAL(dp) :: visc_air(nproma,nlev)   ! viscosity of air
    REAL(dp) :: zdz(nproma,nlev)        !layerthickness = delta height
    REAL(dp) :: zdp(nproma,nlev)        !layerthickness = delta pressure 
    INTEGER  :: i
    REAL(dp) :: tendency(nproma,nlev)   ! sedimentation tendency

   ! calculate air density in kg/m3
    rho_air(1:nproma,1:nlev) = &
! op_pj_20140128+
!!$         density(press(1:nproma,1:nlev) &
         mass_density(press(1:nproma,1:nlev) &
! op_pj_20140128-
         ,temp(1:nproma,1:nlev),sphum(1:nproma,1:nlev))

   ! calculate mean free path of air 
    lambda_air(1:nproma,1:nlev)= &
         mean_free_path(press(1:nproma,1:nlev),temp(1:nproma,1:nlev))

   ! calculate viscosity of air 
    CALL air_viscosity(nproma, nlev, temp(1:nproma,1:nlev),  &
         visc_air(1:nproma,1:nlev))

    ! calculate box height zdz and zdp
     zdp(1:nproma,1:nlev) = pressi(1:nproma,2:nlev+1) &
          - pressi(1:nproma,1:nlev) 
 
    zdz(1:nproma,1:nlev) = zdp(1:nproma,1:nlev)/rho_air(1:nproma,1:nlev)/g

    CALL calc_vt(vt_term(1:nproma,1:nlev), 1, nlev , radaer(1:nproma,1:nlev) &
        , densaer(1:nproma,1:nlev), lambda_air(1:nproma,1:nlev)  &
        , rho_air(1:nproma,1:nlev) & 
        , visc_air(1:nproma,1:nlev), 1._dp ,0 )
 
   ! CALCULATE TOTAL SEDIMENTATION TENDENCY [mol/mol/s]
    SELECT CASE  (order) 
    CASE(0)
       CALL sedi_0(nproma, nlev, tendency(1:nproma,1:nlev)  &
            , xt(1:nproma,1:nlev), vt_term(1:nproma,1:nlev) &
            , zdz(1:nproma,1:nlev), zdp(1:nproma,1:nlev), REAL(DT,dp)  &
            , sfloss(1:nproma,nlev))
    CASE(1) 
       CALL sedi_1(nproma, nlev, press(1:nproma,1:nlev)       &
            , temp(1:nproma,1:nlev), vt_term(1:nproma,1:nlev) &
            , REAL(DT,dp), xt(1:nproma,1:nlev), pressi(1:nproma,1:nlev) &
            , tendency(1:nproma,1:nlev)   &
            , sfloss(1:nproma,nlev))
    CASE DEFAULT
      ! CALL finish(substr, 'UNKOWN ORDER OF SEDIMENTATION SCHEME')   
       STOP
    END SELECT
    
    xt(1:nproma,1:nlev) = xt(1:nproma,1:nlev) + tendency(1:nproma,1:nlev)*DT

  END SUBROUTINE sedi_column_physc
  
  !*****************************************************************************
       
  SUBROUTINE sedi_column_result(iloop)
    
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: iloop

    REAL(dp) :: model_time
    
    model_time = REAL(iloop*DT,dp)

    write(*,*) 'TIME: ', model_time,' WRITE OUTPUT'
    CALL write_nc_file(ncid_messy, model_time, &
         (/ press(1,:),temp(1,:), sphum(1,:)/),klev=SIZE(temp,2) )
    CALL write_nc_file(ncid_aerosol,model_time,&
         (/vt_term(1,:),densaer(1,:),radaer(1,:),xt(1,:)/) &
         ,klev=SIZE(radaer,2))
    CALL write_nc_file(ncid_loss,model_time,sfloss(1,:),klev=SIZE(SFLOSS,2) )
    
  END SUBROUTINE sedi_column_result

  !*****************************************************************************
  SUBROUTINE sedi_column_free_memory

    IMPLICIT NONE

   write(*,*) ' CLOSE NC-FILES'

    CALL close_nc_file(ncid_messy)
    CALL close_nc_file(ncid_aerosol)
    CALL close_nc_file(ncid_loss)

   write(*,*) ' FREE MEMORY'

    IF (ASSOCIATED(press))   DEALLOCATE(press)
    IF (ASSOCIATED(pressi))  DEALLOCATE(pressi)
    IF (ASSOCIATED(sphum))   DEALLOCATE(temp)
    IF (ASSOCIATED(temp))    DEALLOCATE(sphum)
    IF (ASSOCIATED(xt))      DEALLOCATE(xt)
    IF (ASSOCIATED(densaer)) DEALLOCATE(densaer)
    IF (ASSOCIATED(radaer))  DEALLOCATE(radaer)
    IF (ASSOCIATED(vt_term)) DEALLOCATE(vt_term)
    IF (ASSOCIATED(sfloss))  DEALLOCATE(sfloss)

  END SUBROUTINE sedi_column_free_memory
  
  !*****************************************************************************
  !*****************************************************************************
  ! PRIVATE ROUTINES
  !*****************************************************************************
  !*****************************************************************************
  SUBROUTINE sedi_read_nml_cpl_column(status, iou)

    ! SEDI MODULE ROUTINE (CORE)
    !
    ! READ SEDI_COLUMN NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Astrid Kerkweg, MPICH, AUG 2006

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='sedi_read_nml_cpl_column'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    NAMELIST /CPL_COLUMN/  DT,NT,profile,prof_meteo

    status = 1 ! ERROR ON RETURN

    CALL read_nml_open(lex, substr, iou, 'CPL_COLUMN', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_COLUMN, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_COLUMN', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE sedi_read_nml_cpl_column

  !========================================================================
  SUBROUTINE sedi_init_profiles_meteo

    IMPLICIT NONE
    
    INTEGER :: IOU, STATUS
    INTEGER :: i, jk
    LOGICAL :: lex
    CHARACTER(LEN=40) :: inputvar
    IOU = 103

    INQUIRE(file=TRIM(prof_meteo), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) '*** WARNING: FILE '''//TRIM(prof_meteo)//'''  NOT FOUND !'
       WRITE(*,*) ' '//TRIM(modstr)//' SWITCHED OFF !'
       RETURN
    END IF
   

    OPEN(IOU,FILE=TRIM(prof_meteo))

    READ(IOU,'(I2)') nlev
    write(*,*) 'NLEV: ',nlev
    IF (nlev < 10) THEN
       write(*,*) 'Number of vertical layers should be at least 10'
       STOP
    ENDIF

    ALLOCATE(temp(nproma,nlev))
    ALLOCATE(sphum(nproma,nlev))
    ALLOCATE(press(nproma,nlev))
    ALLOCATE(pressi(nproma,nlev+1))
    
    temp(:,:)=0._dp
    sphum(:,:)=0._dp
    press(:,:)=0._dp
    pressi(:,:)=0._dp

    DO i=1, 4
       READ(IOU,*) inputvar
       write(*,*) inputvar
       IF (TRIM(ADJUSTL(inputvar))=='ECHAM5_tm1') THEN
          READ (IOU,*) temp
          write(*,*) 'temperature: ',temp
       ELSEIF (TRIM(ADJUSTL(inputvar))=='ECHAM5_qm1') THEN
          READ (IOU,*) sphum
          write(*,*) 'specific humidity: ',sphum
       ELSEIF (TRIM(ADJUSTL(inputvar))=='ECHAM5_press') THEN
          READ (IOU,*) press
          write(*,*) 'pressure: ',press
       ELSEIF (TRIM(ADJUSTL(inputvar))=='ECHAM5_pressi') THEN
          READ (IOU,*) pressi
          write(*,*) 'interface pressure: ',pressi
       ELSE
          write(*,*) 'READ VARIABLE: ',inputvar
          write(*,*) 'VARIABLE: ',inputvar,' is unkown'
          STOP
       ENDIF
    ENDDO

    CLOSE(IOU)
 
  END SUBROUTINE sedi_init_profiles_meteo
  !========================================================================
  SUBROUTINE sedi_init_aerosol_profiles

    IMPLICIT NONE

    INTEGER :: i,nlevh

    ! ALLOCATE aerosol number density to number of available levels 
    ALLOCATE(xt(nproma,nlev))
    ALLOCATE(densaer(nproma,nlev))
    ALLOCATE(radaer(nproma,nlev))
    xt(:,:)      = 0._dp
    densaer(:,:) = 0._dp
    radaer(:,:)  = 0._dp
    nlevh=INT(nlev/2)

    SELECT CASE (profile)
    CASE(1)
       write (*,*) 'INITIALIZE AEROSOL PROFILE:       ',profile
       write (*,*) '    -> one sharp peak in layer 4  '
       write (*,*) '    -> equal aerosol density in all layers: 1000 kg/m^3'
       write (*,*) '    -> equal aerosol radii   in all layers: 1.e-6 m'
       ! PROFILE No. 1:
       ! Build a sharp peak in the highest 5 layers 
       radaer(:,:) =1.e-6
       densaer(:,:) =1000  
       xt(1:nproma,1) = 0._dp
       xt(1:nproma,2) = 10._dp
       xt(1:nproma,3) = 100._dp
       xt(1:nproma,4) = 1000._dp
       xt(1:nproma,5) = 100._dp
       xt(1:nproma,6) = 10._dp
       xt(1:nproma,7:nlev) = 0._dp
    CASE(2)
       write (*,*) 'INITIALIZE AEROSOL PROFILE:       ',profile
       write (*,*) '    -> one sharp peak in layer 4  '
       write (*,*) '    -> aerosol density is highest in peaklayer:'
       write (*,*) '         --> MAX = 2000 kg/m^3'
       write (*,*) '         --> MIN = 500 kg/m^3'
       write (*,*) '    -> equal aerosol radii highest in peaklayer:'
       write (*,*) '         --> MAX = 1.e-6 m'
       write (*,*) '         --> MIN = 1.e-8 m'
       ! PROFILE No. 1:
       ! Build a sharp peak in the highest 5 layers 
       radaer(1:nproma,1) = 1.e-8
       radaer(1:nproma,2) = 0.5e-7
       radaer(1:nproma,3) = 1.e-7
       radaer(1:nproma,4) = 1.e-6
       radaer(1:nproma,5) = 1.e-7
       radaer(1:nproma,6) = 0.5e-7
       radaer(1:nproma,7:nlev) =1.e-8

       densaer(1:nproma,1) = 500
       densaer(1:nproma,2) = 750
       densaer(1:nproma,3) = 1000
       densaer(1:nproma,4) = 2000
       densaer(1:nproma,5) = 1000
       densaer(1:nproma,6) = 750
       densaer(1:nproma,7:nlev) =500
       
       xt(1:nproma,1) = 0._dp
       xt(1:nproma,2) = 10._dp
       xt(1:nproma,3) = 100._dp
       xt(1:nproma,4) = 1000._dp
       xt(1:nproma,5) = 100._dp
       xt(1:nproma,6) = 10._dp
       xt(1:nproma,7:nlev) = 0._dp
    CASE(3)
       ! PROFILE No. 3:
       ! GAUSSIAN DISTRIBUTION
       write (*,*) 'INITIALIZE AEROSOL PROFILE:       ',profile
       write (*,*) '    -> one large gaussion distribution  '
       write (*,*) '    -> equal aerosol density in all layers: 1000 kg/m^3'
       write (*,*) '    -> equal aerosol radii   in all layers: 1.e-6 m'
       densaer(:,:) = 1000
       radaer(:,:)=1.e-6
       
       DO i=nlevh,1,-1
          xt(1:nproma,i)=100*i/nlevh
       ENDDO
       DO i=nlevh+1, NLEV
          xt(1:nproma,i)=100*(NLEV-i)/nlevh
       ENDDO

    CASE(4)
       ! PROFILE No. 4:
       ! GAUSSIAN DISTRIBUTION
       write (*,*) 'PROFILE 4 NOT YET IMPLEMENTED'
       STOP
    CASE(5)
       ! PROFILE No. 5:
       ! DOUBLE PEAK
       write (*,*) 'INITIALIZE AEROSOL PROFILE:       ',profile
       write (*,*) '    -> two peaks in layer 3 and 9  '
       write (*,*) '    -> equal aerosol density in all layers: 1000 kg/m^3'
       write (*,*) '    -> equal aerosol radii   in all layers: 1.e-6 m'
       densaer(:,:) = 1000
       radaer(:,:)=1.e-6

       xt(1:nproma,1) = 0._dp
       xt(1:nproma,2) = 10._dp
       xt(1:nproma,3) = 100._dp
       xt(1:nproma,4) = 10._dp
       xt(1:nproma,5) = 0._dp
       xt(1:nproma,6) = 0._dp
       xt(1:nproma,7) = 0._dp
       xt(1:nproma,8) = 10._dp
       xt(1:nproma,9) = 100._dp
       xt(1:nproma,10) = 10._dp
       xt(1:nproma,11:nlev) = 0._dp
    CASE DEFAULT
       write (*,*) 'PROFILE No. ', profile,' NOT YET IMPLEMENTED'
       STOP
    ENDSELECT

  END SUBROUTINE sedi_init_aerosol_profiles
  !========================================================================

END MODULE messy_sedi_column

!*******************************************************************************
