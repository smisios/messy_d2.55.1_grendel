!+ Program to generate simple 2D ASCII test data sets for the
!------------------------------------------------------------------------------

MODULE params

!------------------------------------------------------------------------------
!
! Description:
! Program to generate simple 2D ASCII test data sets for the
! soil parameters of the COSMO model. These data sets are
! suitable to be used for idealized runs in case of
!
!    itype_soil_c = 2    and / or   itype_soil_tw = 2  .
!
! Note that these 2D ASCII data sets can be larger than
! the model grid. The dimensions ie and je of these data
! sets can be defined in the below module "params".
!
! The files are written to the subdirectory "TEST_ASCII".
! If it does not exist, it will be created automatically.
!
! Current Code Owner: DWD, Ulrich Blahak
!  phone:  +49  69  8062 2393
!  fax:    +49  69  8062 3721
!  email:  ulrich.blahak@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! V4_20        2011/08/31 <Your name>
!  Initial release
! @VERSION@    @DATE@     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

  IMPLICIT NONE

!==============================================================================

  INTEGER, PARAMETER :: ireals = SELECTED_REAL_KIND (12,200)
  INTEGER, PARAMETER :: iintegers = KIND(1)

  !.. Dimension of the output fields:
  INTEGER(kind=iintegers), PARAMETER :: ie = 561
  INTEGER(kind=iintegers), PARAMETER :: je = 501

  !.. Path of ASCII files to write (will be created if not present):
  CHARACTER(len=*), PARAMETER :: outpath = 'TEST_ASCII'

  PUBLIC ::  ireals, iintegers, ie, je, outpath

END MODULE params



PROGRAM gen_soildata_ascii

  USE params

  IMPLICIT NONE

  REAL(kind=ireals), DIMENSION(ie,je) :: field
  REAL(kind=ireals) :: pi
  INTEGER(kind=iintegers) :: i, j

  pi = 4.0_ireals * ATAN(1.0_ireals)

  !.. Create the output path:
  CALL system('mkdir -p '//TRIM(outpath))

  !====================================================================
  !====================================================================
  !
  ! Write the files for more or less constant parameters one by one:
  !
  !====================================================================
  !====================================================================

  !      hsurf      ,    & ! geometrical heigt of surface topography        (  m  )
  !====================================================================
  DO j=1, je
    DO i=1, ie
      field(i,j) = 500.0 * (1.0 + SIN(0.25*pi*(i-1.0)/ie) * COS(0.125*pi*(j-1.0)/je))
    END DO
  END DO
  CALL out_field(field, 'orofile_test.dat', ' geometical heigt of surface topography [m]')

  !      z0         ,    & ! surface roughness                             (  m  )
  !====================================================================
  field(:,:) = 0.1_ireals * 9.816_ireals
  CALL out_field(field, 'z0file_test.dat', 'roughness length [m]')

  !      fr_land    ,    & ! fraction of land in a grid element            ( --  )
  !====================================================================
  field(:,:) = 1.0_ireals
  CALL out_field(field, 'frlandfile_test.dat', 'fraction of land in a grid element [-]')

  !      plcov      ,    & ! fraction of plant cover                         --
  !====================================================================
  field(:,:) = 0.4_ireals
  CALL out_field(field, 'plcovfile_test.dat', 'fraction of plant cover [-]')

  !      lai        ,    & ! leaf area index of plants                       --
  !====================================================================
  field(:,:) = 2.5_ireals
  CALL out_field(field, 'laifile_test.dat', 'leaf area index of plants [-]')


  !
  !  IF (lsso)
  !
  !      sso_stdh   ,    & ! standard deviation of sub-grid scale orography ( m   )
  !====================================================================
  field(:,:) = 20.0_ireals
  CALL out_field(field, 'ssostdhfile_test.dat', 'standard deviation of sub-grid scale orography [ m ]')

  !      sso_gamma  ,    & ! anisotropy of sub-grid scale orography          --
  !====================================================================
  field(:,:) = 0.0_ireals
  CALL out_field(field, 'ssogammafile_test.dat', 'anisotropy of sub-grid scale orography [ - ]')

  !      sso_theta  ,    & ! angle betw. principal axis of orography and E ( rad )
  !====================================================================
  field(:,:) = 0.25_ireals * 4.0*ATAN(1.0_ireals)
  CALL out_field(field, 'ssothetafile_test.dat', 'angle betw. principal axis of orography and E [ rad ]')

  !      sso_sigma  ,    & ! mean slope of sub-grid scale orography          --
  !====================================================================
  field(:,:) = 0.05_ireals
  CALL out_field(field, 'ssosigmafile_test.dat', 'mean slope of sub-grid scale orography [ - ]')


  !
  !  IF (lsoil)
  !      soiltyp    ,    & ! type of the soil (keys 0-9)                   ( --  )
  !====================================================================
  field(:,:) = 5.0_ireals
  CALL out_field(field, 'soiltypefile_test.dat', 'type of soil [-]')

  !      rootdp     ,    & ! depth of the roots                            (  m  )
  !====================================================================
  field(:,:) = 0.76_ireals
  CALL out_field(field, 'rootdpfile_test.dat', 'depth of the roots [m]')

!!$ STILL MISSING:
!!$    IF (lstomata)
!!$ *    rsmin2d    ,    & ! minimum stomata resistance                    ( s/m )
!!$
!!$    IF (lrad)
!!$      vio3       ,    & ! vertical integrated ozone contents            (pa O3)
!!$      hmo3       ,    & ! ozone maximum                                 ( pa  )
!!$
!!$    IF (lrad .and. itype_aerosol == 2)
!!$ *    aer_su     ,    & ! monthly aerosol climatology sulfate drops     (0 - 1)
!!$ *    aer_du     ,    & ! monthly aerosol climatology total dust        (0 - 1)
!!$ *    aer_or     ,    & ! monthly aerosol climatology organic (water sol.)(0-1)
!!$ *    aer_bc     ,    & ! monthly aerosol climatology black carbon      (0 - 1)
!!$ *    aer_ss     ,    & ! monthly aerosol climatology sea salt          (0 - 1)
!!$
!!$    IF (lrad .and. lemiss)
!!$ *    emis_rad   ,    & ! external thermal emissivity                   (0 - 1)
!!$
!!$    IF (lrad .and. lradtopo)
!!$ *    skyview    ,    & ! sky view
!!$ *    slo_asp    ,    & ! slope aspect
!!$ *    slo_ang    ,    & ! slope angle
!!$ *    horizon    ,    & ! horizon
!!$

  !  IF (lforest)
  !      for_e      ,    & ! ground fraction covered by evergreen forest     --
  !====================================================================
  field(:,:) = 0.25_ireals
  CALL out_field(field, 'forefile_test.dat', 'ground fraction covered by evergreen forest [-]')

  !      for_d      ,    & ! ground fraction covered by deciduous forest     --
  !====================================================================
  field(:,:) = 0.15_ireals
  CALL out_field(field, 'fordfile_test.dat', 'ground fraction covered by deciduous forest [-]')


!!$ STILL MISSING:
!!$    IF (llake)
!!$ *    t_mnw_lk  ,     & ! mean temperature of the water column          (  K  )
!!$ *    t_wml_lk  ,     & ! mixed-layer temperature                       (  K  )
!!$ *    t_bot_lk  ,     & ! temperature at the water-bottom sediment
!!$                        ! interface                                     (  K  )
!!$ *    t_b1_lk   ,     & ! temperature at the bottom of the upper layer
!!$                        ! of the sediments                              (  K  )
!!$ *    c_t_lk    ,     & ! shape factor with respect to the
!!$                        ! temperature profile in lake thermocline       (  -  )
!!$ *    h_ml_lk   ,     & ! thickness of the mixed-layer                  (  m  )
!!$ *    h_b1_lk           ! thickness of the upper layer
!!$                        ! of bottom sediments                           (  m  )
!!$ *    fr_lake    ,    & ! lake fraction in a grid element [0,1]         (  -  )
!!$ *    depth_lk   ,    & ! lake depth                                    (  m  )
!!$ *    fetch_lk   ,    & ! wind fetch over lake                          (  m  )
!!$ *    dp_bs_lk   ,    & ! thickness of the thermally active layer
!!$                        ! of bottom sediments                           (  m  )
!!$ *    t_bs_lk    ,    & ! climatological temperature at the bottom of
!!$                        ! the thermally active layer of sediments       (  K  )
!!$ *    gamso_lk          ! attenuation coefficient for
!!$                        ! solar radiation in lake water                 ( 1/m )
!!$

  !  IF (lseaice .or. llake)
  !      h_ice     ,     & ! ice thickness                                 (  m  )
  !====================================================================
  field(:,:) = 0.1_ireals
  CALL out_field(field, 'hicefile_test.dat', 'ice thickness [ m ]')

  !      t_ice     ,     & ! temperature at the snow-ice or air-ice interface (  K  )
  !====================================================================
  field(:,:) = 270.0_ireals
  CALL out_field(field, 'ticefile_test.dat', 'temperature at the snow-ice or air-ice interface [ K ]')




  !====================================================================
  !====================================================================
  !
  ! Write the initial files for time-varying parameters one by one:
  !
  !====================================================================
  !====================================================================

  
  !       t_soil    ,     & !             soil temperature                  (  k  )
  !====================================================================
  field(:,:) = 276.0_ireals
  CALL out_field(field, 'tsoilfile_test.dat', 'temperature of the soil [K]')

  !       wf_soil   ,     & ! soil water saturation                         (  -  )
  !====================================================================
  field(:,:) = 0.67_ireals
  CALL out_field(field, 'wfsoilfile_test.dat', 'temperature of the snow-surface [K]')

  !       t_snow    ,     & ! temperature of the snow-surface               (  k  )
  !====================================================================
  field(:,:) = 265.0_ireals
  CALL out_field(field, 'tsnowfile_test.dat', 'temperature of the snow-surface [K]')

  !       w_snow    ,     & ! snow water equivalent                         (m H2O)
  !====================================================================
  field(:,:) = 0.16_ireals
  CALL out_field(field, 'wsnowfile_test.dat', 'snow water equivalent [m H2O]')

  !       w_i       ,     & ! water content of interception water           (m H2O)
  !====================================================================
  field(:,:) = 0.01234_ireals
  CALL out_field(field, 'wifile_test.dat', 'water content of interception water [m H2O]')
  


END PROGRAM gen_soildata_ascii


!.. Subroutine to write a 2D-field to an ASCII output file:
SUBROUTINE out_field( field, outfilename, headerline)

  USE params

  IMPLICIT NONE

  REAL(kind=ireals), INTENT(in) :: field(ie,je)
  CHARACTER(len=*), INTENT(in) :: outfilename, headerline

  INTEGER(kind=iintegers) :: i, j, err


  OPEN(10, file=TRIM(ADJUSTL(outpath))//'/'//TRIM(ADJUSTL(outfilename)), &
       status='replace', form='formatted', iostat=err)
  IF (err /= 0 ) THEN
    WRITE (*,*) 'Error opening '//TRIM(ADJUSTL(outfilename))//' !'
    STOP
  END IF
  WRITE(10,'(a)') '# '//TRIM(ADJUSTL(headerline))
  WRITE(10,'(i4,x,i4)') ie, je
  DO j=1, je
    DO i=1, ie
      IF ( ABS(field(i,j)) >= 1e-30_ireals ) THEN
        WRITE (10,'(es12.5)') field(i,j)
      ELSE
        WRITE (10,'(i1)') 0_iintegers
      END IF
    END DO
  END DO
  CLOSE(10)

END SUBROUTINE out_field
