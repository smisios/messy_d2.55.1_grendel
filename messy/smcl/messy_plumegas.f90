! **********************************************************************
MODULE messy_plumegas
! **********************************************************************

  ! MESSy submodel core layer (SMCL) for PLUMEGAS
  !
  ! PLUMEGAS
  !      This MESSy (Joeckel et al. 2005) submodel for the ECHAM5/MESSy
  !      Atmospheric Chemistry (EMAC) model (Joeckel et al. 2006,
  !      Roeckner et al. 2006) aims at applying the Cariolle et al.
  !      (2009) parameterisation of plume chemistry to NOx emissions
  !      from e.g. combustion engines employed in aviation and shipping
  !      (Huszar et al. 2009). See the interface layer header for
  !      further information. This SMCL module includes the
  !      computational mashinery (in particular the equations of the
  !      parameterisation) to be called by either of the MESSy submodel
  !      interface layer (SMIL) for global numerical simulations or a
  !      potential box model interface layer module for box model
  !      simulations. Subroutine and function names at the SMCL begin
  !      with plg_ .
  !
  ! Authors
  !      Mauro Dall'Amico, DLR Oberpfaffenhofen, Germany, 2009-2010
  !      Patrick Joeckel, DLR Oberpfaffenhofen, Germany, 2009-2010
  !
  ! Bug reports
  !      Deutsches Zentrum fuer Luft- und Raumfahrt, Institut fuer
  !      Physik der Atmosphaere, ECHAM-Gruppe, Oberpfaffenhofen,
  !      Germany
  !
  ! Acknowledgements
  !      Financial support for developing submodel PLUMEGAS was
  !      provided by the Helmholtz Young Investigators Group SeaKLIM
  !      led by Veronika Eyring at DLR, Oberpfaffenhofen, Germany.
  !
  ! References
  !      Cariolle et al. (2009) "Parametrisation of plume chemistry into
  !           large-scale atmospheric models: Application to aircraft
  !           NOx emissions", J. Geophys. Res., 114, D19302,
  !           doi:10.1029/2009JD011873.
  !      Huszar et al. (2009) "Modeling the regional impact of ship
  !           emissions on NOx and ozone levels over the Eastern
  !           Atlantic and Western Europe using ship plume
  !           parameterization", Atmos. Chem. Phys. Discuss., 9,
  !           26735-26776.
  !      Joeckel et al. (2006) "The atmospheric chemistry general
  !           circulation model ECHAM5/MESSy1: consistent simulation of
  !           ozone from the surface to the mesosphere", Atmos. Chem.
  !           Phys., 6, 5067-5104.
  !      Roeckner et al. (2006) "Sensitivity of simulated climate to
  !           horizontal and vertical resolution in the ECHAM5
  !           atmosphere model", J. Climate, 19, 3771-3791
  !      Joeckel et al. (2005) Technical Note: The Modular Earth
  !           Submodel System (MESSy) - a new approach towards Earth
  !           System Modeling, Atmos. Chem. Phys., 5, 433-444
  !
  ! Important note
  !      This version is only for internal use and test at DLR,
  !      Oberpfaffenhofen, Germany
  !
  ! Version
  !      1.0 - DRAFT (ENTWURF)


  ! Base model interface layer (BMIL):

  USE messy_main_constants_mem,ONLY: dp

  USE messy_main_tools,        ONLY: ptr_3d_array,ptr_2d_array


  IMPLICIT NONE

  INTRINSIC NULL,TINY

  PRIVATE

  SAVE


  ! As a convention,any new submodel provides a public identification
  ! string and version number (modstr and modver,respectively):
  CHARACTER(LEN=*),PARAMETER,PUBLIC :: modstr='plumegas'
  CHARACTER(LEN=*),PARAMETER,PUBLIC :: modver='0.9'

  ! Maximum number of plume types allowed:
  INTEGER,PARAMETER,PUBLIC :: max_ty=20

  TYPE(ptr_3d_array),DIMENSION(max_ty),PUBLIC :: in_plNOx
  TYPE(ptr_3d_array),DIMENSION(max_ty),PUBLIC :: in_NOx
  TYPE(ptr_3d_array),DIMENSION(max_ty),PUBLIC :: in_HNO3
  TYPE(ptr_3d_array),DIMENSION(max_ty),PUBLIC :: in_O3_titr
  TYPE(ptr_3d_array),DIMENSION(max_ty),PUBLIC :: in_O3_keff

  ! The start values of the tracers used in the submodel are
  ! declared as follows:
  ! Start value of plume NOx:
  TYPE(ptr_2d_array),DIMENSION(max_ty),PUBLIC :: plNOx_0
  ! The above is an array of pointers to 2d arrays (1:nproma, 1:nlev)
  ! to be allocated in init_memory, while external tracers are declared
  ! as module-specific non-channel memory:
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: NO_0
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: NO2_0
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: O3_0

  ! Also the values of the molecular density of air are needed:
  REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC :: molecules_m3

  ! The following channel object includes the 2d field day=1 night=0 :
  REAL(dp),DIMENSION(:,:),POINTER,PUBLIC :: ptr_day => NULL()
  ! (to be retrieved in init_coupling).

  ! The following logical is .FALSE. for void plume types:
  LOGICAL,DIMENSION(max_ty),PUBLIC :: l_pl_ty=.FALSE.

  ! A number which is almost zero for comparing floating point
  ! numbers with zero (epsilon is and intrinsic function):
  REAL(dp),PARAMETER,PUBLIC :: al_0=10._dp*TINY(0.0_dp)


  ! Global CTRL namelist parameters:

  ! Names of plume types (a sequence of blank spaces can
  ! interpreted as a non-used type):
  CHARACTER(LEN=8), DIMENSION(max_ty), PUBLIC :: name_pl_ty='        '

  ! Time constants for plume parameterization (depending
  ! on plume type and phase):
  REAL(dp),DIMENSION(max_ty),PUBLIC :: tau=-1._dp
  REAL(dp),DIMENSION(max_ty),PUBLIC :: beta_1=-1._dp
  REAL(dp),DIMENSION(max_ty),PUBLIC :: beta_0=-1._dp
  REAL(dp),DIMENSION(max_ty),PUBLIC :: emi_ra=-1._dp
  REAL(dp),DIMENSION(max_ty),PUBLIC :: keff=-1._dp


  ! Subroutines at the SMCL must be independent of resolution,base
  ! model grid geometry,parallel environment and other submodels.

  ! Public subroutines:
  PUBLIC :: plg_read_nml_ctrl ! Read control namelist
  PUBLIC :: plg_set_logical
  PUBLIC :: plg_par

  ! The following lines are commented since PRIVATE has been set above:
  !PRIVATE :: plg_in_plNOx
  !PRIVATE :: plg_in_NOx
  !PRIVATE :: plg_in_HNO3
  !PRIVATE :: plg_in_O3_titr
  !PRIVATE :: plg_in_O3_Keff


CONTAINS


! #####################################################################
! Public subroutines
! #####################################################################

! ---------------------------------------------------------------------
  SUBROUTINE plg_set_logical()

    ! plg_set_logical
    !    Set a logical variable to .TRUE. if the plume type is defined.


    IMPLICIT NONE

    INTRINSIC :: TRIM


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_set_logical'


    INTEGER :: ity=0


    all_ty: DO ity=1,max_ty

        l_pl_ty(ity)=(TRIM(name_pl_ty(ity)) /= '')

    ENDDO all_ty


  END SUBROUTINE
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_par( & !status,txt_err &
                    kproma  &
                    ,nlev    &
                    ,jrow    &
                    )

    ! plg_par
    !    Parameterisation of subgrid plume NOx chemistry


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)           :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_par'


    INTEGER :: ity=0


    !status=1
    !txt_err=''


    all_ty: DO ity=1,max_ty

      IF (.NOT. l_pl_ty(ity)) CYCLE

      CALL plg_in_plNOx( & !status,txt_err &
                       ity     &
                       ,kproma  &
                       ,nlev    &
                       ,jrow    &
                       )
      !IF (status /= 0) RETURN ! status=2 ! Commented as there is no
      !                                     RETURN yet in the
      !                                     subroutines.

      CALL plg_in_NOx( & !status,txt_err &
                     ity     &
                     ,kproma  &
                     ,nlev    &
                     ,jrow    &
                     )
      !IF (status /= 0) RETURN ! status=3


      CALL plg_in_HNO3( & !status,txt_err &
                      ity     &
                      ,kproma  &
                      ,nlev    &
                      ,jrow    &
                      )
      !IF (status /= 0) RETURN ! status=4

      CALL plg_in_O3_titr( & !status,txt_err &
                         ity     &
                         ,kproma  &
                         ,nlev    &
                         ,jrow    &
                         )
      !IF (status /= 0) RETURN ! status=5

      CALL plg_in_O3_keff( & !status,txt_err &
                         ity     &
                         ,kproma  &
                         ,nlev    &
                         ,jrow    &
                         )
      !IF (status /= 0) RETURN ! status=6

    ENDDO all_ty

    !status=0
    !txt_err=''

  END SUBROUTINE plg_par
! ---------------------------------------------------------------------


! #####################################################################
! Private subroutines
! #####################################################################

! ---------------------------------------------------------------------
  SUBROUTINE plg_in_plNOx( & !status,txt_err &
                         ity     &
                         ,kproma  &
                         ,nlev    &
                         ,jrow    &
                         )

    ! plmg_in_plNOx
    !    Calculation of increment to plume NOx tendency due to
    !    dilution. The emission of NOx into the plume NOx tracers is
    !    specified through submodel OFFLEM (see the note in
    !    plumegas.nml).


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)          :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: ity
    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_in_plNOx'


    INTEGER :: ikp=0, ile=0


    !status=2
    !txt_err=' (in '//substr//')'


    DO ile=1,nlev
      DO ikp=1,kproma

        in_plNOx(ity)%PTR(ikp,ile,jrow)= &
          -1._dp/tau(ity)                &
          *plNOx_0(ity)%PTR(ikp,ile)

      ENDDO
    ENDDO

    !status=0
    !txt_err=''

  END SUBROUTINE plg_in_plNOx
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_in_NOx( & !status,txt_err &
                       ity     &
                       ,kproma  &
                       ,nlev    &
                       ,jrow    &
                       )

    ! plg_in_NOx
    !    Calculation of increment to the NOx tendency


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)         :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: ity
    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_in_NOx'


    INTEGER :: ikp=0, ile=0


    !status=3
    !txt_err=' (in '//substr//')' 


    DO ile=1,nlev
      DO ikp=1,kproma

        in_NOx(ity)%PTR(ikp,ile,jrow)=                           &
           1._dp/tau(ity)                                        &
          *plNOx_0(ity)%PTR(ikp,ile)                             &
          *(1._dp-beta_1(ity)*ptr_day(ikp,jrow)          &
                 -beta_0(ity)*(1._dp-ptr_day(ikp,jrow)))

      ENDDO
    ENDDO

    !status=0
    !txt_err=''

  END SUBROUTINE plg_in_NOx
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_in_HNO3( & !status,txt_err &
                        ity     &
                        ,kproma  &
                        ,nlev    &
                        ,jrow    &
                        )

    ! plg_in_HNO3
    !    Calculation of increment to the HNO3 tendency


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)          :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: ity
    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_in_HNO3'


    INTEGER :: ikp=0, ile=0


    !status=4
    !txt_err=' (in '//substr//')' 


    DO ile=1,nlev
      DO ikp=1,kproma

        in_HNO3(ity)%PTR(ikp,ile,jrow)=                    &
           1._dp/tau(ity)                                  &
          *plNOx_0(ity)%PTR(ikp,ile)                       &
          *(beta_1(ity)*ptr_day(ikp,jrow)          &
           +beta_0(ity)*(1._dp-ptr_day(ikp,jrow)))

      ENDDO
    ENDDO

    !status=0
    !txt_err=''

  END SUBROUTINE plg_in_HNO3
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_in_O3_titr( & !status,txt_err &
                           ity     &
                           ,kproma  &
                           ,nlev    &
                           ,jrow    &
                           )

    ! plg_in_O3_titr
    !    Calculation of increment to the O3 tendency associated with
    !    titration


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)          :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: ity
    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_in_O3_titr'


    INTEGER :: ikp=0, ile=0


    !status=5
    !txt_err=' (in '//substr//')' 


    DO ile=1,nlev
      DO ikp=1,kproma

        IF (NO_0(ikp,ile)+NO2_0(ikp,ile) > al_0) THEN

          in_O3_titr(ity)%PTR(ikp,ile,jrow)= &
            -1._dp/tau(ity)                  &
            *plNOx_0(ity)%PTR(ikp,ile)       &
            *(NO2_0(ikp,ile)                 &
             /(NO_0(ikp,ile)+NO2_0(ikp,ile)) &
             -emi_ra(ity))                   &
            *ptr_day(ikp,jrow)

        ELSE ! No increment if there are no NOx in the grid box:

          in_O3_titr(ity)%PTR(ikp,ile,jrow)=0._dp

        ENDIF

      ENDDO
    ENDDO

    !status=0
    !txt_err=''

  END SUBROUTINE plg_in_O3_titr
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_in_O3_keff( & !status,txt_err &
                           ity     &
                           ,kproma  &
                           ,nlev    &
                           ,jrow    &
                           )

    ! plg_in_O3_keff
    !    Calculation of increment to O3 tendency associated with the
    !    effective reaction rate (keff)


    IMPLICIT NONE


    ! Input / output:

    !INTEGER,INTENT(OUT)          :: status
    !CHARACTER(LEN=*),INTENT(OUT) :: txt_err

    INTEGER,INTENT(IN) :: ity
    INTEGER,INTENT(IN) :: kproma
    INTEGER,INTENT(IN) :: nlev
    INTEGER,INTENT(IN) :: jrow


    !CHARACTER(LEN=*),PARAMETER :: substr='plg_in_O3_keff'


    INTEGER :: ikp=0, ile=0


    !status=6
    !txt_err=' (in '//substr//')' 


    DO ile=1,nlev
      DO ikp=1,kproma

        in_O3_keff(ity)%PTR(ikp,ile,jrow)= &
          -keff(ity)                       &
          *1.E-6_dp                        & ! cm-3 m3
          *O3_0(ikp,ile)                   &
          *plNOx_0(ity)%PTR(ikp,ile)       &
          *molecules_m3(ikp,ile)           &
          *ptr_day(ikp,jrow)

      ENDDO
    ENDDO

    !status=0
    !txt_err=''

  END SUBROUTINE plg_in_O3_keff
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE plg_read_nml_ctrl(status,iou)

    ! plg_read_nml_ctrl
    !    Read control namelist


    ! BMIL:
    USE messy_main_tools,ONLY: read_nml_open  &
                              ,read_nml_check &
                              ,read_nml_close


    IMPLICIT NONE

    INTRINSIC :: TRIM


    ! Input / output:

    INTEGER,INTENT(OUT) :: status
    INTEGER,INTENT(IN)  :: iou      ! Input / output unit.


    NAMELIST /CTRL/ name_pl_ty &
                   ,tau        &
                   ,beta_1     &
                   ,beta_0     &
                   ,emi_ra     &
                   ,keff


    CHARACTER(LEN=*),PARAMETER :: substr='plg_read_nml_ctrl'


    LOGICAL                    :: lex   ! File exists
    INTEGER                    :: fstat ! File status

    INTEGER          :: ity=0
    CHARACTER(LEN=2) :: ch
    LOGICAL          :: err_tau =.FALSE.
    LOGICAL          :: err_ctrl=.FALSE.


    status=1


    ! Read and check namelist file:

    CALL read_nml_open(lex    &
                      ,substr &
                      ,iou    &
                      ,'CTRL' &
                      ,modstr &
                      )
    IF (.NOT. lex) RETURN ! Namelist file does not exist.

    READ (iou          &
         ,NML=CTRL     &
         ,IOSTAT=fstat &
         )

    CALL read_nml_check(fstat  &
                       ,substr &
                       ,iou    &
                       ,'CTRL' &
                       ,modstr &
                       )
    IF (fstat /= 0) RETURN ! Namelist could not be read.


    ! Diagnose namelist settings:

    WRITE(*,*)
    WRITE(*,*) 'Namelist settings for ',substr

    all_ty: DO ity=1,max_ty
      WRITE(ch,'(i2.2)') ity

      WRITE(*,*)

      IF (TRIM(name_pl_ty(ity)) /= '') THEN
        WRITE(*,*) 'Plume type '//ch//' refers to: '//name_pl_ty(ity)
      ELSE
        WRITE(*,*) 'Plume type '//ch//' is void.'
        CYCLE
      ENDIF


      ! Control that the parameter values are meaningful:

      WRITE(*,*) 'Parameters for plume type '//ch//':'


      WRITE(*,*) 'tau('//ch//')   =', tau(ity)
      IF (tau(ity)    .LE. 0.) err_tau= .TRUE.

      IF (err_tau) WRITE(*,*) "Error: tau is less or equal to zero."


      WRITE(*,*) 'beta_1('//ch//')=', beta_1(ity)
      IF (beta_1(ity) .LT. 0.) err_ctrl=.TRUE.

      WRITE(*,*) 'beta_0('//ch//')=', beta_0(ity)
      IF (beta_0(ity) .LT. 0.) err_ctrl=.TRUE.

      WRITE(*,*) 'emi_ra('//ch//')=', emi_ra(ity)
      IF (emi_ra(ity) .LT. 0.) err_ctrl=.TRUE.

      WRITE(*,*) 'keff('//ch//')  =', keff(ity)
      IF (keff(ity)   .LT. 0.) err_ctrl=.TRUE.

      IF (err_ctrl) WRITE(*,*) "Error: one (or more) parameter(s) "// &
                               "is (are) negative."


      IF (err_tau .OR. err_ctrl) RETURN

    ENDDO all_ty

    WRITE(*,*)

    ! Close namelist file:
    CALL read_nml_close(substr &
                      ,iou     &
                      ,modstr  &
                       )

    status=0 ! If there is no error, the status is 0 - if there was
             ! an error before, this point is not reached.

  END SUBROUTINE plg_read_nml_ctrl
! ---------------------------------------------------------------------

! **********************************************************************
END MODULE messy_plumegas
! **********************************************************************
