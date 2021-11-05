!*****************************************************************************
!                Time-stamp: <2018-02-21 17:50:05 sander>
!*****************************************************************************

PROGRAM oic

  USE messy_main_constants_mem,   ONLY: DP, STRLEN_KPPSPECIES, STRLEN_VLONG
  USE messy_mecca_kpp_function,   ONLY: CalcStoichNum
  USE messy_mecca_kpp_parameters, ONLY: NREACT, NVAR
  USE messy_mecca_kpp_monitor,    ONLY: SPC_NAMES, EQN_TAGS, EQN_NAMES
  USE drgep,                      ONLY: CalcDIC
  USE graph_search,               ONLY: dijkstra_adj
  USE m_mrgref,                   ONLY: mrgref

  IMPLICIT NONE

  INTEGER :: N_targets = 0
  CHARACTER(LEN=STRLEN_KPPSPECIES), DIMENSION(NVAR) :: target_names 
  INTEGER, DIMENSION(NVAR) :: targets
  LOGICAL, PARAMETER :: l_verbose = .FALSE.
  REAL(DP) :: epsilon_ep = 0.001 ! initial value of MaxOIC threshold
  REAL(DP) :: StoichNum(NVAR,NREACT)
  REAL(DP) :: A(NREACT) ! rates
  REAL(DP), DIMENSION(NREACT) :: rates ! overall rates for all reactions
  REAL(DP), DIMENSION(NVAR,NVAR)          :: DIC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: OICdata ! (NVAR,N_targets,N_samplepoints)
  REAL(DP), DIMENSION(NVAR)               :: MaxOIC
  INTEGER,  DIMENSION(NVAR)               :: position
  INTEGER,  DIMENSION(NVAR,NVAR)          :: neighbor
  INTEGER,  DIMENSION(NVAR)               :: n_neigh
  LOGICAL,  DIMENSION(NREACT) :: l_rxn
  INTEGER, PARAMETER :: unit = 10
  INTEGER :: iostatus
  INTEGER :: N_samplepoints, js, js2, jsc, jr, jt
  LOGICAL :: l_ex ! file exists ?
  CHARACTER(LEN=4)  :: str_jsc
  CHARACTER(LEN=STRLEN_VLONG) :: infile
  CHARACTER(LEN=STRLEN_VLONG) :: oneline
  CHARACTER(LEN=*), PARAMETER :: RATESFILE_PATH = '../output/skeleton/fullmech/multirun/runs/'
  CHARACTER(LEN=*), PARAMETER :: RATESFILE_NAME = '/caaba_mecca_a_end.dat'

  CALL CalcStoichNum (StoichNum) ! get stoichiometric numbers
  ! print stoichiometric numbers:
  PRINT *
  OPEN(unit, FILE='StoichNum.dat', status='UNKNOWN')
  DO jr = 1, NREACT
    PRINT *, jr, "<", TRIM(EQN_TAGS(jr)), "> ", TRIM(ADJUSTL(EQN_NAMES(jr)))
    DO js = 1, NVAR
      WRITE (unit, "(1X,F10.5)", ADVANCE='NO') StoichNum(js,jr)
      IF (l_verbose) THEN
        IF (StoichNum(js,jr) /= 0.) THEN
          PRINT *, StoichNum(js,jr), SPC_NAMES(js)
        ENDIF
      ENDIF
    ENDDO
    WRITE (unit,*)
  ENDDO
  CLOSE(unit)

  PRINT *
  PRINT *, "NVAR   = ", NVAR
  PRINT *, "NREACT = ", NREACT

  ! define target array;
  OPEN(unit, FILE=TRIM("targets.txt"),  STATUS='OLD')
  DO
    READ(unit, '(A)', IOSTAT=iostatus) oneline
    IF (iostatus < 0) EXIT ! exit do loop at end of file
    IF (oneline(1:1)=='#') THEN
      CYCLE ! ignore comment, read next line
    ELSE
      N_targets = N_targets + 1
      ! delete leading spaces:
      oneline = ADJUSTL(oneline) 
      ! get substring until first space:
      target_names(N_targets) = oneline(1:INDEX(oneline,' ')-1)
    ENDIF
  ENDDO

  PRINT *
  PRINT *, "Number of target species: N_targets = ", N_targets
  targets(:) = 0
  DO jt = 1, N_targets
    DO js = 1, NVAR
      IF ( TRIM(target_names(jt)) == TRIM(SPC_NAMES(js)) ) THEN
        targets(jt) = js
        EXIT
      ENDIF
    ENDDO
    IF ( targets(jt) == 0 ) THEN
      PRINT *, "ERROR: ", jt, TRIM(target_names(jt)), " not found"
      STOP 1
    ELSE
      PRINT *, "Target", jt, "= species number", targets(jt), &
        " = ", TRIM(target_names(jt))
      ! print *, TRIM(target_names(jt)), TRIM(SPC_NAMES(targets(jt)))
    ENDIF
  ENDDO

  ! find number of sample points:
  jsc = 1
  PRINT *
  PRINT *, "List of reaction rates files from different sample points:"
  DO
    ! check if file exists:
    WRITE(str_jsc,'(I4.4)') jsc ! string of jsc with 4 digits and leading zeroes
    infile  = RATESFILE_PATH//str_jsc//RATESFILE_NAME
    INQUIRE(FILE=TRIM(infile), EXIST=l_ex)
    IF (.NOT.l_ex) EXIT ! exit DO loop
    PRINT *, "Reaction rates file: ", TRIM(infile)
    jsc = jsc + 1
  ENDDO
  N_samplepoints = jsc-1
  IF (N_samplepoints == 0) THEN
    PRINT *, "ERROR, missing input file:"
    PRINT *, TRIM(infile)
    STOP 1
  ELSE
    PRINT *, "N_samplepoints = ", N_samplepoints
  ENDIF
  ALLOCATE(OICdata(NVAR,N_targets,N_samplepoints))

  samplepoint_loop: DO jsc = 1, N_samplepoints ! loop over sample points
    WRITE(str_jsc,'(I4.4)') jsc ! string of jsc with 4 digits and leading zeroes
    infile  = RATESFILE_PATH//str_jsc//RATESFILE_NAME
    PRINT *
    PRINT *, "Working on sample point", jsc
    PRINT *, "Reaction rates file: ", TRIM(infile)
    OPEN(unit, FILE=TRIM(infile), status='OLD')
    DO jr = 1,NREACT
      READ (unit,*) A(jr)
      IF (l_verbose) THEN
        WRITE(*,'(I4,ES15.7,A12,A)') &
          jr, A(jr), " "//EQN_TAGS(jr), TRIM(ADJUSTL(EQN_NAMES(jr)))
      ENDIF
    ENDDO
    CLOSE(unit)

    ! calculate all DICs for current sample point:
    CALL CalcDIC(StoichNum(:,:), A(:), DIC(:,:), neighbor(:,:), n_neigh(:))
    ! IN:  StoichNum, reaction rates A
    ! OUT: DIC, neighbor, n_neigh

    DO jt = 1, N_targets
      ! Dijkstra's algorithm with adjacency list:
      CALL dijkstra_adj (NVAR, DIC(:,:), neighbor(:,:), n_neigh(:), &
        targets(jt), OICdata(:,jt,jsc))
    ENDDO

  ENDDO samplepoint_loop

  PRINT *
  PRINT *, "Maximum Overall Interaction Coefficients"
  PRINT *, "of all targets and all sample points:"
  DO js = 1, NVAR
    MaxOIC(js) = MAXVAL(OICdata(js,:,:))
  ENDDO
  
  ! sort and show results:
  CALL mrgref(MaxOIC, position)
  DO js = NVAR, 1, -1
    js2 = position(js)
    WRITE (*, "(I4,2X,A16,2X,ES14.7)") NVAR+1-js, SPC_NAMES(js2), MaxOIC(js2)
  ENDDO

  ! save OIC to file:
  OPEN(unit, FILE='OIC.dat', status='UNKNOWN')
  DO js = 1, NVAR
    WRITE (unit, '(ES14.7,A16)') MaxOIC(js), TRIM(SPC_NAMES(js))
  ENDDO
  CLOSE(unit)

  ! save EQN_TAGS to file:
  OPEN(unit, FILE='EQN_TAGS.dat', status='UNKNOWN')
  DO jr = 1, NREACT
    WRITE (unit, '(A)') TRIM(EQN_TAGS(jr)) 
  ENDDO
  CLOSE(unit)

  ! save EQN_NAMES to file:
  OPEN(unit, FILE='EQN_NAMES.dat', status='UNKNOWN')
  DO jr = 1, NREACT
    WRITE (unit, '(A)') TRIM(EQN_NAMES(jr)) 
  ENDDO
  CLOSE(unit)

  DEALLOCATE(OICdata)

END PROGRAM oic

!*****************************************************************************
