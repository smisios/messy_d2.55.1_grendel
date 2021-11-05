PROGRAM PLOT_GRID
  USE mo_kind, ONLY: i8!, i4, sp, dp
  USE iso_varying_string
  IMPLICIT NONE
  INTEGER :: ie, je, ito, jto
  CHARACTER*12 CC
  INTEGER :: igg, i, icount, ii1, ii2, ii3, ii4, ioerrstat, &
       it, j, verbose
  INTEGER, PARAMETER :: iunit=100, io_in_anta=80, io_out_anta_diag_l=50, &
       io_out_anta_diag_s=51, io_in_topo=81, io_out_topo=82, &
       io_out_ribek=68, io_out_topo_new=66
  INTEGER(i8) :: ext8_hdr(4)
  REAL, ALLOCATABLE :: ribek(:,:)
  REAL, ALLOCATABLE :: GILA(:,:),GIPH(:,:)
  REAL, ALLOCATABLE :: DEPTO(:,:),  &
       alon(:,:),ALAT(:,:)
  REAL, PARAMETER :: GRARAD=180./3.1415927
  REAL :: ala1, ala2, ala3, ala4, alo1, alo2, alo3, alo4, tt
  TYPE(varying_string) :: grid, depto_edesc, depto_format, anta_fname, &
       topo_fname
  NAMELIST /gridparams/ ie, je, verbose

  ie = -1
  je = -1
  verbose = 0
  READ (*, nml=gridparams)
  IF (verbose .GT. 0) PRINT *, 'ie=', ie, 'je=', je
  IF (ie < 1 .OR. je < 1) THEN
    PRINT *, 'Only positive integral values may be used for ie and je!'
    STOP
  END IF
  CALL get(grid, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit
  CALL get(depto_edesc, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit
  CALL get(anta_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. anta_fname == '') anta_fname=CHAR(grid)//'_anta'
  CALL get(topo_fname, iostat=ioerrstat)
  IF (ioerrstat .GT. 0 .OR. topo_fname == '') topo_fname=CHAR(grid)//'_topo'

  ito = 2 * ie
  jto = 2 * je

  OPEN(io_in_anta, file=CHAR(anta_fname), &
       access='sequential', form='unformatted', action='read')

  OPEN(io_out_anta_diag_l, file='anta.'//CHAR(grid)//'.l.txt', &
       access='sequential', form='formatted', action='write')
  OPEN(io_out_anta_diag_s, file='anta.'//CHAR(grid)//'.s.txt', &
       access='sequential', form='formatted', action='write')

  !      OPEN(io_in_topo,FILE=char(grid)//'_topo.ext'
  !     &   ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

  OPEN(io_in_topo,FILE=CHAR(topo_fname), access='sequential', &
       form='formatted', action='read')

  OPEN(io_out_topo,FILE=CHAR(grid)//'_topo2', access='sequential', &
       form='formatted', action='write')

  IF (verbose .GT. 0) PRINT *, 'read anta'

  ALLOCATE(GILA(ito,jto))


  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) gila
  IF (verbose .GT. 0) PRINT *, 'read anta field gila'

  ALLOCATE(GIPH(ito,jto))
  READ(io_in_anta) ext8_hdr
  READ(io_in_anta) giph
  IF (verbose .GT. 0) PRINT *, 'read anta field giph'

  ALLOCATE(depto(ie,je),alon(ie,je),alat(ie,je),ribek(ie,je))

  depto(:,:)=0.
  depto_format='(I5,20'//CHAR(depto_edesc)//')'
  DO iI1=2,IE-1,20
    iI2=MIN(iI1+19,IE-1)
    READ(io_in_topo,*,END=99)cc
    DO J=1,JE
      READ(io_in_topo,fmt=CHAR(depto_format))igg,(DEPTO(I,J),I=iI1,iI2)
      IF (verbose .GT. 2) PRINT *,'read topo', igg
    END DO
  END DO
99 CONTINUE
  !      read(io_in_topo)il1
  !      read(io_in_topo)((depto(i,j),i=1,362),j=3,180)
  !      read(io_in_topo) depto



  !      DO i=1,ie
  !         il=i
  !         ir=ie+1-i
  !         depto(il,2) = depto(ir,3)                            ! syncronise line 2 with line 3
  !         depto(il,1) = depto(ir,4)                            ! syncronise line 1 with line 4
  !      END DO



  do i=1,ie
    do j=1,je
      if (depto(i,j).le. 0.) depto(i,j)=0.
    enddo
  enddo

  DO iI1=2,IE-1,20
    iI2=MIN(iI1+19,IE-1)
    WRITE(io_out_topo,*)ii1,ii2
    DO J=1,JE
      WRITE(io_out_topo, fmt=CHAR(depto_format)) j, (DEPTO(I,J),I=iI1,iI2)
    ENDDO
  ENDDO

  do i=1,ie
    do j=1,je
      ribek(i,j)=0.
    enddo
  enddo

  ribek(1,:)=1.


  DO it=1,2000

    do i=2,ie-1
      do j=2,je-1

        if (ribek(i,j).lt.0.5)then

          if (ribek(i-1,j).gt.0.5 &
               .or.ribek(i+1,j).gt.0.5 &
               .or.ribek(i,j+1).gt.0.5 &
               .or.ribek(i,j-1).gt.0.5)then
            !      mm=0
            !      if (depto(i-1,j).gt.0.5) mm=mm+1
            !      if (depto(i+1,j).gt.0.5) mm=mm+1
            !      if (depto(i,j-1).gt.0.5) mm=mm+1
            !      if (depto(i,j+1).gt.0.5) mm=mm+1
            !      if (mm.eq.4) ribek(i,j)=1

            if (depto(i,j).gt.0.5) ribek(i,j)=1.

          endif
        endif

      enddo
    enddo

    icount=0
    do i=2,ie-1
      do j=2,je-1
        if (ribek(i,j).gt.0.5) icount=icount+1
      enddo
    enddo

    IF (verbose .GT. 2) PRINT *, it, icount

  enddo


  ii3=0
  ii4=IE*JE
  ext8_hdr(1) = INT(ii1, i8)
  ext8_hdr(2) = INT(ii2, i8)
  ext8_hdr(3) = 0_i8
  ext8_hdr(4) = INT(ie * je, i8)
  OPEN(io_out_ribek, file='ribek.ext', form='unformatted', action='write')
  WRITE(io_out_ribek) ext8_hdr
  WRITE(io_out_ribek) ribek
  CLOSE(io_out_ribek)

  do i=1,ie
    do j=3,je
      if (ribek(i,j).lt.0.5) depto(i,j)=0.
    enddo
  enddo

  DO  iI1=2,IE-1,20
    iI2=MIN(iI1+19,IE-1)
    OPEN(io_out_topo_new,file=CHAR(grid)//'_topo.NEW',form='formatted', &
         action='write')
    WRITE (io_out_topo_new,*) 'STREIFEN ',iI1,iI2
    DO J=1,JE
      WRITE(io_out_topo_new, fmt=CHAR(depto_format)) j, (DEPTO(I,J),I=iI1,iI2)
    ENDDO
  ENDDO

  !      OPEN(81,FILE='weto.'//char(grid)//'.ext4'
  !     &   ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
  !      READ(81) I1,i2,i3,i4
  !      READ(81) DEPTO


  DO I=1,IE
    DO j=1,jE
      !       II(i,j)=i
      !       JJ(i,j)=j
      alat(I,j)=GRARAD*GIPH(2*I,2*j)
      alon(I,j)=GRARAD*GILA(2*I,2*j)
    enddo
  enddo

  !       WRITE(10)I1,1,i3,i4
  !       WRITE(10)II
  !       WRITE(10)I1,2,i3,i4
  !       WRITE(10)JJ





  DO I=2,IE-1
    DO j=1,jE

      if (DEPTO(I,J).le.0.5) then
        WRITE(io_out_anta_diag_l, '(1X,F8.2,1X,F8.2,1X,F6.0)') &
             alon(I,j), alat(I,j), DEPTO(I,j)
      else

        WRITE(io_out_anta_diag_s, '(1X,F8.2,1X,F8.2,1X,F6.0)') &
             alon(I,j), alat(I,j), DEPTO(I,j)
      endif

    enddo
  enddo


  DO I=2,IE-1,20
    WRITE(iunit,'(1a,1X,F8.2,1X,F8.2,1X,F6.0)')'>',alon(I,1) &
         ,alat(I,1) &
         ,DEPTO(I,1)
    DO j=1,jE,20
      WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')alon(I,j),alat(I,j) &
           ,DEPTO(I,j)
    enddo
  enddo

  DO j=1,jE,20
    WRITE(iunit,'(1A,1X,F8.2,1X,F8.2,1X,F6.0)')'>',alon(1,j) &
         ,alat(1,j) &
         ,DEPTO(1,j)
    DO I=1,IE,20
      WRITE(iunit,'(1X,F8.2,1X,F8.2,1X,F6.0)')alon(I,j),alat(I,j) &
           ,DEPTO(I,j)
    enddo
  enddo



  DO I=2,IE-1
    DO j=2,jE-1
      tt=DEPTO(I,J)
      if (tt.le.0.5) then

        ALA1=GRARAD*GIPH(2*I-1,2*j-1)
        ALA2=GRARAD*GIPH(2*I+1,2*j-1)
        ALA3=GRARAD*GIPH(2*I+1,2*j+1)
        ALA4=GRARAD*GIPH(2*I-1,2*j+1)

        ALo1=GRARAD*GIla(2*I-1,2*j-1)
        ALo2=GRARAD*GIla(2*I+1,2*j-1)
        ALo3=GRARAD*GIla(2*I+1,2*j+1)
        ALo4=GRARAD*GIla(2*I-1,2*j+1)



        WRITE(8,'(1A,1X,F8.2,1X,F8.2)') '>',ALO1,ALA1
        WRITE(8,'(1X,F8.2,1X,F8.2)') ALO2,ALA2
        WRITE(8,'(1X,F8.2,1X,F8.2)') ALO3,ALA3
        WRITE(8,'(1X,F8.2,1X,F8.2)') ALO4,ALA4
        WRITE(8,'(1X,F8.2,1X,F8.2)') ALO1,ALA1
      endif
    enddo
  enddo







END PROGRAM PLOT_GRID

SUBROUTINE err_exit
  PRINT *, 'io error occurred!, aborting'
  STOP 1
END SUBROUTINE err_exit
