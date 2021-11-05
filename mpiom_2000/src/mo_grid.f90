MODULE mo_grid
  USE mo_kind, ONLY: i8, i4, dp, wp
  USE mo_param1, ONLY: ie, ie1, je, je1, ke, ie_g, je_g
  USE mo_parallel, ONLY: p_bcast, have_g_je, p_io, p_pe, p_ioff, p_joff, &
       gather, global_max, global_min, scatter, write_slice_sp, &
       have_g_js, read_slice, stop_all
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY: alat, alat_g, alatu, alatv, alatpsi_g &
       ,alon, alon_g, alonu, alonv, alonpsi_g &
       ,amsue, amsuo, area, areain, curvav, &
       ddpo, ddpsio, dduo, ddue, depto, deutie, deutio, deute, deuto, &
       dlxp, dlxp_g, dlxpsi, dlxu, dlxu_g, dlxv, dlxv_g, dlyp, dlyp_g, &
       dlypsi, dlyu, dlyu_g, dlyv, dlyv_g, dpio, dzw, ftwou, ftwov, &
       dlxpi, dlypi, dlxui, dlyui, dlxvi, dlyvi, &
       dtdxuo, dtdxpo, dtdyo, dpyo, dtdxue, dtdxpe, dpye, dt, &
       gila, giph, icontro, istart, istart_new_topo_update, &
       kcondep, lbounds_exch_tp, &
#if !defined WETO_STORE_LOGICAL
       set_wetol1_g, &
#endif
       one, tiestu, tiestw, weto, wetol1_g, zero
  USE mo_constants, ONLY: aradtogra
  USE mo_io_config, ONLY: next_free_unit
  USE MO_UNITS, ONLY: io_ou_alat, io_ou_alatu, io_ou_alatv, &
       io_ou_alon, io_ou_alonu, io_ou_alonv, io_ou_amsue, io_ou_amsuo, &
       io_ou_bek, io_ou_ddpo, io_ou_ddue, io_ou_dduo, io_ou_dept, &
       io_ou_deute, io_ou_deuto, io_ou_dlxp, io_ou_dlxu, io_ou_dlxv, &
       io_ou_dlyp, io_ou_dlyu, io_ou_dlyv, io_ou_weto,  io_stdout
  USE mo_grid_elementals, ONLY: suchij_2d, grid_dist_2d
  USE mo_mpi, ONLY: p_all_comm, generate_mpi_struct_type, create_mpi_op
  IMPLICIT NONE
  INCLUDE 'mpif.h'

  PRIVATE

  REAL(wp), ALLOCATABLE,TARGET :: thkcello(:,:,:) !< cell_thickness

  LOGICAL :: lwith_one_layer_shelfs = .FALSE.
  REAL(wp), ALLOCATABLE :: SHASWO(:,:),SHASWE(:,:), & ! masks for velocity points with only one wet level (zero if true )
                           shamsuo(:,:), shamsue(:,:) ! land sea masks on vector points merged with shaswo/shaswe

  LOGICAL :: write_gridinfo_output = .FALSE.
  INTEGER :: p_grid_dist_2d, p_op_grid_dist_min_2d

  PUBLIC :: coriol, boden, wrte_gridinfo, setup_grid, p_suchij, &
       get_level_index_by_depth,lwith_one_layer_shelfs,shaswo,shaswe, &
       shamsuo, shamsue, thkcello, cell_thickness

#ifdef MESSY
   REAL, ALLOCATABLE :: GILA_G(:,:),GIPH_G(:,:) 
   REAL, ALLOCATABLE :: WETO_G(:,:,:) 
   PUBLIC :: GILA_G,GIPH_G,WETO_G
#endif

CONTAINS

  SUBROUTINE CORIOL
    !********************************************************************
    !
    !
    !      CCCCC   OOO   RRRRR   II   OOO   L
    !     C       O   O  R    R  II  O   O  L
    !     C       O   O  RRRRR   II  O   O  L
    !     C       O   O  R  RR   II  O   O  L
    !      CCCCC   OOO   R   RR  II   OOO   LLLLLL
    !
    !
    !********************************************************************
    !
    ! CORIOLIS-PARAMETER
    ! TIEFE*FAKTOR(CORIOLIS, ZEITSCHRITT)
    !------------------------------------------------------------
    INTEGER i,j,ioff,joff
    INTEGER(kind=i8) :: ibla(4)
    INTEGER :: io_in_anta
    !      write(io_stdout,*)'in coriol'

    REAL(wp), ALLOCATABLE :: gi_g(:,:)

    ALLOCATE(gi_g(2*ie_g,2*je_g))

#ifdef MESSY
    ALLOCATE(GILA_G(2*ie_g,2*je_g))
    ALLOCATE(GIPH_G(2*ie_g,2*je_g))
#endif

    IF(p_pe==p_io) THEN
      io_in_anta = next_free_unit()
#ifdef LITTLE_ENDIAN
#ifndef NOENDIANCONVERT
      OPEN(io_in_anta, file='anta', action='read', &
           access='sequential',form='unformatted',convert='big_endian')
#else
      ! ERROR: compiler does not support convert='big_endian'
#endif
#else
      OPEN(io_in_anta, file='anta', action='read', &
           access='sequential',form='unformatted')
#endif
      READ(io_in_anta) ibla
      READ(io_in_anta) gi_g
#ifdef MESSY
      GILA_G(:,:) = gi_g(:,:)
#endif
    ENDIF
 

#ifdef MESSY
    CALL p_bcast(GILA_G,p_io)
#endif
    CALL p_bcast(gi_g,p_io)

    ioff = 2*p_ioff
    joff = 2*p_joff
    gila(:,:) = gi_g(ioff+1:ioff+2*ie,joff+1:joff+2*je)

    IF(p_pe==p_io) THEN
       do j=0,je_g-1
          DO i=1,ie_g-1
             alonpsi_g(i,j)=gi_g(2*i+1,2*j+1)*aradtogra
          ENDDO
       enddo

       DO i=1,ie_g-1
          alonpsi_g(i,je_g)=gi_g(2*i+1,2*je_g)*aradtogra
       ENDDO

       alonpsi_g(ie_g,:)=alonpsi_g(1,:)
    ENDIF

    IF(p_pe==p_io) THEN
      READ(io_in_anta) ibla
      READ(io_in_anta) gi_g
      CLOSE(io_in_anta)
#ifdef MESSY
      GIPH_G(:,:) = gi_g(:,:) 
#endif
    ENDIF

#ifdef MESSY
      CALL p_bcast(GIPH_G,p_io)
#endif
    CALL p_bcast(gi_g,p_io)

    ioff = 2*p_ioff
    joff = 2*p_joff
    giph(:,:) = gi_g(ioff+1:ioff+2*ie,joff+1:joff+2*je)


    IF(p_pe==p_io) THEN
       do j=0,je_g-1
          DO i=1,ie_g-1
             alatpsi_g(i,j)=gi_g(2*i+1,2*j+1)*aradtogra
          ENDDO
       enddo

       DO i=1,ie_g-1
          alatpsi_g(i,je_g)=gi_g(2*i+1,2*je_g)*aradtogra
       ENDDO

       alatpsi_g(ie_g,:)=alatpsi_g(1,:)
    ENDIF

    DEALLOCATE(gi_g)

    alonv(:,1) = zero
    alonv(1,:) = zero
    alonv(ie,:) = zero
    alonv(:,je) = zero
    DO j=1,je
      DO i=2,ie-1
        alat(i,j)=giph(2*i,2*j)*aradtogra
        alon(i,j)=gila(2*i,2*j)*aradtogra
        alatu(i,j)=giph(2*i+1,2*j)*aradtogra
        alonu(i,j)=gila(2*i+1,2*j)*aradtogra
      ENDDO
    ENDDO
    do j=1,je-1
      do i=2,ie-1
        alatv(i,j)=giph(2*i,2*j+1)*aradtogra
        alonv(i,j)=gila(2*i,2*j+1)*aradtogra
      enddo
    enddo


    call bounds_exch(1,'p',alat,'coriol 1')
    call bounds_exch(1,'p',alon,'coriol 1')
    call bounds_exch(1,'u+',alatu,'coriol 1')
    call bounds_exch(1,'u+',alonu,'coriol 1')


    if (have_g_je) then
      alatv(:,je)=alat(:,je)
      alonv(:,je)=alon(:,je)
    endif


    call bounds_exch(1,'v+',alatv,'coriol 1')
    call bounds_exch(1,'v+',alonv,'coriol 1')

    call gather(alon,alon_g,p_io)
    call gather(alat,alat_g,p_io)


  END SUBROUTINE CORIOL




  SUBROUTINE BODEN
    !*************************************************************
    !
    !     BBBBB    OOO   DDDDD   EEEEEE  NN   N
    !     B    B  O   O  D    D  E       NNN  N
    !     BBBBB   O   O  D    D  EEEEE   N NN N
    !     B    B  O   O  D    D  E       N  NNN
    !     BBBBB    OOO   DDDDD   EEEEEE  N   NN
    !
    !
    !--------------------------------------------------------------
    CHARACTER*6 BLORD

    REAL(wp) DEPTO_G(IE_G,JE_G)
    REAL(wp) DEPALT,DEPMAX,DEPMIN
    REAL(wp) SUGGO, SUPPO, TMAX
    INTEGER I,I1,I2,J,JEV,JJJ,K
    INTEGER ir, io_in_gigi
#ifndef SOR
    INTEGER LAUF, III, JJ, NDIF1
#endif
    ! diagnostic mask (only active with lbounds_exch_tp == .true.
    ! .and. icontro > 0 )
    REAL(wp), ALLOCATABLE :: rmus(:,:)

    ALLOCATE(SHASWO(IE,JE),SHASWE(IE,JE))
    ALLOCATE(shamsuo(ie,je), shamsue(ie,je))

    !------------------------------------------------------------------------
    !
    !   HERE   :  READ TOPOGRAPHY FROM EXTERNAL FILE ON ARRAY DEPTH
    !
    !      (BE CAREFUL WITH THE SIGN OF THE TOPOGR. DATA ! , THIS VERSION
    !       USES POSITIVE DEPTH DATA BECAUSE SYNBAPS COMES WITH POSITIVE
    !       VALUES)
    !
    DEPMAX=SUM(dzw(1:ke))

    IF (lwith_one_layer_shelfs) THEN
      DEPMIN=DZW(1)
    ELSE
      depmin = dzw(1) + dzw(2) + 1._wp
    ENDIF

    !
    WRITE(IO_STDOUT,*)' --------------------------------------------- &
         &                    --------'
    WRITE(IO_STDOUT,*)' MINIMUM WATER DEPTH IS SET TO ',DEPMIN
    WRITE(IO_STDOUT,*)' MAXIMUM WATER DEPTH IS SET TO ',DEPMAX
    WRITE(IO_STDOUT,*)' --------------------------------------------- &
         &                    --------'
    !
    depto_g(:,:) = 0.0_wp

    IF (istart .EQ. istart_new_topo_update) THEN
      j_loop: DO J=1,JE
        i_loop: DO I=1,IE
          !
          IF (depto(i, j) .LE. 1._wp) depto(i, j) = 0._wp
          IF (depto(i, j) .GT. 1._wp .AND. depto(i, j) .LT. depmin)         &
               DEPTO(I,J) = DEPMIN
        END DO i_loop
      END DO j_loop
    ENDIF

    IF (lwith_one_layer_shelfs) THEN
      WHERE (depto(:,:).GT.depmin .AND.                               &
           depto(:,:).LE.tiestu(2)) depto(:,:)=depmin
    ENDIF

    !
    tmax = 0._wp
    DO J=1,JE
      DO I=1,IE
        TMAX=MAX(TMAX,DEPTO(I,J))
      ENDDO
    ENDDO
    CALL global_max(TMAX)
    WRITE(IO_STDOUT,*)'TOPOGRAPHIE CHECK: ',TMAX
    !
!    IF (ISTART .EQ. istart_new_topo_update) CALL gather_arr(DEPTO,DEPTO_G,p_io)
    IF (istart .EQ. istart_new_topo_update) CALL gather(depto, depto_g, p_io)
    !
    IF (p_pe==p_io) THEN
      io_in_gigi = next_free_unit()
#ifdef LITTLE_ENDIAN
#ifndef NOENDIANCONVERT
      OPEN(IO_IN_GIGI,FILE='topo'                                    &
           &               ,ACCESS='SEQUENTIAL',FORM='FORMATTED', convert='big_endian')
#else
      ! ERROR: compiler does not support convert='big_endian'
#endif
#else
         OPEN(IO_IN_GIGI,FILE='topo'                                    &
     &               ,ACCESS='SEQUENTIAL',FORM='FORMATTED')
#endif
      DO I1=2,IE_G-1,20
        I2=MIN(I1+19,IE_G-1)
        WRITE(IO_STDOUT,*) 'LESSTREIFEN ',I1,I2
        IF (istart .EQ. istart_new_topo_update) THEN
          WRITE(IO_IN_GIGI,*)'STREIFEN ',I1,I2
          DO J=1,JE_G
            WRITE(IO_IN_GIGI,63885)J,(NINT(DEPTO_G(I,J)),I=I1,I2)
          ENDDO
        ELSE

          READ(IO_IN_GIGI,*) BLORD

          DO J=1,JE_G
            IF ( lbounds_exch_tp ) THEN
              READ(IO_IN_GIGI,63886)JJJ,(DEPTO_G(I,J),I=I1,I2)
            ELSE
              READ(IO_IN_GIGI,63885)JJJ,(DEPTO_G(I,J),I=I1,I2)
            ENDIF
          ENDDO

          IF ( lbounds_exch_tp ) THEN
            depto_g(:,je_g) = 0._wp
            depto_g(:,je_g-1) = 0._wp
          ENDIF

        ENDIF

63886   FORMAT(I5,20F6.0)    !format for tp10/tp04/tp01
63885   FORMAT(I5,20F5.0)

      ENDDO

      CLOSE(IO_IN_GIGI)

    ENDIF


    IF ( lbounds_exch_tp .and. icontro > 0 ) THEN
      IF ( p_pe == p_io ) THEN
        ALLOCATE(rmus(ie_g,je_g))
        rmus = 0.0_wp
      ELSE
        ALLOCATE(rmus(0,0))
      ENDIF

      IF ( p_pe == p_io ) THEN
        WHERE (depto_g > 0.5_wp) rmus = 13.0_wp
        DO i=2,ie_g/2
          ir=ie_g+1-i
          WRITE(IO_STDOUT,6637)i,ir,(depto_g(i,j),j=12,3,-1),(depto_g(ir,j),j=3,12)
6637      FORMAT(2i4,10f5.0,'|',10f5.0)
        ENDDO
        IF (je_g .GE. 32 ) THEN
          DO i=2,ie_g/2
            ir=ie_g+1-i
            WRITE(IO_STDOUT,6638)i,ir,(rmus(i,j),j=32,3,-1),(rmus(ir,j),j=3,32)
6638        FORMAT(2i3,32f1.0,'|',32f1.0)
          ENDDO
        ENDIF
      ENDIF

      DEALLOCATE(rmus)

    ENDIF


    IF (istart .NE. istart_new_topo_update) CALL scatter(DEPTO_G,DEPTO,p_io)

    !     FLACHE TEILE DER TOPOGRAPHIE KORRIGIEREN
    !

    DO J=1,JE
      DO I=1,IE
        IF (depto(i, j) .GT. 0.5_wp) THEN
          DEPALT=DEPTO(I,J)

          IF (lwith_one_layer_shelfs) THEN

            depto(i,j)=MIN(depto(i,j),depmax)
            IF (depto(i, j) .GT. 1._wp .AND. depto(i, j) .LT. depmin) THEN
              depto(i, j) = depmin
            ENDIF
            IF(depto(i, j) .GT. depmin .AND. depto(i, j) .LE. tiestu(2)) THEN
              depto(i, j) = depmin
            ENDIF

          ELSE

            DEPTO(I,J)=MAX(DEPTO(I,J),DZW(1)+DZW(2))
            depto(i,j)=MIN(depto(i,j),depmax)
          ENDIF

          IF ( icontro > 0 ) THEN
            IF (depalt.NE.depto(i,j)) THEN
              WRITE(IO_STDOUT,*)'TOPOGRAPHY CHANGED AT: ' &
                   ,I+p_ioff,J+p_joff,DEPALT,'==>',DEPTO(I,J)
            ENDIF
          ENDIF

        ELSE
          depto(i, j) = 0._wp
        ENDIF
      ENDDO
    ENDDO

    DO J=1,JE
      DO I=1,IE
        DO K=1,KE
          IF (ABS(depto(i, j) - tiestu(k)) .LT. 0.5_wp) &
               depto(i, j) = tiestu(k) + 1._wp
        ENDDO
      ENDDO
    ENDDO

    CALL bounds_exch(1,'p',DEPTO,'boden 10')
    !
    !     Now the calculations on DEPTO are finished, we gather and
    !     distribute global DEPTO_G (which is also needed for WETO_G)
    !
!    CALL gather_arr(DEPTO,DEPTO_G,0)
    CALL gather(DEPTO,DEPTO_G,0)
    CALL p_bcast(DEPTO_G,0)

    DO J=1,JE1
      DO I=1,IE1
        DEUTE(I,J)=MIN(DEPTO(I,J),DEPTO(I,J+1))
        DEUTO(I,J)=MIN(DEPTO(I,J),DEPTO(I+1,J))
      ENDDO
    ENDDO

    CALL bounds_exch(1,'v+',DEUTE,'boden 11')
    CALL bounds_exch(1,'u+',DEUTO,'boden 12')


    ! option lwith_one_layer_shelfs
    ! mask fields for grid points
    ! with only 1 wet cell in the vertical

    shaswo(:,:) = 1._wp
    shaswe(:,:) = 1._wp

    IF (lwith_one_layer_shelfs) THEN
    WHERE(deuto(:,:) .LT. depmin) shaswo(:,:) = 0._wp
    WHERE(deute(:,:) .LT. depmin) shaswe(:,:) = 0._wp

    WHERE(1.0_wp .LT. deuto(:,:) .AND.                             &
         deuto(:,:) .LT. tiestu(2)) shaswo(:,:) = 0._wp
    WHERE(1.0_wp .LT. deute(:,:) .AND.                             &
         deute(:,:) .LT. tiestu(2)) shaswe(:,:) = 0._wp
    endif


    !----------------------------------------------------------------------
    !     COMPUTE BOUNDARY TOPOGRAPHY , I.E. 3 ROWS
    !     FROM THE NORTHERN AND FROM THE SOUTHERN BOUNDARY ARE SET
    !     TO ZERO DEPTH
    !
    !     CALCULATION OF DEPTH AT SCALAR POINTS (HP : DEPTH AT PRESSURE PT.)
    !        HP IS THE MAXIMUM DEPTH OF THE 4 SURROUNDING U-POINT DEPTHS (H)

#ifndef SOR
    LAUF=0
    JJ=0
    NUM_G(:,:) = 0
    DO III=2,IE_G-1
      JJ=JJ+1
      IF(MOD(III,2).EQ.0)THEN
        I=III/2+1
      ELSE
        I=IE_G-(III/2)
      ENDIF
      DO J=1,JE_G
        IF(DEPTO_G(I,J).GE.1.) THEN
          LAUF=LAUF+1
          NUM_G(I,J)=LAUF
        ENDIF
      ENDDO
    ENDDO
    IF(ICYCLI.GE.1)THEN
      DO J=1,JE_G
        NUM_G(1,J)=NUM_G(IE_G-1,J)
        NUM_G(IE_G,J)=NUM_G(2,J)
      ENDDO
    ENDIF

    !     Set local NUM

    NUM(1:ie,1:je) = NUM_G(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)
    !
    WRITE(IO_STDOUT,*)' FEUCHTPUNKTE',LAUF
    !
    !============WERTE FUER TRIANGULARISIERUNG
    MATR=LAUF
    NMAX=0
    !
    DO J=2,JE_G-1
      DO I=2,IE_G-1
        IF(NUM_G(I,J).EQ.0 .OR. NUM_G(I-1,J).EQ.0) CYCLE
        NDIF1=NUM_G(I,J)-NUM_G(I-1,J)
        IF(NDIF1.GT.NMAX) NMAX=NDIF1
      ENDDO
    ENDDO
    !
    WRITE(IO_STDOUT,*)' MAXDIF ', NMAX

#endif
    !======LAND/OZEAN-STEUERFELDER UND SCHICHTDICKEN
    !
    kloop: DO K=1,KE
      jloop: DO J=1,JE
        JEV=2*J
        iloop: DO I=1,IE
          !------------------------------  EVEN LINES 2,4,6,...,2*JE ---------
          !
          IF(TIESTU(K+1).LT.DEUTE(I,J)) DDUE(I,J,K)=DZW(K)
          IF(TIESTU(K).LT.DEUTE(I,J) .AND. TIESTU(K+1).GE. DEUTE(I,J))      &
               &    DDUE(I,J,K)=DEUTE(I,J)-TIESTW(K)
          !
          IF (ddue(i, j, k) .GT. zero) amsue(i, j, k) = 1._wp
          !====================================================================
          !
          !------------------------------   ODD LINES 1,3,5,...,2*JE-1 -------
          !
          IF(TIESTU(K+1).LT.DEUTO(I,J)) DDUO(I,J,K)=DZW(K)
          IF(TIESTU(K).LT.DEUTO(I,J) .AND. TIESTU(K+1).GE. DEUTO(I,J))      &
               &    DDUO(I,J,K)=DEUTO(I,J)-TIESTW(K)
          IF (dduo(i, j, k) .GT. zero) amsuo(i, j, k) = 1._wp
          IF (tiestu(k) .LT. depto(i, j)) weto(i, j, k) = 1._wp
          IF(TIESTU(K+1).LT.DEPTO(I,J)) DDPO(I,J,K)=DZW(K)
          IF(TIESTU(K).LT.DEPTO(I,J) .AND. TIESTU(K+1).GE. DEPTO(I,J))      &
               &    DDPO(I,J,K)=DEPTO(I,J)-TIESTW(K)
          !
          IF(DDPO(I,J,K) .EQ. ZERO) cycle iloop
          dpio(i, j, k) = 1._wp / ddpo(i, j, k)
        END DO iloop
      END DO jloop
    END DO kloop
#ifdef WETO_STORE_LOGICAL
    lwetol1_g(:,:) = TIESTU(1) .LT. DEPTO_G(:,:)
#elif defined WETO_STORE_BITVECTOR
    bvwetol1_g(:) = 0
#else
    WETOL1_G(:,:) = 0._wp
#endif
#if !defined WETO_STORE_LOGICAL
    K=1
    DO J=1,JE_G
      DO I=1,IE_G
        IF(TIESTU(K).LT.DEPTO_G(I,J)) CALL set_wetol1_g(i, j, .TRUE.)
      ENDDO
    ENDDO
#endif

#ifdef MESSY
    ALLOCATE(WETO_G(IE_G,JE_G,KE))   
    WETO_G(:,:,:) = 0.
    DO K=1,KE
    DO J=1,JE_G
      DO I=1,IE_G
        IF(TIESTU(K).LT.DEPTO_G(I,J)) WETO_G(I,J,K)=1.
      ENDDO
    ENDDO
    ENDDO
#endif

    ddpsio(:,:,:) = 0._wp


    DO K=1,KE
      DO J=1,JE1
        DO I=1,IE1

          if ( icontro > 0 ) then
            IF (ddue(i, j, k) .GT. 0._wp &
                 .AND. ddpo(i, j, k) * ddpo(i, j+1, k) .LT. 1._wp) THEN
              WRITE(IO_STDOUT,*) ' DDUE ',I+p_ioff,J+p_joff,K   &
                   ,DDUE(I,J,K),DDPO(I,J,K),DDPO(I,J+1,K)
            ENDIF

            IF (dduo(i, j, k) .GT. 0._wp &
                 .AND. ddpo(i, j, k) * ddpo(i+1, j, k) .LT. 1._wp) THEN
              WRITE(IO_STDOUT,*) ' DDUO ',I+p_ioff,J+p_joff,k &
                   ,DDUO(I,J,K),DDPO(I,J,K),DDPO(I+1,J,K)
            ENDIF
          endif

          suppo = 0._wp
          suggo = 0._wp
          IF (ddpo(i, j, k) .GT. 1._wp) THEN
            SUPPO=SUPPO+DDPO(I,J,K)
            suggo = suggo + 1._wp
          ENDIF
          IF (ddpo(i+1, j, k) .GT. 1._wp) THEN
            SUPPO=SUPPO+DDPO(I+1,J,K)
            suggo = suggo + 1._wp
          ENDIF
          IF (ddpo(i, j+1, k) .GT. 1._wp) THEN
            SUPPO=SUPPO+DDPO(I,J+1,K)
            suggo = suggo + 1._wp
          ENDIF
          IF (ddpo(i+1, j+1, k) .GT. 1._wp) THEN
            SUPPO=SUPPO+DDPO(I+1,J+1,K)
            suggo = suggo + 1._wp
          ENDIF
          IF (suggo .GE. 1._wp) ddpsio(i, j, k) = suppo/suggo

        enddo
      enddo
    enddo


    CALL bounds_exch(1,'s',DDPSIO,'boden 13')
    !
    weto_j_loop: DO J=1,JE
      weto_i_loop: DO I=1,IE
        KCONDEP(I,J)=NINT(WETO(I,J,1))-99*NINT(1._wp - weto(i, j, 1))
        !
        deutie(i, j) = 0._wp
        deutio(i, j) = 0._wp
        !
        IF (deute(i, j) .GT. one) deutie(i, j) = 1._wp / deute(i, j)
        IF (deuto(i, j) .GT. one) deutio(i, j) = 1._wp / deuto(i, j)
        !
      END DO weto_i_loop
    END DO weto_j_loop

    ! merge amsuo/amsue and shaswo/shaswe masks to use them in occlit
    shamsuo(:,:) = MERGE(amsuo(:,:,1), shaswo(:,:), amsuo(:,:,1) .LT. shaswo(:,:))
    shamsue(:,:) = MERGE(amsue(:,:,1), shaswe(:,:), amsue(:,:,1) .LT. shaswe(:,:))

    CALL bounds_exch(1,'v+',DEUTIE,'boden 14')
    CALL bounds_exch(1,'u+',DEUTIO,'boden 15')
    CALL bounds_exch(1,'v+',DEUTE,'boden 16')
    CALL bounds_exch(1,'u+',DEUTO,'boden 17')
    CALL bounds_exch(1,'u+',DDUO,'boden 18')
    CALL bounds_exch(1,'v+',DDUE,'boden 19')
    CALL bounds_exch(1,'v+',AMSUE,'boden 20')
    CALL bounds_exch(1,'u+',AMSUO,'boden 21')
    CALL bounds_exch(1,'v+',SHASWE,'boden 26')
    CALL bounds_exch(1,'u+',SHASWO,'boden 27')
    CALL bounds_exch(1,'v+',SHAMSUE,'boden 28')
    CALL bounds_exch(1,'u+',SHAMSUO,'boden 29')
    !#ifdef bounds_exch_save
    CALL bounds_exch(1,'p',WETO,'boden 22')
    CALL bounds_exch(1,'p',DDPO,'boden 23')
    !#endif
    !=====================================================================
    WRITE(IO_STDOUT,*)'SR BODEN UEBERLEBT!!!'
    !#ifdef bounds_exch_save
    CALL bounds_exch(1,'p',DLXP,'boden 24')
    CALL bounds_exch(1,'p',DLYP,'boden 25')
    !#endif
    !

    CALL ini_cell_thickness ! initialise thkcello

  END SUBROUTINE BODEN



  SUBROUTINE setup_arcgri
    INTEGER :: i, j
    INTEGER :: io_in_arcg
#if defined(__xlC__) || defined(POOL_DATA_ARCGRI_IS_FIXED_TO_8BYTE_WORD_HEADER)
    INTEGER(kind=i8) :: ext_hdr(4)
#else
    INTEGER(kind=i4) :: ext_hdr(4)
#endif

    REAL(wp) :: dmino, dmaxo

!mz_ap_20100927+
#ifndef LITTLE_ENDIAN
    IF (p_pe==p_io) THEN
      io_in_arcg = next_free_unit()
      OPEN(IO_IN_ARCG,FILE='arcgri', &
           ACCESS='SEQUENTIAL',FORM='UNFORMATTED', ACTION='read')
    write (*,*) "file opened -> big endian"
    END IF
#else
#ifndef NOENDIANCONVERT
    IF(p_pe==p_io) THEN
      io_in_arcg = next_free_unit()
      OPEN(IO_IN_ARCG,FILE='arcgri', &
           ACCESS='SEQUENTIAL',FORM='UNFORMATTED', ACTION='read', CONVERT='BIG_ENDIAN')
    write (*,*) "file opened -> little endian"
    END IF
#else
    ! ERROR: compiler does not support convert='big_endian'
#endif
#endif

    deuto(:,:) = 0._wp
    deute(:,:) = 0._wp
    depto(:,:) = 0._wp

    IF (p_pe==p_io) READ(IO_IN_ARCG) EXT_HDR
    CALL read_slice(IO_IN_ARCG,DEUTO)

    IF (istart .EQ. istart_new_topo_update) THEN

      CALL bounds_exch(1,'p',DEUTO,'mpiom 10')


      DO j=1,je
        DO i=1,ie
          depto(i,j)=deuto(i,j)
        END DO
      END DO


      DO I=1,IE
        IF (have_g_js .AND. .NOT. lbounds_exch_tp ) THEN
          depto(i, 1) = 0._wp
          depto(i, 2) = 0._wp
        END IF
        IF (have_g_je) THEN
          depto(i, je) = 0._wp
          depto(i, je1) = 0._wp
        END IF
      END DO


      yloop: DO J=2,JE1
        xloop: DO I=2,IE1
          DEUTO(I,J)=DEPTO(I,J)
          IF (depto(i, j) .LT. 1._wp) CYCLE xloop
          IF (depto(i-1, j) .LT. 1._wp .AND. depto(i+1, j) .LT. 1._wp  &
               .AND. depto(i, j-1) .LT. 1._wp .AND. depto(i, j+1) .LT. 1._wp) &
               deuto(i, j) = 0._wp
        END DO xloop
      END DO yloop
      CALL bounds_exch(1,'u+',DEUTO,'mpiom 11')


      DO j=1,je
        DO i=1,ie
          depto(i,j) = deuto(i,j)
        END DO
      END DO

    END IF ! istart

    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlxp)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlxu)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlxv)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlyp)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlyu)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, dlyv)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, ftwou)
    IF (p_pe == p_io) READ(io_in_arcg) ext_hdr
    CALL read_slice(io_in_arcg, ftwov)
    IF (p_pe == p_io) CLOSE(io_in_arcg)


    dlxp(:,:) = MAX(1._wp, dlxp(:,:))
    dlxu(:,:) = MAX(1._wp, dlxu(:,:))
    dlxv(:,:) = MAX(1._wp, dlxv(:,:))
    dlyp(:,:) = MAX(1._wp, dlyp(:,:))
    dlyu(:,:) = MAX(1._wp, dlyu(:,:))
    dlyv(:,:) = MAX(1._wp, dlyv(:,:))

    DO j=1,je
      DO i=2,ie1
        dlxpsi(i, j) = 0.5_wp * (dlxv(i, j) + dlxv(i+1, j))
        dlypsi(i, j) = 0.5_wp * (dlyv(i, j) + dlyv(i+1, j))
      END DO
    END DO

    CALL bounds_exch(1,'v+',dlxv,'mpiom 12')
    CALL bounds_exch(1,'p',dlxp,'mpiom 13')
    CALL bounds_exch(1,'u+',dlxu,'mpiom 14')
    CALL bounds_exch(1,'p',dlyp,'mpiom 15')
    CALL bounds_exch(1,'u+',dlyu,'mpiom 16')
    CALL bounds_exch(1,'v+',dlyv,'mpiom 17')
    CALL bounds_exch(1,'s',dlypsi,'mpiom 17a')
    CALL bounds_exch(1,'s',dlxpsi,'mpiom 17b')

!SL
!> Grid deformation
!SL
    DO j=2,je1
      DO i=2,ie1
        curvav(i,j) = (dlxp(i,j+1)-dlxp(i,j)) / (dlyv(i,j)*dlxv(i,j))
      END DO
    END DO

    IF ( lbounds_exch_tp ) curvav(:,:) = 0._wp

    CALL bounds_exch(1,'v+',curvav,'mpiom 18')

    DO I=1,IE
      IF (have_g_je) THEN
        CURVAV(I,JE)=CURVAV(I,JE-1)
        DLYP(I,JE)=DLYP(I,JE-1)
        DLYV(I,JE)=DLYV(I,JE-1)
        DLXV(I,JE)=DLXV(I,JE-1)
        DLXP(I,JE)=DLXP(I,JE-1)
        DLXU(I,JE)=DLXU(I,JE-1)
        DLYU(I,JE)=DLYU(I,JE-1)
      END IF
    END DO

    CALL GATHER(DLXP,DLXP_G,p_io)
    CALL GATHER(DLYP,DLYP_G,p_io)
    CALL GATHER(DLXU,DLXU_G,p_io)
    CALL GATHER(DLYU,DLYU_G,p_io)
    CALL GATHER(DLXV,DLXV_G,p_io)
    CALL GATHER(DLYV,DLYV_G,p_io)

    area(:,:) = dlxp(:,:)*dlyp(:,:)
    areain(:,:) = 1.0_wp / area(:,:)
!
!> Find grid distance maximum and minimum
!
    dmino = MINVAL(dlxv)
    dmaxo = MAXVAL(dlxv)

    CALL global_min(dmino)
    CALL global_max(dmaxo)

#ifdef MESSY
    IF (p_pe==p_io) THEN 
       WRITE(io_stdout,'(a,f10.3,a,f10.3,a)')' RESOLUTION KM  MIN. : ', &
             dmino/1000._wp, '  MAX. : ', dmaxo/1000._wp, ' ODD '
       WRITE(io_stdout,*) (dlxp(5,j),j=1,je)
    ENDIF
#else
    WRITE(io_stdout,'(a,f10.3,a,f10.3,a)')' RESOLUTION KM  MIN. : ', &
          dmino/1000._wp, '  MAX. : ', dmaxo/1000._wp, ' ODD '
    WRITE(io_stdout,*) (dlxp(5,j),j=1,je)
#endif
!
!> Calculate inverse grid distances
!
    dlxpi = 1.0_wp / dlxp
    dlypi = 1.0_wp / dlyp
    dlxui = 1.0_wp / dlxu
    dlyui = 1.0_wp / dlyu
    dlxvi = 1.0_wp / dlxv
    dlyvi = 1.0_wp / dlyv

    !dtdxuo = dt * dlxui
    !dtdxpo = dt * dlxpi
    !dtdyo  = dt * dlypi
    !dpyo   = dt * dlypi
    !dtdxue = dt * dlxvi
    !dtdxpe = dt * dlxvi
    !dpye   = dt * dlyvi

    dtdxuo = dt / dlxu
    dtdxpo = dt / dlxp
    dtdyo  = dt / dlyp
    dpyo   = dt / dlyp
    dtdxue = dt / dlxv
    dtdxpe = dt / dlxv
    dpye   = dt / dlyv

  END SUBROUTINE setup_arcgri



  SUBROUTINE wrte_gridinfo(kky,kkm,kkd)
    !
    !c**** *wrte_gridinfo* - save grid information.
    !c
    !c     ch,    *mpi-met, hh*    10.04.01
    !
    !     modified
    !     --------
    !     s.legutke,        *mpi-mad, hh*    01.10.01
    !     - separate routine extracted from ollie (main)
    !     - netcdf version possible (with cond.comp. pnetcdfo)
    !
    !     purpose
    !     -------
    !
    !
    !     method
    !     -------
    !
    !
    !**   interface.
    !     ----------
    !
    !     *call*       *wrte_gridinfo(kky,kkm,kkd)*
    !
    !     *parameter*  *param1.h*     - grid size parameters for ocean model.
    !     *common*     *commo1.h*     - ocean/sediment tracer arrays.
    !     *common*     *units.h*      - std i/o logical units.
    !
    !**   interface to calling routine (parameter list):
    !     ----------------------------------------------
    !
    !     *integer* *kky*   - actual year.
    !     *integer* *kkm*   - actual month.
    !     *integer* *kkd*    - actual day.
    !
    !
    !     externals
    !     ---------
    !     none.
    !
    !**************************************************************************
    INTEGER(kind=i4) i4i1,i4i2,i4i3,i4i4
    INTEGER kky,kkm,kkd, k

    ! unless user requested gridinfo output...
    IF (.NOT. write_gridinfo_output) RETURN
    ! write to disk (extra format)
    !
    i4i1=((kky*10000)+(kkm*100)+kkd)
    i4i3=0
    i4i4=(ie_g*je_g)


    i4i2=190
    If (p_pe == p_io) Then
      Write(io_ou_alat)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alat,alat)

    i4i2=191
    If (p_pe == p_io) Then
      Write(io_ou_alon)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alon,alon)

    i4i2=154
    If (p_pe == p_io) Then
      Write(io_ou_alatu)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alatu,alatu)

    i4i2=155
    If (p_pe == p_io) Then
      Write(io_ou_alonu)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alonu,alonu)

    i4i2=56
    If (p_pe == p_io) Then
      Write(io_ou_alatv)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alatv,alatv)

    i4i2=57
    If (p_pe == p_io) Then
      Write(io_ou_alonv)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_alonv,alonv)

    i4i2=85
    If (p_pe == p_io) Then
      Write(io_ou_dlxp)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlxp,dlxp)

    i4i2=86
    If (p_pe == p_io) Then
      Write(io_ou_dlyp)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlyp,dlyp)

    i4i2=185
    If (p_pe == p_io) Then
      Write(io_ou_dlxu)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlxu,dlxu)

    i4i2=188
    If (p_pe == p_io) Then
      Write(io_ou_dlxv)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlxv,dlxv)

    i4i2=186
    If (p_pe == p_io) Then
      Write(io_ou_dlyu)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlyu,dlyu)

    i4i2=197
    If (p_pe == p_io) Then
      Write(io_ou_bek)i4i1,i4i2,i4i3,i4i4
    Endif

    i4i2=189
    If (p_pe == p_io) Then
      Write(io_ou_dlyv)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dlyv,dlyv)

    i4i2=193
    If (p_pe == p_io) Then
      Write(io_ou_deuto)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_deuto,deuto)

    i4i2=196
    If (p_pe == p_io) Then
      Write(io_ou_deute)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_deute,deute)

    i4i2=84
    If (p_pe == p_io) Then
      Write(io_ou_dept)i4i1,i4i2,i4i3,i4i4
    Endif
    Call write_slice_sp(io_ou_dept,depto)

    Do k=1,ke
      i4i2=172
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_weto)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_weto,weto(:,:,k))
    Enddo

    Do k=1,ke
      i4i2=194
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_amsue)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_amsue,amsue(:,:,k))
    Enddo

    Do k=1,ke
      i4i2=195
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_amsuo)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_amsuo,amsuo(:,:,k))
    Enddo

    Do k=1,ke
      i4i2=184
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_dduo)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_dduo,dduo(:,:,k))
    Enddo

    Do k=1,ke
      i4i2=187
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_ddue)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_ddue,ddue(:,:,k))
    Enddo

    Do k=1,ke
      i4i2=192
      i4i3=k
      If (p_pe == p_io) Then
        Write(io_ou_ddpo)i4i1,i4i2,i4i3,i4i4
      Endif
      Call write_slice_sp(io_ou_ddpo,ddpo(:,:,k))
    Enddo



    Return
  END SUBROUTINE wrte_gridinfo


  !
  !     subroutine to determine i and j indices of nearest point
  !     from lat and long, parallelized version over multiple tasks
  !
  ! RJ: i and j are GLOBAL indices!!!
  !
  !     input:
  !     alat      latitude (deg.)
  !     alon      longitude (deg.)
  !     k          level (integer)
  !     wetref     if land shall be considered too : 0.
  !                ocean points only                 1.
  !
  !     output:
  !     ipos       i index
  !     jpos       j index
  !     dist       distance from point in m
  !
  SUBROUTINE p_suchij_serial(alat,alon,k,ipos,jpos,dist,wetref)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: alat,alon,wetref
    REAL(wp), INTENT(out) :: dist
    INTEGER, INTENT(in) :: k
    INTEGER, INTENT(out) :: ipos,jpos
    REAL(wp), ALLOCATABLE :: w_g(:,:)
    LOGICAL :: lwetref
    lwetref = wetref .EQ. 1._wp

    if (p_pe == p_io) then
      allocate(w_g(ie_g,je_g))
    else
      allocate(w_g(0,0))

    endif

    call gather(weto(:,:,k),w_g,p_io)

    IF (p_pe == p_io) THEN

      CALL suchij_2d(ie_g, je_g, alat, alon, ipos, jpos, dist, lwetref, w_g, &
           alat_g, alon_g, lbounds_exch_tp)

    ENDIF

    CALL p_bcast(ipos,p_io)
    CALL p_bcast(jpos,p_io)

    DEALLOCATE(w_g)

    RETURN
  END SUBROUTINE p_suchij_serial

  !>
  !! Find level index of layer containing given depth.
  !!
  !! If depth is exactly at layer boundary, take layer above given depth.
  !! Return top layer for zero depth. Use bottom layer for too large depths.
  !!
  function get_level_index_by_depth(depth) result(level_index)

     real(dp), intent(in) :: depth
     integer :: level_index

     level_index = 1
     do while(tiestw(level_index+1) < depth .and. level_index < ke)
        level_index = level_index + 1
     end do

  end function get_level_index_by_depth

  !
  !     subroutine to determine i and j indices of nearest point
  !     from lat and long, parallelized version over multiple tasks
  !
  ! RJ: i and j are GLOBAL indices!!!
  !
  !     input:
  !     lat      latitude (deg.)
  !     lon      longitude (deg.)
  !     k          level (integer)
  !     wetref     if land shall be considered too : 0.
  !                ocean points only                 1.
  !
  !     output:
  !     ipos       i index
  !     jpos       j index
  !     dist       distance from point in m
  !
  SUBROUTINE p_suchij(lat,lon,k,ipos,jpos,dist,wetref,on_pe,llocal_idx)
    IMPLICIT NONE

    REAL(dp), INTENT(IN)  :: lat,lon,wetref
    INTEGER,  INTENT(IN)  :: k
    LOGICAL, INTENT(IN), OPTIONAL  :: llocal_idx
    REAL(dp), INTENT(OUT) :: dist
    INTEGER, INTENT(OUT)  :: ipos, jpos
    INTEGER, INTENT(OUT), OPTIONAL :: on_pe
    LOGICAL :: lwetref
    TYPE(grid_dist_2d) :: my_min, global_min
    INTEGER :: ierr
    LOGICAL :: return_local_index

    lwetref = wetref .EQ. 1._wp

    CALL suchij_2d(ie, je, lat, lon, my_min%i, my_min%j, my_min%dist, lwetref, weto(:,:,k), &
         alat, alon, lbounds_exch_tp .AND. p_joff == 0)

    my_min%pe = p_pe
    my_min%i = my_min%i + p_ioff
    my_min%j = my_min%j + p_joff

    CALL MPI_ALLREDUCE(my_min,global_min,1,p_grid_dist_2d,p_op_grid_dist_min_2d, &
                       p_all_comm, ierr)

    IF (ierr /= MPI_SUCCESS) THEN
      CALL stop_all ('MPI_ALLREDUCE in p_suchij failed.')
    ENDIF

    IF (PRESENT(llocal_idx)) THEN
      return_local_index = llocal_idx
    ELSE
      return_local_index = .FALSE.
    ENDIF

    IF (return_local_index) THEN
      ipos = MERGE(global_min%i-p_ioff, -1, p_pe == global_min%pe)
      jpos = MERGE(global_min%j-p_joff, -1, p_pe == global_min%pe)
    ELSE
      ipos = global_min%i
      jpos = global_min%j
    END IF

    dist = global_min%dist
    IF (PRESENT(on_pe)) THEN
      on_pe=global_min%pe
    ENDIF
  END SUBROUTINE p_suchij

  SUBROUTINE setup_grid

    CALL setup_arcgri
    CALL generate_MPI_struct_type(2, (/1, 3/), &
         (/MPI_DOUBLE_PRECISION, MPI_INTEGER/), &
         p_grid_dist_2d)
    CALL create_MPI_op(op_grid_dist_min_2d, .TRUE., p_op_grid_dist_min_2d)


  END SUBROUTINE setup_grid


  SUBROUTINE op_grid_dist_min_2d(a, b, n, type)

    INTEGER, INTENT(IN) :: n, type
    TYPE(grid_dist_2d), INTENT(IN)   :: a(n)
    TYPE(grid_dist_2d), INTENT(INOUT)  :: b(n)
    INTEGER :: i

    DO i=1,n
      IF (b(i)%dist < a(i)%dist) THEN
         CYCLE
      ELSEIF (b(i)%dist == a(i)%dist) THEN
         IF (b(i)%j > a(i)%j) THEN
           b(i) = a(i)
         ELSEIF (b(i)%j == a(i)%j .AND. b(i)%i > a(i)%i ) THEN
           b(i) = a(i)
         ENDIF
      ELSE
         b(i) = a(i)
      ENDIF
    ENDDO

  END SUBROUTINE op_grid_dist_min_2d

  SUBROUTINE ini_cell_thickness

    ALLOCATE(thkcello(ie,je,ke))

    ! initialise the cell thickness with the layer thickness
    thkcello(:,:,:)=ddpo(:,:,:)


  END SUBROUTINE ini_cell_thickness

  SUBROUTINE cell_thickness

    USE mo_commo1, ONLY : ddpo,zo,sicsno,sictho
    USE mo_planetary_constants, ONLY : rhosnwa,rhoicwa

    thkcello(:,:,1) = ddpo(:,:,1) + zo(:,:) &
         - (rhosnwa * sicsno(:,:) + rhoicwa * sictho(:,:))

  END SUBROUTINE cell_thickness

END MODULE mo_grid
