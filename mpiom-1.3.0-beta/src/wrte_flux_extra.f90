      SUBROUTINE WRTE_FLUX_EXTRA
#ifdef __coupled
!C****************************************************************
!C
!C**** *WRTE_FLUX_EXTRA* - save atmosph. fluxes at coupled timestep
!
!C     JJ,    *MPI-Met, HH*    03.08.2003
!C     SL,    *MPI-Met, M&D*   08.08.2003
!C      - field name correction
!C     NSK,   *IFM-GEOMAR*     28.09.2004  
!C      - MPI version
!C
!C     Modified
!C     --------
!C
!C     Purpose
!C     -------
!C     sace atmosph flux fileds at each coupled timestep
!C
!C     Method
!C     -------
!C**   Interface.
!C     ----------
!C!
!C     *CALL*       *WRTE_MEAN(kdtday,kdays,kmonts,kmean)*
!C
!C     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!C     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!C     *COMMON*     *UNITS.h*      - std I/O logical units.
!C
!C     Externals
!C     ---------
!C     none.
!C
!C**************************************************************************
      USE MO_PARAM1
      USE MO_COMMO1

      USE MO_UNITS
      USE MO_FLUXES1
      USE MO_PARALLEL
      USE MO_KIND

      INTEGER(KIND=i4) I4I1,I4I2,I4I3,I4I4

!      REAL(wp), POINTER ::         glfld(:,:)
      REAL(wp), ALLOCATABLE ::         glfld(:,:)

      IDATE= (LYEARS*10000)+(LMONTS*100)+LDAYS
      i4i1=idate
      i4i2=270
      i4i3=0
      i4i4=ie_g*je_g
!     WRITE(IO_STDOUT,*) 'in sbr flux_extra' 
!     WRITE(IO_STDOUT,*) 'i4i1= ',i4i1
!     WRITE(IO_STDOUT,*) 'i4i2= ',i4i2
!     WRITE(IO_STDOUT,*) 'i4i3= ',i4i3
!     WRITE(IO_STDOUT,*) 'i4i4= ',i4i4
!

      IF ( p_parallel_io ) THEN
         ALLOCATE (glfld(ie_g,je_g))
!!$      ELSE
!!$         glfld => NULL()
      ENDIF

      call gather_arr(aoflnhwo,glfld,p_io)
      IF (p_parallel_io) THEN
         OPEN(IO_OU_AONHW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AONHW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AONHW) real(glfld,sp)
      CLOSE(IO_OU_AONHW)
      ENDIF
!
      call gather_arr(aoflshwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=271
         OPEN(IO_OU_AOSHW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOSHW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOSHW) real(glfld,sp)
      CLOSE(IO_OU_AOSHW)
      ENDIF
!
      call gather_arr(aoflrhio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=272
         OPEN(IO_OU_AORHI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AORHI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AORHI) real(glfld,sp)
      CLOSE(IO_OU_AORHI)
      ENDIF
!
      call gather_arr(aoflchio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=273
         OPEN(IO_OU_AOCHI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOCHI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOCHI) real(glfld,sp)
      CLOSE(IO_OU_AOCHI)
      ENDIF
!
      call gather_arr(aoflfrwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=274
         OPEN(IO_OU_AOFRW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRW) real(glfld,sp)
      CLOSE(IO_OU_AOFRW)
      ENDIF
!
      call gather_arr(aoflfrio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=275
         OPEN(IO_OU_AOFRI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOFRI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOFRI) real(glfld,sp)
      CLOSE(IO_OU_AOFRI)
      ENDIF
!
      call gather_arr(aofltxwo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=276
         OPEN(IO_OU_AOTXW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTXW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTXW) real(glfld,sp)
      CLOSE(IO_OU_AOTXW)
      ENDIF
!
      call gather_arr(aofltywe,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=277
         OPEN(IO_OU_AOTYW, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTYW) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTYW) real(glfld,sp)
      CLOSE(IO_OU_AOTYW)
      ENDIF
!
      call gather_arr(aofltxio,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=278
         OPEN(IO_OU_AOTXI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTXI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTXI) real(glfld,sp)
      CLOSE(IO_OU_AOTXI)
      ENDIF
!
      call gather_arr(aofltyie,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=279
         OPEN(IO_OU_AOTYI, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOTYI) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOTYI) real(glfld,sp)
      CLOSE(IO_OU_AOTYI)
      ENDIF
!
      call gather_arr(aoflwsvo,glfld,p_io)
      IF (p_parallel_io) THEN
      i4i2=280
         OPEN(IO_OU_AOWSV, STATUS='UNKNOWN',                            &
     &       ACCESS='SEQUENTIAL',                                       &
     &       POSITION='APPEND',                                         &
     &       FORM='UNFORMATTED')
      WRITE(IO_OU_AOWSV) i4i1,i4i2,i4i3,i4i4
      WRITE(IO_OU_AOWSV) real(glfld,sp)
      CLOSE(IO_OU_AOWSV)
      ENDIF
!      

      IF (p_parallel_io) THEN      
         DEALLOCATE(glfld)
      ENDIF
!
      RETURN
#endif
      END
