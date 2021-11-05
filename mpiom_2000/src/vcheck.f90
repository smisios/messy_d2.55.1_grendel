SUBROUTINE vcheck(l,feld)

  USE mo_mpi
  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_units

  IMPLICIT NONE

  REAL(wp), POINTER :: helfie_g(:,:),vq(:)
  REAL(wp) :: feld(:,:,:)
  INTEGER i,k,ir,l,id1(1),id2(1)

  IF (p_pe==p_io) THEN
     ALLOCATE(helfie_g(ie_g,je_g),vq(ie_g))
  ELSE
     helfie_g => NULL()
     vq       => NULL()
  ENDIF

  IF(l.LE.50)RETURN

  CALL gather(feld(:,:,1),helfie_g,p_io)

  IF (p_pe==p_io) THEN
     vq = 0._wp
     DO i=2,ie_g/2 - 1
        ir=ie_g+1-i
         vq(i)=helfie_g(i,2)+helfie_g(ir,2)
      ENDDO
      !    write(0,*)'linie',(i,vq(i),i=1,ie_g/2)
      id1=MINLOC(vq)
      id2=MAXLOC(vq)
      i=id1(1)
      k=id2(1)
      WRITE(0,*)'vcheck',l,i,k,MINVAL(vq),MAXVAL(vq),helfie_g(i,2)&
           &, helfie_g(k,2)
   ENDIF

   DEALLOCATE(helfie_g,vq)

   RETURN
 END SUBROUTINE vcheck
