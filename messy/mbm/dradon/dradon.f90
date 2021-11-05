PROGRAM DRADON

  USE messy_dradon_box, ONLY: initialize, integrate, finalize

  CALL initialize

  WRITE(*,*) 'START INTEGRATION ...'
  CALL integrate
  WRITE(*,*) '... INTEGRATION FINISHED!'

  CALL finalize

END PROGRAM DRADON
