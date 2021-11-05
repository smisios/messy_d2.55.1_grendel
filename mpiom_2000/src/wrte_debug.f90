      SUBROUTINE wrte_debug(III)

      USE MO_PARAM1
      USE MO_MPI
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_COMMOAU2
      USE mo_restart, ONLY: dump_restart_data_extra
      USE MO_UNITS
#ifndef MESSY                                                                                                                           
  USE MO_ISO_C_KINDS                                                                                                                    
#else                                                                                                                                   
  USE mo_kind, ONLY : c_int64_t => i8, c_float => sp                                                                                    
#endif  

      IMPLICIT NONE

      INTEGER(KIND=C_INT64_T) IDATE
      INTEGER :: III, iunit
!
      IUNIT=777
      IF(p_pe==p_io) THEN
        OPEN(IUNIT,FILE='D37000',FORM='UNFORMATTED')
      ENDIF

      IDATE=INT(III, C_INT64_T)


      CALL dump_restart_data_extra(iunit, idate)
     END SUBROUTINE wrte_debug
