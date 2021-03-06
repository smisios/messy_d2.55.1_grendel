MODULE messy_main_import_bi

! *********************************************************
! Author: Astrid Kerkweg, UNI-MZ, Okt 2009
! *********************************************************

 USE messy_main_import_grid, ONLY: import_grid_submodstr => submodstr

 PUBLIC :: import_grid_submodstr

 PUBLIC :: main_import_initialize
 PUBLIC :: main_import_init_memory
 PUBLIC :: main_import_init_tracer
 PUBLIC :: main_import_global_start
 PUBLIC :: main_import_free_memory

! IMPLICIT NONE

 CONTAINS

   SUBROUTINE main_import_initialize

     USE messy_main_import_grid_bi, ONLY: import_grid_initialize
     USE messy_main_import_ts_bi,   ONLY: import_ts_initialize
     USE messy_main_import_lt_bi,   ONLY: import_lt_initialize

     IMPLICIT NONE

     CALL import_grid_initialize
     CALL import_ts_initialize
     CALL import_lt_initialize

   END SUBROUTINE main_import_initialize

! ----------------------------------------------

   SUBROUTINE main_import_init_memory

     USE messy_main_import_grid_bi, ONLY: import_grid_init_memory
     USE messy_main_import_ts_bi,   ONLY: import_ts_init_memory

     IMPLICIT NONE

     CALL import_grid_init_memory
     CALL import_ts_init_memory
   END SUBROUTINE main_import_init_memory

! ----------------------------------------------

   SUBROUTINE main_import_init_tracer

     USE messy_main_import_grid_bi, ONLY: import_grid_init_tracer

     IMPLICIT NONE

     CALL import_grid_init_tracer

   END SUBROUTINE main_import_init_tracer

! ----------------------------------------------

   SUBROUTINE main_import_global_start

     USE messy_main_import_grid_bi, ONLY: import_grid_global_start
     USE messy_main_import_lt_bi,   ONLY: import_lt_global_start
     USE messy_main_import_ts_bi,   ONLY: import_ts_global_start

     IMPLICIT NONE

     CALL import_grid_global_start
     CALL import_lt_global_start
     CALL import_ts_global_start

   END SUBROUTINE main_import_global_start

! ----------------------------------------------

   SUBROUTINE main_import_free_memory

     USE messy_main_import_grid_bi, ONLY: import_grid_free_memory
     USE messy_main_import_ts_bi,   ONLY: import_ts_free_memory
     USE messy_main_import_lt_bi,   ONLY: import_lt_free_memory

     IMPLICIT NONE

     CALL import_grid_free_memory
     CALL import_ts_free_memory
     CALL import_lt_free_memory

   END SUBROUTINE main_import_free_memory

! ----------------------------------------------

END MODULE messy_main_import_bi
