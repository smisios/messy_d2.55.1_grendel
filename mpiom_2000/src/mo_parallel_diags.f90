MODULE mo_parallel_diags
#ifndef MESSY
  USE mo_kind, ONLY: dp
  USE mo_mpi, ONLY: p_tagged_minmeanmax, min_mean_max_dp, p_pe, p_io
  USE mo_units, ONLY: io_stdout
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timing_diagnostics
CONTAINS
  SUBROUTINE timing_diagnostics(time_diffs, time_aggs, region_tags, &
       per_pe_diag_arg, ldt)
    REAL(dp), INTENT(in) :: time_diffs(:)
    TYPE(min_mean_max_dp) :: time_aggs(:)
    CHARACTER(len=*), INTENT(in) :: region_tags(:)
    LOGICAL, OPTIONAL, INTENT(in) :: per_pe_diag_arg
    INTEGER, INTENT(in) :: ldt

    LOGICAL :: per_pe_diag
    INTEGER :: i, n

    IF (PRESENT(per_pe_diag_arg)) THEN
      per_pe_diag  = per_pe_diag_arg
    ELSE
      per_pe_diag = .FALSE.
    END IF

    n = SIZE(time_diffs)

    DO i = 1, n
      time_aggs(i)%min = time_diffs(i)
      time_aggs(i)%pe_min = p_pe
      time_aggs(i)%max = time_diffs(i)
      time_aggs(i)%pe_max = p_pe
      time_aggs(i)%mean = time_diffs(i)
      time_aggs(i)%rcount = 1
    END DO
    CALL p_tagged_minmeanmax(time_aggs, p_io)
    IF (p_pe == p_io) THEN
      DO i = 1, n
        WRITE(io_stdout,'(a,a,a,i7,a,f8.3,a,i4,a,f8.3,a,f8.3,a,i4,a)') &
             'Time for ', TRIM(region_tags(i)), ' timestep: ', ldt, &
             '; min:', time_aggs(i)%min, 's (on rank ', &
             time_aggs(i)%pe_min, '), mean:', &
             time_aggs(i)%mean, &
             's, max:', time_aggs(i)%max, &
             's (on rank ', time_aggs(i)%pe_max, ')'
      END DO
    END IF
    IF (per_pe_diag) THEN
      DO i = 1, n
        WRITE(io_stdout,'(i4,a,a,a,i7,a,f8.3,a)') &
             p_pe, ': Time for ', region_tags(i), ', timestep', ldt, &
             ':', time_diffs(i), ' s'
      END DO
    END IF
  END SUBROUTINE timing_diagnostics
#endif
END MODULE mo_parallel_diags
