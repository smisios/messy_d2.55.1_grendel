MODULE MESSY_GMXE_OC_AGING

  USE MESSY_GMXE_MEM

  IMPLICIT NONE


  INTEGER,PARAMETER :: num_wsoc   = 5
  INTEGER,PARAMETER :: decay_rate = 2

  INTEGER           :: spec_idx_wsoc(0:num_wsoc)
  REAL(dp)          :: oc_ratio(0:num_wsoc)
  REAL(dp)          :: kappa(0:num_wsoc)
  REAL(dp)          :: exposure_time(0:num_wsoc)

CONTAINS 

!=============================================================================
SUBROUTINE INIT_OC_AGING

  INTEGER  :: jt, ji
  INTEGER, PARAMETER :: looklimit = 500
  REAL(dp) :: oc_rat_func(looklimit)
  REAL(dp) :: val1

  DO ji=1,looklimit
    val1 = REAL(ji,dp)
    oc_rat_func(ji) = LOG10(1.34_dp * val1) * 0.122_dp + &
                      0.0002_dp * val1 + 0.3_dp
  ENDDO

  oc_ratio(0)      = 0.3_dp
  kappa(0)         = 0.01_dp
  exposure_time(0) = 0._dp
  DO jt=1,num_wsoc
     oc_ratio(jt) = 0.3_dp + REAL(jt,dp)/REAL(num_wsoc,dp) * 0.4_dp
     kappa(jt)    = 0.01_dp + (oc_ratio(jt) - 0.3_dp) * 0.5_dp
!     kappa(jt)    = 0.06_dp + (oc_ratio(jt) - 0.3_dp) * 0.5_dp
     DO ji=1,looklimit
        IF ( oc_rat_func(ji) > oc_ratio(jt) ) THEN
           exposure_time(jt) = REAL(ji,dp) * 1.e9_dp
           EXIT
        END IF
     END DO
  END DO
  

END SUBROUTINE INIT_OC_AGING
!=============================================================================
SUBROUTINE gmxe_oc_aging(kproma, nlev, naertot, nmod, nsoluble, ztmst, &
                         paerml, kappa_oc)

  INTEGER,  INTENT(IN)    :: kproma, nlev, naertot, nmod, nsoluble
  REAL(dp), INTENT(IN)    :: ztmst
  REAL(dp), INTENT(INOUT) :: paerml(kproma,nlev,0:naertot)
  REAL(dp), INTENT(INOUT) :: kappa_oc(kproma,nlev,1:nmod)

  INTEGER :: jt, jm, jl, jk, jc, idx1, idx2, km
  REAL(dp), DIMENSION(kproma,nlev) :: reaction_rate
  REAL(dp), DIMENSION(kproma,nlev) :: sum_tot_oc
  REAL(dp)                         :: val1, val2, val0, new_val, ohval

  DO jk = 1,nlev
    DO jl=1,kproma
      SELECT CASE (decay_rate)
      CASE (1) 
        ! constant loss rate 
        ! equivalent OH 1e5 molec/cm^3
        ohval = 1.e5_dp 
      CASE (2)
        ! using real OH value
        ohval = paerml(jl,jk,species(spec_idx_oh)%aermlidx(0))
      END SELECT
      
      reaction_rate(jl,jk) = ohval * ztmst
      
    ENDDO
  ENDDO
!  DO jm=1,nmod    
  DO jm=1,nsoluble
    km = jm
    IF (jm > nsoluble) km = jm - nsoluble + 1
    DO jt=num_wsoc-1,0,-1
      idx1 = species(spec_idx_wsoc(jt))%aermlidx(jm)
      IF (idx1 == 0) CYCLE
      DO jc=num_wsoc,jt+1,-1
        idx2 = species(spec_idx_wsoc(jc))%aermlidx(km)
        IF (idx2 == 0) CYCLE
        DO jk=1,nlev
          DO jl=1,kproma
            val0 = exposure_time(jt) + reaction_rate(jl,jk)
            val1 = val0 - exposure_time(jc-1)
            IF (val1 > 0._dp) THEN
              val2 = MIN (0.9999_dp, reaction_rate(jl,jk) &
                       / ( exposure_time(jc) - exposure_time(jc-1) ) )
              new_val = paerml(jl,jk,idx1) * (1._dp - val2)
!!$              print*, "val2", val2, jl, jk, jc, jt, val1, val0, &
!!$                "before: ",paerml(jl,jk,idx1), paerml(jl,jk,idx2), &
!!$                paerml(jl,jk,idx1)+ paerml(jl,jk,idx2), &
!!$                "after: ",new_val+paerml(jl,jk,idx2)+paerml(jl,jk,idx1)*val2, &
!!$                new_val, paerml(jl,jk,idx2) + paerml(jl,jk,idx1) * val2
              paerml(jl,jk,idx2) = paerml(jl,jk,idx2) &
                                 + paerml(jl,jk,idx1) * val2
              paerml(jl,jk,idx1) = new_val
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  kappa_oc(:,:,:) = 0._dp
  DO jm=1,nmod
    sum_tot_oc(:,:) = 0._dp
    DO jt=0,num_wsoc
      idx1 = species(spec_idx_wsoc(jt))%aermlidx(jm)
      IF (idx1 == 0) CYCLE
      DO jk=1,nlev
        DO jl=1,kproma
          sum_tot_oc(jl,jk) = sum_tot_oc(jl,jk) + paerml(jl,jk,idx1)
        ENDDO
      ENDDO
    END DO
    DO jt=0,num_wsoc
      idx1 = species(spec_idx_wsoc(jt))%aermlidx(jm)
      IF (idx1 == 0) CYCLE
      DO jk=1,nlev
        DO jl=1,kproma
          IF (sum_tot_oc(jl,jk) > 0._dp) &
            kappa_oc(jl,jk,jm) = kappa_oc(jl,jk,jm) + &
                                 kappa(jt)          * &
                                 paerml(jl,jk,idx1) / sum_tot_oc(jl,jk)
        END DO
      END DO
    END DO
  ENDDO


END SUBROUTINE gmxe_oc_aging
!============================================================================

END MODULE MESSY_GMXE_OC_AGING
