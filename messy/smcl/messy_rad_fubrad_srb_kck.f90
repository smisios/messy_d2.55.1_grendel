! *********************************************************************************
MODULE  messy_rad_fubrad_srb_kck
  !
  ! PURPOSE:
  !    Provide a parameterization of the short wave heating rates due to 
  !    absorption of UV radiation in the Schumann-Runge bands (175 - 205 nm)
  !    by oxygen.
  !
  ! REFERENCE: 
  !    G.Kockarts, Penetration of solar radiation in the Schumann-Runge bands of 
  !    molecular oxygen:  A robust approximation,
  !    Annales Geophysicae, 12, 12, 1207-2017, doi:10.1007/BF03191317, 1994.
  !
  ! CODE HISTORY:
  !    Original code is from the Model for OZone and Related chemical Tracers 
  !    (MOZART), provided by NCAR: MOZART-4 module mo_schu, subroutine schu, 
  !    original code to calculate dto2, and xscho2
  !
  !    Markus Kunze, FU-Berlin, 09/2015: converted to f90 and implemented in FUBRAD
  !
  !********************************************************************************
  
  USE messy_main_constants_mem, ONLY: dp
  USE messy_rad_fubrad_mem
  
  IMPLICIT NONE
  PRIVATE
  SAVE
  
  PUBLIC :: fubrad_srb_kck_schu             &
          , fubrad_srb_kck_schu_init        &
          , fubrad_srb_kck_mem_ini          &
          , fubrad_srb_kck_mem_clean
  
  !
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: d_table
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: x_table
  REAL(dp), DIMENSION(:),   ALLOCATABLE :: o2_table
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: ac, bc

  ! Interface of public subroutines
  
  INTERFACE fubrad_srb_kck_schu
     MODULE PROCEDURE fubrad_srb_kck_schu
  END INTERFACE
  
  INTERFACE fubrad_srb_kck_schu_init
     MODULE PROCEDURE fubrad_srb_kck_schu_init
  END INTERFACE
  
  INTERFACE fubrad_srb_kck_mem_ini
     MODULE PROCEDURE fubrad_srb_kck_mem_ini
  END INTERFACE
  
CONTAINS
  !
  !=============================================================================*
  !
  SUBROUTINE fubrad_srb_kck_schu( nbdim, nproma, nswlev, nswlevp1, nsrb, flux_srb)
    !-----------------------------------------------------------------------------
    !   Purpose:
    ! ----------
    !   Calculate the equivalent absorption cross section of o2 in the sr bands.
    !
    !   Reference:
    ! ------------
    !   The algorithm is based on:
    !   G.Kockarts, Penetration of solar radiation in the Schumann-Runge bands of 
    !               molecular oxygen:  A robust approximation,
    !               Annales Geophysicae, v12, n12, pp. 1207ff, dec 1994.
    !
    !   Method:
    ! ---------
    !   Calculation is done on the wavelength grid used by Kockarts (1994).
    !   final values do include effects from the herzberg continuum.
    !
    !   History:
    ! ----------
    !   The original code is taken from the Model for OZone and Related 
    !   chemical Tracers (MOZART), provided by NCAR.
    !   MOZART-4 module mo_schu, subroutine schu, original code to calculate 
    !                            dto2, and xscho2
    !   Markus Kunze,  implemented in FUBRAD, September 2015.
    !-----------------------------------------------------------------------------
    !
    ! INPUT
    !-----------------------------------------------------------------------------
    INTEGER, INTENT(in) :: nbdim
    INTEGER, INTENT(in) :: nproma
    INTEGER, INTENT(in) :: nswlev
    INTEGER, INTENT(in) :: nswlevp1
    INTEGER, INTENT(in) :: nsrb
    ! OUTPUT:
    ! -------
    ! flux_srb  - REAL, flux in SR bands at each specified vertical layer 
    REAL(dp), DIMENSION(nbdim,nswlevp1), INTENT(out) :: flux_srb
    !
    !  ... local variables
    !-----------------------------------------------------------------------------
    !   dto2    - real, optical depth due to O2 absorption at each specified 
    !             vertical layer at each specified wavelength
    !   xscho2  - real, molecular absorption cross section in sr bands at  
    !             each specified wavelength.  includes herzberg continuum
    !   o2col   - real, slant overhead O2 column (molec/cc) at each specified
    !             altitude
    !   secchi  - ratio of slant to vertical O2 columns
    !
    INTEGER                                  :: j, k, iw, ki
    REAL(dp)                                 :: lo2col
    INTEGER,  DIMENSION(nbdim,nswlevp1)      :: index
    REAL(dp), DIMENSION(nbdim,nswlevp1,nsrb) :: dto2
    REAL(dp), DIMENSION(nbdim,nswlevp1,nsrb) :: xscho2
    REAL(dp), DIMENSION(nbdim,nswlevp1)      :: o2col
    REAL(dp), DIMENSION(nbdim,nswlevp1)      :: secchi
    REAL(dp), DIMENSION(nbdim,nswlevp1)      :: rjm, rjo2
    REAL(dp), DIMENSION(nbdim,nswlevp1)      :: dels
    REAL(dp), DIMENSION(nswlevp1)            :: trans

    !-----------------------------------------------------------------------------
    !  ... initialize 
    !-----------------------------------------------------------------------------
    flux_srb(:,1:nswlevp1) = 0._dp
    rjm (:,1:nswlev)   = 0._dp
    rjo2(:,1:nswlevp1) = 0._dp

    !-----------------------------------------------------------------------------
    ! ... initialize the table interpolation
    !-- ---------------------------------------------------------------------------
    DO j = 1, nproma
       DO ki = 1, nswlevp1
          o2col (j,ki) = po2c (j,ki) * przsec(j)
          IF ( o2col(j,ki) /= 0._dp ) THEN
             secchi(j,ki) = o2col(j,ki) / po2c(j,ki)
             lo2col = LOG10( o2col(j,ki) )
             IF ( lo2col <= o2_table(1) ) THEN
                dels(j,ki)  = 0._dp
                index(j,ki) = 1
             ELSE IF ( lo2col >= o2_table(tdim) ) THEN
                dels(j,ki)  = 1._dp
                index(j,ki) = tdim-1
             ELSE
                DO k = 2, tdim
                   IF ( lo2col <= o2_table(k) ) THEN
                      dels(j,ki)  = t_fac*(lo2col - o2_table(k-1))
                      index(j,ki) = k-1
                      EXIT
                   END IF
                END DO
             END IF
          ELSE
             index(j,ki) = 0
             dels(j,ki)  = 0._dp
          END IF
       END DO
       ! at the uppermost level po2c is zero: take the same value as for k=2
       secchi(j,1) = secchi(j,2)
       !-----------------------------------------------------------------------------
       ! ... calculate sum of exponentials (eqs 7 and 8 of kockarts 1994)
       !-----------------------------------------------------------------------------
       DO iw = 1, nsrb
          DO k = 1, nswlevp1
             ki = index(j,k)
             rjm(j,k)  = x_table(ki,iw) + dels(j,k)*(x_table(ki+1,iw) - x_table(ki,iw))
             rjo2(j,k) = d_table(ki,iw) + dels(j,k)*(d_table(ki+1,iw) - d_table(ki,iw))
          END DO
          trans(1) = 1._dp
          DO k = 2, nswlevp1
             IF ( rjm(j,k) > 1.E-100_dp ) THEN
                IF ( rjm(j,k-1) > 0._dp ) THEN
                   dto2(j,k,iw) = LOG( rjm(j,k-1) ) / secchi(j,k-1) - &
                                  LOG( rjm(j,k) )   * secchi(j,k)
                ELSE
                   dto2(j,k,iw) = 1000._dp
                END IF
             ELSE
                dto2(j,k,iw)   = 1000._dp
             END IF
             trans(k) = trans(k-1) * EXP( -dto2(j,k,iw) )
             IF (ldb_first) PRINT *,'rjm(j,',k,') = ',rjm(j,k),' rjo2(j,',k,') = ',rjo2(j,k), &
                                    ' secchi(j,',k-1,') = ',secchi(j,k-1), &
                                    ' secchi(j,',k,') = ',secchi(j,k)
             IF (ldb_first) PRINT *,'trans(',k-1,')=',trans(k-1), &
                                    ' EXP( -dto2(',j,',',k,',',iw,') ) = ',EXP( -dto2(j,k,iw) ), &
                                    ' dto2(',j,',',k,',',iw,') = ',dto2(j,k,iw)
          END DO
          DO k = 1, nswlevp1
             IF ( rjm(j,k) > 1.E-100_dp ) THEN
                IF ( rjo2(j,k) > 1.E-100_dp ) THEN
                   xscho2(j,k,iw) = rjo2(j,k)/rjm(j,k)
                ELSE
                   xscho2(j,k,iw) = 0._dp
                END IF
             ELSE
                xscho2(j,k,iw) = 0._dp
             END IF
             IF (ldb_first) PRINT *,' xscho2(',j,',',k,',',iw,') = ',xscho2(j,k,iw), &
                                    ' o2col(',j,',',k,') = ',o2col(j,k), &
                                    ' EXP(-xscho2(j,k,iw)*o2col(j,k)) = ',EXP(-xscho2(j,k,iw)*o2col(j,k))
             IF (ldb_first) PRINT *,' trans(',k,')=',trans(k)
             IF (ldb_first) PRINT *,' flux_srb_diff1(j,',k,')=',zFsrb(iw) * EXP(-xscho2(j,k,iw)*o2col(j,k))
             IF (ldb_first) PRINT *,' flux_srb_diff2(j,',k,')=',zFsrb(iw) * trans(k)
             flux_srb(j,k) = flux_srb(j,k) + zFsrb(iw) * trans(k)
          END DO
       END DO
    END DO
    ldb_first = .FALSE.
    RETURN
  END SUBROUTINE fubrad_srb_kck_schu

  !
  !-----------------------------------------------------------------------------   
  SUBROUTINE fubrad_srb_kck_schu_init(status)
    !
    INTEGER, INTENT(out) :: status
    !
    ! local variables
    INTEGER                :: iw, k
    REAL(dp)               :: col
    REAL(dp), DIMENSION(6) :: a0, a1, b0, b1
 
    status = -1
    !-----------------------------------------------------------------------------
    !  a(16,12)             coefficients for rj(m) (table 1 in kockarts 1994)
    !  b(16,12)                              rj(o2)(table 2 in kockarts 1994)
    !  rjm                  attenuation coefficients rj(m)
    !  rjo2                 rj(o2)
    !-----------------------------------------------------------------------------
    ac(:,:) = RESHAPE( (/  &
    !a 57000-56500.5 cm-1  (175.44 - 176.99 nm)
       1.13402e-01,1.00088e-20,3.48747e-01,2.76282e-20,3.47322e-01,1.01267e-19, &
       1.67351e-01,5.63588e-19,2.31433e-02,1.68267e-18,0.00000e+00,0.00000e+00, &
    !a 56500-56000.5 cm-1  (176.99 - 178.57 nm)
       2.55268e-03,1.64489e-21,1.85483e-01,2.03591e-21,2.60603e-01,4.62276e-21, &
       2.50337e-01,1.45106e-20,1.92340e-01,7.57381e-20,1.06363e-01,7.89634e-19, &
    !a 56000-55500.5 cm-1  (178.57 - 180.18 nm)
       4.21594e-03,8.46639e-22,8.91886e-02,1.12935e-21,2.21334e-01,1.67868e-21, &
       2.84446e-01,3.94782e-21,2.33442e-01,1.91554e-20,1.63433e-01,2.25346e-19, &
    !a 55500-55000.5 cm-1  (180.18 - 181.82 nm)
       3.93529e-03,6.79660e-22,4.46906e-02,9.00358e-22,1.33060e-01,1.55952e-21, &
       3.25506e-01,3.43763e-21,2.79405e-01,1.62086e-20,2.10316e-01,1.53883e-19, &
    !a 55000-54500.5 cm-1  (181.82 - 183.48 nm)
       2.60939e-03,2.33791e-22,2.08101e-02,3.21734e-22,1.67186e-01,5.77191e-22, &
       2.80694e-01,1.33362e-21,3.26867e-01,6.10533e-21,1.96539e-01,7.83142e-20, &
    !a 54500-54000.5 cm-1  (183.48 - 185.18 nm)
       9.33711e-03,1.32897e-22,3.63980e-02,1.78786e-22,1.46182e-01,3.38285e-22, &
       3.81762e-01,8.93773e-22,2.58549e-01,4.28115e-21,1.64773e-01,4.67537e-20, &
    !a 54000-53500.5 cm-1  (185.18 - 186.91 nm)
       9.51799e-03,1.00252e-22,3.26320e-02,1.33766e-22,1.45962e-01,2.64831e-22, &
       4.49823e-01,6.42879e-22,2.14207e-01,3.19594e-21,1.45616e-01,2.77182e-20, &
    !a 53500-53000.5 cm-1  (186.91 - 188.68 nm)
       7.87331e-03,3.38291e-23,6.91451e-02,4.77708e-23,1.29786e-01,8.30805e-23, &
       3.05103e-01,2.36167e-22,3.35007e-01,8.59109e-22,1.49766e-01,9.63516e-21, &
    !a 53000-52500.5 cm-1  (188.68 - 190.47 nm)
       6.92175e-02,1.56323e-23,1.44403e-01,3.03795e-23,2.94489e-01,1.13219e-22, &
       3.34773e-01,3.48121e-22,9.73632e-02,2.10693e-21,5.94308e-02,1.26195e-20, &
    !a 52500-52000.5 cm-1  (190.47 - 192.31 nm)
       1.47873e-01,8.62033e-24,3.15881e-01,3.51859e-23,4.08077e-01,1.90524e-22, &
       8.08029e-02,9.93062e-22,3.90399e-02,6.38738e-21,8.13330e-03,9.93644e-22, &
    !a 52000-51500.5 cm-1  (192.31 - 194.17 nm)
       1.50269e-01,1.02621e-23,2.39823e-01,3.48120e-23,3.56408e-01,1.69494e-22, &
       1.61277e-01,6.59294e-22,8.89713e-02,2.94571e-21,3.25063e-03,1.25548e-20, &
    !a 51500-51000.5 cm-1  (194.17 - 196.08 nm)
       2.55746e-01,8.49877e-24,2.94733e-01,2.06878e-23,2.86382e-01,9.30992e-23, &
       1.21011e-01,3.66239e-22,4.21105e-02,1.75700e-21,0.00000e+00,0.00000e+00, &
    !a 51000-50500.5 cm-1  (196.08 - 198.02 nm)
       5.40111e-01,7.36085e-24,2.93263e-01,2.46742e-23,1.63417e-01,1.37832e-22, &
       3.23781e-03,2.15052e-21,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &
    !a 50500-50000.5 cm-1  (198.02 - 200.00 nm)
       8.18514e-01,7.17937e-24,1.82262e-01,4.17496e-23,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &
    !a 50000-49500.5 cm-1  (200.00 - 202.02 nm)
       8.73680e-01,7.13444e-24,1.25583e-01,2.77819e-23,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &
    !a 49500-49000.5 cm-1  (202.02 - 204.08 nm)
       3.32476e-04,7.00362e-24,9.89000e-01,6.99600e-24,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00  &
    /), (/ ntr_max, nsrb /) )
    !
    bc(:,:) = RESHAPE( (/  &
    !  57000-56500.5 cm-1  (175.44 - 176.99 nm)
       1.07382e-21,9.95029e-21,7.19430e-21,2.48960e-20,2.53735e-20,7.54467e-20, &
       4.48987e-20,2.79981e-19,9.72535e-20,9.29745e-19,2.30892e-20,4.08009e-17, &
    !  56500-56000.5 cm-1  (176.99 - 178.57 nm)
       3.16903e-22,1.98251e-21,5.87326e-22,3.44057e-21,2.53094e-21,8.81484e-21, &
       8.82299e-21,4.17179e-20,2.64703e-20,2.43792e-19,8.73831e-20,1.46371e-18, &
    !  56000-55500.5 cm-1  (178.57 - 180.18 nm)
       1.64421e-23,9.26011e-22,2.73137e-22,1.33640e-21,9.79188e-22,2.99706e-21, &
       3.37768e-21,1.39438e-20,1.47898e-20,1.04322e-19,4.08014e-20,6.31023e-19, &
    !  55500-55000.5 cm-1  (180.18 - 181.82 nm)
       8.68729e-24,7.31056e-22,8.78313e-23,1.07173e-21,8.28170e-22,2.54986e-21, &
       2.57643e-21,9.42698e-21,9.92377e-21,5.21402e-20,3.34301e-20,2.91785e-19, &
    !  55000-54500.5 cm-1  (181.82 - 183.48 nm)
       1.20679e-24,2.44092e-22,2.64326e-23,4.03998e-22,2.53514e-22,8.53166e-22, &
       1.29834e-21,3.74482e-21,5.12103e-21,2.65798e-20,2.10948e-20,2.35315e-19, &
    !  54500-54000.5 cm-1  (183.48 - 185.18 nm)
       2.79656e-24,1.40820e-22,3.60824e-23,2.69510e-22,4.02850e-22,8.83735e-22, &
       1.77198e-21,6.60221e-21,9.60992e-21,8.13558e-20,4.95591e-21,1.22858e-17, &
    !  54000-53500.5 cm-1  (185.18 - 186.91 nm)
       2.36959e-24,1.07535e-22,2.83333e-23,2.16789e-22,3.35242e-22,6.42753e-22, &
       1.26395e-21,5.43183e-21,4.88083e-21,5.42670e-20,3.27481e-21,1.58264e-17, &
    !  53500-53000.5 cm-1  (186.91 - 188.68 nm)
       8.65018e-25,3.70310e-23,1.04351e-23,6.43574e-23,1.17431e-22,2.70904e-22, &
       4.88705e-22,1.65505e-21,2.19776e-21,2.71172e-20,2.65257e-21,2.13945e-17, &
    !  53000-52500.5 cm-1  (188.68 - 190.47 nm)
       9.63263e-25,1.54249e-23,4.78065e-24,2.97642e-23,6.40637e-23,1.46464e-22, &
       1.82634e-22,7.12786e-22,1.64805e-21,2.37376e-17,9.33059e-22,1.13741e-20, &
    !  52500-52000.5 cm-1  (190.47 - 192.31 nm)
       1.08414e-24,8.37560e-24,9.15550e-24,2.99295e-23,9.38405e-23,1.95845e-22, &
       2.84356e-22,3.39699e-21,1.94524e-22,2.72227e-19,1.18924e-21,3.20246e-17, &
    !  52000-51500.5 cm-1  (192.31 - 194.17 nm)
       1.52817e-24,1.01885e-23,1.22946e-23,4.16517e-23,9.01287e-23,2.34869e-22, &
       1.93510e-22,1.44956e-21,1.81051e-22,5.17773e-21,9.82059e-22,6.22768e-17, &
    !  51500-51000.5 cm-1  (194.17 - 196.08 nm)
       2.12813e-24,8.48035e-24,5.23338e-24,1.93052e-23,1.99464e-23,7.48997e-23, &
       4.96642e-22,6.15691e-17,4.47504e-23,2.76004e-22,8.26788e-23,1.65278e-21, &
    !  51000-50500.5 cm-1  (196.08 - 198.02 nm)
       3.81336e-24,7.32307e-24,5.60549e-24,2.04651e-23,3.36883e-22,6.15708e-17, &
       2.09877e-23,1.07474e-22,9.13562e-24,8.41252e-22,0.00000e+00,0.00000e+00, &
    !  50500-50000.5 cm-1  (198.02 - 200.00 nm)
       5.75373e-24,7.15986e-24,5.90031e-24,3.05375e-23,2.97196e-22,8.92000e-17, &
       8.55920e-24,1.66709e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &
    !  50000-49500.5 cm-1  (200.00 - 202.02 nm)
       6.21281e-24,7.13108e-24,3.30780e-24,2.61196e-23,1.30783e-22,9.42550e-17, &
       2.69241e-24,1.46500e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &
    !  49500-49000.5 cm-1  (202.02 - 204.08 nm)
       6.81118e-24,6.98767e-24,7.55667e-25,2.75124e-23,1.94044e-22,1.45019e-16, &
       1.92236e-24,3.73223e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00  &
    /), (/ ntr_max, nsrb /) )
       
    DO iw = 1, nsrb
       x_table(0,iw) = SUM( ac(1:11:2,iw) )
       d_table(0,iw) = SUM( bc(1:11:2,iw) )
       DO k = 1, tdim
          col = 22. + t_del*REAL(k-1,dp)
          o2_table(k) = col
          col = 10.**col
          a1(:) = ac(2:12:2,iw) * col
          b1(:) = bc(2:12:2,iw) * col
          WHERE( a1(:) < 500._dp )
             a0(:) = EXP( -a1(:) )
          ELSE WHERE
             a0(:) = 0._dp
          END WHERE
          WHERE( b1(:) < 500._dp )
             b0(:) = EXP( -b1(:) )
          ELSE WHERE
             b0(:) = 0._dp
          END WHERE
          x_table(k,iw) = DOT_PRODUCT( ac(1:11:2,iw),a0(:) )
          d_table(k,iw) = DOT_PRODUCT( bc(1:11:2,iw),b0(:) )
       END DO
    END DO
    status = 0
    RETURN
  END SUBROUTINE fubrad_srb_kck_schu_init
  !
  ! =============================================================================
  !
  SUBROUTINE fubrad_srb_kck_mem_ini (nbdim,klev)
    !
    INTEGER, INTENT(in) :: nbdim,klev
    !
    ALLOCATE (ac(ntr_max,nsrb), bc(ntr_max,nsrb))
    ALLOCATE (d_table(0:tdim,nsrb), x_table(0:tdim,nsrb))
    ALLOCATE (o2_table(tdim))
    !
    RETURN
  END SUBROUTINE fubrad_srb_kck_mem_ini
  !
  ! =============================================================================
  !
  SUBROUTINE fubrad_srb_kck_mem_clean
    !
    DEALLOCATE (ac, bc)
    DEALLOCATE (d_table, x_table)
    DEALLOCATE (o2_table)
    !
  END SUBROUTINE fubrad_srb_kck_mem_clean
  
END MODULE  messy_rad_fubrad_srb_kck
