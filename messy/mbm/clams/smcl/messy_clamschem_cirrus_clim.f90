Module messy_clamschem_cirrus_clim

Contains


!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function densice_cirrus (ch2o,parth2o)

  USE messy_clams_global, ONLY: prec

  implicit none

  real(PREC),intent(in) :: ch2o,parth2o
  real(PREC) ::   densice_cirrus


  real(PREC) :: iwc_ppm,p_fit,iwc_g_m3

! Constants
  real(PREC),parameter:: m_h2o = 18.       ! g/mol
  real(PREC),parameter:: R     = 8.314     ! J/Kmol 
  real(PREC),parameter:: N_A = 6.022d23    ! molec/mol
  real(PREC),parameter:: rho_ice = 0.92    ! g/cm3
  real(PREC),parameter:: pi = 3.14159265358979324


! Particle size assumptions for cirrus:

! real(PREC),parameter:: ri = 5.  !min (um)
  real(PREC),parameter:: ri = 10.E-4 !mean (cm)
! real(PREC),parameter:: ri = 20. !max (um)

  densice_cirrus =   ch2o * (1.0 -parth2o)* m_h2o/( N_A * rho_ice) * 0.75/(pi*ri**3)

  ! Info bei Initialisierung...
  if (parth2o < -99.)densice_cirrus=ri

end function densice_cirrus


!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function rhice_freeze_clim(temp)
  ! Tempeartur-abhaengige Gefrierschwellen
  ! laut Martina Kraemer abgeleitet aus AIDA Messungen

  USE messy_clams_global, ONLY: prec

  implicit none

  real(PREC),intent(in) :: temp
  real(PREC) ::  rhice_freeze_clim 

! real(PREC),parameter:: a0= 242.3,  a1=0.41 ! sulfuric acid (Koop)
! character(16),parameter::part_type='sulfuric acid (Koop)'

! real(PREC),parameter:: a0= 305.375,a1=0.717! sulfuric acid
! character(16),parameter::part_type=' sulfuric acid'

  real(PREC),parameter:: a0= 230., a1=0.43333! soot+coating
  character(16),  parameter:: part_type='soot+coating'

! real(PREC),parameter:: a0= 206., a1=0.4    ! soot
! character(16),parameter::part_type='soot'

! real(PREC),parameter:: a0= 224., a1=0.48333! ammonium sulfate
! character(16),parameter::part_type='ammonium sulfate'

! real(PREC),parameter:: a0= 134., a1=0.1    ! mineral dust
! character(16),parameter::part_type='mineral dust'
  
  rhice_freeze_clim=max(1.1d0,(a0-a1*temp)/100.)
  rhice_freeze_clim=min(rhice_freeze_clim,1.8d0)
  if (temp==0.d0) write(*,*)'Ice nucleation on '//part_type
end function rhice_freeze_clim

End Module messy_clamschem_cirrus_clim
