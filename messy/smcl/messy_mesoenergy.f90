MODULE messy_mesoenergy

  USE messy_main_constants_mem, ONLY: &
    dp, &                                 ! kind parameter for real
    Pi

  IMPLICIT NONE

!-----------------------------------------------------------------
! Everything is PRIVATE, except when explicitely stated otherwise:
!-----------------------------------------------------------------
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modver = '0.1'
  CHARACTER(LEN=*), PARAMETER :: modstr = 'mesoenergy'
 
  PUBLIC :: mesoenergy_read_nml_ctrl 
  PUBLIC :: oh_calc_prod_reak
!!$  PUBLIC :: vmr_to_nd ! op_pj_20180713: moved to messy_main_tools to be 
!!$  PUBLIC :: nd_to_vmr !                 shard between submodels
  PUBLIC :: oh_calc_relaxation
  PUBLIC :: oh_calc_radlifetime
  PUBLIC :: oh_calc_quenching

  PUBLIC :: dp, modver, modstr

  LOGICAL, SAVE :: chemheat        ! switch to turn on/off chemical heating
  REAL(dp), SAVE :: rea_ent_ho3    ! reaction enthalpy for reaction H + O3 -> OH + O2
  REAL(dp), DIMENSION(9), SAVE :: fv  ! Quasi-nascent fraction of molecules produced in vibrational state v
  REAL(dp), DIMENSION(10), SAVE :: a8, a9, a10 ! Quenching coefficients

  
  PUBLIC :: chemheat, rea_ent_ho3, fv, a8, a9, a10

CONTAINS

  SUBROUTINE mesoenergy_read_nml_ctrl(status, iou)
  !--------------------------------------------------------------------
  ! This routine reads and checks the ageofair CTRL namelist. It is designed
  ! according to the MESSy (Mainz Earth Submodel System) standard. 
  !--------------------------------------------------------------------

  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  USE messy_main_blather, ONLY: start_message, end_message
  
  IMPLICIT NONE
  
  !-----
  ! input/output variables
  !-----
  INTEGER, INTENT(out) :: status   ! error status
  INTEGER, INTENT(in) :: iou       ! logical I/O unit
  
  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'mesoenergy_read_nml_ctrl'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status

  NAMELIST /CTRL/ chemheat, rea_ent_ho3, fv, a8, a9, a10
  
  CALL start_message(TRIM(modstr),'INITIALISATION', substr)
  !-----
  ! initialisation
  !-----
  status = 1

  CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
  IF (.not.lex) RETURN   ! mesoenergy.nml does not exist

  READ(iou, NML=CTRL, IOSTAT=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  IF (fstat /= 0) RETURN   ! error while reading namelist

  a8=a8*1.d-11
  a9=a9*1.d-13
  a10=a10*1.d-14

  write(*,*) ' /---------------------\ '
  write(*,*) ' | M E S O E N E R G Y | '
  write(*,*) ' \---------------------/ '
  write(*,*) ' '
  write(*,*) ' ----------------------------------'
  write(*,*) '    C A U T I O N ! ! ! ! ! ! !'
  write(*,*) '      Still under development'
  write(*,*) ' ----------------------------------'
  write(*,*) ' '
  write(*,*) ' Calculate energies in the mesosphere'
  write(*,*) '  Chemical heating for'
  if (chemheat) then
     write(*,*) '   - H + O3 -> OH + O2  (turned on)'
     write(*,*) '       Reaction Enthalpy: ',rea_ent_ho3
  ELSE
     write(*,*) '   - H + O3 -> OH + O2  (turned off)'
  ENDIF
  write(*,*) ' Quasi-nascent fraction of molecules produced in vibrational state v:', fv
  write(*,*) ' Chemical loss rate for OH(v)+O->H+O2: ', a8
  write(*,*) ' Quenching coefficients for OH(v)+O2->OH(v-1)+O2: ',a9
  write(*,*) ' Quenching coefficients for OH(v)+N2->OH(v-1)+N2: ',a10

  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error

  CALL end_message(TRIM(modstr),'INITIALISATION', substr)

END SUBROUTINE mesoenergy_read_nml_ctrl

 SUBROUTINE oh_calc_prod_reak(prod,fv,Pc)

 IMPLICIT NONE
 REAL(dp), INTENT(in) :: fv, Pc
 REAL(dp), INTENT(out) :: prod

  prod = fv*Pc

END SUBROUTINE oh_calc_prod_reak

! op_pj_20180713+: moved to messy_main_tools to be shared between submodels
!!$SUBROUTINE vmr_to_nd(nd,vmr,pressure,temp)
!!$! nd : number density in #/m-3
!!$! vmr : mixing ratio
!!$! pressure : pressure in Pa
!!$! temp : temperature in K
!!$
!!$  USE messy_main_constants_mem, ONLY: R_gas, N_A
!!$
!!$  IMPLICIT NONE
!!$  REAL(dp), INTENT(out) :: nd
!!$  REAL(dp), INTENT(in) :: vmr, pressure, temp
!!$
!!$  nd = vmr * (N_A) * pressure / (R_gas*temp)
!!$
!!$END SUBROUTINE vmr_to_nd
!!$
!!$SUBROUTINE nd_to_vmr(vmr,nd,pressure,temp)
!!$! nd : number density in #/m-3
!!$! vmr : mixing ratio
!!$! pressure : pressure in Pa
!!$! temp : temperature in K
!!$
!!$  USE messy_main_constants_mem, ONLY: R_gas, N_A
!!$
!!$  IMPLICIT NONE
!!$  REAL(dp), INTENT(out) :: vmr
!!$  REAL(dp), INTENT(in) :: nd, pressure, temp
!!$
!!$  vmr = nd / (N_A) / pressure * (R_gas*temp)
!!$
!!$END SUBROUTINE nd_to_vmr
! op_pj_20180713-

SUBROUTINE oh_calc_relaxation(relaxation_matrix, OH, Avv, time_step_len)
  IMPLICIT NONE
   REAL(dp), DIMENSION(9,6), INTENT(out) :: relaxation_matrix
   REAL(dp), DIMENSION(9,6), INTENT(in) :: Avv
   REAL(dp), DIMENSION(9), INTENT(in) :: OH
   REAL(dp), INTENT(in) :: time_step_len
   INTEGER :: i,j

   DO i=1,9
    DO j=1,6
     relaxation_matrix(i,j)=OH(i)*Avv(i,j)
    END DO
   END DO

END SUBROUTINE oh_calc_relaxation


SUBROUTINE oh_calc_radlifetime(radlifetime,Avv)
  IMPLICIT NONE
   REAL(dp), DIMENSION(9,6), INTENT(in) :: Avv
   REAL(dp), DIMENSION(9), INTENT(out) :: radlifetime
   INTEGER :: i

!   DO i=1,9
     radlifetime(:)=sum(Avv,2)
!   END DO
END SUBROUTINE oh_calc_radlifetime

SUBROUTINE oh_calc_quenching(out_numb, in_numb, kvv)
  IMPLICIT NONE
  REAL(dp), INTENT(out) :: out_numb
  REAL(dp), INTENT(in) :: in_numb
  REAL(dp), INTENT(in) :: kvv
   out_numb=in_numb*kvv
END SUBROUTINE oh_calc_quenching

END MODULE messy_mesoenergy
