MODULE messy_clamschem_data

USE messy_clamschem_defs_mod,   ONLY: chch_t, ratb_t, rath_t, ratj_t, ratt_t,  &
                                      chch_defs, ratb_defs, rath_defs,         &
                                      ratj_defs, ratt_defs

USE messy_clamschem_data_clim3, ONLY: chch_defs_clim3, ratb_defs_clim3,        &
                                      ratj_defs_clim3, rath_defs_clim3,        &
                                      ratt_defs_clim3, &
                                      clams_chem_init_clim3, clams_chem_clean_clim3

USE messy_clamschem_data_standard, ONLY: chch_defs_std, ratb_defs_std,        &
                                         ratj_defs_std, rath_defs_std,        &
                                         ratt_defs_std, &
                                         clams_chem_init_standard, &
                                         clams_chem_clean_standard

USE messy_clamschem_data_ust, ONLY: chch_defs_ust, ratb_defs_ust,            &
                                 ratj_defs_ust, rath_defs_ust,          &
                                 ratt_defs_ust,                         &
                                 clams_chem_init_ust, clams_chem_clean_ust

USE messy_clamschem_data_cfc, ONLY: chch_defs_cfc, ratb_defs_cfc,            &
                                 ratj_defs_cfc, rath_defs_cfc,          &
                                 ratt_defs_cfc,                         &
                                 clams_chem_init_cfc, clams_chem_clean_cfc

USE messy_clamschem_data_vsls, ONLY: chch_defs_vsls, ratb_defs_vsls,     &
                                 ratj_defs_vsls, rath_defs_vsls,         &
                                 ratt_defs_vsls,                         &
                                 clams_chem_init_vsls, clams_chem_clean_vsls

USE messy_clamschem_asad_mod_clams, ONLY: jpctr, jpspec, jpbk, jptk, jphk,  &
                                          jppj, jpdd, jpdw, &
                                          jpnr, jpbkp1, jptkp1, jppjp1, jphkp1, &
                                          chemdata_type, mype

IMPLICIT NONE

CONTAINS

  SUBROUTINE clams_chem_init

    USE messy_clamschem_asad_mod,   ONLY: method
    USE messy_clamschem_asad_dummy, ONLY: ereport
    
    IMPLICIT NONE

    INTEGER :: errcode                ! Variable passed to ereport
    INTEGER :: i
    INTEGER :: isizes(5)
    CHARACTER(LEN=72) :: cmessage

    chemdata_type = uppercase(chemdata_type)

    SELECT CASE (chemdata_type)
    CASE ('CLIM3') 
       jpctr = 10
       jpspec = 15
       jpbk = 7
       jptk = 2
       jppj = 4
       jphk = 11
    CASE ('STANDARD') 
       jpctr = 27
       jpspec = 40
       jpbk = 64
       jptk = 13
       jppj = 27
       jphk = 11
    CASE ('UST') 
       jpctr = 30
       jpspec = 43
       jpbk = 76
       jptk = 14
       jppj = 30
       jphk = 11
    CASE ('CFC') 
       jpctr = 36
       jpspec = 49
       jpbk = 83
       jptk = 14
       jppj = 36
       jphk = 11
    CASE ('VSLS') 
       jpctr = 44
       jpspec = 57
       jpbk = 100
       jptk = 14
       jppj = 45
       jphk = 11

    CASE default
       cmessage = ' Unknown chemistry, unable to set chch_defs'
       errcode = 1
       CALL EREPORT('CLAMS_CHEM_INIT',errcode,cmessage)
    END SELECT
 

    jpnr   = jpbk + jptk + jppj + jphk 
    jpbkp1 = jpbk + 1
    jptkp1 = jptk + 1 
    jppjp1 = jppj + 1
    jphkp1 = jphk + 1


    IF (.NOT. ALLOCATED(chch_defs)) ALLOCATE(chch_defs(jpspec))
    IF (.NOT. ALLOCATED(ratb_defs)) ALLOCATE(ratb_defs(jpbk))
    IF (.NOT. ALLOCATED(rath_defs)) ALLOCATE(rath_defs(jphk))
    IF (.NOT. ALLOCATED(ratj_defs)) ALLOCATE(ratj_defs(jppj))
    IF (.NOT. ALLOCATED(ratt_defs)) ALLOCATE(ratt_defs(jptk))

    SELECT CASE (chemdata_type)
    CASE ('CLIM3') 
       if (mype==0) then
          write (*,*)
          write (*,*) 'Chemistry input: ',trim(chemdata_type)
          write (*,*)
       endif
       CALL clams_chem_init_clim3
       isizes(1) = SIZE(chch_defs_clim3)
       isizes(2) = SIZE(ratb_defs_clim3)
       isizes(3) = SIZE(rath_defs_clim3)
       isizes(4) = SIZE(ratj_defs_clim3)
       isizes(5) = SIZE(ratt_defs_clim3)
       CALL CHECK_DIMS
       chch_defs = chch_defs_clim3
       ratb_defs = ratb_defs_clim3
       rath_defs = rath_defs_clim3
       ratj_defs = ratj_defs_clim3
       ratt_defs = ratt_defs_clim3
       CALL clams_chem_clean_clim3
    CASE ('STANDARD') 
       if (mype==0) then
          write (*,*)
          write (*,*) 'Chemistry input: ',trim(chemdata_type)
          write (*,*)
       endif
       CALL clams_chem_init_standard
       isizes(1) = SIZE(chch_defs_std)
       isizes(2) = SIZE(ratb_defs_std)
       isizes(3) = SIZE(rath_defs_std)
       isizes(4) = SIZE(ratj_defs_std)
       isizes(5) = SIZE(ratt_defs_std)
       CALL CHECK_DIMS
       chch_defs = chch_defs_std
       ratb_defs = ratb_defs_std
       rath_defs = rath_defs_std
       ratj_defs = ratj_defs_std
       ratt_defs = ratt_defs_std
       CALL clams_chem_clean_standard
    CASE ('UST') 
       if (mype==0) then
          write (*,*)
          write (*,*) 'Chemistry input: ',trim(chemdata_type)
          write (*,*)
       endif
       CALL clams_chem_init_ust
       isizes(1) = SIZE(chch_defs_ust)
       isizes(2) = SIZE(ratb_defs_ust)
       isizes(3) = SIZE(rath_defs_ust)
       isizes(4) = SIZE(ratj_defs_ust)
       isizes(5) = SIZE(ratt_defs_ust)
       CALL CHECK_DIMS
       chch_defs = chch_defs_ust
       ratb_defs = ratb_defs_ust
       rath_defs = rath_defs_ust
       ratj_defs = ratj_defs_ust
       ratt_defs = ratt_defs_ust
       CALL clams_chem_clean_ust
    CASE ('CFC') 
       if (mype==0) then
          write (*,*)
          write (*,*) 'Chemistry input: ',trim(chemdata_type)
          write (*,*)
       endif
       CALL clams_chem_init_cfc
       isizes(1) = SIZE(chch_defs_cfc)
       isizes(2) = SIZE(ratb_defs_cfc)
       isizes(3) = SIZE(rath_defs_cfc)
       isizes(4) = SIZE(ratj_defs_cfc)
       isizes(5) = SIZE(ratt_defs_cfc)
       CALL CHECK_DIMS
       chch_defs = chch_defs_cfc
       ratb_defs = ratb_defs_cfc
       rath_defs = rath_defs_cfc
       ratj_defs = ratj_defs_cfc
       ratt_defs = ratt_defs_cfc
       CALL clams_chem_clean_cfc
    CASE ('VSLS') 
       if (mype==0) then
          write (*,*)
          write (*,*) 'Chemistry input: ',trim(chemdata_type)
          write (*,*)
       endif
       CALL clams_chem_init_vsls
       isizes(1) = SIZE(chch_defs_vsls)
       isizes(2) = SIZE(ratb_defs_vsls)
       isizes(3) = SIZE(rath_defs_vsls)
       isizes(4) = SIZE(ratj_defs_vsls)
       isizes(5) = SIZE(ratt_defs_vsls)
       CALL CHECK_DIMS
       chch_defs = chch_defs_vsls
       ratb_defs = ratb_defs_vsls
       rath_defs = rath_defs_vsls
       ratj_defs = ratj_defs_vsls
       ratt_defs = ratt_defs_vsls
       CALL clams_chem_clean_vsls

    CASE default
       cmessage = ' Unknown chemistry, unable to set chch_defs'
       errcode = 1
       CALL EREPORT('CLAMS_CHEM_INIT',errcode,cmessage)
    END SELECT

    ! N-R solver / SVODE
    ! Ignore Family and steady-state settings in the case of 
    ! N-R solver (3) and SVODE (11) and count new dimension jpctr
    IF (method==3) THEN
       jpctr = 0
       do i = 1, jpspec
          if (chch_defs(i)%ctype=='FM' .or. chch_defs(i)%ctype=='SS') then
             chch_defs(i)%ctype = 'TR'
          endif
          if (chch_defs(i)%speci=='O(1D)' .and. chemdata_type /= 'CLIM3') then
             chch_defs(i)%ctype = 'SS'
          endif
          if (chch_defs(i)%ctype=='TR') jpctr = jpctr + 1
       enddo
    ENDIF

    IF (method==11) THEN
       jpctr = 0
       do i = 1, jpspec
          if (chch_defs(i)%ctype=='FM' .or. chch_defs(i)%ctype=='SS') then
             chch_defs(i)%ctype = 'TR'
          endif
          if (chch_defs(i)%ctype=='TR') jpctr = jpctr + 1
       enddo
    ENDIF

  CONTAINS

    SUBROUTINE CHECK_DIMS

      implicit none

      integer :: ierr

      IF (isizes(1) /= jpspec) THEN
         cmessage= ' chch_defs: size not equal to jpspec'
         write(6,*) 'SIZES: ',isizes(1),' AND ',jpspec,' ARE NOT EQUAL!'
         ierr=91
      ELSE IF (isizes(2) /= jpbk) THEN
         cmessage= ' ratb_defs: size not equal to jpbk'
         write(6,*) 'SIZES: ',isizes(2),' AND ',jpbk,' ARE NOT EQUAL!'
         ierr=92
      ELSE IF (isizes(3) /= jphk) THEN
         cmessage= ' rath_defs: size not equal to jphk'
         write(6,*) 'SIZES: ',isizes(3),' AND ',jphk,' ARE NOT EQUAL!'
         ierr=93
      ELSE IF (isizes(4) /= jppj) THEN
         cmessage= ' ratj_defs: size not equal to jppj'
         write(6,*) 'SIZES: ',isizes(4),' AND ',jppj,' ARE NOT EQUAL!'
         ierr=94
      ELSE IF (isizes(5) /= jptk) THEN
         cmessage= ' ratt_defs: size not equal to jptk'
         write(6,*) 'SIZES: ',isizes(5),' AND ',jptk,' ARE NOT EQUAL!'
         ierr=95
      ELSE
         ierr=0
      ENDIF
      IF (ierr > 1) THEN
         CALL EREPORT('CLAMS_CHEM_INIT.check_dims',ierr,cmessage)
      ENDIF
      
    END SUBROUTINE CHECK_DIMS

    !**************************************************************************
    ! convert string to uppercase
    !**************************************************************************
    function uppercase (str)

      implicit none
      
      character(*) :: str
      character(len(str)) :: uppercase
      
      character(26), parameter :: char_up='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(26), parameter :: char_low='abcdefghijklmnopqrstuvwxyz'
      
      integer :: i, pos
      
      uppercase = str
      
      do i = 1, len_trim(str)
         
         pos = index(char_low, str(i:i))
         if (pos /= 0) uppercase(i:i) = char_up(pos:pos)
         
      enddo
      
    end function uppercase
    
  END SUBROUTINE clams_chem_init

END MODULE messy_clamschem_data
