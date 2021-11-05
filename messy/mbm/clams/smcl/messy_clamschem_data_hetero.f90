Module messy_clamschem_data_hetero


CONTAINS

  Subroutine clams_chem_init_hetero

    USE messy_clamschem_global, ONLY: hetvar, nhetspec

    implicit none

    integer :: i

    ! hetvar(1:nhetspec) => added to SPECARR 
    ! nhetspec=9 !

    hetvar(1)%name  = 'HCl(c)'
    hetvar(2)%name  = 'HNO3(c)'
    hetvar(3)%name  = 'H2O(c)'
    hetvar(4)%name  = 'HBr(c)'
    do i = 1,4
       hetvar(i)%units = 'm^3/m^3'
       hetvar(i)%longname = hetvar(i)%name(1:len_trim(hetvar(i)%name)-3)//'(cond.)'
   enddo

    hetvar(5)%name     = 'N(aero)'
    hetvar(5)%longname = 'Aerosol Number Density'
    hetvar(6)%name     = 'N(SAT)'
    hetvar(6)%longname = 'SAT Number Density'
    hetvar(7)%name     = 'N(NAT)'
    hetvar(7)%longname = 'NAT Number Density'
    hetvar(8)%name     = 'N(ice)'
    hetvar(8)%longname = 'Ice Number Density'
    do i = 5,8
       hetvar(i)%units = 'cm^-3' 
    enddo

    hetvar(9)%name     = 'aer_H2SO4'
    hetvar(9)%longname = 'Particle H2SO4 as equivalent gasphase mixing ratio'
    hetvar(9)%units    = 'ppb' 


    ! hetvar (nhetspec+1:nhetspec+nhetpar) only for control output (channel clamschem_HETPAR)
    ! nhetspec=9, nhetpar=11: hetvar (10:20) 

    hetvar(nhetspec+1)%name  = 'WT_S'
    hetvar(nhetspec+2)%name  = 'WT_N'
    hetvar(nhetspec+3)%name  = 'WT_Cl'
    hetvar(nhetspec+4)%name  = 'WT_Br'
    hetvar(nhetspec+5)%name  = 'WT_HOCl'
    hetvar(nhetspec+6)%name  = 'WT_HOBr'
    do i = nhetspec+1,nhetspec+6
       hetvar(i)%units = '%'
       hetvar(i)%longname = 'Weight Percent '// hetvar(i)%name(4:len_trim(hetvar(i)%name))
    enddo

    hetvar(nhetspec+7)%name = 'A(ice)'
    hetvar(nhetspec+8)%name = 'A(NAT)'
    hetvar(nhetspec+9)%name = 'A(liq)'
    hetvar(nhetspec+10)%name = 'A(SAT)'
    do i = nhetspec+7,nhetspec+10
       hetvar(i)%units = 'cm^2/cm^3' 
       hetvar(i)%longname = 'Surface Area  '// hetvar(i)%name(3:5)
    enddo

    hetvar(nhetspec+11)%name     = 'V(liq)'
    hetvar(nhetspec+11)%longname = 'Liquid Aerosol Volume Density'
    hetvar(nhetspec+11)%units     = 'cm^3 cm^-3' 


  End Subroutine clams_chem_init_hetero


End Module messy_clamschem_data_hetero
