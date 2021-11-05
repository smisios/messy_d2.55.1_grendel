!****************************************************************************
!                Time-stamp: <2011-02-17 13:00:30 joec_pa>
!****************************************************************************

MODULE messy_main_compilerinfo_mem

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! THIS WILL BE REPLACED BY configure and/or gmake
  CHARACTER(LEN=*), PARAMETER :: compiler_version  = 'ifort (IFORT) 19.1.1.217 20200306'
  CHARACTER(LEN=*), PARAMETER :: compiler_call     = 'mpifort'
  CHARACTER(LEN=*), PARAMETER :: compiler_flags    = '-O3 -fpp -heap-arrays -fp-model strict -lpthread -save-temps -fno-alias -alig&
&n all'
  CHARACTER(LEN=*), PARAMETER :: compiler_cppdefs  = '-DMESSY -DLITTLE_ENDIAN -D_LINUX64 -DPNCREGRID -DMESSYMMD -DMPIOM_13B'
  CHARACTER(LEN=*), PARAMETER :: compiler_includes = '-I/home/stergios/sw/netcdf-3.6.3-intel/include -I/home/stergios/sw/netcdf-3.6&
&.3-intel/include       '
  
END MODULE messy_main_compilerinfo_mem
