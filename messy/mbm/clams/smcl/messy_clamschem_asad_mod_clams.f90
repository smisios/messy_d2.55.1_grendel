Module messy_clamschem_asad_mod_clams

  !    Set up for CLaMS chemistry 
  !
  ! Parameter    Description
  ! ---------    -----------
  ! jpctr        Number of model tracers/families in the model.
  ! jpspec       Number of species in the chemistry.
  ! jpbk         Number of bimolecular reactions. May be 0.
  ! jptk         Number of trimolecular reactions. ""
  ! jppj         Number of photolysis reactions. ""
  ! jphk         Number of heterogeneous reactions. ""
  ! jpnr         Total number of reactions.
  ! jpbkp1       jpbk + 1. Used to size reaction rate arrays since
  !              the user may set these parameters to zero. ie. we
  !              oversize the arrays by one.
  ! jptkp1       jptk + 1.
  ! jppjp1       jppj + 1.
  ! jphkp1       jphk + 1.
  ! lhet         .true. if heterogeneous reactions are to be included
  !              in the chemistry.
  ! nfhet        Frequency of calls to routine hetero to compute
  !              heterogeneous rates. If set to zero, hetero is only called
  !              once. If nonzero and postive hetero is called once
  !              every nfhet chemical substeps. If negative, hetero is
  !              called every -nfphot seconds. Note that hetero is 
  !              ALWAYS called for the first chemical substep.
  ! hetdt        Timestep for heterogeneous chemistry (GG, 30.8.2000)
       
  ! bisher in include/paramcdef.h:
  integer :: jpctr, jpspec, jpbk, jptk, jppj, jphk
  ! => Initialisierung in clams_chem_init (abhaengig von chemdata_type!)

  ! bisher in include/paramc.h
  integer :: jpnr, jpbkp1, jptkp1, jppjp1, jphkp1
  ! => Initialisieren in clams_chem_init (abhaengig von chemdata_type!)

  ! jpspb, jpspt, jpspj, jpsph (bisher in include/paramc.h) 
  !  => Konstanten in messy_clamschem_asad_mod 

  ! bisher in ASAD-2.0/src/icgx/cmycn.h:
  integer :: jl_current
 
  ! bisher in ASAD-2.0/src/cmcctl.h
  logical, dimension(:), allocatable :: failed 
  logical                            :: lhet, lphotol
  integer                            :: nfhet
  real(kind(0d0))                    :: hetdt

!!!!! NEU: 
  integer, parameter :: jpdw = 0   ! wet dep rates
  integer, parameter :: jpdd = 0   ! Dry dep rates
  integer, parameter :: model_levels = 1
  integer            :: mype
!!!!! auf welche Werte setzen ???
  integer, parameter :: prstatus_oper = 0
  integer, parameter :: prstatus_normal = 0
  integer, parameter :: prstatus_diag = 0
  integer            :: printstatus = 0
  
  ! max. number of non-advected tracers ???
  integer, parameter :: n_chem_diags = 22 !???
  ! => wird in asad_inrats genutzt !

  ! Fuer Dimensionierungen in messy_clamschem_asad_mod: 
  integer :: theta_field_size  ! bisher jpnl !!!

  ! Welcher Chemie-Input wird genutzt?
  ! -> wird in clams_chem_init abgefragt
  ! -> Chemie-Input-Typ kann aus chem.inp eingelesen werden
  character(10) :: chemdata_type = 'clim3'

End Module messy_clamschem_asad_mod_clams
