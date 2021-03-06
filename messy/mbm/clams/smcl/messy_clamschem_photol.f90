Module messy_clamschem_photol

contains

subroutine photol 

  !
  !     Glenn D. Carver,             Centre for Atmospheric Science,
  !                                  University of Cambridge.
  !
  !     ASAD: photol                 Version: photol.f 4.1 01/15/97
  !     -------
  !     ASAD (V2) dummy routine merged with photol routine from
  !     ASAD(1) version for calling dissoc (Lary photolysis)
  !     Rolf Mueller 20.3.98
  !
  !     Purpose
  !     -------
  !      The role of
  !     this driver routine is to set photolysis rates to the
  !     reactions stored in the common block cmrate. The names of
  !     the species involved in each reaction are read from the
  !     ratefile and kept in common block cmreac; their indeces
  !     are kept in cmrate.
  !
  !     The rates are calculated at the frequency determined by variable
  !     nfphot in common cmcctl. This subroutine is ALWAYS called on the
  !     first chemical timestep and may be called every chemical substep.
  !
  !     The pressure, temperature and number density are kept in the
  !     common block cmphys.
  !
  !     If the calculation of the photolysis rates involves computing
  !     the zenith angle, then position information (e.g. latitude and
  !     longitude) and time of year should be passed via common from the
  !     calling model so that this routine may compute the zenith angle.
  !     Then use the variables jsubs and cdt from common
  !     cmcctl to compute the time into the chemical substepping.
  !
  !     Interface.
  !     ----------
  !     Called from cdrive.
  !
  !---------------------------------------------------------------------

  !-----------------------------------------------
  !   I n t e r f a c e   B l o c k s
  !-----------------------------------------------
  USE messy_main_constants_mem, ONLY: DP

  use messy_clamschem_asad_mod,       only: t, p, rk, nprkx
  use messy_clamschem_asad_mod_clams, only: jppj

  !USE messy_clamschem_phys_const,ONLY: pi, torad
  USE messy_cmn_photol_mem,     ONLY: IP_MAX

  USE messy_clamschem_global,   ONLY: ntraj, pi, torad, dissoc_rate, &
                                      jpnl_count, jpnl_offset, ipart, &
                                      iodump, iphot, cindex, zangle
  use messy_clamschem_read_roeth, only: read_roeth

  implicit none

  !-----------------------------------------------
  !   L o c a l   P a r a m e t e r s
  !-----------------------------------------------
  integer, parameter :: njsr = 110 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: ireac, itr, ij, i, itr_dissoc
  real(kind=DP) :: sza_rad
  real(kind=DP), dimension(njsr),save :: js1r, a, b, c 

  ! iphot=1  means call of DISSOC, 
  ! the "Roeth photolysis" is currently disabled
!!!!!
!  integer, parameter:: iphot=1

!-----------------------------------------------
!
!
!     include any commons here for getting position and time info
!     from the calling model
!
!     njs = (maximum) number of photolysis rates to be assigned in subroutine
!     dissoc; njs must be >= jpdim-1 that is defined in subroutine inphot
!     (jug, 01/2014) njs replaced by IT_MAX      
!
!
!
!     ---------------------------------------------------------------
!          1. Calculate photolysis rates
!             --------- ---------- -----
!
!
 
  if (iodump) write (6, *) 'photol called' 
  if (iphot == 2) then 
     
     write(*,*) 'Roeth parametrisation needs to be adapted'
     stop
     call read_roeth (njsr, a, b, c) 
  endif


  do itr = 1, ntraj 

     itr_dissoc = jpnl_offset(ipart) + (itr - 1)
     !if (itr==1) write (*,*) 'itr, itr_dissoc=',itr,itr_dissoc

     sza_rad = zangle(itr)*torad
     !
     !     calculate photolysis rate
     if (iphot /= 2) then 
        !     Lary - Photyolysis
        !-----------------------------------------------------------------------
        !

        ! call to photolysis routine 
        ! call dissoc (t(itr), p(itr), lats(itr), sza_rad, photrates)
        ! now in messy_dissoc_si
        ! photrates now transferred through module messy_clamschem


        do ij = 1, jppj 
 
!JUG     next line changed to .ge. since cindex = 0 for first photolysis rate
           if (cindex(ij) >= 0) then 
!!$              rk(itr,nprkx(ij)) = photrates(itr, cindex(ij)) 
              rk(itr,nprkx(ij)) = dissoc_rate(cindex(ij))%values(itr_dissoc) 
           else 
              rk(itr,nprkx(ij)) = 0.0D0 
           endif
           !if (iodump) write (6,*) 'PHOTO_INFO', ij, cindex(ij), rk(itr,nprkx(ij)) 
           !if (itr==1) write (6,*) 'PHOTO_INFO', ij, cindex(ij), rk(itr,nprkx(ij)) 
        end do
        
 
     else if (iphot == 2) then 
        !           Roeth - Photyolysis
        !-----------------------------------------------------------------------
        if (zangle(itr) >= 100.0) then 
           js1r(:njsr) = 0.0 
        else 
           do ireac = 1, njsr 
              if (c(ireac)*sza_rad >= pi/2.) then 
                 js1r(ireac) = 0.0 
              else 
                 js1r(ireac) = a(ireac)*exp(b(ireac)*(1.-1./cos(c(ireac)*&
                      sza_rad))) 
              endif
           end do
        endif
 
        do ij = 1, jppj 
           ! JUG  next line changed to .ge. since cindex = 0 for first
           ! photolysis rate
           if (cindex(ij) >= 0) then 
              rk(itr,nprkx(ij)) = js1r(cindex(ij)) 
           else 
              rk(itr,nprkx(ij)) = 0.0D0 
           endif
           if(iodump) write (6, *) 'PHOTO INFO', ij, cindex(ij), rk(1,nprkx(ij)) 
        end do
     endif 
!-------------------------------------------------------------------------------
  end do 
! op_pj_20160616+
!!$ if (iodump) write(*,'(15h photol JO2info: p10e12.4)')  rk(1:10,3)
  if (iodump) write(*,*) 'photol JO2info:', rk(1:10,3)
! op_pj_20160616-
  if (iodump) write (6, *) 'photol finished!' 
  return  
end subroutine photol 
 
 
end Module messy_clamschem_photol


