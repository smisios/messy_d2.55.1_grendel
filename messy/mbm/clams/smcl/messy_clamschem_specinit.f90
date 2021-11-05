Module messy_clamschem_specinit

Contains

! For MESSY: large parts commented out for simplification

      subroutine specinit
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!      use netcdf

      ! MESSy Main
      use messy_main_constants_mem, only: DP
      USE messy_main_constants_mem, ONLY: DP

      ! ASAD
      use messy_clamschem_asad_mod,       only: speci, ctype, family, jpfm, jpif, moffam, advt
      use messy_clamschem_asad_mod_clams, only: jpspec

      ! CLaMS
      USE messy_clams_global,        ONLY: rank
      USE messy_clamschem_global,    ONLY: iodump, nfamily, ftrindex, fnames

      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer       :: ip, if
      character(10) :: lastfam
!-----------------------------------------------
!

!     Expand the dataset with  spaces for the chemical species not at this stage
!     we want space for all of the individual species and not just the sum over
!     familiy values.  this means that we need to read CHCH.d alternatively
!     we can access this internal information directly from
!     common block CMCHCH which gives all the info that we want.
!
!     From now on species will be denote the resolved species
!     and tracer will denote the values that are passed to the dynamic code
!     the chemical code unbundles the families internally and so we must access
!     these values at the end of a chemical timestep by reading from the
!     array tnd in the common block cmspec that must be included the species
!     number densities are stored in y and to convert to mixing ratio must be
!     divided by TND in common block CMPHYS.
!     i.e. both these arrays must be included if the complete list of species
!     is to be accessed.  We will also store the TND value as a new parameter
!     as it will be useful to have this available.
!     values
!
!     18.7.2001 added output of species emissions to NetCDF file  (jug)
!     17.10.2008 changed the way of duplicate reactands in output file 
!                from possible 3 to up to 26 (jug)
 
 
!  Now identify the families and generate
!  ftrindex
   nfamily = 0 
   lastfam = ''

!  First check to see if the family already exists
!  If RCODE = -1 then it probably does not and should be created
!  Otherwise the creation step should be skipped.
 
      do ip = 1, jpspec 
         if (iodump .and. rank==0) write (*,*) 'nfamily = ', nfamily 
         if (ctype(ip)/=jpfm .and. ctype(ip)/=jpif) cycle 
         if (family(ip) == lastfam) cycle  
         lastfam=family(ip)
         nfamily = nfamily + 1 
         ftrindex(nfamily) = moffam(ip)

      end do
     
      if (iodump) then 
         do if = 1, nfamily 
            if (rank==0) write (*,*) 'family ', advt(ftrindex(if)), ' identified' 
            fnames(if) = advt(ftrindex(if)) 
         end do 
      else 
         fnames(:nfamily) = advt(ftrindex(:nfamily)) 
      endif 
 
      if(rank==0) then
 
   endif
   return  
 
 end subroutine specinit
 
End Module messy_clamschem_specinit
