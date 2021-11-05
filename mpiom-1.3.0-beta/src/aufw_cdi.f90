SUBROUTINE aufw_cdi
#ifndef NOCDI
  !****************************************************************
  !
  !     aaaaaa  u    u  ffffff  w    w
  !     a    a  u    u  f       w    w
  !     aaaaaa  u    u  fffff   w ww w
  !     a    a  u    u  f       wwwwww
  !     a    a  uuuuuu  f       ww  ww
  !
  !
  !*****************************************************************
  !     this sbr writes routinely at certain time intervals
  !     a restart file on z37000 or z38000
  !
  !-----------------------------------------------------------------

  USE mo_kind
  USE mo_param1
  USE mo_mpi
  USE mo_parallel
  USE mo_commo1

  USE mo_units


  USE mo_cdi

  IMPLICIT NONE

  INTEGER  :: idate,itime,n,ne,k
  INTEGER  :: streamid,gridid,status,zaxis2id         ! cdi stream handler
  INTEGER  :: vlistid,varid,varid1,gridid1     ! cdi variable list handler
  INTEGER  :: zaxisid,zaxis1id,taxisid   ! cdi variable list handler
  real(wp) :: zlevels(ke),zlevels2(kep),rrdt(2)


  CALL create_restart_varlist

  IF(p_pe==p_io) THEN

    gridid = gridcreate(grid_curvilinear, (ie_g*je_g))
    CALL griddefxsize(gridid, ie_g)
    CALL griddefysize(gridid, je_g)

    gridid1 = gridcreate(grid_generic,2)

    zaxisid = zaxiscreate(zaxis_depth_below_sea, ke)
    zaxis2id = zaxiscreate(zaxis_depth_below_sea, ke+1)
    zaxis1id = zaxiscreate(zaxis_surface, 1)

!    DO k = 1, ke
!      zlevels(k) = tiestu(k)
!    END DO

    CALL zaxisdeflevels(zaxisid, tiestu)

!    DO k = 1, kep
!      zlevels(k) = tiestw(k)
!    END DO
    CALL zaxisdeflevels(zaxis2id, tiestw)

    CALL zaxisdeflevels(zaxis1id, 0)
    
    vlistid = vlistcreate()

    taxisid = taxiscreate(taxis_relative)
    CALL vlistdeftaxis(vlistid, taxisid)

    idate=(lyears*10000)+(lmonts*100)+ldays
    itime=1200
    CALL taxisdefvdate(taxisid, idate)
    CALL taxisdefvtime(taxisid, itime)


    rrdt(1)=real(ldt)
    rrdt(2)=real(ldays)
 
    ne = SIZE(variables)
      
    varid1 = vlistdefvar(vlistid,gridid1,zaxis1id, time_constant)
    call vlistdefvardatatype(vlistid,varid1,datatype_flt64)    

    DO n = 1, ne

      if     ( variables(n)%dim .eq. ke  ) then 

         varid  = vlistdefvar(vlistid, gridid, zaxisid, time_variable)

      elseif ( variables(n)%dim .eq. kep ) then

         varid  = vlistdefvar(vlistid, gridid, zaxis2id, time_variable)

      elseif  ( variables(n)%dim .eq. 1 ) then

         varid  = vlistdefvar(vlistid, gridid, zaxis1id, time_variable)
         write(0,*)'in aufw_cdi: write: ',n,vlistid,variables(n)%name,variables(n)%dim,variables(n)%varid,variables(n)%array

      endif

      CALL vlistdefvarname(vlistid, varid, variables(n)%name )

      variables(n)%varid = varid

      call vlistdefvarcode(vlistid, varid, variables(n)%code )
      call vlistdefvardatatype(vlistid,varid,datatype_flt64)  


    ENDDO

    streamid = streamopenwrite('z37000.ext',filetype_ext)

    CALL streamdefvlist(streamid, vlistid)

    status = streamdeftimestep(streamid,0)       

    call streamwritevar(streamid,varid1,rrdt,0)
    
    DO n=1,ne
!       write(io_stdout,*)'in aufw_cdi: write: ',n,streamid,variables(n)%name,variables(n)%dim,variables(n)%varid
!      write(io_stdout,*)'in aufw_cdi: write: ',n,variables(n)%name,size(variables(n)%array)
       CALL streamwritevar(streamid,variables(n)%varid,variables(n)%array,0)
!      write(io_stdout,*)'written',streamid,variables(n)%varid,variables(n)%name 
    ENDDO
  
    CALL streamclose(streamid)

  ENDIF
#endif
  END SUBROUTINE aufw_cdi

