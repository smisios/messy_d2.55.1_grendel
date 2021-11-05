Module messy_clamschem_asad_dummy



Contains

 
  subroutine ereport (str,errcode,errstr)
    character(*) :: str, errstr
    integer      :: errcode

    write (*,*) 'ERROR: ', errcode
    write (*,*) str,': ',errstr 

!!$#ifdef USE_MPI
#ifndef NOMPI
   !CALL MPI_ABORT (MPI_COMM_WORLD, 0, p_error)
   CALL MPI_FINALIZE(ierr)
#endif    

    stop

  end subroutine ereport

  subroutine asad_jac_dummy
    return
  end subroutine asad_jac_dummy


  subroutine inwdep 
    return  
  end subroutine inwdep
  
  subroutine inddep 
    return  
  end subroutine inddep
 
 
  subroutine wetdep 
    return  
  end subroutine wetdep
 
  subroutine drydep 
    return  
  end subroutine drydep


End Module messy_clamschem_asad_dummy
