MODULE mo_bgc_iolist

#ifndef NO_NEW_IO

  USE mo_file_list, ONLY: file_list
  USE mo_iolist, ONLY: iolist_read_config, iolist_create_file_list

  USE mo_control_bgc,ONLY : io_stdi_bgc
  USE mo_bgc_varlist, ONLY: bgc_varlist

  IMPLICIT NONE

  TYPE (file_list), SAVE :: bgc_iolist

CONTAINS


  SUBROUTINE read_namelist_bgc(ierror)

    INTEGER, INTENT(out) :: ierror

    CALL iolist_read_config(io_stdi_bgc, ierror)

  END SUBROUTINE read_namelist_bgc


  SUBROUTINE build_iolist_bgc

    CALL iolist_create_file_list(bgc_iolist, bgc_varlist)

  END SUBROUTINE build_iolist_bgc

#endif/*ndef NO_NEW_IO */

END MODULE mo_bgc_iolist
