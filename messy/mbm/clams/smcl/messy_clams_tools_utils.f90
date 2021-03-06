!******************************************************************************
! File utils/src/utils.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! N. Thomas, J. Ankenbrand, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Thu Apr 11 13:30:41 2019
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version. This program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more
! details. You should have received a copy of the GNU General Public
! License along with this program; if not, see <https://www.gnu.org/licenses/>.
!
!******************************************************************************
!
! Module utils
!
! This module contains the following subroutines:
!   function lowercase (str)
!   function uppercase (str)
!   function str_found (nstr, str_array, str)
!   function str_pos (nstr, str_array, str)
!   function delete_control_chars (str)
!   subroutine isnumber(str,valid,number)
!   subroutine n_numbers(str,n,valid,numbers)
!   subroutine bubble_sort (x,dim)
!   subroutine bubble_sort_index (x,ind,dim)
!   subroutine quick_sort(array)
!   subroutine quick_sort_index(array,ind)
!   function get_metfilename (prefix, dirname, iyear, imon, iday, ihr)
!   subroutine pack_values (array, mask, mdi)
!
!**************************************************************************

MODULE messy_clams_tools_utils

PRIVATE ::  quicksort_index 

INTERFACE bubble_sort_index
   MODULE PROCEDURE bubble_sort_index_real 
   MODULE PROCEDURE bubble_sort_index_int
END INTERFACE

INTERFACE quick_sort_index
   MODULE PROCEDURE quick_sort_index_real
   MODULE PROCEDURE quick_sort_index_int
END INTERFACE

INTERFACE quicksort_index
   MODULE PROCEDURE quicksort_index_real
   MODULE PROCEDURE quicksort_index_int
END INTERFACE

CONTAINS

!**************************************************************************
! convert string to lowercase 
!**************************************************************************
function lowercase (str)

  implicit none

  character(*) :: str
  character(len(str)) :: lowercase

  character(26), parameter :: char_up='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: char_low='abcdefghijklmnopqrstuvwxyz'

  integer :: i, pos

  lowercase = str

  do i = 1, len_trim(str)
     
     pos = index(char_up, str(i:i))
     if (pos /= 0) lowercase(i:i) = char_low(pos:pos)
     
  enddo

end function lowercase

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

!**************************************************************************
! 
!**************************************************************************
function str_found (nstr, str_array, str)
  
  implicit none
  
  logical      :: str_found
  integer      :: nstr
  character(*) :: str
  character(*) :: str_array(nstr)

  integer :: i

  str_found = .false.

  do i = 1, nstr
     if (str_array(i) == str) str_found = .true.
  enddo

end function str_found

!**************************************************************************
! 
!**************************************************************************
function str_pos (nstr, str_array, str)
  
  implicit none
  
  integer      :: str_pos
  integer      :: nstr
  character(*) :: str
  character(*) :: str_array(nstr)

  integer :: i

  str_pos = -1

  do i = 1, nstr
     if (str_array(i) == str) then
        str_pos = i
        exit
     endif
  enddo

end function str_pos

!**************************************************************************
! Delete control characters in string "str":
! All not printable ASCII characters (1-31 and 127) will be removed.
!**************************************************************************
function delete_control_chars (str)

  implicit none

  character(*)   :: str
  character(len(str)) :: delete_control_chars

  integer ::  i, k

  delete_control_chars = ""

  ! delete control characters (\0 etc.)
  k = 1
  do i = 1, len_trim(str)
     if ((IACHAR(str(i:i))>31) .AND. (IACHAR(str(i:i))<127)) then
        delete_control_chars(k:k) = str(i:i)
        k = k+1
     endif
  enddo 

end function delete_control_chars


!**************************************************************************
!* Prozedur:  isnumber                                                    *
!* Aufruf: call isnumber(str,valid,number)                                *
!* Funktion: Ueberprueft, ob auf einem String eine korrekte Integerzahl   *
!*           steht. Dabei duerfen sich vor der Zahl Leerzeichen und       *
!*           fuehrende Nullen befinden und hinter der Zahl, durch         *
!*           mindestens ein Komma oder Leerzeichen abgetrennt, beliebige  *
!*           andere Zeichen.                                              *
!* Parameter: str: (CHARACTER) Eingabeparameter, der den zu               *
!*                 ueberpruefenden String enthaelt.                       *  
!*            valid: (LOGICAL) Ist TRUE , falls der String eine gueltige  *
!*                   Integerzahl enthaelt.                                *
!*            number: (INTEGER) Enthaelt, die ermittelte Integerzahl      *
!* Autor: J. Ankenbrand                                                   *
!************************************************************************** 
SUBROUTINE isnumber(str,valid,number)
IMPLICIT NONE
character(*):: str
character(len(str))::hstr
integer::number,laenge,i,zahllaenge
logical::valid,vorzeichen

hstr=str           !damit der Eingabeparameter nicht veraendert wird
hstr=adjustl(hstr) !entfernen von Leerzeichen am Anfang des Strings

!Vorzeichen ist 'true', falls ein Minuszeichen vor der Zahl steht, sonst
!'false'. Steht ein '+' oder ein '-' vor der Zahl wird es entfernt  
vorzeichen=.false.
if ((hstr(1:1) == '+') .or. (hstr(1:1) == '-')) then 
  if (hstr(1:1) == '-') then
    vorzeichen=.true.
  endif
  hstr=hstr(2:)
endif


laenge=len_trim(hstr) ! Maximales Ende der zu ueberpruefenden Zahl
valid=.true.          ! Voreinstellung: es ist eine Zahl
zahllaenge=0          ! Laenge der Zahl
i=1
do while (i<=laenge)
  !Ueberpruefung ,ob das zu untersuchende Zeichen eine ziffer ist
  if ((hstr(i:i)=='0') .or. (hstr(i:i)=='1') .or. (hstr(i:i)=='2') .or. & 
  (hstr(i:i)=='3') .or. (hstr(i:i)=='4') .or. (hstr(i:i)=='5') .or. &
  (hstr(i:i)=='6') .or. (hstr(i:i)=='7') .or. (hstr(i:i)=='8') .or. &
  (hstr(i:i)=='9') ) then
    zahllaenge=i
    i=i+1
  else
    if((hstr(i:i) == ',') .or. (hstr(i:i) == ' ')) then
      i=laenge+1 !Sprung aus while-Schleife, da das Ende der Zahl erreicht ist
    else
      valid=.false.!Sprung aus while-Schleife, da es sich um keine Zahl mehr 
      i=laenge+1   !handeln kann
    endif
  endif
end do

if (zahllaenge == 0) then
  valid=.false.    !das erste Zeichen war also keine Ziffer
endif

!Es liegt eine Zahl vor:
if (valid) then
  hstr=hstr(1:zahllaenge)
  read(hstr,*) number   ! Schreiben des Strings auf eine Integer Variable
  if (vorzeichen) then  ! Falls negatives Vorzeichen 
    number=0-number
  endif
endif

end subroutine isnumber


!**************************************************************************
!* Prozedur:  n_numbers                                                   *
!* Aufruf: call n_numbers(str,n,valid,numbers                             *
!* Funktion: Ueberprueft, ob auf einem String n  korrekte Integerzahlen   *
!*           stehen. Dabei duerfen sich vor den Zahlen Leerzeichen        *
!*           befinden und hinter den Zahlen, durch mindestens ein Komma   *
!*           oder Leerzeichen abgetrennt , beliebige andere Zeichen.      *                 
!*           Die einzelnen Zahlen sind durch ein Komma und oder ein oder  *
!*           mehrere Leerzeichen voneinander getrennt                     *
!* Parameter: str: (CHARACTER) Eingabeparameter, der den zu               *
!*                 ueberpruefenden String  enthaelt.                      *  
!*            n: Anzahl der gewuenschten Integerzahlen aus dem String     *
!*            valid: (LOGICAL) Ist TRUE , falls der String eine gueltige  *
!*                   Integerzahl enthaelt.                                *
!*            number: (INTEGER,DIMENSION(n)) Enthaelt, die ermittelte     *
!*                    Integerzahl.                                        *
!* Autor: J. Ankenbrand                                                   *
!**************************************************************************

SUBROUTINE n_numbers(str,n,valid,numbers)
IMPLICIT NONE
character(*):: str
character(len(str)):: hstr
integer::n,i,break,a,b,zaehler,komma
integer,dimension(n):: numbers,hnumbers
logical::valid

hstr=str
hstr=adjustl(hstr)
valid=.true.
!Falls keine Zahl gewuenscht wird
if (n==0) then
  valid =.true.
else
  i=1
  do while ((i<=n) .and. valid)
    call isnumber(hstr,valid,hnumbers(i))  !Suchen nach EINER Integerzahl
    if (valid) then ! Falls kein Fehler auftrat
      zaehler=i !zum Zaehlen der gefundenen Integerzahlen
      a=index(hstr,' ') !suchen nach Komma und Leerzeichen 
      b=index(hstr,',')
      if (a==0 .and. b==0) then  !keine weiteren Integerzahlen vorhanden
        i=n+1                    !deshalb Abbruch
      else
        if (((a<b).and.(a.ne.0)).or.(b==0)) then !der kleinere Wert ungleich Null
          break=a                                !gilt
        else
          break=b
        endif
        
        !Suchen nach einem Komma und/oder mehreren Leerzeichen 
        komma=0    
        do while ((hstr(break:break)==' ').or.(hstr(break:break)==','))
           if (hstr(break:break)==',') then
             komma=komma+1
           endif
           break=break+1
        enddo
        if ((komma>=2).and.(i.ne.n)) then !es darf nur ein Komma 
          valid=.false.                 !zwischen zwei Zahlen stehen
        endif
        hstr=hstr(break:) ! Substring ab Zeichen nach letztem Komma oder Leerzeichen 
        i=i+1
      endif
    endif  
  enddo
  
  if ((zaehler==n) .and. valid) then ! es stehen n Zahlen auf dem String
    valid=.true.
    numbers=hnumbers
    write(*,*)numbers
  endif
  
endif

end subroutine n_numbers


!-----------------------------------------------------------------------    
! This subroutine sorts the elements of x by size.
!----------------------------------------------- N THOMAS 7 AUG 1996 ---    
SUBROUTINE bubble_sort (x,dim)

  use messy_clams_global, only: prec

  IMPLICIT NONE

  INTEGER    :: dim
  REAL(PREC) :: x(dim)
  
  REAL(PREC) :: help
  INTEGER    :: i, k
  LOGICAL    :: sorted

  sorted = .FALSE.
  k = 0

  DO WHILE (.NOT. sorted)
     sorted = .TRUE.
     k = k+1
     DO i = 1, dim-k
        IF (x(i) > x(i+1)) THEN
           help = x(i)
           x(i) = x(i+1)
           x(i+1) = help
           sorted = .FALSE.
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE bubble_sort

!-----------------------------------------------------------------------    
! This subroutine returns an index array  of the order of elements in 
! array x 
!----------------------------------------------- N THOMAS 7 AUG 1996 ---    
SUBROUTINE bubble_sort_index_real (x,ind,dim)

  use messy_clams_global, only: prec

  IMPLICIT NONE

  INTEGER,     intent(in)  :: dim
  REAL(PREC),  intent(in)  :: x(dim)
  INTEGER,     intent(out) :: ind(dim)
  
  INTEGER     :: help
  INTEGER     :: i, k
  LOGICAL     :: sorted

  ! initialize index array
  DO i=1,dim
     ind(i)=i
  END DO

  sorted = .FALSE.
  k = 0

  DO WHILE (.NOT. sorted)
     sorted = .TRUE.
     k = k+1
     DO i = 1, dim-k
        IF (x(ind(i)) > x(ind(i+1))) THEN
           help = ind(i)
           ind(i) = ind(i+1)
           ind(i+1) = help
           sorted = .FALSE.
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE bubble_sort_index_real

SUBROUTINE bubble_sort_index_int (x,ind,dim)

  IMPLICIT NONE

  INTEGER,  intent(in)  :: dim
  INTEGER,  intent(in)  :: x(dim)
  INTEGER,  intent(out) :: ind(dim)
  
  INTEGER     :: help
  INTEGER     :: i, k
  LOGICAL     :: sorted

  ! initialize index array
  DO i=1,dim
     ind(i)=i
  END DO

  sorted = .FALSE.
  k = 0

  DO WHILE (.NOT. sorted)
     sorted = .TRUE.
     k = k+1
     DO i = 1, dim-k
        IF (x(ind(i)) > x(ind(i+1))) THEN
           help = ind(i)
           ind(i) = ind(i+1)
           ind(i+1) = help
           sorted = .FALSE.
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE bubble_sort_index_int

!-----------------------------------------------------------------------    
! This subroutine sorts the elements of x by size with 
! a quick sort algorithm
!-----------------------------------------------------------------------    
RECURSIVE SUBROUTINE quick_sort(array)

  use messy_clams_global, only: prec

  IMPLICIT NONE

  real(prec), dimension(:), intent(inout) :: array
  integer                                 :: x, low, high
  
  low = 1
  high = size(array)
  if(low < high) then
     x = partition(low,high)
     call quick_sort(array(low:x-1))
     call quick_sort(array(x+1:high))
  endif
  
CONTAINS
      
    FUNCTION partition(low,high)

      integer :: partition
      integer :: low, high
      
      integer :: i, j
      real    :: pivot, temp
      
      i = low
      j = high-1
      pivot = array(high)
      
      do
         do while (array(i)<=pivot .and. i<high)
            i = i+1
         end do
         do while (array(j)>=pivot .and. j>low)
            j = j-1
         end do
         if (i<j) then
            temp     = array(i)
            array(i) = array(j)
            array(j) = temp
         endif
         if (i>=j) exit
      enddo
      
      temp        = array(i)
      array(i)    = array(high)
      array(high) = temp
      
      partition = i
      
    END FUNCTION partition
    
END SUBROUTINE quick_sort

!-----------------------------------------------------------------------    
! This subroutine returns an index array  of the order of elements in 
! array x 
! (quicksort is used)
!-----------------------------------------------------------------------    
SUBROUTINE quick_sort_index_real(array,ind)

  use messy_clams_global, only: prec

  implicit none 
  real(prec), dimension(:), intent(in)  :: array
  integer,    dimension(:), intent(out) :: ind
  
  real(prec), dimension(:), allocatable :: helparr
  integer :: i

  allocate (helparr(size(array)))
  helparr = array
  
  ! initialize index array
  DO i=1,size(ind)
     ind(i)=i
  END DO

  call quicksort_index (helparr, ind) 
   
END SUBROUTINE quick_sort_index_real

SUBROUTINE quick_sort_index_int(array,ind)

  implicit none 
  integer, dimension(:), intent(in)  :: array
  integer, dimension(:), intent(out) :: ind
  
  integer, dimension(:), allocatable :: helparr
  integer :: i
  
  allocate (helparr(size(array)))
  helparr = array
  
  ! initialize index array
  DO i=1,size(ind)
     ind(i)=i
  END DO

  call quicksort_index (helparr, ind) 
   
END SUBROUTINE quick_sort_index_int



RECURSIVE SUBROUTINE quicksort_index_real(array,ind)

  use messy_clams_global, only: prec

  IMPLICIT NONE
  
  real(prec), dimension(:), intent(inout) :: array
  integer,    dimension(:), intent(inout) :: ind
  integer                                 :: x, low, high
  
  integer :: itemp
  
  low = 1
  high = size(array)
  if(low < high) then
     x = partition(low,high)
     call quicksort_index_real(array(low:x-1),ind(low:x-1))
     call quicksort_index_real(array(x+1:high),ind(x+1:high))
  endif
  
  CONTAINS
  
    FUNCTION partition(low,high)
    
      integer :: partition
      integer :: low, high
      
      integer   :: i, j
      real(prec):: pivot, temp
      
      i = low
      j = high-1
      pivot = array(high)
      
      do
         do while (array(i)<=pivot .and. i<high)
            i = i+1
         end do
         do while (array(j)>=pivot .and. j>low)
            j = j-1
         end do
         if (i<j) then
            temp     = array(i)
            array(i) = array(j)
            array(j) = temp
            itemp    = ind(i)
            ind(i)   = ind(j)
            ind(j)   = itemp
         endif
         if (i>=j) exit
      enddo
      
      temp        = array(i)
      array(i)    = array(high)
      array(high) = temp
      itemp       = ind(i)
      ind(i)      = ind(high)
      ind(high)   = itemp
       
      partition = i
      
    END FUNCTION partition

END SUBROUTINE quicksort_index_real

RECURSIVE SUBROUTINE quicksort_index_int(array,ind)

  IMPLICIT NONE
  
  integer, dimension(:), intent(inout) :: array
  integer, dimension(:), intent(inout) :: ind
  integer                              :: x, low, high
  
  integer :: itemp
  
  low = 1
  high = size(array)
  if(low < high) then
     x = partition(low,high)
     call quicksort_index_int(array(low:x-1),ind(low:x-1))
     call quicksort_index_int(array(x+1:high),ind(x+1:high))
  endif
  
  CONTAINS
  
    FUNCTION partition(low,high)
    
      integer :: partition
      integer :: low, high
      
      integer :: i, j
      integer :: pivot, temp
      
      i = low
      j = high-1
      pivot = array(high)
      
      do
         do while (array(i)<=pivot .and. i<high)
            i = i+1
         end do
         do while (array(j)>=pivot .and. j>low)
            j = j-1
         end do
         if (i<j) then
            temp     = array(i)
            array(i) = array(j)
            array(j) = temp
            itemp    = ind(i)
            ind(i)   = ind(j)
            ind(j)   = itemp
         endif
         if (i>=j) exit
      enddo
      
      temp        = array(i)
      array(i)    = array(high)
      array(high) = temp
      itemp       = ind(i)
      ind(i)      = ind(high)
      ind(high)   = itemp
       
      partition = i
      
    END FUNCTION partition

END SUBROUTINE quicksort_index_int

  !*********************************************************************
  !
  !*********************************************************************                  
  function get_metfilename (prefix, dirname, iyear, imon, iday, ihr)

    use messy_clams_global,  only: filenamelen
   
    implicit none

    character(filenamelen) :: get_metfilename
    character(*)           :: prefix, dirname
    integer                :: iyear, imon, iday, ihr

    character(filenamelen) :: dir
    character(10)  :: subdir
    integer        :: pos

    
    write(get_metfilename,'(A,4I2.2,A)') TRIM(prefix)//'_',MOD(iyear,100), &
         imon,iday,ihr,'.nc'
    
    IF (dirname /= ' ')  THEN
       
       dir = ADJUSTL(dirname)
       if (dir(LEN_TRIM(dir):LEN_TRIM(dir))=='/')  dir = dir(1:LEN_TRIM(dir)-1)
       
       if (index (uppercase(dir),"YYYY/MM",back=.true.)>0) then
          write(subdir,'(I4.4,A,I2.2)') iyear,'/',imon
          pos = index (uppercase(dir),"YYYY/MM",back=.true.)
          dir = trim(dir(1:pos-1))//trim(subdir)//trim(dir(pos+7:len_trim(dir)))
       elseif (index (uppercase(dir),"YYYYMM",back=.true.)>0) then
          write(subdir,'(I4.4,I2.2)') iyear,imon 
          pos = index (uppercase(dir),"YYYYMM",back=.true.)
          dir = trim(dir(1:pos-1))//trim(subdir)//trim(dir(pos+6:len_trim(dir)))
       elseif (index (uppercase(dir),"YYMM",back=.true.)>0) then
          write(subdir,'(2I2.2)') MOD(iyear,100),imon 
          pos = index (uppercase(dir),"YYMM",back=.true.)
          dir = trim(dir(1:pos-1))//trim(subdir)//trim(dir(pos+4:len_trim(dir)))
       endif
       get_metfilename = TRIM(dir)//'/'//TRIM(ADJUSTL(get_metfilename))
       
    ENDIF

  end function get_metfilename

  !*********************************************************************
  !
  !*********************************************************************                  
  subroutine pack_values (array, mask, mdi)

    use messy_clams_global, only: prec

    implicit none

    real(prec), dimension(:), pointer :: array
    logical,    dimension(:), pointer :: mask
    real(prec)                        :: mdi
    
    real(prec), dimension(:), pointer :: helparr
    
    ! allocate help array with size of valid values
    allocate (helparr(count(mask)))

    ! write valid values to helparr
    helparr = pack(array,mask)  

    ! write valid values to array without gaps and fill with missing value
    array(1:count(mask)) = helparr
    array(count(mask)+1:) = mdi

    ! deallocate help arrray
    deallocate (helparr)
    
  end subroutine pack_values

END MODULE messy_clams_tools_utils
