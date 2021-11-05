!******************************************************************************
! File utils/src/dateconv.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! N. Thomas, M. Fisher, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Tue Jul 26 12:40:37 2011
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
! Module dateconv
!
! This module contains subroutines for working with dates and time:
!
! subroutine js2ymds_interface (js, y, m, d, s, dates30)
! function ymds2js_interface (y, m, d, s, dates30) result(erg)
! subroutine incdat_interface (jsec,jday,jmon,jyer,isec,iday,imon,iyer,dates30)          
! Function ymd2jd (y, m, d) RESULT(erg)
! Subroutine jd2ymd (jd, y, m, d)
! Function js2jd (js)
! Subroutine js2ymds (js, y, m, d, s)
! Function ymds2js (y, m, d, s) RESULT(erg)
! Logical function leapyear(year) 
! Integer function  yrday (day, month, year)
! Subroutine yd2md (year, yd, month, day)
! Subroutine check_date (hour, day, month, year, ok)
! Subroutine incdat (JSEC,JDAY,JMON,JYER,ISEC,IDAY,IMON,IYER)   
! Subroutine nnorm (JSEC,JDAY,JMON,JYER)                                
! Function d1gtd2 (JSEC,JDAY,JMON,JYER,ISEC,IDAY,IMON,IYER) RESULT(erg) 
! Function datsec (ISEC,IDAY,IMON,IYEAR) RESULT(erg)
! subroutine js30_to_ymds30 (js, y, m, d, s)
! function ymds30_to_js30 (y, m, d, s) result(erg)
! subroutine check_date30 (hour, day, month, year, ok)
! subroutine incdat30 (jsec,jday,jmon,jyer,isec,iday,imon,iyer)          
! subroutine nnorm30 (JSEC,JDAY,JMON,JYER)                                
!
!******************************************************************************
 
MODULE messy_clams_tools_dateconv

CONTAINS

!-------------------------------------------------------------
! Conversion from julian seconds to date (year, month, day, sec)
! => call js2ymds or js30_to_ymds30 (if dates30=.true.)
!-------------------------------------------------------------
subroutine js2ymds_interface (js, y, m, d, s, dates30)
  
  implicit none

  real(kind(0d0)) :: js
  integer :: y, m, d, s
  logical :: dates30

  if (dates30) then
     call js30_to_ymds30 (js,y,m,d,s)
  else
     call js2ymds (js,y,m,d,s)
  endif

end subroutine js2ymds_interface


!-------------------------------------------------------------
! Conversion from date (year, month, day, seconds) to
! julian seconds
! => call ymds2js or ymds30_to_js30 (if dates30=.true.)
!-------------------------------------------------------------
function ymds2js_interface (y, m, d, s, dates30) result(erg)

  implicit none
 
  integer :: y, m, d, s
  logical :: dates30
  real(kind(0d0)) :: erg

  if (dates30) then
     erg = ymds30_to_js30 (y, m, d, s)
  else
     erg = ymds2js (y, m, d, s)
  endif

end function ymds2js_interface

!-----------------------------------------------------------------------    
! This subroutine increments the date jsec,jday,jmon,jyer by              
! isec seconds iday days imon months and iyer years   
! => call incdat or incdat30 (if dates30=.true.)                      
!-----------------------------------------------------------------------    
subroutine incdat_interface (jsec,jday,jmon,jyer,isec,iday,imon,iyer,dates30)          

  implicit none

  integer :: jsec,jday,jmon,jyer,isec,iday,imon,iyer
  logical :: dates30

  if (dates30) then
     call incdat30 (jsec,jday,jmon,jyer,isec,iday,imon,iyer)          
  else
     call incdat (jsec,jday,jmon,jyer,isec,iday,imon,iyer)          
  endif

end subroutine incdat_interface


!-------------------------------------------------------------
! Conversion of date to julian day 
!
! input parameters:  
!   y : INTEGER; year 
!   m : INTEGER; month
!   d : INTEGER; day
!
! output parameters:
!   ymd2jd : INTEGER; julian day
!
! Reference: Fliegel, H.F. and van Flandern, T.C. (1968)
!-------------------------------------------------------------
function ymd2jd (y, m, d) result(erg)

  implicit none
  integer :: y, m, d ,erg

  erg = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +        &
        ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -  &
        ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 + &
        d - 32075

end function ymd2jd

!-------------------------------------------------------------
! Conversion from julian day to date
!
! input parameter:
!    jd : INTEGER; julian day
! 
! output parameters:
!    y : INTEGER; year
!    m : INTEGER; month
!    d : INTEGER; day
!
! Reference: Fliegel, H.F. and van Flandern, T.C. (1968)
!-------------------------------------------------------------
subroutine jd2ymd (jd, y, m, d)

  implicit none
  integer :: jd, y, m, d
  integer :: l, n, i, j
    
  l = jd + 68569
  n = ( 4 * l ) / 146097
  l = l - ( 146097 * n + 3 ) / 4
  i = ( 4000 * ( l + 1 ) ) / 1461001
  l = l - ( 1461 * i ) / 4 + 31
  j = ( 80 * l ) / 2447
  d = l - ( 2447 * j ) / 80
  l = j / 11
  m = j + 2 - ( 12 * l )
  y = 100 * ( n - 49 ) + i + l

end subroutine jd2ymd

!-------------------------------------------------------------
! Conversion from julian seconds to julian day
!
! input parameter:
!   js : DOUBLE; julian seconds
!
! output parameters:
!   jd : DOUBLE; julian day
! 
!-------------------------------------------------------------
function js2jd (js)

  implicit none

  real(kind(0d0)) :: js2jd, js  
 
  js2jd = ymd2jd(2000,1,1) - 0.5d0 + js/86400d0
 
end function js2jd

!-------------------------------------------------------------
! Conversion from julian seconds to date (year, month, day, sec)
!
! input parameter:
!    js : DOUBLE; julian seconds
!
! output parameters:
!    y : INTEGER; year
!    m : INTEGER; month
!    d : INTEGER; day
!    s : INTEGER; seconds of day 
!
!-------------------------------------------------------------
subroutine js2ymds (js, y, m, d, s)

  implicit none

  real(kind(0d0)) :: js
  integer :: y, m, d, s

  integer, parameter :: jd2000=2451545 ! 1.1.2000
  integer :: days

  days = floor(js/86400)
  s = js - 86400d0*days
  call jd2ymd (days+jd2000, y, m, d)
        
end subroutine js2ymds

!-------------------------------------------------------------
! Conversion from date (year, month, day, seconds) to
! julian seconds
!
! input parameters: 
!   y : INTEGER; year
!   m : INTEGER; month
!   d : INTEGER; day
!   s : INTEGER; seconds of day
!
! output parameter:
!   ymds2js : DOUBLE; julian seconds
!
!-------------------------------------------------------------
function ymds2js (y, m, d, s) result(erg)

  implicit none
 
  integer :: y, m, d, s
  real(kind(0d0)) :: erg

  integer, parameter :: jd2000=2451545 ! 1.1.2000

  erg = s + (ymd2jd(y,m,d)-jd2000)*86400d0
  
end function ymds2js

!---------------------------------------------------------------    
! This function returns the value 'TRUE' if year is a leap-year.
!---------------------------------------------------------------    
logical function leapyear(year) 

  implicit none
  integer  :: year

  leapyear = (MOD(year,4)==0 .AND. MOD(year,100)/=0)  &
       .OR. (MOD(year,400)==0) 

end function leapyear

!---------------------------------------------------------------    
! Convert from date to day of year
!---------------------------------------------------------------    
integer function  yrday (day, month, year)

  implicit none
     
  integer :: day, month, year
  integer :: days(12)
 
  data days / 0,31,59,90,120,151,181,212,243,273,304,334 /
  
  yrday = days(month) + day
  if (leapyear(year) .and.  &
       (month>2 .or. (month==2 .and. day==29))) then
     yrday = yrday + 1
  endif
  
end function yrday

!---------------------------------------------------------------    
! Convert from year and day of year to month and day of month.
!--------------------------------------------------------------- 
subroutine yd2md (year, yd, month, day)

  implicit none

  integer :: year, yd, month, day

  integer :: i
  integer :: ydays(13)
 
  data ydays / 0,31,59,90,120,151,181,212,243,273,304,334,365 /

  if (leapyear(year)) ydays = (/0,31,60,91,121,152,182,213,244,274,305,335,366 /)

  if (yd < 1 .or. yd > ydays(13)) then
     month = -1
     day = -1
     return
  endif

  month = 1
  do while (month <=12 .and. yd > ydays(month+1))
     month = month + 1
  enddo

  day = yd - ydays(month)

end subroutine yd2md

!---------------------------------------------------------------    
! This subroutine checks if the date is correct.
!---------------------------------------------------------------    
subroutine check_date (hour, day, month, year, ok)

  IMPLICIT NONE

  INTEGER :: hour, day, month, year
  LOGICAL :: ok

  INTEGER :: days(12)

  ok = .false.

  days(1) = 31
  days(2) = 29
  days(3) = 31
  days(4) = 30
  days(5) = 31
  days(6) = 30
  days(7) = 31
  days(8) = 31
  days(9) = 30
  days(10) = 31
  days(11) = 30
  days(12) = 31
   
  if (year < 1960) then
     write (*,*) 'Year is incorrect !!!'
     write (*,*) '    year >= 1960'
  elseif ((month<0) .or. (month>12)) then
     write (*,*) 'Month is incorrect !!!'
  elseif (day<0) then
     write (*,*) 'Day is incorrect !!!'
  elseif ((month==2) .and. &
       (((leapyear(year)) .and. (day>29)) .or.  &
       ((.not. leapyear(year)) .and. (day>28)))) then
     write (*,*) 'leapyear=', leapyear(year)
     write (*,*) 'Day is incorrect !!!'
  elseif (day>days(month)) then  
     write (*,*) 'Day is incorrect !!!'
  elseif ((hour<0) .or. (hour>=24)) then
     write (*,*) 'Hour is incorrect !!!'
  else
     ok = .true.
  endif
  
end subroutine check_date

!-----------------------------------------------------------------------    
!   This subroutine increments the date jsec,jday,jmon,jyer by              
!   isec seconds iday days imon months and iyer years                       
!----------------------------------------- AUTHOR: M FISHER 1984--------    
subroutine incdat (jsec,jday,jmon,jyer,isec,iday,imon,iyer)          

  implicit none

  integer :: jsec,jday,jmon,jyer,isec,iday,imon,iyer
   
   jyer = jyer + iyer                                                    
   jmon = jmon + imon                                                    
   call nnorm (jsec,jday,jmon,jyer)                                      
   jday = jday + iday                                                    
   call nnorm (jsec,jday,jmon,jyer)                                      
   jsec = jsec + isec                                                    
   call nnorm (jsec,jday,jmon,jyer)                                      
   
end subroutine incdat

!-----------------------------------------------------------------------    
!   This subroutine normalises the date jsec,jday,jmon,jyer so that         
!   -1<jsec<86400 0<jday<1month and 0<jmon<13                               
!----------------------------------------- AUTHOR: M FISHER 1984--------    
subroutine nnorm (JSEC,JDAY,JMON,JYER)                                

   implicit none   
   integer :: jsec,jday,jmon,jyer
   integer :: lenmon(12)                                                  
   
   data lenmon /31,0,31,30,31,30,31,31,30,31,30,31/                      
   
   jday = jday + jsec/86400                                              
   jsec = mod (jsec,86400)                                               
   
   if (jsec < 0) then                                                 
      jsec = jsec + 86400                                             
      jday = jday - 1                                                 
   endif
   
   jyer = jyer + (jmon-1)/12                                             
   jmon = 1 + mod (jmon-1,12)                                            
   
   if (jmon <= 0) then                                                 
      jmon = jmon + 12                                                
      jyer = jyer - 1                                                 
   endif
   
   lenmon(2) = 28                                                        
   if (leapyear(jyer)) lenmon(2)=29

   do while ((jday > lenmon(jmon)) .or. (jday <= 0))
   
      if (jday > lenmon(jmon)) then                                      
         jday = jday - lenmon(jmon)                                      
         jmon = jmon + 1                                                 
         if (jmon > 12) then                                          
            jmon = 1                                                  
            jyer = jyer + 1                                           
         end if
      else                              
         jmon = jmon - 1                                                 
         if (jmon <= 0) then                                           
            jmon = 12                                                 

            jyer = jyer - 1                                           
         end if
         jday = jday + lenmon(jmon)                                      
      endif     
                                                            
      lenmon(2) = 28                                                        
      if (leapyear(jyer)) lenmon(2)=29
   
   enddo
   
end subroutine nnorm
   
!-----------------------------------------------------------------------    
!   This function checks the date 'jsec',...,'jyer' against the             
!   date 'isec',...,'iyer' and returns with the value '.true.' if        
!   'j' > 'i'                                                               
!-----------------------------------------------------------------------    
function d1gtd2 (jsec,jday,jmon,jyer,isec,iday,imon,iyer) result(erg) 

   implicit none
   integer :: jsec,jday,jmon,jyer,isec,iday,imon,iyer
   logical :: erg
   
   erg = .false.                                                      
   
   if (jyer > iyer) then                                         
       erg = .true.                                                 
   else if (jyer == iyer) then                                         
      if (jmon > imon) then                                   
          erg = .true.                                           
      else if (jmon == imon) then                                   
          if (jday > iday) then                             
              erg = .true.                                     
          else if (jday == iday) then                             
              if (jsec > isec) then                       
                  erg = .true.                               
              end if                                              
          end if                                                    
      end if                                                          
   end if     
  
end function d1gtd2
   

!-----------------------------------------------------------------------    
!   This function returns the date in seconds from 0,1,1,1960               
!-----------------------------------------------------------------------    
function datsec (isec,iday,imon,iyear) result(erg)
   
   implicit none
   integer :: isec,iday,imon,iyear
   integer :: mnday(12)                                            
   double precision :: erg      
   
   data mnday /0,31,52,90,120,151,181,212,243,273,304,334/               
   
   erg = (365*(iyear-1960)+(iyear-1960)/4   &                          
        &          +mnday(imon)+iday-1)*86400d0 +isec                          
   
   ! if ((mod(iyear,4) == 0).and.(imon >= 2)) erg=erg+86400d0    
   
   ! Hector Maclean Sept 1994 -- correction
   if ((mod(iyear,4) == 0).and.(imon < 2)) erg=erg-86400d0    
   
end function datsec                                                            

!-------------------------------------------------------------
! Conversion from "julian seconds 30" to "date 30"
!
! input parameter:
!    js : DOUBLE; julian seconds 30 (all months with 30 days!)
!
! output parameters:
!    y : INTEGER; year
!    m : INTEGER; month
!    d : INTEGER; day
!    s : INTEGER; seconds of day 
!
!-------------------------------------------------------------
subroutine js30_to_ymds30 (js, y, m, d, s)

  implicit none

  real(kind(0d0)) :: js
  integer :: y, m, d, s

  integer :: days

  days = floor(js/86400)
  s = js - 86400d0*days

  if (js<0)  days = days + 1
  y = days/360 
  days = days - y*360
  m = days/30  
  days = days - m*30
  if (js<0) then
     y = y + 1999
     m = 12+m 
     d = 30+days 
  else
     y = y + 2000
     m = m+1 
     d = days+1 
  endif
        
end subroutine js30_to_ymds30

!-------------------------------------------------------------
! Conversion from "date 30" (year, month, day, seconds) to
! "julian seconds 30"
!
! input parameters: 
!   y : INTEGER; year
!   m : INTEGER; month
!   d : INTEGER; day
!   s : INTEGER; seconds of day
!
! output parameter:
!   ymds2js : DOUBLE; julian seconds 30 (all months with 30 days!)
!
!-------------------------------------------------------------
function ymds30_to_js30 (y, m, d, s) result(erg)

  implicit none
 
  integer :: y, m, d, s
  real(kind(0d0)) :: erg

  if (y<2000) then
     erg = (y-1999)*360*86400+(m-12)*30*86400+(d-30)*86400+(s-86400)
  else
     erg = (y-2000)*360*86400+(m-1)*30*86400+(d-1)*86400+s
  endif
  
end function ymds30_to_js30

!---------------------------------------------------------------    
! This subroutine checks if the date (all months with 30 days)
! is correct.
!---------------------------------------------------------------    
subroutine check_date30 (hour, day, month, year, ok)

  IMPLICIT NONE

  INTEGER :: hour, day, month, year
  LOGICAL :: ok

  ok = .false.

  if (year < 1960) then
     write (*,*) 'Year is incorrect !!!'
     write (*,*) '    year >= 1960'
  elseif (month<0 .or. month>12) then
     write (*,*) 'Month is incorrect !!!'
  elseif (day<0 .or. day>30) then
     write (*,*) 'Day is incorrect !!!'
  elseif (hour<0 .or. hour>=24) then
     write (*,*) 'Hour is incorrect !!!'
  else
     ok = .true.
  endif
  
end subroutine check_date30

!-----------------------------------------------------------------------    
!   This subroutine increments the date jsec,jday,jmon,jyer by              
!   isec seconds iday days imon months and iyer years   
!   All months with 30 days !!!                    
!-----------------------------------------------------------------------    
subroutine incdat30 (jsec,jday,jmon,jyer,isec,iday,imon,iyer)          

  implicit none

  integer :: jsec,jday,jmon,jyer,isec,iday,imon,iyer
   
   jyer = jyer + iyer                                                    
   jmon = jmon + imon                                                    
   call nnorm30 (jsec,jday,jmon,jyer)                                      
   jday = jday + iday                                                    
   call nnorm30 (jsec,jday,jmon,jyer)                                      
   jsec = jsec + isec                                                    
   call nnorm30 (jsec,jday,jmon,jyer)                                      
   
end subroutine incdat30

!-----------------------------------------------------------------------    
!   This subroutine normalises the date jsec,jday,jmon,jyer so that         
!   -1<jsec<86400 0<jday<1month and 0<jmon<13                   
!   All months with 30 days !!!            
!-----------------------------------------------------------------------    
subroutine nnorm30 (JSEC,JDAY,JMON,JYER)                                

   implicit none   
   integer :: jsec,jday,jmon,jyer
   
   jday = jday + jsec/86400                                              
   jsec = mod (jsec,86400)                                               
   
   if (jsec < 0) then                                                 
      jsec = jsec + 86400                                             
      jday = jday - 1                                                 
   endif
   
   jyer = jyer + (jmon-1)/12                                             
   jmon = 1 + mod (jmon-1,12)                                            
   
   if (jmon <= 0) then                                                 
      jmon = jmon + 12                                                
      jyer = jyer - 1                                                 
   endif
   
   do while ((jday > 30) .or. (jday <= 0))
   
      if (jday > 30) then                                      
         jday = jday - 30                                      
         jmon = jmon + 1                                                 
         if (jmon > 12) then                                          
            jmon = 1                                                  
            jyer = jyer + 1                                           
         end if
      else                              
         jmon = jmon - 1                                                 
         if (jmon <= 0) then                                           
            jmon = 12                                                 

            jyer = jyer - 1                                           
         end if
         jday = jday + 30                                      
      endif     
                                                            
   enddo
   
end subroutine nnorm30
   
end module messy_clams_tools_dateconv
