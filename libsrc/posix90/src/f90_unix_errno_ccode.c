/*
  f90_unix_errno_ccode.c : c stubs for f90_unix_errno
  
   This file is part of Posix90

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA. 

*/

#include <stdio.h>
#include <errno.h>
#include <string.h>

void c_strerror_(char * str, int * err, int *errc, int lstr){
    char * estr;
    int errno_save;

    errno_save = errno ; errno = 0;
    estr = strerror(*err);
    strncpy(str, estr, (size_t) lstr);
    *errc = errno; errno = errno_save;
}

void c_perror_(char * str, int * errc){

    errno = *errc;
    perror(str);
}

void c_set_errno_(int * errc){
    
    errno = *errc;
}

int c_get_errno_(int * errc){
    
    return errno;
}
