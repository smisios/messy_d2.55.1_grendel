/*
  f90_unix_proc_const.c : Determine values of constants for f90_unix_proc

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
#include <sys/types.h>
#include <sys/wait.h>

main(){

    printf("! DO NOT EDIT THIS FILE. \n");
    printf("! It is created from f90_unix_proc_const.c\n");
    printf("integer, parameter :: pid_kind = %d\n", sizeof(pid_t));
    printf("integer, parameter :: WNOHANG = %d\n", WNOHANG);
    printf("integer, parameter :: WUNTRACED = %d\n", WUNTRACED);
    printf("integer, parameter :: WCONTINUED = %d\n", WCONTINUED);
    return 0;
}


