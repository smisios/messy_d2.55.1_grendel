/*
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

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

void c_chdir_(char *path, int pathlen){
    chdir(path);
}

void c_getcwd_(char *path, int pathlen){
    getcwd(path, (size_t) pathlen);
}

void c_link_(char *existing, char *new, int elen, int nlen){
    link(existing, new);
}

void c_mkdir_(char *path, mode_t *mode, int plen){
    mkdir(path, *mode);
}

void c_mkfifo_(char *path, mode_t *mode, int plen){
    mkfifo(path, *mode);
}

void c_rmdir_(char *path, int plen){
    rmdir(path);
}

void c_unlink_(char *path, int plen){
    unlink(path);
}

