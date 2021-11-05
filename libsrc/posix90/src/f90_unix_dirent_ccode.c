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

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>

void c_opendir_(char *dirname, DIR **dirp, int dirnamelen){
    *dirp = opendir(dirname);
}

void c_readdir_(DIR **dirp, char *name, int *namelen, int namelenmax){
    struct dirent *ent;

    ent = readdir(*dirp);
    if(ent) {
	strncpy(name, ent->d_name, (size_t) namelenmax);
	*namelen = strlen(ent->d_name);
    } else {
	*name=0 ; *namelen = -1;
    }
}

void c_rewinddir_(DIR **dirp){
    rewinddir(*dirp);
}

void c_closedir_(DIR **dirp){
    closedir(*dirp);
}


