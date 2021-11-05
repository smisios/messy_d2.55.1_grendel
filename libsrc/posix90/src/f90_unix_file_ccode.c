/*
  f90_unix_file_ccode: c stub code for f90_unix_file

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

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

void c_access_(char *path, int *mode, int lpath){
    access(path, *mode);
}

void c_chmod_(char *path, mode_t *mode, int lpath){
    chmod(path, *mode);
}

void c_chown_(char *path, uid_t *owner, gid_t *group, int lpath){
    chown(path, *owner, *group);
}

void c_stat_(char *path, dev_t *st_dev, ino_t *st_ino, mode_t *st_mode, 
	     nlink_t *st_nlink, uid_t *st_uid, gid_t *st_gid, dev_t *st_rdev,
	     off_t *st_size, time_t *st_atimes, time_t *st_mtimes,
	     time_t *st_ctimes, int lpath){
    struct stat buf;
    stat(path, &buf);
    *st_dev = buf.st_dev;
    *st_ino = buf.st_ino;
    *st_mode = buf.st_mode;
    *st_nlink = buf.st_nlink;
    *st_uid = buf.st_uid;
    *st_gid = buf.st_gid;
    *st_rdev = buf.st_rdev;
    *st_size = buf.st_size;
    *st_atimes = buf.st_atime;
    *st_mtimes = buf.st_mtime;
    *st_ctimes = buf.st_ctime;
}
