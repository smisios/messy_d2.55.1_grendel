/*
  f90_unix_env_ccode : c stubs for f90_unix_env

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
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/utsname.h>

int c_clk_tck_(){
    /* return CLK_TCK; */
    return _SC_CLK_TCK;
}

void c_ctermid_(char *s,  int ls){
    ctermid(s);
}

gid_t c_getegid_(){
    return getegid();
}

gid_t c_getgid_(){
    return getgid();
}

uid_t c_geteuid_(){
    return geteuid();
}

void c_getgroups_(gid_t * glist, int * gmax, int * ngroups, int * err){
    
    * ngroups = getgroups(*gmax, glist);   
    * err = errno;
}

void c_gethostname_(char * name, int maxlen){
    gethostname(name, (size_t) maxlen);
}

void c_getlogin_(char * name, int maxlen){
    getlogin_r(name, (size_t) maxlen);
    /* kludge..
       getlogin and getlogin_err raise error 13-permission denied
       on my system. Why, oh why?
     */
    errno = 0;
}

pid_t c_getpgrp_(){
    return getpgrp();
}

pid_t c_getpid_(){
    return getpid();
}

pid_t c_getppid_(){
    return getppid();
}

uid_t c_getuid_(){
    return getuid();
}

void c_setgid_(gid_t *gid){
    setgid(*gid);
}

void c_setpgid_(gid_t *gid, gid_t *pgid){
    setpgid(*gid, *pgid);
}

void c_setsid_(){
    setsid();
}

void c_setuid_(uid_t *uid){
    setuid(*uid);
}

void c_sysconf_(int *name, int *val, int *err){
    *val = sysconf(*name);
    *err = errno;
}

time_t c_time_(){
    return time(NULL);
}

clock_t c_times_(clock_t *utime, clock_t *stime,
		 clock_t *cutime, clock_t *cstime){
    struct tms buf;
    return times(&buf);
    *utime = buf.tms_utime ; *stime = buf.tms_stime;
    *cutime = buf.tms_cutime ; *cstime = buf.tms_cstime;
}

void c_uname_(char *sysname, char *nodename, char *release, char *version, 
	      char *machine, int sl, int nl, int rl, int vl, int ml)
{
    struct utsname buf;

    uname(&buf);
    strncpy(sysname, buf.sysname, (size_t) sl);
    strncpy(nodename, buf.nodename, (size_t) nl);
    strncpy(release, buf.release, (size_t) rl);
    strncpy(version, buf.version, (size_t) vl);
    strncpy(machine, buf.machine, (size_t) ml);
}
