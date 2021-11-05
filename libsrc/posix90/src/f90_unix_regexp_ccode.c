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
#include <sys/stat.h>
#include <sys/types.h>
#include <regex.h>

void c_regcomp_(regex_t **preg, char *regex, int *cflags, int * errc, int lregex){
    if((*preg = (regex_t *) malloc(sizeof(regex_t)))==NULL){
	fprintf(stderr, "Memory allocation failed in posix90::c_regcomp_\n");
    }
    *errc = regcomp(*preg, regex, *cflags);
}

void c_regexec_(regex_t  **preg, char *string, int *nmatch, 
		regmatch_t *pmatch, int *eflags, int * errc, int lstring){
    *errc = regexec(*preg, string, (size_t) *nmatch, pmatch, *eflags);
}

void c_regerror_(int *errcode, regex_t  **preg, char *errbuf, int lerrbuf){
    regerror(*errcode, *preg, errbuf, (size_t) lerrbuf);
}

void c_regfree_(regex_t **preg){
    regfree(*preg);
    free(*preg);
}
