/*
  f90_unix_io_ccode.c : c stub code for f90_unix_io

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

int c_feof_(FILE **stream){
    return feof(*stream);
}

void c_fseek_(FILE **stream, long *offset, int *whence){
    fseek(*stream, *offset, *whence);
}

long c_ftell_(FILE **stream){
    return ftell(*stream);
}

void c_fopen_(FILE **fp, char *path, char *mode){
    *fp = fopen(path, mode);
}

void c_fclose_(FILE **fp){
    fclose(*fp);
}

void c_popen_(FILE **fp, char *cmd, char *mode){
    *fp = popen(cmd, mode);
}

void c_pclose_(FILE **fp){
    pclose(*fp);
}

int c_fread_str_(char *str, int *size, FILE **fp, int maxsize){
    return fread(str, (size_t) *size, (size_t) 1, *fp);
}

int c_fread_str_array_(char *str, int *size, int *nelem, FILE **fp, int maxsize){
    if(*size == maxsize) {
	return fread(str, (size_t) *size, (size_t) *nelem, *fp);
    } else {
	int i,j=0,k;
	for(i=0; i<*nelem; i++){
	    k = fread(str+i*maxsize, (size_t) *size, (size_t) 1, *fp);
	    if(k==0) break;
	    j+=k;
	}
	return j;
    }
}

int c_fread_real_(float *f, int *nf, FILE **fp){
    return fread(f, (size_t)sizeof(float), (size_t)*nf, *fp);
}

int c_fread_double_(double *f, int *nf, FILE **fp){
    return fread(f, (size_t)sizeof(double), (size_t)*nf, *fp);
}

int c_fread_int_(int *i, int *ni,  FILE **fp){
    return fread(i, (size_t)sizeof(int), (size_t)*ni, *fp);
}

int c_fwrite_str_(char *str, int *size, FILE **fp, int maxsize){
    return fwrite(str, (size_t) *size, (size_t) 1, *fp);
}

int c_fwrite_str_array_(char *str, int *size, int *nelem, FILE **fp, int maxsize){
    if(*size == maxsize) {
	return fwrite(str, (size_t) *size, (size_t) *nelem, *fp);
    } else {
	int i,j=0,k;
	for(i=0; i<*nelem; i++){
	    k = fwrite(str+i*maxsize, (size_t) *size, (size_t) 1, *fp);
	    if(k==0) break;
	    j+=k;
	}
	return j;
    }
}

int c_fwrite_real_(float *f, int *nf, FILE **fp){
    return fwrite(f, (size_t)sizeof(float), (size_t)*nf, *fp);
}

int c_fwrite_double_(double *f, int *nf, FILE **fp){
    return fwrite(f, (size_t)sizeof(double), (size_t)*nf, *fp);
}

int c_fwrite_int_(int *i, int *ni,  FILE **fp){
    return fwrite(i, (size_t)sizeof(int), (size_t)*ni, *fp);
}
void c_fgets_(char *str, FILE **fp, int len){
    fgets(str, len, *fp);
}

void  c_fputs_(char *str, FILE **fp, int len){
    fputs(str, *fp);
}

void c_stdin_(FILE **fp){
    *fp = stdin;
}

void c_stdout_(FILE **fp){
    *fp = stdout;
}

void c_stderr_(FILE **fp){
    *fp = stderr;
}
