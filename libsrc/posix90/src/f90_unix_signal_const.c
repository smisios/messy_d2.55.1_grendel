/*

 f90_unix_signal_const.c : Determine values of constants for f90_unix_signal

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
#include <signal.h>

main(){

    printf("! DO NOT EDIT THIS FILE. \n");
    printf("! It is created from f90_unix_signal_const.c\n");
    printf("integer, parameter :: funcp_kind = %d\n", 
	   sizeof( int *));
    printf("type sigset_type ; character(len=%d) :: sigset ; end type\n", 
	   sizeof(sigset_t));
    printf("integer, parameter :: SIG_IGN = %d\n", SIG_IGN);
    printf("integer, parameter :: SIG_DFL = %d\n", SIG_DFL);
    printf("integer, parameter :: SIGHUP = %d\n", SIGHUP);
    printf("integer, parameter :: SIGINT = %d\n", SIGINT);
    printf("integer, parameter :: SIGQUIT = %d\n", SIGQUIT);
    printf("integer, parameter :: SIGILL = %d\n", SIGILL);
    printf("integer, parameter :: SIGTRAP = %d\n", SIGTRAP);
    printf("integer, parameter :: SIGABRT = %d\n", SIGABRT);
    printf("integer, parameter :: SIGIOT = %d\n", SIGIOT);
    printf("integer, parameter :: SIGBUS = %d\n", SIGBUS);
    printf("integer, parameter :: SIGFPE = %d\n", SIGFPE);
    printf("integer, parameter :: SIGKILL = %d\n", SIGKILL);
    printf("integer, parameter :: SIGUSR1 = %d\n", SIGUSR1);
    printf("integer, parameter :: SIGSEGV = %d\n", SIGSEGV);
    printf("integer, parameter :: SIGUSR2 = %d\n", SIGUSR2);
    printf("integer, parameter :: SIGPIPE = %d\n", SIGPIPE);
    printf("integer, parameter :: SIGALRM = %d\n", SIGALRM);
    printf("integer, parameter :: SIGTERM = %d\n", SIGTERM);
    printf("integer, parameter :: SIGSTKFLT = %d\n", SIGSTKFLT);
    printf("integer, parameter :: SIGCLD = %d\n", SIGCLD);
    printf("integer, parameter :: SIGCHLD = %d\n", SIGCHLD);
    printf("integer, parameter :: SIGCONT = %d\n", SIGCONT);
    printf("integer, parameter :: SIGSTOP = %d\n", SIGSTOP);
    printf("integer, parameter :: SIGTSTP = %d\n", SIGTSTP);
    printf("integer, parameter :: SIGTTIN = %d\n", SIGTTIN);
    printf("integer, parameter :: SIGTTOU = %d\n", SIGTTOU);
    printf("integer, parameter :: SIGURG = %d\n", SIGURG);
    printf("integer, parameter :: SIGXCPU = %d\n", SIGXCPU);
    printf("integer, parameter :: SIGXFSZ = %d\n", SIGXFSZ);
    printf("integer, parameter :: SIGVTALRM = %d\n", SIGVTALRM);
    printf("integer, parameter :: SIGPROF = %d\n", SIGPROF);
    printf("integer, parameter :: SIGWINCH = %d\n", SIGWINCH);
    printf("integer, parameter :: SIGPOLL = %d\n", SIGPOLL);
    printf("integer, parameter :: SIGIO = %d\n", SIGIO);
    printf("integer, parameter :: SIGPWR = %d\n", SIGPWR);
    printf("integer, parameter :: SIGSYS = %d\n", SIGSYS);
    /* printf("integer, parameter :: SIGUNUSED = %d\n", SIGUNUSED); */
    printf("integer, parameter :: NSIG = %d\n", _NSIG);
    printf("integer, parameter :: SA_NOCLDSTOP = %d\n", SA_NOCLDSTOP);
    printf("integer, parameter :: SA_NOCLDWAIT = %d\n", SA_NOCLDWAIT);
    printf("integer, parameter :: SA_RESETHAND = %d\n", SA_RESETHAND);
    printf("integer, parameter :: SA_ONSTACK = %d\n", SA_ONSTACK);
    printf("integer, parameter :: SA_RESTART = %d\n", SA_RESTART);
    printf("integer, parameter :: SA_NODEFER = %d\n", SA_NODEFER);
    printf("!integer, parameter :: SA_SIGINFO = %d\n", SA_SIGINFO); /* Do not pass this flag */
    printf("integer, parameter :: SA_NOWRAPPER = %d\n", 128); /* see also f90_unix_signal_ccode.c */
    printf("integer, parameter :: SIG_BLOCK = %d\n",SIG_BLOCK);
    printf("integer, parameter :: SIG_UNBLOCK = %d\n",SIG_UNBLOCK);
    printf("integer, parameter :: SIG_SETMASK = %d\n",SIG_SETMASK);

    return 0;
}
 
