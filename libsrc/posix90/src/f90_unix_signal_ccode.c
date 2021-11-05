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

#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <errno.h>

#define SA_NOWRAPPER 128 /* see also f90_unix_signal_const.c */

typedef void (*sighandler_t)(int);
typedef void (*fsighandler_t)(int *);

static fsighandler_t sigtable[_NSIG];

void c_sigaction_compile_( fsighandler_t *s_handler,
			   sigset_t *s_mask, int *s_flags, 
			   fsighandler_t handler,
			   sigset_t *mask, int *flags,
			   int dmy1, int dmy2){
    *s_handler = handler;
    *s_mask = *mask;
    *s_flags = *flags;
    return;
}

void c_sa_handler_wrapper_(int sig){
    fsighandler_t handler;
    handler = sigtable[sig];
    (*handler)(&sig);
}

void c_sigaction_(int *sig,
		 fsighandler_t *handler, sigset_t *mask, int *flags,
		 fsighandler_t *oldhandler, sigset_t *oldmask, int *oldflags,
		  int dmy1, int dmy2){
    struct sigaction new_sig, old_sig;

    if( (*sig <0) || (*sig>=_NSIG)) {
	errno = EINVAL;
	return;
    }

    if(dmy1>0) {
	/* We got a handler. But we install the wrapper instead and store
	   it in sigtable[], except the SA_NOWRAPPER flag is set
	*/
	if( ((int)*handler)!=SIG_IGN && 
	    ((int)*handler)!= SIG_DFL && 
	    !( *flags&SA_NOWRAPPER)) {
	    /* We have a true handler (aka not SIG_IGN and not SIG_DFL
	       and the SA_NOWRAPPER flag is not set. We install the wrapper.
	    */
	    new_sig.sa_handler = &(c_sa_handler_wrapper_);
	}else{
	    /* We install *handler as is: It is either SIG_IGN, SIG_DFL or this was
	       requested by setting SA_NOWRAPPER
	    */
	    new_sig.sa_handler = *handler;
	}
	sigtable[*sig] = * handler;
	new_sig.sa_mask = *mask;
	new_sig.sa_flags = (*flags & (!SA_NOWRAPPER)); /* Keep the SA_NOWRAPPER bit 
							  to ourselves */
	sigaction(*sig, &new_sig, &old_sig);
    } else {
	/* We only want to read the old sig */
	sigaction(*sig, NULL, &old_sig);
    }
    /* Copy old_sig to my arguments */
    if (old_sig.sa_handler== &(c_sa_handler_wrapper_)){
	/* We found the wrapper, give the fortran function instead */
	*oldhandler = sigtable[*sig];
    } else {
	/* We found SIG_IGN, SIG_DFL or another routine. return this and set the
	   SA_NOWRAPPER flag.
	 */
	*oldhandler = old_sig.sa_handler;
	old_sig.sa_flags |= SA_NOWRAPPER;
    }

    *oldflags = old_sig.sa_flags;
    *oldmask = old_sig.sa_mask;
}

void c_sigemptyset_(sigset_t *set, int dmy){
     sigemptyset(set);
}

void c_sigfillset_(sigset_t *set, int dmy){
     sigfillset(set);
}

void c_sigaddset_(sigset_t *set, int *sig, int dmy){
     sigaddset(set, *sig);
}

void c_sigdelset_(sigset_t *set, int *sig, int dmy){
     sigdelset(set, *sig);
}

void c_sigismember_(sigset_t *set, int *sig, int dmy){
     sigismember(set, *sig);
}

void c_kill_(pid_t *pid, int *sig){
    kill(*pid, *sig);
}

void c_raise_(int *sig){
    raise(*sig);
}
