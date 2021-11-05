
/* If WALLCLOCK_TIME is defined,  measure wall clock time, 
else measure CPU time */

#define WALLCLOCK_TIME 1


#include <time.h>
#include "mmdc_util.h"

#ifdef WALLCLOCK_TIME
#include <sys/timeb.h>
  static struct timeb  tp_first;
#endif

/* The following function is used by client and server */

double MMDc_U_Time () {

#ifdef WALLCLOCK_TIME
  double               ft,ft_sec,ft_ms;
  struct timeb         tp;
  static int           first=1;
   
  if(first == 1) {
    (void) ftime (&tp_first);
    first = 0;
  }

  (void) ftime (&tp);


  ft_sec = (double) (tp.time - tp_first.time);
  ft_ms  = (double) tp.millitm - (double) tp_first.millitm;
  ft     = ft_sec + (ft_ms/1000.);

  /*  fprintf (stderr,"Time %f %f %f %d %d  \n", ft, ft_sec, ft_ms, tp.time, tp.millitm);*/
#else

  clock_t        t;
  double         ft;
  struct timeb   *tp;


  t = clock();
  ft = (double) t / (double) CLOCKS_PER_SEC;

#endif

  return ft;
}
