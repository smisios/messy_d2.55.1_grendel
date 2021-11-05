/*
 * drop-in replacement for fortran 2003 procedure
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cfortran.h"

static void
die(char *msg)
{
  perror(msg);
  exit(EXIT_FAILURE);
}

static void
getenvWrapper(const char *name, char *val, int *length, int *status,
                     int trim_name, int trailing_blanks)
{
  char *name_copy;
  const char *env_value;
  /* without trimming we need to re-add the trailing blanks stripped
   * by cfortran */
  if (!trim_name)
  {
    int name_len = strlen(name);
    name_copy = malloc(name_len + trailing_blanks + 1);
    if (!name_copy)
      die("Failed to copy input string");
    strcpy(name_copy, name);
    memset(name_copy + name_len, ' ', trailing_blanks);
    name_copy[name_len + trailing_blanks] = '\0';
    name = name_copy;
  }
  env_value = getenv(name);
  strncpy(val, env_value, *length);
  *length = strlen(env_value);
  *status = env_value != NULL;
  if (!trim_name)
    free(name_copy);
}


FCALLSCSUB6(getenvWrapper,GETENV_WRAPPER,getenv_wrapper,\
            STRING,PSTRING,PINT,PINT,LOGICAL,INT)

