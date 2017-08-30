/*$Id: util_message.c */

#include <sys/types.h>
#include <time.h>


char *util_date()
{
  time_t t = time((time_t *) 0);
  char *tmp = ctime(&t);
  return tmp;
}

char *util_message(const int n)
{
   
}
