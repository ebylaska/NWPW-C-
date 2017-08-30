
#include        <time.h>

#define Clock   clock()


void    seconds(double *tt)
{
   *tt = Clock/((double) CLOCKS_PER_SEC);
}

