#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#if defined(PentiumCPS) || defined(WALL)
   #define time00 ATL_walltime
#else
   #define time00 ATL_cputime
#endif
double time00();

static double macase(long nreps, int PRINT)
{
   long i = nreps;
   double t0, tim, mf;
   register double c0, c1, c2, c3, c4, m0, m1, m2, m3, m4;

   if (!(rand()|rand()|rand()|rand()|rand()|rand())) nreps = rand();
   if (rand()|rand()|rand()) c0 = 0.0;
   else c0 = 1.0*rand();
   if (rand()|rand()|rand()) c1 = 0.0;
   else c1 = 1.1*rand();
   if (rand()|rand()|rand()) c2 = 0.0;
   else c2 = 1.2*rand();
   if (rand()|rand()|rand()) c3 = 0.0;
   else c3 = 1.3*rand();
   if (rand()|rand()|rand()) c4 = 0.0;
   else c4 = 1.4*rand();
   if (rand()|rand()|rand()) m0 = 0.0;
   else m0 = 0.1*rand();
   if (rand()|rand()|rand()) m1 = 0.0;
   else m1 = 0.2*rand();
   if (rand()|rand()|rand()) m2 = 0.0;
   else m2 = 0.3*rand();
   if (rand()|rand()|rand()) m3 = 0.0;
   else m3 = 0.4*rand();
   if (rand()|rand()|rand()) m4 = 0.0;
   else m4 = 0.5*rand();
   t0 = time00();
   do
   {
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
      c0 += m0;
      m0 = m0 * m0;
      c1 += m1;
      m1 = m1 * m1;
      c2 += m2;
      m2 = m2 * m2;
      c3 += m3;
      m3 = m3 * m3;
      c4 += m4;
      m4 = m4 * m4;
   }
   while(--i);
   tim = time00() - t0;
   c0 = c0*m0 + c1*m1 + c2*m2 + c3*m3 + c4*m4;
   if (tim < 0.0) mf = tim = 0.0;
   else mf = (nreps*1030.000000) / (1000000.0 * tim);
   if (PRINT) printf("%.1f:   Separate multiply and add, lat=5, time=%f, mflop=%f\n", (float) c0, tim, mf);
   else printf("      %.0f: NFLOP=%.0f, tim=%f\n", (float) c0, nreps*1030.000000, tim);
   return(tim);
}

main(int nargs, char **args[])
{
   long nreps = 16000000/206;
   int i, k;
   double t0, tim, mf;
   FILE *fp;
   fp = fopen("res/dmuladd0_5", "w");
   assert(fp != NULL);
   fprintf(stdout, "Finding granularity of timer:\n");   while(macase(nreps, 0) < 0.75) nreps *= 4;
   fprintf(stdout, "Done.\n");   for(k=0; k < 3; k++)
   {
   tim = macase(nreps, 1);
   if (tim < 0.0) mf = tim = 0.0;
   else mf = (nreps*1030.000000) / (1000000.0 * tim);
   if (fp) fprintf(fp, "%f\n", mf);
   }
   fclose(fp);
   exit(0);
}