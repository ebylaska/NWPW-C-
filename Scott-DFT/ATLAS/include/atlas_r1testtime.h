#ifndef ATLAS_R1TESTTIME_H
   #define ATLAS_R1TESTTIME_H

#include "atlas_r1parse.h"
#include "atlas_gentesttime.h"


static char *GetAlignStr(int alignX, int alignY, int alignA)
{
   static char ln[256];
   int i=0;

   ln[0] = '\0';
   if (alignX)
      i += sprintf(ln, "-FAx %d ", alignX);
   if (alignY)
      i += sprintf(ln+i, "-FAy %d ", alignY);
   if (alignA)
      i += sprintf(ln+i, "-FAa %d ", alignA);
   return(ln);
}

static int R1KernelFailsTest
   (int verb, char pre, ATL_INT M, ATL_INT N, ATL_INT lda, ATL_r1node_t *kn)
{
   char ln[1024];
   int i;
/*
 * If the file is generated, call generator to create it
 */
   if (kn->genstr)
   {
      if (system(kn->genstr))
      {
         fprintf(stderr, "ERROR, LINE %d of %s\n", __LINE__, __FILE__);
         fprintf(stderr, "UNABLE TO GENERATE WITH COMMAND: %s\n", kn->genstr);
         exit(-1);
      }
   }
   assert(kn->rout);
   assert (M >= kn->minX);
   assert (N >= kn->minY);
   i = sprintf(ln, "make %cr1ktest r1rout=%s align=\"%s\" ", pre, kn->rout,
               GetAlignStr(kn->alignX, kn->alignY, kn->alignA));
   if (kn->comp)
      i += sprintf(ln+i, "%cR1CC=\"%s\" ", pre, kn->comp);
   if (kn->cflags)
      i += sprintf(ln+i, "%cR1FLAGS=\"%s\" ", pre, kn->cflags);
   i += sprintf(ln+i, "Mt=%d Nt=%d ldat=%d", M, N, lda);
   if (verb < 3)
      i += sprintf(ln+i, " > /dev/null 2>&1\n");
   else
      i += sprintf(ln+i, "\n");
   if (verb > 1)
      fprintf(stdout, "system call:%s\n", ln);
   i = system(ln);
   if (verb)
   {
      if (i)
         fprintf(stderr, "%s(ID=%d) FAILS TESTER!!\n", kn->rout,kn->ID);
      else
         fprintf(stderr, "%s(ID=%d) *PASSES* TESTER!!\n", kn->rout,kn->ID);
   }
   return(i);
}


static char *GetResIdStr(ATL_r1node_t *r1p, ATL_INT M, ATL_INT N,
                         ATL_INT lda, ATL_INT percL1, int mflop)
{
   /* <ID>_<M>x<N>_<lda>_<percL1>_a<alignA>x<aX>x<aY> */
   static char ln[512];
   sprintf(ln, "%d_%dx%d_%d_%d_a%dx%dx%d", r1p->ID, M, N, lda, percL1,
           r1p->alignA, r1p->alignX, r1p->alignY);
   return(ln);
}

static double TimeR1Kernel
(int verb,              /* 0: no output, 1 min ouput, 2: full output */
 int FORCETIME,         /* if nonzero, ignore existing timing file */
 ATL_r1node_t *r1p,     /* ptr to kernel structure */
 char pre,              /* precision prefix */
 ATL_INT M, ATL_INT N,  /* dimensions to time */
 ATL_INT lda,           /* stride between row elements */
 ATL_INT percL1,        /* if 0, time kernel directly wt no blocking */
                        /* if non-zero, block for that % of L1 cache size */
 int mflop,             /* force mflop flops in each timing interval */
 int cflush             /* if >= 0, size of cache flush area, else ignored */
)
{
   char ln[2048], resf[256];
   double *dp, mf;
   int i;
/*
 * If the file is generated, call generator to create it
 */
   if (r1p->genstr)
   {
      if (system(r1p->genstr))
      {
         fprintf(stderr, "ERROR, LINE %d of %s\n", __LINE__, __FILE__);
         fprintf(stderr, "UNABLE TO GENERATE WITH COMMAND: %s\n", r1p->genstr);
         exit(-1);
      }
   }

   if (r1p->minY)
      N = Mmax(N, r1p->minY);
   if (r1p->minX)
      M = Mmax(M, r1p->minX);
   i = r1p->ldamul / pre2size(pre);
   lda = (i) ? ((lda+i-1)/i)*i : lda;

   sprintf(resf, "res/%cr1%s", pre, GetResIdStr(r1p, M, N, lda, percL1, mflop));
   dp = FORCETIME ? NULL : ReadResultsFile(0, resf);
   if (dp)
   {
      if (verb > 0)
         fprintf(stdout, "   %d:%s gets %.2f MFLOPS\n",
                 r1p->ID, r1p->rout, *dp);
      return(*dp);
   }

   if (percL1)
      i = sprintf(ln, "make %cr1time M=%d N=%d lda=%d l1mul=%d r1rout=\"%s\"",
                  pre, M, N, lda, percL1, r1p->rout);
   else
      i = sprintf(ln, "make %cr1ktime M=%d N=%d lda=%d r1rout=\"%s\"",
                  pre, M, N, lda, r1p->rout);
   if (r1p->comp)
      i += sprintf(ln+i, " %dR1CC=\"%s\"", pre, r1p->comp);
   if (r1p->cflags)
      i += sprintf(ln+i, " %dR1CFLAGS=\"%s\"", pre, r1p->cflags);
   if (r1p->alignA || r1p->alignX || r1p->alignY)
   {
      i += sprintf(ln+i, " align=\"");
      if (r1p->alignA)
         i += sprintf(ln+i, " -Fa %d", r1p->alignA);
      if (r1p->alignY)
         i += sprintf(ln+i, " -Fy %d", r1p->alignY);
      if (r1p->alignX)
         i += sprintf(ln+i, " -Fx %d", r1p->alignX);
      if (mflop >=0)
         i += sprintf(ln+i, " -F %d", mflop);
      i += sprintf(ln+i, "\"");
   }
   i += sprintf(ln+i, " tflags=\"-f %s", resf);
   if (mflop >= 0)
      i += sprintf(ln+i, " -F %d", mflop);
   if (cflush >=0)
      i += sprintf(ln+i, " -C %d", cflush);
   i += sprintf(ln+i, "\"");
   if (verb < 3)
      i += sprintf(ln+i, " > /dev/null 2>&1\n");
   else
      i += sprintf(ln+i, "\n");
   if (system(ln))
   {
      fprintf(stderr, "\nERROR, LINE %d OF %s\n", __LINE__, __FILE__);
      fprintf(stderr, "SYSTEM CALL FAILED: %s\n", ln);
      exit(-1);
   }
   if (verb > 1)
   {
      dp = ReadResultsFile(1, resf);
      mf = PrintResultsFromFile(stdout, dp);
      free(dp);
      dp = &mf;
   }
   else
      dp = ReadResultsFile(0, resf);
   assert(dp);
   if (verb == 1)
      fprintf(stdout, "   %d:%s gets %.2f MFLOPS\n", r1p->ID, r1p->rout, *dp);
   return(*dp);
}

#endif  /* end guard around atlas_r1testtime.h */
